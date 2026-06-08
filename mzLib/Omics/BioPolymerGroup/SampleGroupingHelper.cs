using MassSpectrometry;
using Omics.SpectralMatch;

namespace Omics.BioPolymerGroup;

/// <summary>
/// Static helper for building <see cref="SampleGroupResult"/> collections from PSMs,
/// experimental design, and optional intensity data. Encapsulates the logic for grouping
/// label-free samples by (Condition, BiologicalReplicate) and isobaric samples by
/// (File, Channel), as well as the fallback file-only grouping when no experimental
/// design is available.
/// </summary>
internal static class SampleGroupingHelper
{
    /// <summary>
    /// Builds a list of <see cref="SampleGroupResult"/> objects from the provided PSMs,
    /// samples, and intensities, then invokes the caller-supplied occupancy callback for each result.
    /// </summary>
    /// <param name="allPsms">All PSMs to distribute into sample groups.</param>
    /// <param name="samplesForQuantification">
    /// Experimental design samples. May contain <see cref="SpectraFileInfo"/> (label-free)
    /// and/or <see cref="IsobaricQuantSampleInfo"/> (isobaric). Null or empty triggers
    /// file-only grouping.
    /// </param>
    /// <param name="intensitiesBySample">
    /// Optional per-sample intensity values. When non-null, results are created with
    /// <see cref="SampleGroupResult.IntensitiesBySample"/> populated from the matching samples.
    /// </param>
    /// <param name="occupancyCallback">
    /// Action invoked for each created <see cref="SampleGroupResult"/> with the PSMs that
    /// belong to that group. The callback is responsible for populating occupancy data
    /// (e.g., via <see cref="ModificationOccupancyCalculator"/>).
    /// </param>
    /// <returns>List of <see cref="SampleGroupResult"/> objects in deterministic order.</returns>
    public static List<SampleGroupResult> BuildSampleGroupResults(
        HashSet<ISpectralMatch> allPsms,
        List<ISampleInfo>? samplesForQuantification,
        Dictionary<ISampleInfo, double>? intensitiesBySample,
        Action<SampleGroupResult, List<ISpectralMatch>> occupancyCallback)
    {
        var results = new List<SampleGroupResult>();

        var spectraFiles = samplesForQuantification?.OfType<SpectraFileInfo>().ToList() ?? [];
        var isobaricSamples = samplesForQuantification?.OfType<IsobaricQuantSampleInfo>().ToList() ?? [];

        if (spectraFiles.Count > 0)
        {
            bool unfractionated = spectraFiles.Select(p => p.Fraction).Distinct().Count() == 1;
            bool conditionsUndefined = spectraFiles.All(p => string.IsNullOrEmpty(p.Condition));
            bool silacExperimentalDesign = spectraFiles.Any(p => !File.Exists(p.FullFilePathWithExtension));

            foreach (var conditionGroup in spectraFiles.GroupBy(p => p.Condition))
            {
                foreach (var bioRepGroup in conditionGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    var filesInGroup = bioRepGroup.ToList();
                    string label = (conditionsUndefined && unfractionated) || silacExperimentalDesign
                        ? filesInGroup.First().FilenameWithoutExtension
                        : $"{conditionGroup.Key}_{bioRepGroup.Key + 1}";

                    var filePaths = new HashSet<string>(filesInGroup.Select(f => f.FullFilePathWithExtension));
                    var psmsInGroup = allPsms
                        .Where(p => filePaths.Contains(p.FullFilePath))
                        .ToList();

                    // Create SampleGroupResult with per-sample intensities if available.
                    // Otherwise, create with empty intensities (HasIntensityData = false) for spectral counting.
                    var resultIntensities = new Dictionary<string, double>();
                    SampleGroupResult result;
                    if (intensitiesBySample != null)
                    {
                        foreach (var file in filesInGroup)
                        {
                            if (intensitiesBySample.TryGetValue(file, out var fileIntensity))
                                resultIntensities[file.FilenameWithoutExtension] = fileIntensity;
                        }

                        result = new SampleGroupResult(conditionGroup.Key, bioRepGroup.Key)
                        {
                            Label = label,
                            SpectralCount = psmsInGroup.Count,
                            FilesInGroup = filesInGroup.ToDictionary(kvp => kvp.FilenameWithoutExtension, kvp => (ISampleInfo)kvp),
                            IntensitiesBySample = resultIntensities
                        };
                    }
                    else
                    {
                        result = new SampleGroupResult(conditionGroup.Key, bioRepGroup.Key)
                        {
                            Label = label,
                            SpectralCount = psmsInGroup.Count,
                            FilesInGroup = filesInGroup.ToDictionary(kvp => kvp.FilenameWithoutExtension, kvp => (ISampleInfo)kvp)
                            // IntensitiesBySample left empty → HasIntensityData = false
                        };
                    }

                    occupancyCallback(result, psmsInGroup);
                    results.Add(result);
                }
            }
        }
        else if (isobaricSamples.Count > 0)
        {
            foreach (var fileGroup in isobaricSamples.GroupBy(p => p.FullFilePathWithExtension).OrderBy(g => g.Key))
            {
                var psmsInFile = allPsms
                    .Where(p => p.FullFilePath.Equals(fileGroup.Key))
                    .ToList();

                foreach (var sample in fileGroup.OrderBy(p => p.ChannelLabel))
                {
                    string label = $"{Path.GetFileNameWithoutExtension(sample.FullFilePathWithExtension)}_{sample.ChannelLabel}";

                    // Build per-channel intensity lookup for this result
                    SampleGroupResult result;
                    if (intensitiesBySample != null && intensitiesBySample.TryGetValue(sample, out var channelIntensity))
                    {
                        result = new SampleGroupResult(sample.Condition, sample.BiologicalReplicate)
                        {
                            Label = label,
                            SpectralCount = psmsInFile.Count,
                            FilesInGroup = new Dictionary<string, ISampleInfo> { { label, sample } },
                            IntensitiesBySample = new Dictionary<string, double> { { label, channelIntensity } }
                        };
                    }
                    else
                    {
                        result = new SampleGroupResult(sample.Condition, sample.BiologicalReplicate)
                        {
                            Label = label,
                            SpectralCount = psmsInFile.Count,
                            FilesInGroup = new Dictionary<string, ISampleInfo> { { label, sample } }
                            // IntensitiesBySample left empty → HasIntensityData = false
                        };
                    }

                    occupancyCallback(result, psmsInFile);
                    results.Add(result);
                }
            }
        }
        else
        {
            // No experimental design — group PSMs by source file for count-only results
            foreach (var fileGroup in allPsms.GroupBy(p => p.FullFilePath).OrderBy(g => g.Key))
            {
                var psmsInFile = fileGroup.ToList();
                string label = Path.GetFileNameWithoutExtension(fileGroup.Key);

                var result = new SampleGroupResult(string.Empty, 0)
                {
                    Label = label,
                    SpectralCount = psmsInFile.Count
                    // FilesInGroup and IntensitiesByFile left empty → HasIntensityData = false
                };

                occupancyCallback(result, psmsInFile);
                results.Add(result);
            }
        }

        return results;
    }
}
