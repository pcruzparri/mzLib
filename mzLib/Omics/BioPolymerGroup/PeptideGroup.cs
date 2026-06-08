using MassSpectrometry;
using Omics.BioPolymer;
using Omics.SpectralMatch;
using System.Text;

namespace Omics.BioPolymerGroup;

/// <summary>
/// Represents a group of PSMs sharing the same peptide or oligonucleotide sequence
/// (same <see cref="BaseSequence"/>, <see cref="OneBasedStartResidue"/>, <see cref="OneBasedEndResidue"/>)
/// from one or more parent biopolymers. Used for peptide-level modification occupancy calculation
/// and future quantified peptide group reporting.
/// </summary>
public class PeptideGroup : IEquatable<PeptideGroup>
{
    #region Identity

    /// <summary>The base sequence of this peptide (e.g., "ACDEF").</summary>
    public string BaseSequence { get; }

    /// <summary>One-based start residue on the parent biopolymer.</summary>
    public int OneBasedStartResidue { get; }

    /// <summary>One-based end residue on the parent biopolymer.</summary>
    public int OneBasedEndResidue { get; }

    /// <summary>
    /// Parent biopolymers this peptide originates from.
    /// A peptide can map to multiple parents (shared peptides).
    /// </summary>
    public List<IBioPolymer> Parents { get; }

    /// <summary>
    /// Display name for the group: "BaseSequence_Start_End_Parent1|Parent2|...".
    /// </summary>
    public string PeptideGroupName { get; }

    #endregion

    #region PSMs

    /// <summary>
    /// All PSMs for this peptide that pass the 1% FDR threshold.
    /// </summary>
    public HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

    #endregion

    #region Quantification

    /// <summary>
    /// Samples that contribute quantification data for this peptide.
    /// </summary>
    private List<ISampleInfo>? _samplesForQuantification;
    public List<ISampleInfo>? SamplesForQuantification
    {
        get => _samplesForQuantification;
        set
        {
            _samplesForQuantification = value;
            SampleGroupResults = null;
        }
    }

    /// <summary>
    /// Measured intensity values for this peptide, keyed by sample.
    /// </summary>
    private Dictionary<ISampleInfo, double>? _intensitiesBySample;
    public Dictionary<ISampleInfo, double>? IntensitiesBySample
    {
        get => _intensitiesBySample;
        set
        {
            _intensitiesBySample = value;
            SampleGroupResults = null;
        }
    }

    /// <summary>
    /// Per-sample-group quantification and modification occupancy results.
    /// </summary>
    public List<SampleGroupResult>? SampleGroupResults { get; set; }

    #endregion

    #region Flags (derived from Parents)

    public bool IsDecoy => Parents.Any(p => p.IsDecoy);
    public bool IsContaminant => Parents.Any(p => p.IsContaminant);
    public bool IsEntrapment => Parents.Any(p => p.IsEntrapment);

    #endregion

    #region Occupancy

    /// <summary>
    /// Intensity rollup strategy for modification occupancy calculations.
    /// Defaults to <see cref="IntensityRollupStrategy.Sum"/> for backward compatibility.
    /// Changing the value invalidates <see cref="SampleGroupResults"/>.
    /// </summary>
    public IntensityRollupStrategy OccupancyRollupStrategy
    {
        get => _occupancyRollupStrategy;
        set
        {
            _occupancyRollupStrategy = value;
            SampleGroupResults = null;
        }
    }
    private IntensityRollupStrategy _occupancyRollupStrategy = IntensityRollupStrategy.Sum;

    /// <summary>
    /// Populates <see cref="SampleGroupResults"/> from PSMs and quantification data.
    /// </summary>
    public void PopulateSampleGroupResults()
    {
        SampleGroupResults = SampleGroupingHelper.BuildSampleGroupResults(
            AllPsmsBelowOnePercentFDR,
            SamplesForQuantification,
            IntensitiesBySample,
            PopulateOccupancy);
    }

    private void PopulateOccupancy(SampleGroupResult result, List<ISpectralMatch> psms)
    {
        var occupancy = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy(
            psms, OccupancyRollupStrategy);

        if (occupancy.Count > 0)
            result.ParentOccupancy[BaseSequence] = occupancy;
    }

    #endregion

    #region FDR (future — add when needed)

    // public double QValue { get; set; }
    // public double Score { get; set; }
    // public int CumulativeTarget { get; set; }
    // public int CumulativeDecoy { get; set; }

    #endregion

    #region Output

    public string GetTabSeparatedHeader()
    {
        if (SampleGroupResults is null) PopulateSampleGroupResults();

        var fields = new List<string>
        {
            "Peptide Sequence",
            "Start Residue",
            "End Residue",
            "Parent Accessions"
        };

        foreach (var group in SampleGroupResults!)
        {
            fields.Add($"SpectralCount_{group.Label}");
            fields.Add($"CountOccupancy_{group.Label}");
            if (group.HasIntensityData)
            {
                string prefix = OccupancyRollupStrategy switch
                {
                    IntensityRollupStrategy.Mean => "MeanIntensityOccupancy",
                    IntensityRollupStrategy.Median => "MedianIntensityOccupancy",
                    IntensityRollupStrategy.Sum => "IntensityOccupancy",
                    _ => throw new ArgumentOutOfRangeException(nameof(OccupancyRollupStrategy), OccupancyRollupStrategy, $"Unsupported rollup strategy: {OccupancyRollupStrategy}")
                };
                fields.Add($"{prefix}_{group.Label}");
            }
        }

        return string.Join("\t", fields);
    }

    public override string ToString()
    {
        if (SampleGroupResults is null) PopulateSampleGroupResults();

        var fields = new List<string>
        {
            BaseSequence,
            OneBasedStartResidue.ToString(),
            OneBasedEndResidue.ToString(),
            string.Join("|", Parents.Select(p => p.Accession))
        };

        foreach (var group in SampleGroupResults!)
        {
            fields.Add(group.SpectralCount.ToString());

            var occupancy = group.ParentOccupancy.TryGetValue(BaseSequence, out var pos) ? pos : null;
            if (occupancy != null)
            {
                fields.Add(FormatOccupancy(occupancy, intensityBased: false));
                if (group.HasIntensityData)
                    fields.Add(FormatOccupancy(occupancy, intensityBased: true));
            }
            else
            {
                fields.Add(string.Empty);
                if (group.HasIntensityData) fields.Add(string.Empty);
            }
        }

        return string.Join("\t", fields);
    }

    private static string FormatOccupancy(
        Dictionary<int, List<SiteSpecificModificationOccupancy>> occupancy,
        bool intensityBased)
    {
        var parts = new List<string>();
        foreach (var kvp in occupancy.OrderBy(k => k.Key))
        {
            foreach (var site in kvp.Value)
            {
                parts.Add(site.ToModInfoString(intensityBased));
            }
        }
        return string.Join(";", parts);
    }

    #endregion

    #region Subsetting

    /// <summary>
    /// Creates a new peptide group containing only PSMs and quantification data from a specific spectra file.
    /// </summary>
    /// <param name="fullFilePath">The full path to the spectra file to filter by.</param>
    /// <returns>A new <see cref="PeptideGroup"/> containing only data from the specified file.</returns>
    public PeptideGroup ConstructSubsetPeptideGroup(string fullFilePath)
    {
        var psmsInFile = AllPsmsBelowOnePercentFDR
            .Where(p => p.FullFilePath == fullFilePath)
            .ToHashSet();

        var subset = new PeptideGroup(BaseSequence, OneBasedStartResidue, OneBasedEndResidue, Parents, psmsInFile)
        {
            OccupancyRollupStrategy = OccupancyRollupStrategy
        };

        if (SamplesForQuantification != null)
        {
            subset.SamplesForQuantification = SamplesForQuantification
                .Where(s => s.FullFilePathWithExtension == fullFilePath)
                .ToList();
        }

        if (IntensitiesBySample != null && subset.SamplesForQuantification != null)
        {
            subset.IntensitiesBySample = IntensitiesBySample
                .Where(kvp => subset.SamplesForQuantification.Contains(kvp.Key))
                .ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
        }

        return subset;
    }

    #endregion

    #region Constructor

    public PeptideGroup(
        string baseSequence,
        int oneBasedStartResidue,
        int oneBasedEndResidue,
        IEnumerable<IBioPolymer> parents,
        IEnumerable<ISpectralMatch> psms)
    {
        BaseSequence = baseSequence ?? throw new ArgumentNullException(nameof(baseSequence));
        OneBasedStartResidue = oneBasedStartResidue;
        OneBasedEndResidue = oneBasedEndResidue;
        Parents = parents?.ToList() ?? throw new ArgumentNullException(nameof(parents));
        AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>(psms ?? throw new ArgumentNullException(nameof(psms)));
        PeptideGroupName = $"{baseSequence}_{oneBasedStartResidue}_{oneBasedEndResidue}_{string.Join("|", Parents.OrderBy(p => p.Accession).Select(p => p.Accession))}";
    }

    #endregion

    #region Equality

    public bool Equals(PeptideGroup? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        return PeptideGroupName == other.PeptideGroupName;
    }

    public override bool Equals(object? obj) => Equals(obj as PeptideGroup);

    public override int GetHashCode() => PeptideGroupName.GetHashCode();

    #endregion
}
