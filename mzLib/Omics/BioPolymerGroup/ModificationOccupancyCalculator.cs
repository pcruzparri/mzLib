using MzLibUtil;
using Omics.BioPolymer;
using Omics.Modifications;
using Omics.SpectralMatch;

namespace Omics.BioPolymerGroup;

/// <summary>
/// Calculates modification occupancy/stoichiometry from identified peptides.
/// Supports both count-based and intensity-based metrics, at the protein or peptide level.
/// </summary>
public static class ModificationOccupancyCalculator
{
    /// <summary>
    /// Mod types to exclude from occupancy calculations.
    /// </summary>
    private static readonly string[] ExcludedModTypes = ["Common Variable", "Common Fixed"];

    /// <summary>
    /// Location restrictions to exclude (peptide-terminal, not protein-terminal).
    /// </summary>
    private static readonly string[] ExcludedLocations = ["NPep", "PepC"];

    /// <summary>
    /// Calculates per-site modification occupancy mapped to protein coordinates directly from PSMs.
    /// PSM grouping, form filtering, TotalCount derivation, and intensity lookup are all handled internally.
    /// </summary>
    /// <param name="bioPolymer">The parent biopolymer whose length defines the coordinate space.</param>
    /// <param name="psms">
    /// All PSMs to consider. Forms are filtered to <paramref name="bioPolymer"/> internally.
    /// PSMs whose <see cref="ISpectralMatch.Intensities"/> is a single-element array contribute
    /// to intensity-based stoichiometry; others contribute only to count-based metrics.
    /// </param>
    /// <param name="strategy">Intensity rollup strategy for computing intensity-based stoichiometry.
    /// Defaults to <see cref="IntensityRollupStrategy.Sum"/> for backward compatibility.</param>
    public static Dictionary<int, List<SiteSpecificModificationOccupancy>> CalculateParentLevelOccupancy(
        IBioPolymer bioPolymer,
        IEnumerable<ISpectralMatch> psms,
        IntensityRollupStrategy strategy = IntensityRollupStrategy.Sum)
    {
        var psmList = psms as IList<ISpectralMatch> ?? psms.ToList();

        // Pre-compute the matching form for each PSM by position.
        var psmForms = psmList
            .Select(p => p.GetIdentifiedBioPolymersWithSetMods()
                .FirstOrDefault(s => s.FullSequence != null
                    && s.BaseSequence == p.BaseSequence
                    && s.FullSequence == p.FullSequence
                    && s.Parent.Accession == bioPolymer.Accession))
            .ToArray();

        // Pass 1: accumulate totals per position.
        var positionTotals = new Dictionary<int, (int totalCount, double totalIntensitySum)>();
        var positionIntensities = strategy == IntensityRollupStrategy.Sum
            ? null
            : new Dictionary<int, List<double>>();

        for (int j = 0; j < psmList.Count; j++)
        {
            var psm = psmList[j];
            var sequence = psmForms[j];
            if (sequence is null) // PSM for this protein might be ambiguous (e.g. missing full sequence)
            {
                try
                {
                    // Still want to count it toward TotalCount/TotalIntensity for any positions it covers,
                    // so find the best-matching form without the full sequence requirement.
                    sequence = psm.GetIdentifiedBioPolymersWithSetMods()
                        .FirstOrDefault(s => s.BaseSequence == psm.BaseSequence
                            && s.Parent.Accession == bioPolymer.Accession);
                }
                catch (Exception)
                {
                    continue; // If we can't find any form for this PSM, skip it entirely.
                }
            }

            if (sequence is null) // No form found for this PSM, skip it entirely.
                continue;

            int rangeStart = sequence.OneBasedStartResidue + (sequence.OneBasedStartResidue == 1 ? 0 : 1); // Include position 1 if sequence starts at the protein N-terminus
            int rangeEnd = sequence.OneBasedEndResidue + (sequence.OneBasedEndResidue == bioPolymer.Length ? 2 : 1); // Include last position if sequence ends at the protein C-terminus
            for (int i = rangeStart; i <= rangeEnd; i++)
            {
                if (!positionTotals.ContainsKey(i))
                    positionTotals[i] = (0, 0.0);
                var totals = positionTotals[i];
                totals.totalCount++;
                if (psm.Intensities is { Length: 1 })
                {
                    totals.totalIntensitySum += psm.Intensities[0];
                    positionIntensities?.AddToList(i, psm.Intensities[0]);
                }
                positionTotals[i] = totals;
            }
        }

        // Pass 2: build occupancy entries per modification site.
        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();
        var modifiedIntensities = strategy == IntensityRollupStrategy.Sum
            ? null
            : new Dictionary<(int position, string modId), List<double>>();

        for (int j = 0; j < psmList.Count; j++)
        {
            var psm = psmList[j];
            var sequence = psmForms[j];
            if (sequence is null)  // PSM has no form for this protein, skip
                continue;

            foreach (var mod in sequence.AllModsOneIsNterminus)
            {
                if (IsExcludedMod(mod.Value))
                    continue;

                if (!TryGetProteinPosition(mod, sequence, bioPolymer.Length, out int indexInProtein))
                    continue;

                if (!working.TryGetValue(indexInProtein, out var modsAtPosition))
                {
                    modsAtPosition = new Dictionary<string, SiteSpecificModificationOccupancy>();
                    working[indexInProtein] = modsAtPosition;
                }

                if (!modsAtPosition.ContainsKey(mod.Value.IdWithMotif))
                {
                    if (!positionTotals.TryGetValue(indexInProtein, out var posTotals))
                        continue;

                    modsAtPosition[mod.Value.IdWithMotif] = new SiteSpecificModificationOccupancy(indexInProtein, mod.Value.IdWithMotif)
                    {
                        TotalCount = posTotals.totalCount,
                        TotalIntensity = posTotals.totalIntensitySum
                    };
                }

                var siteOcc = modsAtPosition[mod.Value.IdWithMotif];
                siteOcc.ModifiedCount++;
                if (psm.Intensities is { Length: 1 })
                {
                    siteOcc.ModifiedIntensity += psm.Intensities[0];
                    modifiedIntensities?.AddToList((indexInProtein, mod.Value.IdWithMotif), psm.Intensities[0]);
                }
            }
        }

        // Post-process: apply rollup strategy for non-Sum strategies.
        if (strategy != IntensityRollupStrategy.Sum && working.Count > 0)
        {
            foreach (var kvp in working)
            {
                int position = kvp.Key;
                var allIntensities = positionIntensities?.GetValueOrDefault(position) ?? new List<double>();

                foreach (var siteKvp in kvp.Value)
                {
                    var site = siteKvp.Value;
                    var modIntensities = modifiedIntensities?.GetValueOrDefault((position, site.ModificationIdWithMotif)) ?? new List<double>();

                    site.TotalIntensity = ApplyStrategy(allIntensities, strategy);
                    site.ModifiedIntensity = ApplyStrategy(modIntensities, strategy);
                }
            }
        }

        return working.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.Values.ToList());
    }

    /// <summary>
    /// Calculates per-site modification occupancy in peptide-local coordinates directly from PSMs,
    /// returning results for all observed base sequences in a single call.
    /// PSM grouping, intensity derivation, and base-sequence bucketing are all handled internally.
    /// </summary>
    /// <param name="psms">
    /// All PSMs to consider. PSMs are grouped internally by <see cref="ISpectralMatch.BaseSequence"/>.
    /// PSMs whose <see cref="ISpectralMatch.Intensities"/> is a single-element array contribute
    /// to intensity-based stoichiometry; others contribute only to count-based metrics.
    /// </param>
    /// <param name="strategy">Intensity rollup strategy for computing intensity-based stoichiometry.
    /// Defaults to <see cref="IntensityRollupStrategy.Sum"/> for backward compatibility.</param>
    /// <returns>
    /// Dictionary keyed by peptide-local position (AllModsOneIsNterminus convention) containing 
    /// <see cref="SiteSpecificModificationOccupancy"/> entries.
    /// </returns>
    public static Dictionary<int, List<SiteSpecificModificationOccupancy>> CalculateDigestionProductLevelOccupancy(
        IEnumerable<ISpectralMatch> psms,
        IntensityRollupStrategy strategy = IntensityRollupStrategy.Sum)
    {
        var psmList = psms as IList<ISpectralMatch> ?? psms.ToList();

        var psmsWithBaseSeq = psmList.Where(p => p.BaseSequence != null).ToList();

        if (psmsWithBaseSeq.Count == 0)
            return new Dictionary<int, List<SiteSpecificModificationOccupancy>>();

        if (!psmsWithBaseSeq.Select(p => p.BaseSequence).AllSame())
        {
            throw new ArgumentException("All PSMs must have the same BaseSequence for peptide-level occupancy calculation.");
        }

        // Map each PSM to the single form that owns its intensity: matching FullSequence + Accession.
        // Ambiguous forms (psms without a full sequence match) are filtered out and do not contribute to occupancy.
        var psmToForm = psmsWithBaseSeq
            .ToDictionary(
                p => p,
                p => p.GetIdentifiedBioPolymersWithSetMods()
                    .FirstOrDefault(s => s.FullSequence == p.FullSequence));

        var totalCount = psmsWithBaseSeq.Count;
        var totalIntensitySum = psmsWithBaseSeq
            .Where(p => p.Intensities is { Length: 1 })
            .Sum(p => p.Intensities[0]);

        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();

        foreach (var psm in psmsWithBaseSeq)
        {
            var form = psmToForm[psm];
            if (form is null)
                continue;

            foreach (var mod in form.AllModsOneIsNterminus)
            {
                if (IsExcludedMod(mod.Value, ignoreLocation: true))
                    continue;

                if (!working.TryGetValue(mod.Key, out var modsAtPosition))
                {
                    modsAtPosition = new Dictionary<string, SiteSpecificModificationOccupancy>();
                    working[mod.Key] = modsAtPosition;
                }

                if (!modsAtPosition.ContainsKey(mod.Value.IdWithMotif))
                {
                    modsAtPosition[mod.Value.IdWithMotif] = new SiteSpecificModificationOccupancy(mod.Key, mod.Value.IdWithMotif)
                    {
                        TotalCount = totalCount,
                        TotalIntensity = totalIntensitySum
                    };
                }

                var siteOcc = modsAtPosition[mod.Value.IdWithMotif];
                siteOcc.ModifiedCount++;
                if (psm.Intensities is { Length: 1 })
                {
                    siteOcc.ModifiedIntensity += psm.Intensities[0];
                }
            }
        }

        if (working.Count == 0)
            return new Dictionary<int, List<SiteSpecificModificationOccupancy>>(); // Return empty if no mods passed filtering

        return working.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.Values.ToList());
    }

    /// <summary>
    /// Calculates per-site modification occupancy mapped to protein coordinates by first computing
    /// per-peptide occupancy ratios, then aggregating those ratios across peptides using the specified
    /// strategy. This produces bounded [0, 1] stoichiometry values that are not susceptible to the
    /// ratio-of-aggregates bias that affects <see cref="CalculateParentLevelOccupancy"/> with
    /// non-Sum strategies.
    /// </summary>
    /// <param name="bioPolymer">The parent biopolymer whose length defines the coordinate space.</param>
    /// <param name="psms">All PSMs to consider. Forms are filtered to <paramref name="bioPolymer"/> internally.</param>
    /// <param name="strategy">Aggregation strategy for combining per-peptide occupancies.
    /// Must be <see cref="IntensityRollupStrategy.Mean"/> or <see cref="IntensityRollupStrategy.Median"/>.
    /// For <see cref="IntensityRollupStrategy.Sum"/>, use <see cref="CalculateParentLevelOccupancy"/> instead.</param>
    /// <returns>
    /// Dictionary keyed by one-based protein position (AllModsOneIsNterminus convention) containing
    /// <see cref="SiteSpecificModificationOccupancy"/> entries with aggregated stoichiometry.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="strategy"/> is <see cref="IntensityRollupStrategy.Sum"/>.
    /// </exception>
    public static Dictionary<int, List<SiteSpecificModificationOccupancy>> CalculateParentLevelOccupancyByPeptide(
        IBioPolymer bioPolymer,
        IEnumerable<ISpectralMatch> psms,
        IntensityRollupStrategy strategy)
    {
        if (strategy == IntensityRollupStrategy.Sum)
            throw new ArgumentOutOfRangeException(nameof(strategy), strategy,
                "Sum strategy should use CalculateParentLevelOccupancy directly.");

        var psmList = psms as IList<ISpectralMatch> ?? psms.ToList();

        // Resolve forms via parallel arrays (avoids ToDictionary duplicate-key crash).
        var psmForms = new IBioPolymerWithSetMods?[psmList.Count];
        for (int j = 0; j < psmList.Count; j++)
        {
            var psm = psmList[j];
            var form = psm.GetIdentifiedBioPolymersWithSetMods()
                .FirstOrDefault(s => s.FullSequence != null
                    && s.BaseSequence == psm.BaseSequence
                    && s.FullSequence == psm.FullSequence
                    && s.Parent.Accession == bioPolymer.Accession);

            if (form is null)
            {
                try
                {
                    form = psm.GetIdentifiedBioPolymersWithSetMods()
                        .FirstOrDefault(s => s.BaseSequence == psm.BaseSequence
                            && s.Parent.Accession == bioPolymer.Accession);
                }
                catch (Exception)
                {
                    // Leave as null — PSM will be skipped.
                }
            }
            psmForms[j] = form;
        }

        // Group PSMs by peptide: unique combination of BaseSequence, Start, End on this protein.
        var peptideGroups = new List<List<int>>(); // indices into psmList
        var peptideKeys = new HashSet<(string baseSeq, int start, int end)>();
        for (int j = 0; j < psmList.Count; j++)
        {
            if (psmForms[j] is null) continue;
            var key = (psmForms[j]!.BaseSequence, psmForms[j]!.OneBasedStartResidue, psmForms[j]!.OneBasedEndResidue);
            if (!peptideKeys.TryGetValue(key, out _))
            {
                peptideKeys.Add(key);
                var indices = new List<int>();
                for (int k = 0; k < psmList.Count; k++)
                {
                if (psmForms[k] is not null
                    && psmForms[k]!.BaseSequence == key.BaseSequence
                    && psmForms[k]!.OneBasedStartResidue == key.OneBasedStartResidue
                    && psmForms[k]!.OneBasedEndResidue == key.OneBasedEndResidue)
                    {
                        indices.Add(k);
                    }
                }
                peptideGroups.Add(indices);
            }
        }

        // Per-peptide data: occupancy dict, covered positions, total intensity, total count
        var peptideData = new List<(
            Dictionary<int, List<SiteSpecificModificationOccupancy>> occupancy,
            int rangeStart,
            int rangeEnd,
            double totalIntensity,
            int totalCount)>();
        var allKeys = new HashSet<(int position, string modId)>();

        foreach (var indices in peptideGroups)
        {
            var peptidePsms = indices.Select(i => psmList[i]).ToList();
            var occupancy = CalculateParentLevelOccupancy(bioPolymer, peptidePsms, IntensityRollupStrategy.Sum);

            var form = psmForms[indices[0]]!;
            int rangeStart = form.OneBasedStartResidue + (form.OneBasedStartResidue == 1 ? 0 : 1);
            int rangeEnd = form.OneBasedEndResidue + (form.OneBasedEndResidue == bioPolymer.Length ? 2 : 1);

            double totalIntensity = peptidePsms
                .Where(p => p.Intensities is { Length: 1 })
                .Sum(p => p.Intensities[0]);
            int totalCount = peptidePsms.Count;

            peptideData.Add((occupancy, rangeStart, rangeEnd, totalIntensity, totalCount));

            foreach (var posKvp in occupancy)
            {
                foreach (var site in posKvp.Value)
                {
                    allKeys.Add((posKvp.Key, site.ModificationIdWithMotif));
                }
            }
        }

        // Collect per-peptide occupancy values, including 0.0 for unmodified peptides that cover the site.
        var occupancyByPeptide = new Dictionary<(int position, string modId), List<double>>();
        var countByPeptide = new Dictionary<(int position, string modId), List<int>>();
        var totalCountByPeptide = new Dictionary<(int position, string modId), List<int>>();
        var modifiedIntensityByPeptide = new Dictionary<(int position, string modId), List<double>>();
        var totalIntensityByPeptide = new Dictionary<(int position, string modId), List<double>>();

        foreach (var pd in peptideData)
        {
            foreach (var key in allKeys)
            {
                if (key.position < pd.rangeStart || key.position > pd.rangeEnd)
                    continue; // Peptide does not cover this position

                double stoichiometry = 0.0;
                int modifiedCount = 0;
                double modifiedIntensity = 0.0;

                if (pd.occupancy.TryGetValue(key.position, out var mods)
                    && mods.FirstOrDefault(m => m.ModificationIdWithMotif == key.modId) is { } site)
                {
                    stoichiometry = site.IntensityBasedStoichiometry;
                    modifiedCount = site.ModifiedCount;
                    modifiedIntensity = site.ModifiedIntensity;
                }

                occupancyByPeptide.AddToList(key, stoichiometry);
                if (!countByPeptide.TryGetValue(key, out var countList))
                {
                    countList = new List<int>();
                    countByPeptide[key] = countList;
                }
                countList.Add(modifiedCount);
                if (!totalCountByPeptide.TryGetValue(key, out var totalCountList))
                {
                    totalCountList = new List<int>();
                    totalCountByPeptide[key] = totalCountList;
                }
                totalCountList.Add(pd.totalCount);
                modifiedIntensityByPeptide.AddToList(key, modifiedIntensity);
                totalIntensityByPeptide.AddToList(key, pd.totalIntensity);
            }
        }

        // Aggregate per-peptide occupancies using the specified strategy.
        var result = new Dictionary<int, List<SiteSpecificModificationOccupancy>>();
        foreach (var kvp in occupancyByPeptide)
        {
            var (position, modId) = kvp.Key;
            var stoichiometries = kvp.Value;

            double aggregatedStoichiometry = ApplyStrategy(stoichiometries, strategy);

            // Sum counts and intensities across all peptides for this site.
            int totalModifiedCount = countByPeptide[kvp.Key].Sum();
            int totalCount = totalCountByPeptide[kvp.Key].Sum();
            double totalModifiedIntensity = modifiedIntensityByPeptide[kvp.Key].Sum();
            double totalTotalIntensity = totalIntensityByPeptide[kvp.Key].Sum();

            var site = new SiteSpecificModificationOccupancy(position, modId)
            {
                ModifiedCount = totalModifiedCount,
                TotalCount = totalCount,
                ModifiedIntensity = totalModifiedIntensity,
                TotalIntensity = totalTotalIntensity,
                IntensityBasedStoichiometry = aggregatedStoichiometry,
                PeptideCount = totalCountByPeptide[kvp.Key].Count,
                ModifiedPeptideCount = countByPeptide[kvp.Key].Count(c => c > 0)
            };

            if (!result.TryGetValue(position, out var modsAtPosition))
            {
                modsAtPosition = new List<SiteSpecificModificationOccupancy>();
                result[position] = modsAtPosition;
            }
            modsAtPosition.Add(site);
        }

        return result;
    }

    /// <summary>
    /// Applies the specified rollup strategy to a list of intensity values.
    /// </summary>
    private static double ApplyStrategy(List<double> values, IntensityRollupStrategy strategy)
    {
        if (values.Count == 0)
            return 0;

        return strategy switch
        {
            IntensityRollupStrategy.Mean => values.Average(),
            IntensityRollupStrategy.Median => Median(values),
            IntensityRollupStrategy.Sum => values.Sum(),
            _ => throw new ArgumentOutOfRangeException(nameof(strategy), strategy, $"Unsupported rollup strategy: {strategy}")
        };
    }

    /// <summary>
    /// Computes the median of a list of values.
    /// For an even number of elements, returns the average of the two middle values.
    /// </summary>
    private static double Median(List<double> values)
    {
        if (values.Count == 0)
            return 0;

        var sorted = values.OrderBy(v => v).ToList();
        int n = sorted.Count;
        if (n % 2 == 1)
            return sorted[n / 2];
        return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
    }

    private static void AddToList<TKey>(this Dictionary<TKey, List<double>> dict, TKey key, double value)
        where TKey : notnull
    {
        if (!dict.TryGetValue(key, out var list))
        {
            list = new List<double>();
            dict[key] = list;
        }
        list.Add(value);
    }

    private static bool TryGetProteinPosition(
        KeyValuePair<int, Modification> mod,
        IBioPolymerWithSetMods sequence,
        int bioPolymerLength,
        out int indexInProtein)
    {
        indexInProtein = 0;

        if (IsExcludedMod(mod.Value))
            return false;

        if (mod.Value.LocationRestriction.Equals("N-terminal."))
        {
            if (sequence.OneBasedStartResidue != 1)
                return false;

            indexInProtein = 1;
        }
        else if (mod.Value.LocationRestriction.Equals("Anywhere."))
        {
            indexInProtein = sequence.OneBasedStartResidue + mod.Key - 1;
        }
        else if (mod.Value.LocationRestriction.Equals("C-terminal."))
        {
            if (sequence.OneBasedEndResidue != bioPolymerLength)
                return false;

            indexInProtein = bioPolymerLength + 2;
        }
        else
        {
            return false;
        }

        return true;
    }

    private static bool IsExcludedMod(Modification mod, bool ignoreLocation = false)
    {
        if (ExcludedLocations.Contains(mod.LocationRestriction) && !ignoreLocation)
            return true;

        if (ExcludedModTypes.Contains(mod.ModificationType))
            return true;

        return false;
    }
}
