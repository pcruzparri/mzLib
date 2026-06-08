namespace Omics.BioPolymerGroup;

/// <summary>
/// Strategy for aggregating PSM intensities when computing intensity-based
/// modification occupancy/stoichiometry.
/// </summary>
public enum IntensityRollupStrategy
{
    /// <summary>
    /// Sum of all PSM intensities. The default and historical behavior.
    /// Stoichiometry = Sum(modified PSM intensities) / Sum(all PSM intensities covering the site).
    /// </summary>
    Sum,

    /// <summary>
    /// Arithmetic mean of per-peptide occupancy ratios. For each peptide covering a site,
    /// occupancy is computed as Sum(modified PSM intensities) / Sum(all PSM intensities),
    /// then the mean of those per-peptide ratios is taken. This produces bounded [0, 1]
    /// stoichiometry values that are not susceptible to the ratio-of-aggregates bias.
    /// </summary>
    Mean,

    /// <summary>
    /// Median of per-peptide occupancy ratios. For each peptide covering a site,
    /// occupancy is computed as Sum(modified PSM intensities) / Sum(all PSM intensities),
    /// then the median of those per-peptide ratios is taken. This produces bounded [0, 1]
    /// stoichiometry values that are robust against outlier peptides.
    /// </summary>
    Median
}
