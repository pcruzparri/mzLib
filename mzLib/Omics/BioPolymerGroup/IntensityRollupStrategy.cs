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
    /// Arithmetic mean of PSM intensities. Reduces influence of high-abundance outlier PSMs.
    /// Stoichiometry = Mean(modified PSM intensities) / Mean(all PSM intensities covering the site).
    /// </summary>
    Mean,

    /// <summary>
    /// Median of PSM intensities. Robust against outlier intensity values.
    /// Stoichiometry = Median(modified PSM intensities) / Median(all PSM intensities covering the site).
    /// Note: The median ratio can exceed 1.0 when the median modified intensity is greater
    /// than the median total intensity (e.g., when modified PSMs are systematically higher
    /// intensity than unmodified PSMs).
    /// </summary>
    Median
}
