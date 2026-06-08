using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Omics.SpectralMatch;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics.Occupancy;

[TestFixture]
[ExcludeFromCodeCoverage]
public class ModificationOccupancyCalculatorTests
{
    #region CalculateProteinLevelOccupancy Tests

    [Test]
    public void ProteinLevelWithSingleModOnSinglePeptide()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var psm = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [peptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm]);

        Assert.That(result.ContainsKey(4), Is.True);
        Assert.That(result[4].Count, Is.EqualTo(1));
        Assert.That(result[4][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[4][0].TotalCount, Is.EqualTo(1));
        Assert.That(result[4][0].CountBasedOccupancy, Is.EqualTo(1.0));
    }

    [Test]
    public void ProteinLevelWithModifiedAndUnmodifiedPeptides()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm1, psm2]);

        Assert.That(result.ContainsKey(4), Is.True);
        Assert.That(result[4][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[4][0].TotalCount, Is.EqualTo(2));
        Assert.That(result[4][0].CountBasedOccupancy, Is.EqualTo(0.5));
    }

    [Test]
    public void ProteinLevelModIsExcluded()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var commonMod = new Modification("Oxidation", null, "Common Variable", null, motif, "Anywhere.", null, 15.995);

        var mods = new Dictionary<int, Modification> { { 4, commonMod } };
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Oxidation]EF", protein, 1, 5, mods);
        var psm = new MockSpectralMatch("test.raw", "ACD[Oxidation]EF", "ACDEF", 1.0, 1, [peptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm]);

        Assert.That(result, Is.Empty);
    }

    [Test]
    public void ProteinLevelPeptideTerminalModIsExcluded()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("A", out var motif);
        var pepNMod = new Modification("Acetylation", null, "Biological", null, motif, "NPep", null, 42.011);

        var mods = new Dictionary<int, Modification> { { 1, pepNMod } };
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "[Acetylation]ACDEF", protein, 1, 5, mods);
        var psm = new MockSpectralMatch("test.raw", "[Acetylation]ACDEF", "ACDEF", 1.0, 1, [peptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm]);

        Assert.That(result, Is.Empty);
    }

    [Test]
    public void ProteinLevelWithIntensities()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [3_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm1, psm2]);

        var site = result[4][0];
        Assert.That(site.ModifiedIntensity, Is.EqualTo(1_000_000));
        Assert.That(site.TotalIntensity, Is.EqualTo(4_000_000));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0.25));
    }

    [Test]
    public void ProteinLevelWithOverlappingPeptidesCoveringPosition()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        // Peptide 1: ACDEF (positions 1-5), modified at D (position 3 in protein)
        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);

        // Peptide 2: CDEFG (positions 2-6), unmodified but covers position 3
        var overlappingPeptide = new MockBioPolymerWithSetMods("CDEFG", "CDEFG", protein, 2, 6);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        var psm2 = new MockSpectralMatch("test.raw", "CDEFG", "CDEFG", 1.0, 2, [overlappingPeptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm1, psm2]);

        Assert.That(result[4][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[4][0].TotalCount, Is.EqualTo(2));
        Assert.That(result[4][0].CountBasedOccupancy, Is.EqualTo(0.5));
    }

    [Test]
    public void ProteinLevelWithNoPeptides()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, Enumerable.Empty<ISpectralMatch>());

        Assert.That(result, Is.Empty);
    }

    /// <summary>
    /// Regression test: a single PSM whose <see cref="ISpectralMatch.GetIdentifiedBioPolymersWithSetMods"/>
    /// returns two ambiguous forms must not inflate TotalCount.
    /// Only the form whose FullSequence matches <see cref="ISpectralMatch.FullSequence"/> is counted;
    /// the alternative form is ignored entirely.
    /// </summary>
    [Test]
    public void ProteinLevel_SinglePsmTwoAmbiguousInterpretations_OccupancyIsNotInflated()
    {
        var protein = new MockBioPolymer("IVENGSEQGSYDADK", "Q6PI26");
        ModificationMotif.TryGetMotif("N", out var motif);
        var deamidation = new Modification("Deamidation on N", null, "Biological", null, motif, "Anywhere.", null, 0.984);
        var deamidatedAsp = new Modification("Deamidated asparagine on N", null, "Biological", null, motif, "Anywhere.", null, 0.984);

        var form1 = new MockBioPolymerWithSetMods(
            "IVEN", "IVEN[Deamidation on N]", protein, 1, 4,
            new Dictionary<int, Modification> { { 5, deamidation } });
        var form2 = new MockBioPolymerWithSetMods(
            "IVEN", "IVEN[Deamidated asparagine on N]", protein, 1, 4,
            new Dictionary<int, Modification> { { 5, deamidatedAsp } });

        // Single PSM: FullSequence matches form1. form2 is an alternative returned by
        // GetIdentifiedBioPolymersWithSetMods() but must not be counted.
        var psm = new MockSpectralMatch("test.raw", "IVEN[Deamidation on N]", "IVEN", 1.0, 1, [form1, form2]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm]);

        Assert.That(result.ContainsKey(5), Is.True);
        var modsAtSite = result[5];

        // Only the form matching PSM.FullSequence should be discovered.
        var deamSite = modsAtSite.FirstOrDefault(s => s.ModificationIdWithMotif == "Deamidation on N");
        Assert.That(deamSite, Is.Not.Null);
        Assert.That(deamSite!.TotalCount, Is.EqualTo(1), "TotalCount must be 1 (one PSM), not 2 (two forms)");
        Assert.That(deamSite.ModifiedCount, Is.EqualTo(1));
        Assert.That(deamSite.CountBasedOccupancy, Is.EqualTo(1.0));

        // The unmatched alternative form's mod must not appear.
        var deamAspSite = modsAtSite.FirstOrDefault(s => s.ModificationIdWithMotif == "Deamidated asparagine on N");
        Assert.That(deamAspSite, Is.Null, "Alternative form not matching PSM.FullSequence must be excluded");
    }

    /// <summary>
    /// Regression test: when a PSM maps to two proteins and protein B's form appears first in
    /// <see cref="ISpectralMatch.GetIdentifiedBioPolymersWithSetMods"/>, TotalCount for protein A
    /// must still be 1 — the Accession filter inside the calculator ensures the correct form is found.
    /// </summary>
    [Test]
    public void ProteinLevel_SharedPeptideTwoProteins_TotalCountIsNotUnderCounted()
    {
        var proteinA = new MockBioPolymer("ACDEFGHIK", "P00001");
        var proteinB = new MockBioPolymer("ACDEFKLMN", "P00002");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var formB = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", proteinB, 1, 5, mods);
        var formA = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", proteinA, 1, 5, mods);

        // PSM returns protein B's form first; the calculator must still resolve protein A's form.
        var psm = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [formB, formA]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(proteinA, [psm]);

        Assert.That(result.ContainsKey(4), Is.True);
        var site = result[4][0];
        Assert.That(site.TotalCount, Is.EqualTo(1), "TotalCount must be 1 — protein-A form must be found even when protein B's form comes first");
        Assert.That(site.ModifiedCount, Is.EqualTo(1));
        Assert.That(site.CountBasedOccupancy, Is.LessThanOrEqualTo(1.0));
    }

    #endregion

    #region CalculatePeptideLevelOccupancy Tests

    [Test]
    public void PeptideLevelOccupancyReturnedPerGroup()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);

        var result = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([psm1, psm2]);

        Assert.That(result.ContainsKey(4), Is.True);
        var site = result[4][0];
        Assert.That(site.ModifiedCount, Is.EqualTo(1));
        Assert.That(site.TotalCount, Is.EqualTo(2));
        Assert.That(site.CountBasedOccupancy, Is.EqualTo(0.5));
    }

    [Test]
    public void PeptideLevelWithIntensities()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [2_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [8_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([psm1, psm2]);

        var site = result[4][0];
        Assert.That(site.ModifiedIntensity, Is.EqualTo(2_000_000));
        Assert.That(site.TotalIntensity, Is.EqualTo(10_000_000));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0.2));
    }

    [Test]
    public void PeptideLevelCommonFixedModIsExcluded()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("C", out var motif);
        var fixedMod = new Modification("Carbamidomethyl", null, "Common Fixed", null, motif, "Anywhere.", null, 57.021);

        var mods = new Dictionary<int, Modification> { { 3, fixedMod } };
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "AC[Carbamidomethyl]DEF", protein, 1, 5, mods);
        var psm = new MockSpectralMatch("test.raw", "AC[Carbamidomethyl]DEF", "ACDEF", 1.0, 1, [peptide]);

        var result = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([psm]);

        Assert.That(result, Is.Empty);
    }

    /// <summary>
    /// The method requires all PSMs to share the same BaseSequence.
    /// Passing PSMs with different base sequences must throw <see cref="ArgumentException"/>.
    /// To calculate occupancy across multiple peptides, call the method separately per base sequence.
    /// </summary>
    [Test]
    public void PeptideLevelWithMultipleBaseSequences()
    {
        var protein = new MockBioPolymer("ACDEFGHIKLMNPQR", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods1 = new Dictionary<int, Modification> { { 4, mod } };
        var peptide1 = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods1);

        var mods2 = new Dictionary<int, Modification> { { 3, mod } };
        var peptide2 = new MockBioPolymerWithSetMods("GHIKLM", "GH[Phosphorylation]IKLM", protein, 6, 11, mods2);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [peptide1]);
        var psm2 = new MockSpectralMatch("test.raw", "GH[Phosphorylation]IKLM", "GHIKLM", 1.0, 2, [peptide2]);

        // PSMs with different base sequences must not be mixed in a single call.
        Assert.That(
            () => ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([psm1, psm2]),
            Throws.ArgumentException);
    }

    #endregion

    #region Intensity Rollup Strategy Tests

    [Test]
    public void ProteinLevel_DefaultRollupIsSum()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [3_000_000.0];

        // Default (no parameter) should match Sum
        var resultDefault = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm1, psm2]);
        var resultExplicit = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm1, psm2], IntensityRollupStrategy.Sum);

        Assert.That(resultDefault[4][0].ModifiedIntensity, Is.EqualTo(resultExplicit[4][0].ModifiedIntensity));
        Assert.That(resultDefault[4][0].TotalIntensity, Is.EqualTo(resultExplicit[4][0].TotalIntensity));
        Assert.That(resultDefault[4][0].IntensityBasedStoichiometry, Is.EqualTo(0.25));
    }

    [Test]
    public void ProteinLevel_MeanStrategy_CorrectAverages()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        // 1 modified PSM at intensity 1e6, 2 unmodified PSMs at 2e6 and 3e6
        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [2_000_000.0];
        var psm3 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 3, [unmodifiedPeptide]);
        psm3.Intensities = [3_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, [psm1, psm2, psm3], IntensityRollupStrategy.Mean);

        var site = result[4][0];
        // Mean(modified) = 1e6 / 1 = 1e6
        // Mean(total) = (1e6 + 2e6 + 3e6) / 3 = 2e6
        // Stoichiometry = 1e6 / 2e6 = 0.5
        Assert.That(site.ModifiedIntensity, Is.EqualTo(1_000_000.0));
        Assert.That(site.TotalIntensity, Is.EqualTo(2_000_000.0));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0.5));
    }

    [Test]
    public void ProteinLevel_MedianStrategy_CorrectMedian()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        // 1 modified PSM at intensity 10, 3 unmodified PSMs at 1, 20, 100
        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [10.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [1.0];
        var psm3 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 3, [unmodifiedPeptide]);
        psm3.Intensities = [20.0];
        var psm4 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 4, [unmodifiedPeptide]);
        psm4.Intensities = [100.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, [psm1, psm2, psm3, psm4], IntensityRollupStrategy.Median);

        var site = result[4][0];
        // Median(modified) = [10] -> 10
        // Median(total) = [1, 10, 20, 100] sorted -> (10 + 20) / 2 = 15
        // Stoichiometry = 10 / 15 = 0.6667
        Assert.That(site.ModifiedIntensity, Is.EqualTo(10.0));
        Assert.That(site.TotalIntensity, Is.EqualTo(15.0));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(10.0 / 15.0).Within(1e-10));
    }

    [Test]
    public void ProteinLevel_MedianWithEvenCount_AveragesMiddleTwo()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        // 2 modified PSMs: [10, 100], 2 unmodified PSMs: [1, 20]
        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [10.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 2, [modifiedPeptide]);
        psm2.Intensities = [100.0];
        var psm3 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 3, [unmodifiedPeptide]);
        psm3.Intensities = [1.0];
        var psm4 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 4, [unmodifiedPeptide]);
        psm4.Intensities = [20.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, [psm1, psm2, psm3, psm4], IntensityRollupStrategy.Median);

        var site = result[4][0];
        // Median(modified) = [10, 100] -> (10 + 100) / 2 = 55
        // Median(total) = [1, 10, 20, 100] -> (10 + 20) / 2 = 15
        Assert.That(site.ModifiedIntensity, Is.EqualTo(55.0));
        Assert.That(site.TotalIntensity, Is.EqualTo(15.0));
    }

    [Test]
    public void ProteinLevel_SinglePsm_AllStrategiesEqual()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);

        var psm = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm.Intensities = [500.0];

        var resultSum = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, [psm], IntensityRollupStrategy.Sum);
        var resultMean = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, [psm], IntensityRollupStrategy.Mean);
        var resultMedian = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, [psm], IntensityRollupStrategy.Median);

        // Single PSM: Sum = Mean = Median = 500
        Assert.That(resultMean[4][0].ModifiedIntensity, Is.EqualTo(resultSum[4][0].ModifiedIntensity));
        Assert.That(resultMedian[4][0].ModifiedIntensity, Is.EqualTo(resultSum[4][0].ModifiedIntensity));
        Assert.That(resultMean[4][0].TotalIntensity, Is.EqualTo(resultSum[4][0].TotalIntensity));
        Assert.That(resultMedian[4][0].TotalIntensity, Is.EqualTo(resultSum[4][0].TotalIntensity));
        Assert.That(resultMean[4][0].IntensityBasedStoichiometry, Is.EqualTo(1.0));
        Assert.That(resultMedian[4][0].IntensityBasedStoichiometry, Is.EqualTo(1.0));
    }

    [Test]
    public void ProteinLevel_NoIntensityData_AllStrategiesReturnZero()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        // No Intensities set

        foreach (var strategy in new[] { IntensityRollupStrategy.Sum, IntensityRollupStrategy.Mean, IntensityRollupStrategy.Median })
        {
            var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
                protein, [psm1, psm2], strategy);
            var site = result[4][0];
            Assert.That(site.ModifiedIntensity, Is.EqualTo(0), $"Strategy {strategy}: ModifiedIntensity should be 0");
            Assert.That(site.TotalIntensity, Is.EqualTo(0), $"Strategy {strategy}: TotalIntensity should be 0");
            Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0), $"Strategy {strategy}: Stoichiometry should be 0");
        }
    }

    [Test]
    public void PeptideLevel_MeanStrategy_FallsBackToSum()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [3_000_000.0];
        var psm3 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 3, [unmodifiedPeptide]);
        psm3.Intensities = [5_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy(
            [psm1, psm2, psm3], IntensityRollupStrategy.Mean);

        var site = result[4][0];
        // At peptide level, Mean/Median fall back to Sum (single peptide, per-PSM aggregation is invalid)
        Assert.That(site.ModifiedIntensity, Is.EqualTo(1_000_000.0));
        Assert.That(site.TotalIntensity, Is.EqualTo(9_000_000.0));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(1.0 / 9.0).Within(1e-10));
    }

    [Test]
    public void PeptideLevel_MedianStrategy_FallsBackToSum()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [10.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [1.0];
        var psm3 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 3, [unmodifiedPeptide]);
        psm3.Intensities = [20.0];
        var psm4 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 4, [unmodifiedPeptide]);
        psm4.Intensities = [100.0];

        var result = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy(
            [psm1, psm2, psm3, psm4], IntensityRollupStrategy.Median);

        var site = result[4][0];
        // At peptide level, Mean/Median fall back to Sum
        Assert.That(site.ModifiedIntensity, Is.EqualTo(10.0));
        Assert.That(site.TotalIntensity, Is.EqualTo(131.0));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(10.0 / 131.0).Within(1e-10));
    }

    [Test]
    public void PeptideLevel_EmptyInput_ReturnsEmpty()
    {
        var result = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy(
            Array.Empty<ISpectralMatch>());

        Assert.That(result, Is.Empty);
    }

    #endregion

    #region CalculateParentLevelOccupancyByPeptide Tests (Option 3)

    [Test]
    public void ParentLevelByPeptide_Mean_NeverExceedsOne()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        // Peptide 1: ACDEF (pos 1-5), modified at D (pos 4)
        var mods1 = new Dictionary<int, Modification> { { 4, mod } };
        var peptide1Modified = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods1);
        var peptide1Unmod = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        // Peptide 2: CDEFG (pos 2-6), covers D (pos 4) but unmodified
        var peptide2 = new MockBioPolymerWithSetMods("CDEFG", "CDEFG", protein, 2, 6);

        // PSMs for peptide 1: 1 modified (intensity 1e6), 1 unmodified (intensity 3e6)
        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [peptide1Modified]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [peptide1Unmod]);
        psm2.Intensities = [3_000_000.0];

        // PSMs for peptide 2: 1 unmodified (intensity 1e6)
        var psm3 = new MockSpectralMatch("test.raw", "CDEFG", "CDEFG", 1.0, 3, [peptide2]);
        psm3.Intensities = [1_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
            protein, [psm1, psm2, psm3], IntensityRollupStrategy.Mean);

        var site = result[4][0];
        // Per-peptide occupancies:
        //   Peptide 1: 1e6 / (1e6 + 3e6) = 0.25
        //   Peptide 2: 0 / 1e6 = 0
        // Mean = (0.25 + 0) / 2 = 0.125
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0.125).Within(1e-10));
        Assert.That(site.IntensityBasedStoichiometry, Is.LessThanOrEqualTo(1.0));
        Assert.That(site.IntensityBasedStoichiometry, Is.GreaterThanOrEqualTo(0.0));
    }

    [Test]
    public void ParentLevelByPeptide_Median_NeverExceedsOne()
    {
        var protein = new MockBioPolymer("ACDEFGHIKLMNPQR", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        // Three peptides covering position 4 (D), each with different occupancy ratios
        var mods = new Dictionary<int, Modification> { { 4, mod } };

        var pep1Mod = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var pep1Unmod = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);
        var pep2 = new MockBioPolymerWithSetMods("CDEFG", "CDEFG", protein, 2, 6);
        var pep3Mod = new MockBioPolymerWithSetMods("BCDEF", "B[Phosphorylation]CDEF", protein, 2, 6,
            new Dictionary<int, Modification> { { 3, mod } });
        var pep3Unmod = new MockBioPolymerWithSetMods("BCDEF", "BCDEF", protein, 2, 6);

        // Peptide 1: 1 modified (1e6), 1 unmodified (3e6) -> ratio = 0.25
        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [pep1Mod]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [pep1Unmod]);
        psm2.Intensities = [3_000_000.0];

        // Peptide 2: unmodified only (1e6) -> ratio = 0
        var psm3 = new MockSpectralMatch("test.raw", "CDEFG", "CDEFG", 1.0, 3, [pep2]);
        psm3.Intensities = [1_000_000.0];

        // Peptide 3: 1 modified (3e6), 1 unmodified (1e6) -> ratio = 0.75
        var psm4 = new MockSpectralMatch("test.raw", "B[Phosphorylation]CDEF", "BCDEF", 1.0, 4, [pep3Mod]);
        psm4.Intensities = [3_000_000.0];
        var psm5 = new MockSpectralMatch("test.raw", "BCDEF", "BCDEF", 1.0, 5, [pep3Unmod]);
        psm5.Intensities = [1_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
            protein, [psm1, psm2, psm3, psm4, psm5], IntensityRollupStrategy.Median);

        var site = result[4][0];
        // Per-peptide ratios: 0.25, 0, 0.75 -> sorted: 0, 0.25, 0.75 -> median = 0.25
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0.25).Within(1e-10));
        Assert.That(site.IntensityBasedStoichiometry, Is.LessThanOrEqualTo(1.0));
        Assert.That(site.IntensityBasedStoichiometry, Is.GreaterThanOrEqualTo(0.0));
    }

    [Test]
    public void ParentLevelByPeptide_PeptideCountCorrect()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods1 = new Dictionary<int, Modification> { { 4, mod } };
        var peptide1 = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods1);
        var peptide2 = new MockBioPolymerWithSetMods("CDEFG", "CDEFG", protein, 2, 6);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [peptide1]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "CDEFG", "CDEFG", 1.0, 2, [peptide2]);
        psm2.Intensities = [1_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
            protein, [psm1, psm2], IntensityRollupStrategy.Mean);

        var site = result[4][0];
        Assert.That(site.PeptideCount, Is.EqualTo(2));
    }

    [Test]
    public void ParentLevelByPeptide_ModifiedPeptideCountCorrect()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods1 = new Dictionary<int, Modification> { { 4, mod } };
        var peptide1 = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods1);
        var peptide2 = new MockBioPolymerWithSetMods("CDEFG", "CDEFG", protein, 2, 6);
        var peptide3 = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods1);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [peptide1]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "CDEFG", "CDEFG", 1.0, 2, [peptide2]);
        psm2.Intensities = [1_000_000.0];
        var psm3 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 3, [peptide3]);
        psm3.Intensities = [2_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
            protein, [psm1, psm2, psm3], IntensityRollupStrategy.Mean);

        var site = result[4][0];
        // Two distinct peptides have the modification (ACDEF from psm1 and psm3 are same peptide)
        Assert.That(site.ModifiedPeptideCount, Is.EqualTo(1));
    }

    [Test]
    public void ParentLevelByPeptide_SumStrategyThrowsArgumentOutOfRange()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        Assert.That(
            () => ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
                protein, [psm], IntensityRollupStrategy.Sum),
            Throws.InstanceOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void ParentLevelByPeptide_SinglePeptide_EqualsSumStrategy()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [3_000_000.0];

        var resultMean = ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
            protein, [psm1, psm2], IntensityRollupStrategy.Mean);
        var resultMedian = ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
            protein, [psm1, psm2], IntensityRollupStrategy.Median);
        var resultSum = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, [psm1, psm2], IntensityRollupStrategy.Sum);

        // With only one peptide, Mean/Median/Sum should all give the same per-peptide ratio
        Assert.That(resultMean[4][0].IntensityBasedStoichiometry, Is.EqualTo(0.25).Within(1e-10));
        Assert.That(resultMedian[4][0].IntensityBasedStoichiometry, Is.EqualTo(0.25).Within(1e-10));
        Assert.That(resultSum[4][0].IntensityBasedStoichiometry, Is.EqualTo(0.25).Within(1e-10));
    }

    [Test]
    public void ParentLevelByPeptide_EmptyInput_ReturnsEmpty()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
            protein, Array.Empty<ISpectralMatch>(), IntensityRollupStrategy.Mean);

        Assert.That(result, Is.Empty);
    }

    [Test]
    public void ParentLevelByPeptide_AllUnmodified_ReturnsEmpty()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1, [peptide]);
        psm.Intensities = [1_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
            protein, [psm], IntensityRollupStrategy.Mean);

        // No modifications observed, so no occupancy entries to return
        Assert.That(result, Is.Empty);
    }

    [Test]
    public void ParentLevelByPeptide_UnmodifiedPeptidesIncludedAsZero()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptideA = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);
        var unmodifiedPeptideB = new MockBioPolymerWithSetMods("CDEFG", "CDEFG", protein, 2, 6);

        // Modified peptide: 1 modified PSM
        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [1_000_000.0];

        // Unmodified peptide A: 1 unmodified PSM (same sequence as modified)
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptideA]);
        psm2.Intensities = [3_000_000.0];

        // Unmodified peptide B: 1 unmodified PSM (different sequence, covers position 4)
        var psm3 = new MockSpectralMatch("test.raw", "CDEFG", "CDEFG", 1.0, 3, [unmodifiedPeptideB]);
        psm3.Intensities = [1_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancyByPeptide(
            protein, [psm1, psm2, psm3], IntensityRollupStrategy.Mean);

        var site = result[4][0];
        // Per-peptide ratios:
        //   ACDEF group (psm1 + psm2, same peptide): 1e6 / (1e6 + 3e6) = 0.25
        //   CDEFG group (psm3, different peptide, unmodified): 0.0
        // Mean = (0.25 + 0.0) / 2 = 0.125
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0.125).Within(1e-10));
        Assert.That(site.PeptideCount, Is.EqualTo(2));
        Assert.That(site.ModifiedPeptideCount, Is.EqualTo(1));
    }

    #endregion

    #region SiteSpecificModificationOccupancy Tests

    [Test]
    public void ToSpectralCountModInfoStringMatchesExpectedFormat()
    {
        var site = new SiteSpecificModificationOccupancy(5, "Phosphorylation on S")
        {
            ModifiedCount = 3,
            TotalCount = 10
        };

        string expected = "pos4[Phosphorylation on S,info:fraction=0.30(3/10)]";
        Assert.That(site.ToModInfoString(intensityBased: false), Is.EqualTo(expected));
    }

    [Test]
    public void IntensityBasedStoichiometryZeroTotalIntensityDoesNotThrowDivByZero()
    {
        var site = new SiteSpecificModificationOccupancy(1, "TestMod")
        {
            ModifiedIntensity = 0,
            TotalIntensity = 0
        };

        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0));
    }

    [Test]
    public void CountBasedOccupancyZeroTotalIntensityDoesNotThrowDivByZero()
    {
        var site = new SiteSpecificModificationOccupancy(1, "TestMod")
        {
            ModifiedCount = 0,
            TotalCount = 0
        };

        Assert.That(site.CountBasedOccupancy, Is.EqualTo(0));
    }

    #endregion
}
