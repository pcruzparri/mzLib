using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics.Occupancy;

[TestFixture]
[ExcludeFromCodeCoverage]
public class PeptideGroupTests
{
    #region Constructor and Identity

    [Test]
    public void Constructor_InitializesCorrectly()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm]);

        Assert.That(group.BaseSequence, Is.EqualTo("ACDEF"));
        Assert.That(group.OneBasedStartResidue, Is.EqualTo(1));
        Assert.That(group.OneBasedEndResidue, Is.EqualTo(5));
        Assert.That(group.Parents.Count, Is.EqualTo(1));
        Assert.That(group.Parents[0].Accession, Is.EqualTo("P00001"));
        Assert.That(group.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(1));
        Assert.That(group.PeptideGroupName, Is.EqualTo("ACDEF_1_5_P00001"));
    }

    [Test]
    public void Constructor_MultipleParents()
    {
        var proteinA = new MockBioPolymer("ACDEFGHIK", "P00001");
        var proteinB = new MockBioPolymer("ACDEFKLMN", "P00002");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group = new PeptideGroup("ACDEF", 1, 5, [proteinA, proteinB], [psm]);

        Assert.That(group.Parents.Count, Is.EqualTo(2));
        Assert.That(group.PeptideGroupName, Is.EqualTo("ACDEF_1_5_P00001|P00002"));
    }

    [Test]
    public void Constructor_NullBaseSequenceThrows()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        Assert.That(
            () => new PeptideGroup(null!, 1, 5, [protein], [psm]),
            Throws.InstanceOf<ArgumentNullException>());
    }

    [Test]
    public void Constructor_NullParentsThrows()
    {
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        Assert.That(
            () => new PeptideGroup("ACDEF", 1, 5, null!, [psm]),
            Throws.InstanceOf<ArgumentNullException>());
    }

    [Test]
    public void Constructor_NullPsmsThrows()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");

        Assert.That(
            () => new PeptideGroup("ACDEF", 1, 5, [protein], null!),
            Throws.InstanceOf<ArgumentNullException>());
    }

    #endregion

    #region Flags

    [Test]
    public void IsDecoy_TrueWhenAnyParentIsDecoy()
    {
        var target = new MockBioPolymer("ACDEFGHIK", "P00001", false);
        var decoy = new MockBioPolymer("ACDEFGHIK", "P00002", true);
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group = new PeptideGroup("ACDEF", 1, 5, [target, decoy], [psm]);

        Assert.That(group.IsDecoy, Is.True);
    }

    [Test]
    public void IsDecoy_FalseWhenNoParentIsDecoy()
    {
        var target = new MockBioPolymer("ACDEFGHIK", "P00001", false);
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group = new PeptideGroup("ACDEF", 1, 5, [target], [psm]);

        Assert.That(group.IsDecoy, Is.False);
    }

    #endregion

    #region Equality

    [Test]
    public void Equals_SameSequenceAndParents_ReturnsTrue()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm1 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2);

        var group1 = new PeptideGroup("ACDEF", 1, 5, [protein], [psm1]);
        var group2 = new PeptideGroup("ACDEF", 1, 5, [protein], [psm2]);

        Assert.That(group1.Equals(group2), Is.True);
        Assert.That(group1.GetHashCode(), Is.EqualTo(group2.GetHashCode()));
    }

    [Test]
    public void Equals_DifferentBaseSequence_ReturnsFalse()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group1 = new PeptideGroup("ACDEF", 1, 5, [protein], [psm]);
        var group2 = new PeptideGroup("CDEFG", 2, 6, [protein], [psm]);

        Assert.That(group1.Equals(group2), Is.False);
    }

    [Test]
    public void Equals_DifferentParents_ReturnsFalse()
    {
        var proteinA = new MockBioPolymer("ACDEFGHIK", "P00001");
        var proteinB = new MockBioPolymer("ACDEFGHIK", "P00002");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group1 = new PeptideGroup("ACDEF", 1, 5, [proteinA], [psm]);
        var group2 = new PeptideGroup("ACDEF", 1, 5, [proteinB], [psm]);

        Assert.That(group1.Equals(group2), Is.False);
    }

    [Test]
    public void Equals_Null_ReturnsFalse()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);
        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm]);

        Assert.That(group.Equals(null), Is.False);
    }

    [Test]
    public void Equals_SameReference_ReturnsTrue()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);
        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm]);

        Assert.That(group.Equals(group), Is.True);
    }

    #endregion

    #region PopulateSampleGroupResults

    [Test]
    public void PopulateSampleGroupResults_NoExperimentalDesign_GroupsByFile()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("fileA.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        var psm2 = new MockSpectralMatch("fileB.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);

        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm1, psm2]);
        group.PopulateSampleGroupResults();

        Assert.That(group.SampleGroupResults, Is.Not.Null);
        Assert.That(group.SampleGroupResults!.Count, Is.EqualTo(2));

        var fileAResult = group.SampleGroupResults.First(r => r.Label == "fileA");
        Assert.That(fileAResult.SpectralCount, Is.EqualTo(1));
        Assert.That(fileAResult.ParentOccupancy.ContainsKey("ACDEF"), Is.True);

        var fileBResult = group.SampleGroupResults.First(r => r.Label == "fileB");
        Assert.That(fileBResult.SpectralCount, Is.EqualTo(1));
    }

    #endregion

    #region Subsetting

    [Test]
    public void ConstructSubsetPeptideGroup_FiltersByFile()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm1 = new MockSpectralMatch("fileA.raw", "ACDEF", "ACDEF", 1.0, 1);
        var psm2 = new MockSpectralMatch("fileB.raw", "ACDEF", "ACDEF", 1.0, 2);

        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm1, psm2])
        {
            OccupancyRollupStrategy = IntensityRollupStrategy.Mean
        };
        var subset = group.ConstructSubsetPeptideGroup("fileA.raw");

        Assert.That(subset.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(1));
        Assert.That(subset.AllPsmsBelowOnePercentFDR.First().FullFilePath, Is.EqualTo("fileA.raw"));
        Assert.That(subset.BaseSequence, Is.EqualTo("ACDEF"));
        Assert.That(subset.Parents, Is.EqualTo(group.Parents));
        Assert.That(subset.OccupancyRollupStrategy, Is.EqualTo(IntensityRollupStrategy.Mean));
    }

    [Test]
    public void HeaderAndToString_ColumnCountsMatch()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);
        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);

        var psm1 = new MockSpectralMatch("fileA.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("fileB.raw", "ACDEF", "ACDEF", 1.0, 2, [modifiedPeptide]);
        psm2.Intensities = [2_000_000.0];

        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm1, psm2]);
        var header = group.GetTabSeparatedHeader();
        var row = group.ToString();

        int headerCols = header.Split('\t').Length;
        int rowCols = row.Split('\t').Length;
        Assert.That(rowCols, Is.EqualTo(headerCols), $"Header has {headerCols} columns but row has {rowCols}");
    }

    #endregion

    #region Output

    [Test]
    public void GetTabSeparatedHeader_ContainsExpectedColumns()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm]);
        var header = group.GetTabSeparatedHeader();

        Assert.That(header, Does.Contain("Peptide Sequence"));
        Assert.That(header, Does.Contain("Start Residue"));
        Assert.That(header, Does.Contain("End Residue"));
        Assert.That(header, Does.Contain("Parent Accessions"));
        Assert.That(header, Does.Contain("SpectralCount_"));
        Assert.That(header, Does.Contain("CountOccupancy_"));
    }

    [Test]
    public void ToString_ContainsExpectedFields()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm]);
        var output = group.ToString();

        var parts = output.Split('\t');
        Assert.That(parts[0], Is.EqualTo("ACDEF"));
        Assert.That(parts[1], Is.EqualTo("1"));
        Assert.That(parts[2], Is.EqualTo("5"));
        Assert.That(parts[3], Is.EqualTo("P00001"));
    }

    #endregion

    #region OccupancyRollupStrategy

    [Test]
    public void OccupancyRollupStrategy_DefaultIsSum()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm]);

        Assert.That(group.OccupancyRollupStrategy, Is.EqualTo(IntensityRollupStrategy.Sum));
    }

    [Test]
    public void OccupancyRollupStrategy_SetInvalidatesSampleGroupResults()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 1);

        var group = new PeptideGroup("ACDEF", 1, 5, [protein], [psm]);
        group.PopulateSampleGroupResults();
        Assert.That(group.SampleGroupResults, Is.Not.Null);

        group.OccupancyRollupStrategy = IntensityRollupStrategy.Mean;
        Assert.That(group.SampleGroupResults, Is.Null);
    }

    #endregion
}
