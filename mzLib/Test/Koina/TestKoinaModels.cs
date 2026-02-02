using NUnit.Framework;
using Predictions.Koina.SupportedModels;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.Koina
{
    /// <summary>
    /// Tests for Koina prediction model infrastructure.
    /// Focuses on critical validation, transformation, and batching functionality.
    /// </summary>
    [TestFixture]
    public class TestKoinaModels
    {
        #region KoinaModelBase Tests

        /// <summary>
        /// TEST: IsValidBaseSequence with valid sequences
        /// WHY: Core validation must accept standard amino acid sequences within length limits
        /// </summary>
        [Test]
        public void IsValidBaseSequence_ValidSequences_ReturnsTrue()
        {
            var model = new PFly2024FineTuned(new List<string> { "PEPTIDE" }, out _);

            Assert.IsTrue(model.IsValidBaseSequence("PEPTIDE"));
            Assert.IsTrue(model.IsValidBaseSequence("ACDEFGHIKLMNPQRSTVWY")); // All valid AAs
            Assert.IsTrue(model.IsValidBaseSequence("A")); // Minimum length
        }

        /// <summary>
        /// TEST: IsValidBaseSequence with invalid sequences
        /// WHY: Core validation must reject sequences with invalid characters or exceeding length
        /// </summary>
        [Test]
        public void IsValidBaseSequence_InvalidSequences_ReturnsFalse()
        {
            var model = new PFly2024FineTuned(new List<string> { "PEPTIDE" }, out _);

            Assert.IsFalse(model.IsValidBaseSequence("")); // Empty
            Assert.IsFalse(model.IsValidBaseSequence("PEPTXDE")); // Invalid AA 'X'
            Assert.IsFalse(model.IsValidBaseSequence("PEPTBDE")); // Invalid AA 'B'
            Assert.IsFalse(model.IsValidBaseSequence("peptide")); // Lowercase
            Assert.IsFalse(model.IsValidBaseSequence("PEP TIDE")); // Space
        }

        /// <summary>
        /// TEST: IsValidBaseSequence respects MaxPeptideLength
        /// WHY: Different models have different length limits that must be enforced
        /// </summary>
        [Test]
        public void IsValidBaseSequence_RespectsMaxPeptideLength()
        {
            // PFly has MaxPeptideLength = 40
            var pfly = new PFly2024FineTuned(new List<string> { "PEPTIDE" }, out _);
            Assert.IsTrue(pfly.IsValidBaseSequence(new string('A', 40)));
            Assert.IsFalse(pfly.IsValidBaseSequence(new string('A', 41)));

            // Prosit models have MaxPeptideLength = 30
            var prosit = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);
            Assert.IsTrue(prosit.IsValidBaseSequence(new string('A', 30)));
            Assert.IsFalse(prosit.IsValidBaseSequence(new string('A', 31)));
        }

        /// <summary>
        /// TEST: IsValidBaseSequence handles sequences with modifications
        /// WHY: Validation must work on base sequence after stripping modification annotations
        /// </summary>
        [Test]
        public void IsValidBaseSequence_WithModifications_ValidatesBaseSequence()
        {
            var model = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            // Should validate the base sequence, ignoring mod annotations
            Assert.IsTrue(model.IsValidBaseSequence("PEPTM[Common Variable:Oxidation on M]IDE"));
            Assert.IsTrue(model.IsValidBaseSequence("[SomeMod]PEPTIDE"));
        }

        /// <summary>
        /// TEST: ToBatchedRequests creates correct batch structure
        /// WHY: Batching is critical for efficient API communication with Koina server
        /// </summary>
        [Test]
        public void ToBatchedRequests_CreatesCorrectStructure()
        {
            var peptides = new List<string> { "PEPTIDE", "SEQUENCE", "ANOTHER" };
            var model = new PFly2024FineTuned(peptides, out _);

            var requests = model.ToBatchedRequests();

            Assert.AreEqual(1, requests.Count); // All fit in one batch
            Assert.IsTrue(requests[0].ContainsKey("id"));
            Assert.IsTrue(requests[0].ContainsKey("inputs"));
        }

        /// <summary>
        /// TEST: ToBatchedRequests respects MaxBatchSize
        /// WHY: Large requests must be split to avoid server limits and timeouts
        /// </summary>
        [Test]
        public void ToBatchedRequests_RespectsMaxBatchSize()
        {
            // PFly has MaxBatchSize = 128
            var peptides = Enumerable.Range(0, 200).Select(i => "PEPTIDE").ToList();
            var model = new PFly2024FineTuned(peptides, out _);

            var requests = model.ToBatchedRequests();

            Assert.AreEqual(2, requests.Count); // 200 peptides / 128 max = 2 batches
        }

        /// <summary>
        /// TEST: Constructor filters invalid sequences and returns warnings
        /// WHY: Invalid inputs must be filtered with clear feedback to the user
        /// </summary>
        [Test]
        public void Constructor_InvalidSequences_FiltersAndWarns()
        {
            var peptides = new List<string>
            {
                "PEPTIDE",      // Valid
                "INVALID123",   // Invalid - numbers
                "ACDEPKR",      // Valid
                ""              // Invalid - empty
            };

            var model = new PFly2024FineTuned(peptides, out var warnings);

            Assert.AreEqual(2, model.PeptideSequences.Count);
            Assert.IsNotNull(warnings);
            Assert.IsTrue(warnings.Message.Contains("invalid"));
        }

        /// <summary>
        /// TEST: Constructor with empty input returns warning
        /// WHY: Empty inputs should not throw but should warn the user
        /// </summary>
        [Test]
        public void Constructor_EmptyInput_ReturnsWarning()
        {
            var model = new PFly2024FineTuned(new List<string>(), out var warnings);

            Assert.AreEqual(0, model.PeptideSequences.Count);
            Assert.IsNotNull(warnings);
            Assert.IsTrue(warnings.Message.Contains("empty"));
        }

        #endregion

        #region PrositModelBase Tests

        /// <summary>
        /// TEST: HasValidModifications accepts valid Prosit modifications
        /// WHY: Only specific modifications are supported by Prosit models
        /// </summary>
        [Test]
        public void HasValidModifications_ValidMods_ReturnsTrue()
        {
            var model = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            Assert.IsTrue(model.HasValidModifications("PEPTIDE")); // No mods
            Assert.IsTrue(model.HasValidModifications("PEPTM[Common Variable:Oxidation on M]IDE"));
            Assert.IsTrue(model.HasValidModifications("PEPTIC[Common Fixed:Carbamidomethyl on C]DE"));
        }

        /// <summary>
        /// TEST: HasValidModifications rejects unsupported modifications
        /// WHY: Unsupported mods would cause prediction errors or incorrect results
        /// </summary>
        [Test]
        public void HasValidModifications_InvalidMods_ReturnsFalse()
        {
            var model = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            Assert.IsFalse(model.HasValidModifications("PEPT[Phospho]IDE"));
            Assert.IsFalse(model.HasValidModifications("PEPT[Unknown Mod]IDE"));
        }

        /// <summary>
        /// TEST: ConvertToPrositModificationFormat converts mzLib format to UNIMOD
        /// WHY: Prosit server requires UNIMOD format for modifications
        /// </summary>
        [Test]
        public void ConvertToPrositModificationFormat_ConvertsToUnimod()
        {
            var model = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            var result = model.ConvertToPrositModificationFormat(
                "PEPTM[Common Variable:Oxidation on M]IDE");

            Assert.IsTrue(result.Contains("[UNIMOD:35]"));
            Assert.IsFalse(result.Contains("Common Variable:Oxidation on M"));
        }

        /// <summary>
        /// TEST: ConvertToPrositModificationFormat carbamidomethylates cysteines
        /// WHY: Prosit expects all cysteines to be carbamidomethylated
        /// </summary>
        [Test]
        public void ConvertToPrositModificationFormat_CarbamidomethylatesCysteines()
        {
            var model = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            var result = model.ConvertToPrositModificationFormat("PEPTCIDE");

            Assert.IsTrue(result.Contains("C[UNIMOD:4]"));
        }

        /// <summary>
        /// TEST: ConvertToPrositModificationFormat doesn't double-modify cysteines
        /// WHY: Already modified cysteines should not be carbamidomethylated again
        /// </summary>
        [Test]
        public void ConvertToPrositModificationFormat_DoesNotDoubleModifyCysteines()
        {
            var model = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            var result = model.ConvertToPrositModificationFormat(
                "PEPTC[Common Fixed:Carbamidomethyl on C]IDE");

            // Should have exactly one [UNIMOD:4], not two
            var count = result.Split("[UNIMOD:4]").Length - 1;
            Assert.AreEqual(1, count);
        }

        /// <summary>
        /// TEST: ConvertToMzLibModificationFormat reverses UNIMOD to mzLib format
        /// WHY: Results need to be converted back for downstream mzLib processing
        /// </summary>
        [Test]
        public void ConvertToMzLibModificationFormat_ReversesUnimod()
        {
            var model = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            var result = model.ConvertToMzLibModificationFormat(
                "PEPTM[UNIMOD:35]IDE");

            Assert.IsTrue(result.Contains("[Common Variable:Oxidation on M]"));
            Assert.IsFalse(result.Contains("UNIMOD:35"));
        }

        /// <summary>
        /// TEST: ConvertToMzLibModificationFormatWithMassesOnly creates mass-based notation
        /// WHY: Peptide objects require numeric mass values, not named modifications
        /// </summary>
        [Test]
        public void ConvertToMzLibModificationFormatWithMassesOnly_CreatesMassNotation()
        {
            var model = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            var result = model.ConvertToMzLibModificationFormatWithMassesOnly(
                "PEPTM[Common Variable:Oxidation on M]IDE");

            Assert.IsTrue(result.Contains("[15.994915]"));
        }

        #endregion

        #region Prosit2020iRTTMT Extended Modification Tests

        /// <summary>
        /// TEST: Prosit2020iRTTMT accepts TMT modifications
        /// WHY: This model specifically supports TMT labeling for quantitative proteomics
        /// </summary>
        [Test]
        public void Prosit2020iRTTMT_AcceptsTMTModifications()
        {
            var model = new Prosit2020iRTTMT(new List<string> { "PEPTIDE" }, out _);

            Assert.IsTrue(model.HasValidModifications(
                "[Common Fixed:TMT6plex on N-terminus]PEPTK[Common Fixed:TMT6plex on K]IDE"));
            Assert.IsTrue(model.HasValidModifications(
                "[Common Fixed:TMTpro on N-terminus]PEPTIDE"));
        }

        /// <summary>
        /// TEST: Prosit2020iRTTMT accepts iTRAQ modifications
        /// WHY: iTRAQ is another common labeling strategy this model supports
        /// </summary>
        [Test]
        public void Prosit2020iRTTMT_AcceptsiTRAQModifications()
        {
            var model = new Prosit2020iRTTMT(new List<string> { "PEPTIDE" }, out _);

            Assert.IsTrue(model.HasValidModifications(
                "[Common Fixed:iTRAQ4plex on N-terminus]PEPTIDE"));
            Assert.IsTrue(model.HasValidModifications(
                "[Common Fixed:iTRAQ8plex on K]PEPTIDE"));
        }

        /// <summary>
        /// TEST: Prosit2020iRTTMT accepts SILAC modifications
        /// WHY: SILAC labeling is supported for metabolic labeling experiments
        /// </summary>
        [Test]
        public void Prosit2020iRTTMT_AcceptsSILACModifications()
        {
            var model = new Prosit2020iRTTMT(new List<string> { "PEPTIDE" }, out _);

            Assert.IsTrue(model.HasValidModifications(
                "PEPTK[Common Variable:Label:13C(6)15N(2) on K]IDE"));
            Assert.IsTrue(model.HasValidModifications(
                "PEPTIDER[Common Variable:Label:13C(6)15N(4) on R]"));
        }

        #endregion

        #region Prosit2020IntensityHCD Tests

        /// <summary>
        /// TEST: Prosit2020IntensityHCD validates charge states
        /// WHY: Only charges 1-6 are supported by the model
        /// </summary>
        [Test]
        public void Prosit2020IntensityHCD_ValidatesChargeStates()
        {
            var peptides = new List<string> { "PEPTIDE", "PEPTIDE", "PEPTIDE" };
            var charges = new List<int> { 2, 7, 3 }; // 7 is invalid
            var energies = new List<int> { 30, 30, 30 };
            var rts = new List<double?> { 10.0, 10.0, 10.0 };

            var model = new Prosit2020IntensityHCD(
                peptides, charges, energies, rts, out var warnings);

            Assert.AreEqual(2, model.PeptideSequences.Count); // Only 2 valid
            Assert.AreEqual(2, model.PrecursorCharges.Count);
            Assert.IsNotNull(warnings);
        }

        /// <summary>
        /// TEST: Prosit2020IntensityHCD validates collision energies
        /// WHY: Collision energy must be positive for meaningful predictions
        /// </summary>
        [Test]
        public void Prosit2020IntensityHCD_ValidatesCollisionEnergies()
        {
            var peptides = new List<string> { "PEPTIDE", "PEPTIDE" };
            var charges = new List<int> { 2, 2 };
            var energies = new List<int> { 30, 0 }; // 0 is invalid
            var rts = new List<double?> { 10.0, 10.0 };

            var model = new Prosit2020IntensityHCD(
                peptides, charges, energies, rts, out var warnings);

            Assert.AreEqual(1, model.PeptideSequences.Count);
            Assert.IsNotNull(warnings);
        }

        /// <summary>
        /// TEST: Prosit2020IntensityHCD requires matching input list lengths
        /// WHY: Parallel lists must be aligned for correct peptide-charge-energy association
        /// </summary>
        [Test]
        public void Prosit2020IntensityHCD_RequiresMatchingListLengths()
        {
            var peptides = new List<string> { "PEPTIDE", "SEQUENCE" };
            var charges = new List<int> { 2 }; // Mismatched length
            var energies = new List<int> { 30, 30 };
            var rts = new List<double?> { 10.0, 10.0 };

            Assert.Throws<ArgumentException>(() =>
                new Prosit2020IntensityHCD(peptides, charges, energies, rts, out _));
        }

        /// <summary>
        /// TEST: Prosit2020IntensityHCD ToBatchedRequests includes charges and energies
        /// WHY: MS2 predictions require charge and collision energy for each peptide
        /// </summary>
        [Test]
        public void Prosit2020IntensityHCD_ToBatchedRequests_IncludesAllInputs()
        {
            var peptides = new List<string> { "PEPTIDE" };
            var charges = new List<int> { 2 };
            var energies = new List<int> { 30 };
            var rts = new List<double?> { 10.0 };

            var model = new Prosit2020IntensityHCD(
                peptides, charges, energies, rts, out _);

            var requests = model.ToBatchedRequests();
            var inputs = requests[0]["inputs"] as List<object>;

            Assert.AreEqual(3, inputs.Count); // peptides, charges, energies
        }

        #endregion

        #region Model Property Tests

        /// <summary>
        /// TEST: Each model has correct ModelName
        /// WHY: ModelName must match the Koina server endpoint exactly
        /// </summary>
        [Test]
        public void Models_HaveCorrectModelNames()
        {
            var pfly = new PFly2024FineTuned(new List<string> { "PEPTIDE" }, out _);
            var prosit2019 = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);
            var prositTmt = new Prosit2020iRTTMT(new List<string> { "PEPTIDE" }, out _);
            var prositHcd = new Prosit2020IntensityHCD(
                new List<string> { "PEPTIDE" },
                new List<int> { 2 },
                new List<int> { 30 },
                new List<double?> { null }, out _);

            Assert.AreEqual("pfly_2024_fine_tuned", pfly.ModelName);
            Assert.AreEqual("Prosit_2019_irt", prosit2019.ModelName);
            Assert.AreEqual("Prosit_2020_irt_TMT", prositTmt.ModelName);
            Assert.AreEqual("Prosit_2020_intensity_HCD", prositHcd.ModelName);
        }

        /// <summary>
        /// TEST: Each model has appropriate MaxBatchSize
        /// WHY: Batch sizes affect performance and must not exceed server limits
        /// </summary>
        [Test]
        public void Models_HaveCorrectMaxBatchSize()
        {
            var pfly = new PFly2024FineTuned(new List<string> { "PEPTIDE" }, out _);
            var prosit = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            Assert.AreEqual(128, pfly.MaxBatchSize);
            Assert.AreEqual(1000, prosit.MaxBatchSize);
        }

        /// <summary>
        /// TEST: Each model has appropriate MaxPeptideLength
        /// WHY: Models are trained on specific peptide length ranges
        /// </summary>
        [Test]
        public void Models_HaveCorrectMaxPeptideLength()
        {
            var pfly = new PFly2024FineTuned(new List<string> { "PEPTIDE" }, out _);
            var prosit = new Prosit2019iRT(new List<string> { "PEPTIDE" }, out _);

            Assert.AreEqual(40, pfly.MaxPeptideLength);
            Assert.AreEqual(30, prosit.MaxPeptideLength);
        }

        #endregion
    }
}