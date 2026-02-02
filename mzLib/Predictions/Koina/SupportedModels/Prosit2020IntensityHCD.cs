using Chemistry;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Proteomics.AminoAcidPolymer;
using Readers.SpectralLibrary;
using System.ComponentModel;

namespace Predictions.Koina.SupportedModels
{
    /// <summary>
    /// Koina model for MS2 fragment intensity prediction using Prosit 2020 HCD.
    /// Predicts fragment ion intensities and creates spectral library entries.
    /// 
    /// Model details: https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_HCD/infer
    /// </summary>
    public class Prosit2020IntensityHCD : PrositModelBase
    {
        #region Model-Specific Properties

        public override string ModelName => "Prosit_2020_intensity_HCD";

        /// <summary>
        /// Valid precursor charge states for this model.
        /// </summary>
        public HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };

        /// <summary>
        /// Number of fragment ions predicted per peptide (b and y ions, charges 1-3, up to length 29).
        /// </summary>
        public int NumberOfPredictedFragmentIons => 174;

        /// <summary>
        /// Minimum intensity threshold for including fragment ions in the library spectrum.
        /// </summary>
        public double MinIntensityFilter { get; }

        #endregion

        #region Additional Input Data

        /// <summary>
        /// Precursor charge states for each peptide (parallel to PeptideSequences).
        /// </summary>
        public List<int> PrecursorCharges { get; } = new();

        /// <summary>
        /// HCD collision energies for each peptide (parallel to PeptideSequences).
        /// Model performs best for collision energies 20, 23, 25, 28, 30, and 35.
        /// </summary>
        public List<int> CollisionEnergies { get; } = new();

        /// <summary>
        /// Retention times for each peptide (parallel to PeptideSequences).
        /// Used for LibrarySpectrum creation.
        /// </summary>
        public List<double?> RetentionTimes { get; } = new();

        #endregion

        #region Output Data

        /// <summary>
        /// Predicted spectral library entries for each valid input.
        /// </summary>
        public List<LibrarySpectrum> PredictedSpectra { get; private set; } = new();

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new Prosit 2020 HCD intensity prediction model.
        /// All input lists must have the same length.
        /// </summary>
        /// <param name="peptideSequences">Peptide sequences with valid modifications.</param>
        /// <param name="precursorCharges">Precursor charge states (1-6).</param>
        /// <param name="collisionEnergies">HCD collision energies (best: 20, 23, 25, 28, 30, 35).</param>
        /// <param name="retentionTimes">Retention times for library spectrum creation.</param>
        /// <param name="warnings">Output parameter containing any validation warnings.</param>
        /// <param name="minIntensityFilter">Minimum intensity threshold for fragment ions (default: 1e-4).</param>
        /// <exception cref="ArgumentException">Thrown when input lists have different lengths.</exception>
        public Prosit2020IntensityHCD(
            List<string> peptideSequences,
            List<int> precursorCharges,
            List<int> collisionEnergies,
            List<double?> retentionTimes,
            out WarningException? warnings,
            double minIntensityFilter = 1e-4)
        {
            MinIntensityFilter = minIntensityFilter;

            // Verify input lists are of the same length
            if (peptideSequences.Count != precursorCharges.Count
                || precursorCharges.Count != collisionEnergies.Count
                || collisionEnergies.Count != retentionTimes.Count)
            {
                throw new ArgumentException("Input lists must have the same length.");
            }

            if (peptideSequences.Count == 0)
            {
                warnings = new WarningException("Inputs were empty. No predictions will be made.");
                return;
            }

            // Validate and add entries
            var invalidArguments = new List<string>();
            for (int i = 0; i < peptideSequences.Count; i++)
            {
                var peptide = peptideSequences[i];
                var charge = precursorCharges[i];
                var energy = collisionEnergies[i];
                var retentionTime = retentionTimes[i];

                if (!IsValidPrositSequence(peptide) ||
                    !AllowedPrecursorCharges.Contains(charge) ||
                    energy <= 0)
                {
                    invalidArguments.Add($"Index {i}: Peptide '{peptide}' (Length: {GetBaseSequence(peptide).Length}), Charge: {charge}, Collision Energy: {energy}");
                }
                else
                {
                    PeptideSequences.Add(ConvertToPrositModificationFormat(peptide));
                    PrecursorCharges.Add(charge);
                    CollisionEnergies.Add(energy);
                    RetentionTimes.Add(retentionTime);
                }
            }

            warnings = null;
            if (invalidArguments.Count > 0)
            {
                warnings = new WarningException(
                    "The following input entries are invalid and will be skipped:\n"
                    + string.Join("\n", invalidArguments)
                    + "\nModel Requirements:\n"
                    + $"- Peptide length <= {MaxPeptideLength}\n"
                    + "- Peptide sequence is not empty\n"
                    + "- Peptide has valid modifications\n"
                    + $"- Precursor charge in [{string.Join(", ", AllowedPrecursorCharges)}]\n"
                    + "- Collision energy > 0");
            }
        }

        /// <summary>
        /// Creates a model from existing library spectra (not yet implemented).
        /// </summary>
        public Prosit2020IntensityHCD(List<LibrarySpectrum> spectralLibrary)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Creates a model from a spectral library file (not yet implemented).
        /// </summary>
        public Prosit2020IntensityHCD(string filePath)
        {
            throw new NotImplementedException();
        }

        #endregion

        #region Overrides

        /// <summary>
        /// Creates batched requests including peptide sequences, charges, and collision energies.
        /// </summary>
        public override List<Dictionary<string, object>> ToBatchedRequests()
        {
            var batchedPeptides = PeptideSequences.Chunk(MaxBatchSize).ToList();
            var batchedCharges = PrecursorCharges.Chunk(MaxBatchSize).ToList();
            var batchedEnergies = CollisionEnergies.Chunk(MaxBatchSize).ToList();

            var batchedRequests = new List<Dictionary<string, object>>();

            for (int i = 0; i < batchedPeptides.Count; i++)
            {
                var request = new Dictionary<string, object>
                {
                    { "id", $"Batch{i}_" + Guid.NewGuid() },
                    { "inputs", new List<object>
                        {
                            new {
                                name = "peptide_sequences",
                                shape = new[] { batchedPeptides[i].Length, 1 },
                                datatype = "BYTES",
                                data = batchedPeptides[i]
                            },
                            new {
                                name = "precursor_charges",
                                shape = new[] { batchedCharges[i].Length, 1 },
                                datatype = "INT32",
                                data = batchedCharges[i]
                            },
                            new {
                                name = "collision_energies",
                                shape = new[] { batchedEnergies[i].Length, 1 },
                                datatype = "FP32",
                                data = batchedEnergies[i]
                            }
                        }
                    }
                };
                batchedRequests.Add(request);
            }

            return batchedRequests;
        }

        /// <inheritdoc/>
        protected override void ProcessResponses(string[] responses)
        {
            if (PeptideSequences.Count == 0)
                return;

            var deserializedResponses = DeserializeResponses(responses);

            for (int batchIndex = 0; batchIndex < deserializedResponses.Count; batchIndex++)
            {
                var responseBatch = deserializedResponses[batchIndex];
                if (responseBatch.Outputs.Count != 3)
                {
                    throw new Exception($"API response is not in the expected format. Expected 3 outputs, got {responseBatch.Outputs.Count}.");
                }

                var currentBatchSize = responseBatch.Outputs[0].Shape[0];
                var outputIonAnnotations = responseBatch.Outputs[0].Data.Select(d => (string)d).ToArray();
                var outputMZs = responseBatch.Outputs[1].Data.Select(d => Convert.ToDouble(d)).ToArray();
                var outputIntensities = responseBatch.Outputs[2].Data.Select(d => Convert.ToDouble(d)).ToArray();

                for (int precursorIndex = 0; precursorIndex < currentBatchSize; precursorIndex++)
                {
                    int linearIndex = batchIndex * MaxBatchSize + precursorIndex;
                    var spectrum = CreateLibrarySpectrum(
                        linearIndex,
                        precursorIndex,
                        outputIonAnnotations,
                        outputMZs,
                        outputIntensities);
                    PredictedSpectra.Add(spectrum);
                }
            }

            // Check for duplicates
            var unique = PredictedSpectra.DistinctBy(p => p.Name).ToList();
            if (unique.Count != PredictedSpectra.Count)
            {
                throw new WarningException(
                    $"Duplicate spectra found in predictions. Reduced from {PredictedSpectra.Count} predicted spectra to {unique.Count} unique spectra.");
            }
        }

        #endregion

        #region Helper Methods

        /// <summary>
        /// Creates a LibrarySpectrum from the model output for a single peptide.
        /// </summary>
        private LibrarySpectrum CreateLibrarySpectrum(
            int linearIndex,
            int precursorIndex,
            string[] ionAnnotations,
            double[] mzValues,
            double[] intensities)
        {
            var peptideSequence = PeptideSequences[linearIndex];
            var peptide = new Peptide(
                ConvertToMzLibModificationFormatWithMassesOnly(
                    ConvertToMzLibModificationFormat(peptideSequence)));

            var fragmentIons = new List<MatchedFragmentIon>();
            int startIndex = precursorIndex * NumberOfPredictedFragmentIons;
            int endIndex = startIndex + NumberOfPredictedFragmentIons;

            for (int fragmentIndex = startIndex; fragmentIndex < endIndex; fragmentIndex++)
            {
                // Skip impossible ions (-1) and low intensity peaks
                if (mzValues[fragmentIndex] == -1 || intensities[fragmentIndex] < MinIntensityFilter)
                    continue;

                var fragmentIon = ParseFragmentIon(
                    ionAnnotations[fragmentIndex],
                    mzValues[fragmentIndex],
                    intensities[fragmentIndex]);

                fragmentIons.Add(fragmentIon);
            }

            return new LibrarySpectrum(
                sequence: peptideSequence,
                precursorMz: peptide.ToMz(PrecursorCharges[linearIndex]),
                chargeState: PrecursorCharges[linearIndex],
                peaks: fragmentIons,
                rt: RetentionTimes[linearIndex]);
        }

        /// <summary>
        /// Parses a Prosit ion annotation (e.g., "b5+1") into a MatchedFragmentIon.
        /// </summary>
        private static MatchedFragmentIon ParseFragmentIon(string annotation, double mz, double intensity)
        {
            // Parse annotation like "b5+1" or "y12+2"
            var ionType = annotation[0].ToString(); // 'b' or 'y'
            var plusIndex = annotation.IndexOf('+');
            var fragmentNumber = int.Parse(annotation.Substring(1, plusIndex - 1));
            var fragmentCharge = int.Parse(annotation.Substring(plusIndex + 1));

            return new MatchedFragmentIon(
                neutralTheoreticalProduct: new Product(
                    productType: Enum.Parse<ProductType>(ionType),
                    terminus: ionType == "b" ? FragmentationTerminus.N : FragmentationTerminus.C,
                    neutralMass: 0.0, // Placeholder - not directly provided by Prosit
                    fragmentNumber: fragmentNumber,
                    residuePosition: fragmentNumber,
                    neutralLoss: 0), // Prosit annotations don't encode neutral losses
                experMz: mz,
                experIntensity: intensity,
                charge: fragmentCharge);
        }

        /// <summary>
        /// Saves the predicted spectra to a spectral library file.
        /// </summary>
        /// <param name="filePath">Output file path.</param>
        public void SavePredictedSpectralLibrary(string filePath)
        {
            var spectralLibrary = new SpectralLibrary
            {
                Results = PredictedSpectra
            };
            spectralLibrary.WriteResults(filePath);
        }

        #endregion
    }
}