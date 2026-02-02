using System.ComponentModel;

namespace Predictions.Koina.SupportedModels
{
    /// <summary>
    /// Koina model for peptide detectability/flyability prediction.
    /// Predicts the probability of a peptide being detectable in mass spectrometry.
    /// </summary>
    public class PFly2024FineTuned : KoinaModelBase
    {
        public override string ModelName => "pfly_2024_fine_tuned";
        public override int MaxBatchSize => 128;
        public override int MaxPeptideLength => 40;

        /// <summary>
        /// Number of detectability classification categories.
        /// </summary>
        public int NumberOfDetectabilityClasses => 4;

        /// <summary>
        /// Labels for each detectability class.
        /// </summary>
        public List<string> DetectabilityClasses => new()
        {
            "Not Detectable",
            "Low Detectability",
            "Intermediate Detectability",
            "High Detectability"
        };

        /// <summary>
        /// Predicted detectability probabilities for each peptide.
        /// Each inner list contains probabilities for each detectability class.
        /// </summary>
        public List<List<double>> DetectabilityProbabilityTable { get; private set; } = new();

        /// <summary>
        /// Creates a new PFly2024 detectability prediction model.
        /// </summary>
        /// <param name="peptideSequences">List of peptide sequences to predict detectability for.</param>
        /// <param name="warnings">Output parameter containing any validation warnings.</param>
        public PFly2024FineTuned(List<string> peptideSequences, out WarningException? warnings)
        {
            warnings = ValidateAndAddSequences(
                peptideSequences,
                isValidSequence: IsValidBaseSequence);
        }

        /// <inheritdoc/>
        protected override void ProcessResponses(string[] responses)
        {
            if (PeptideSequences.Count == 0)
                return;

            var deserialized = DeserializeResponses(responses);

            DetectabilityProbabilityTable = deserialized
                .SelectMany(batch => batch.Outputs[0].Data)
                .Chunk(NumberOfDetectabilityClasses)
                .Select(chunk => chunk.Select(d => Convert.ToDouble(d)).ToList())
                .ToList();
        }
    }
}