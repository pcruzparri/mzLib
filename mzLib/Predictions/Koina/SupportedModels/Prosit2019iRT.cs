using System.ComponentModel;

namespace Predictions.Koina.SupportedModels
{
    /// <summary>
    /// Koina model for indexed retention time (iRT) prediction using Prosit 2019.
    /// Predicts normalized retention times for peptides.
    /// </summary>
    public class Prosit2019iRT : PrositModelBase
    {
        public override string ModelName => "Prosit_2019_irt";

        /// <summary>
        /// Predicted indexed retention times for each peptide.
        /// Values are in the same order as the validated PeptideSequences.
        /// </summary>
        public List<double> PredictedIndexedRetentionTimes { get; private set; } = new();

        /// <summary>
        /// Creates a new Prosit 2019 iRT prediction model.
        /// </summary>
        /// <param name="peptideSequences">List of peptide sequences to predict retention times for.</param>
        /// <param name="warnings">Output parameter containing any validation warnings.</param>
        public Prosit2019iRT(List<string> peptideSequences, out WarningException? warnings)
        {
            warnings = ValidateAndAddSequences(
                peptideSequences,
                isValidSequence: IsValidPrositSequence,
                transformSequence: ConvertToPrositModificationFormat);
        }

        /// <inheritdoc/>
        protected override void ProcessResponses(string[] responses)
        {
            if (PeptideSequences.Count == 0)
                return;

            var deserialized = DeserializeResponses(responses);

            PredictedIndexedRetentionTimes = deserialized
                .SelectMany(batch => batch.Outputs[0].Data)
                .Select(irt => Convert.ToDouble(irt))
                .ToList();
        }
    }
}