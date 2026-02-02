using System.ComponentModel;
using System.Text.RegularExpressions;

namespace Predictions.Koina.SupportedModels
{
    /// <summary>
    /// Koina model for indexed retention time (iRT) prediction using Prosit 2020 with TMT support.
    /// Predicts normalized retention times for peptides with TMT, iTRAQ, and SILAC modifications.
    /// </summary>
    public class Prosit2020iRTTMT : PrositModelBase
    {
        public override string ModelName => "Prosit_2020_irt_TMT";

        /// <summary>
        /// Extended modification mapping including TMT, iTRAQ, and SILAC labels.
        /// </summary>
        public override Dictionary<string, string> ValidModificationUnimodMapping => new()
        {
            // Standard modifications
            { "[Common Variable:Oxidation on M]", "[UNIMOD:35]" },
            { "[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]" },
            // SILAC modifications
            { "[Common Variable:Label:13C(6)15N(2) on K]", "[UNIMOD:259]" },
            { "[Common Variable:Label:13C(6)15N(4) on R]", "[UNIMOD:267]" },
            // TMT modifications
            { "[Common Fixed:TMT6plex on K]", "[UNIMOD:737]" },
            { "[Common Fixed:TMT6plex on N-terminus]", "[UNIMOD:737]-" },
            // TMTpro modifications
            { "[Common Fixed:TMTpro on K]", "[UNIMOD:2016]" },
            { "[Common Fixed:TMTpro on N-terminus]", "[UNIMOD:2016]-" },
            // iTRAQ 4-plex modifications
            { "[Common Fixed:iTRAQ4plex on K]", "[UNIMOD:214]" },
            { "[Common Fixed:iTRAQ4plex on N-terminus]", "[UNIMOD:214]-" },
            // iTRAQ 8-plex modifications
            { "[Common Fixed:iTRAQ8plex on K]", "[UNIMOD:730]" },
            { "[Common Fixed:iTRAQ8plex on N-terminus]", "[UNIMOD:730]-" }
        };

        /// <summary>
        /// Predicted indexed retention times for each peptide.
        /// Values are in the same order as the validated PeptideSequences.
        /// </summary>
        public List<double> PredictedIndexedRetentionTimes { get; private set; } = new();

        /// <summary>
        /// Creates a new Prosit 2020 iRT TMT prediction model.
        /// </summary>
        /// <param name="peptideSequences">List of peptide sequences to predict retention times for.</param>
        /// <param name="warnings">Output parameter containing any validation warnings.</param>
        public Prosit2020iRTTMT(List<string> peptideSequences, out WarningException? warnings)
        {
            warnings = ValidateAndAddSequences(
                peptideSequences,
                isValidSequence: IsValidPrositSequence,
                transformSequence: ConvertToPrositModificationFormat);
        }

        /// <summary>
        /// Extended validation for TMT model - ensures N-terminal modifications are valid.
        /// </summary>
        public override bool HasValidModifications(string sequence)
        {
            var matches = Regex.Matches(sequence, ModificationPattern);
            if (matches.Count == 0)
                return true;

            // Check all modifications are in the valid list
            if (!matches.All(m => ValidModificationUnimodMapping.ContainsKey(m.Value)))
                return false;

            // If there's a modification at position 0, it must be a valid N-terminal mod
            var firstMatch = matches.FirstOrDefault(m => m.Index == 0);
            if (firstMatch != null)
            {
                return ValidModificationUnimodMapping.ContainsKey(firstMatch.Value);
            }

            return true;
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