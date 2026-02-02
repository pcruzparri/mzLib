using System.Text.RegularExpressions;

namespace Predictions.Koina.SupportedModels
{
    /// <summary>
    /// Abstract base class for Prosit models providing UNIMOD modification conversion
    /// and cysteine carbamidomethylation.
    /// </summary>
    public abstract class PrositModelBase : KoinaModelBase
    {
        #region Modification Handling

        /// <summary>
        /// Mapping from mzLib modification format to Prosit UNIMOD format.
        /// Override in derived classes to add model-specific modifications (e.g., TMT, iTRAQ).
        /// </summary>
        public virtual Dictionary<string, string> ValidModificationUnimodMapping => new()
        {
            { "[Common Variable:Oxidation on M]", "[UNIMOD:35]" },
            { "[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]" }
        };

        /// <summary>
        /// Mapping from mzLib modification format to monoisotopic mass values.
        /// Used for creating Peptide objects with mass-based modification notation.
        /// Override in derived classes to add model-specific modifications.
        /// </summary>
        public virtual Dictionary<string, double> ValidModificationsMonoisotopicMasses => new()
        {
            { "[Common Variable:Oxidation on M]", 15.994915 },
            { "[Common Fixed:Carbamidomethyl on C]", 57.021464 }
        };

        /// <summary>
        /// Whether to automatically carbamidomethylate unmodified cysteines.
        /// Default is true for Prosit models.
        /// </summary>
        protected virtual bool CarbamidomethylateCysteines => true;

        /// <summary>
        /// Validates that all modifications in the sequence are supported by this model.
        /// </summary>
        /// <param name="sequence">The peptide sequence to validate.</param>
        /// <returns>True if all modifications are valid; otherwise, false.</returns>
        public virtual bool HasValidModifications(string sequence)
        {
            var matches = Regex.Matches(sequence, ModificationPattern);
            if (matches.Count == 0)
                return true; // No modifications is valid

            return matches.All(m => ValidModificationUnimodMapping.ContainsKey(m.Value));
        }

        /// <summary>
        /// Validates that a sequence meets all Prosit requirements (base sequence + modifications).
        /// </summary>
        /// <param name="sequence">The peptide sequence to validate.</param>
        /// <returns>True if the sequence is valid for this Prosit model; otherwise, false.</returns>
        public virtual bool IsValidPrositSequence(string sequence)
            => IsValidBaseSequence(sequence) && HasValidModifications(sequence);

        /// <summary>
        /// Converts a peptide sequence from the mzLib modification format to the Prosit UNIMOD format.
        /// By default, all unmodified cysteines are carbamidomethylated (UNIMOD:4) to match the 
        /// expectations of Prosit models.
        /// </summary>
        /// <param name="sequence">Peptide sequence in mzLib modification format.</param>
        /// <returns>The sequence converted to Prosit UNIMOD format with cysteines carbamidomethylated.</returns>
        public virtual string ConvertToPrositModificationFormat(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Key, mod.Value);
            }

            if (CarbamidomethylateCysteines)
            {
                // Carbamidomethylate all cysteines that are not already modified
                sequence = Regex.Replace(sequence, @"C(?!\[UNIMOD:4\])", "C[UNIMOD:4]");
            }

            return sequence;
        }

        /// <summary>
        /// Converts a peptide sequence from Prosit UNIMOD format back to mzLib modification format.
        /// </summary>
        /// <param name="sequence">Peptide sequence in Prosit UNIMOD format.</param>
        /// <returns>The sequence converted to mzLib modification format.</returns>
        public virtual string ConvertToMzLibModificationFormat(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Value, mod.Key);
            }
            return sequence;
        }

        /// <summary>
        /// Converts a peptide sequence from mzLib modification format to mass-only notation.
        /// Useful for creating Peptide objects that require numeric mass values.
        /// </summary>
        /// <param name="sequence">Peptide sequence in mzLib modification format.</param>
        /// <returns>The sequence with modifications replaced by their monoisotopic masses.</returns>
        public virtual string ConvertToMzLibModificationFormatWithMassesOnly(string sequence)
        {
            foreach (var mod in ValidModificationsMonoisotopicMasses)
            {
                sequence = sequence.Replace(mod.Key, $"[{mod.Value:F6}]");
            }
            return sequence;
        }

        #endregion
    }
}