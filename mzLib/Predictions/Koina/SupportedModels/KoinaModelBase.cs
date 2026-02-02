using Easy.Common.Extensions;
using MzLibUtil;
using Predictions.Koina.Client;
using Predictions.Koina.Interfaces;
using System.ComponentModel;
using System.Text.RegularExpressions;

namespace Predictions.Koina.SupportedModels
{
    /// <summary>
    /// Abstract base class for all Koina prediction models providing common validation,
    /// batching, and HTTP inference infrastructure.
    /// </summary>
    public abstract class KoinaModelBase : IKoinaModelIO
    {
        #region Abstract/Virtual Properties (Model-Specific)

        /// <summary>
        /// The model name as registered on the Koina server.
        /// </summary>
        public abstract string ModelName { get; }

        /// <summary>
        /// Maximum number of peptides per batch request.
        /// </summary>
        public virtual int MaxBatchSize => 1000;

        /// <summary>
        /// Maximum allowed peptide length for this model.
        /// </summary>
        public virtual int MaxPeptideLength => 30;

        #endregion

        #region Common Properties

        /// <summary>
        /// Regex pattern to match modification annotations in sequences.
        /// </summary>
        public string ModificationPattern => @"\[[^\]]+\]";

        /// <summary>
        /// Regex pattern for valid canonical amino acids.
        /// </summary>
        public string CanonicalAminoAcidPattern => @"^[ACDEFGHIKLMNPQRSTVWY]+$";

        /// <summary>
        /// Validated peptide sequences ready for prediction.
        /// </summary>
        public List<string> PeptideSequences { get; } = new();

        #endregion

        #region Validation Methods

        /// <summary>
        /// Validates that the base sequence (without mods) meets length and amino acid requirements.
        /// </summary>
        /// <param name="sequence">The peptide sequence to validate (may contain modifications).</param>
        /// <returns>True if the base sequence is valid; otherwise, false.</returns>
        public virtual bool IsValidBaseSequence(string sequence)
        {
            var baseSequence = GetBaseSequence(sequence);
            return baseSequence.Length > 0
                && baseSequence.Length <= MaxPeptideLength
                && Regex.IsMatch(baseSequence, CanonicalAminoAcidPattern);
        }

        /// <summary>
        /// Extracts the base sequence by removing all modification annotations.
        /// </summary>
        /// <param name="sequence">The peptide sequence with modifications.</param>
        /// <returns>The base sequence without modifications.</returns>
        protected string GetBaseSequence(string sequence)
            => Regex.Replace(sequence, ModificationPattern, string.Empty);

        #endregion

        #region Batching Infrastructure

        /// <summary>
        /// Creates batched request dictionaries for the Koina API.
        /// Override in derived classes that need additional input parameters (e.g., charge states).
        /// </summary>
        /// <returns>List of request dictionaries ready to send to the Koina server.</returns>
        public virtual List<Dictionary<string, object>> ToBatchedRequests()
        {
            var batchedPeptides = PeptideSequences.Chunk(MaxBatchSize).ToList();
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
                            }
                        }
                    }
                };
                batchedRequests.Add(request);
            }
            return batchedRequests;
        }

        #endregion

        #region HTTP Inference

        /// <summary>
        /// Executes inference requests against the Koina server.
        /// </summary>
        public virtual async Task RunInferenceAsync()
        {
            if (PeptideSequences.Count == 0)
                return;

            // Typically a full batch takes about a minute. Setting timeout to double that for safety.
            int timeoutMinutes = PeptideSequences.Count / MaxBatchSize * 2 + 2;
            var http = new HTTP(timeoutInMinutes: timeoutMinutes);

            try
            {
                var requests = ToBatchedRequests();
                var responses = await Task.WhenAll(
                    requests.Select(request => http.InferenceRequest(ModelName, request)));
                ProcessResponses(responses);
            }
            finally
            {
                http.Dispose();
            }
        }

        /// <summary>
        /// Processes the raw JSON responses from the Koina server.
        /// Must be implemented by derived classes to parse model-specific outputs.
        /// </summary>
        /// <param name="responses">Array of JSON response strings from the Koina server.</param>
        protected abstract void ProcessResponses(string[] responses);

        /// <summary>
        /// Deserializes responses with standard null checking.
        /// </summary>
        /// <param name="responses">Array of JSON response strings.</param>
        /// <returns>List of deserialized response objects.</returns>
        /// <exception cref="Exception">Thrown when deserialization fails.</exception>
        protected List<ResponseJSONStruct> DeserializeResponses(string[] responses)
        {
            var deserialized = responses
                .Select(r => Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(r))
                .ToList();

            if (deserialized.IsNullOrEmpty() || deserialized.Any(r => r == null))
            {
                throw new Exception("Something went wrong during deserialization of responses.");
            }

            return deserialized!;
        }

        #endregion

        #region Constructor Helpers

        /// <summary>
        /// Standard validation logic for constructors. Validates input sequences and populates PeptideSequences.
        /// </summary>
        /// <param name="inputSequences">List of peptide sequences to validate.</param>
        /// <param name="isValidSequence">Function to determine if a sequence is valid.</param>
        /// <param name="transformSequence">Optional function to transform valid sequences (e.g., modification format conversion).</param>
        /// <returns>WarningException if there were issues, null otherwise.</returns>
        protected WarningException? ValidateAndAddSequences(
            List<string> inputSequences,
            Func<string, bool> isValidSequence,
            Func<string, string>? transformSequence = null)
        {
            if (inputSequences.IsNullOrEmpty())
            {
                return new WarningException("Inputs were empty. No predictions will be made.");
            }

            var invalidSequences = new List<string>();

            foreach (var seq in inputSequences)
            {
                if (isValidSequence(seq))
                {
                    PeptideSequences.Add(transformSequence?.Invoke(seq) ?? seq);
                }
                else
                {
                    invalidSequences.Add(seq);
                }
            }

            if (invalidSequences.Count > 0)
            {
                return new WarningException(
                    $"The following peptide sequences were invalid and will be skipped: {string.Join(", ", invalidSequences)}");
            }

            return null;
        }

        #endregion
    }
}