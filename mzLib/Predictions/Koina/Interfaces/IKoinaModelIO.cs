namespace Predictions.Koina.Interfaces
{
    /// <summary>
    /// Interface for Koina prediction model input/output operations.
    /// </summary>
    public interface IKoinaModelIO
    {
        /// <summary>
        /// The model name as registered on the Koina server.
        /// </summary>
        string ModelName { get; }

        /// <summary>
        /// Maximum number of peptides per batch request.
        /// </summary>
        int MaxBatchSize { get; }

        /// <summary>
        /// Maximum allowed peptide length for this model.
        /// </summary>
        int MaxPeptideLength { get; }

        /// <summary>
        /// Creates batched request dictionaries for the Koina API.
        /// </summary>
        List<Dictionary<string, object>> ToBatchedRequests();

        /// <summary>
        /// Executes inference requests against the Koina server.
        /// </summary>
        Task RunInferenceAsync();
    }
}