﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ClassExtensions = Chemistry.ClassExtensions;
using FlashLFQ.PEP;

namespace FlashLFQ
{
    public class ChromatographicPeak : IEquatable<ChromatographicPeak>
    {
        public double Intensity;
        public double ApexRetentionTime => Apex?.IndexedPeak.RetentionTime ?? -1;
        public readonly SpectraFileInfo SpectraFileInfo;
        public List<IsotopicEnvelope> IsotopicEnvelopes;
        public int ScanCount => IsotopicEnvelopes.Count;
        public double SplitRT;
        public readonly bool IsMbrPeak;
        public DetectionType DetectionType { get; set; }
        public double PredictedRetentionTime { get; init; }
        public double MbrScore;
        public double PpmScore { get; set; }
        public double IntensityScore { get; set; }
        public double RtScore { get; set; }
        public double ScanCountScore { get; set; }
        public double IsotopicDistributionScore { get; set; }
        /// <summary>
        /// Stores the pearson correlation between the apex isotopic envelope and the theoretical isotopic distribution
        /// </summary>
        public double IsotopicPearsonCorrelation => Apex?.PearsonCorrelation ?? -1;
        public double RtPredictionError { get; set; }
        public List<int> ChargeList { get; set; }
        internal double MbrQValue { get; set; }
        public ChromatographicPeakData PepPeakData { get; set; }
        public double? MbrPep { get; set; }

        public ChromatographicPeak(Identification id, bool isMbrPeak, SpectraFileInfo fileInfo, bool randomRt = false)
        {
            SplitRT = 0;
            NumChargeStatesObserved = 0;
            MassError = double.NaN;
            NumIdentificationsByBaseSeq = 1;
            NumIdentificationsByFullSeq = 1;
            Identifications = new List<Identification>() { id };
            IsotopicEnvelopes = new List<IsotopicEnvelope>();
            IsMbrPeak = isMbrPeak;
            SpectraFileInfo = fileInfo;
            RandomRt = randomRt;
        }

        /// <summary>
        /// overloaded constructor for Isobaric_ambiguity peaks. In this case, the peak is identified by multiple identifications
        /// </summary>
        /// <param name="ids"></param>
        /// <param name="isMbrPeak"></param>
        /// <param name="fileInfo"></param>
        /// <param name="randomRt"></param>
        public ChromatographicPeak(List<Identification> ids, bool isMbrPeak, SpectraFileInfo fileInfo, double predictedRetentionTime, DetectionType detectionType = DetectionType.Imputed) :
            this(null, isMbrPeak, fileInfo)
        {
            PredictedRetentionTime = predictedRetentionTime;
            DetectionType = detectionType; // default to imputed
            Identifications = ids;
        }

        public ChromatographicPeak(Identification id, bool isMbrPeak, SpectraFileInfo fileInfo, double predictedRetentionTime, DetectionType detectionType = DetectionType.Default) :
            this(id, isMbrPeak, fileInfo)
        {
            PredictedRetentionTime = predictedRetentionTime;
            DetectionType = detectionType; // default to imputed

            if (detectionType == DetectionType.Default && isMbrPeak)
            {
                DetectionType = DetectionType.MBR;
            }
        }

        public bool Equals(ChromatographicPeak peak)
        {
            return SpectraFileInfo.Equals(peak.SpectraFileInfo) 
                && Identifications.First().ModifiedSequence.Equals(peak.Identifications.First().ModifiedSequence)
                && ApexRetentionTime == peak.ApexRetentionTime;
        }

        public IsotopicEnvelope Apex { get; private set; }
        public List<Identification> Identifications { get; private set; }
        public int NumChargeStatesObserved { get; private set; }
        public int NumIdentificationsByBaseSeq { get; private set; }
        public int NumIdentificationsByFullSeq { get; private set; }
        public double MassError { get; private set; }
        /// <summary>
        /// Bool that describes whether the retention time of this peak was randomized
        /// If true, implies that this peak is a decoy peak identified by the MBR algorithm
        /// </summary>
        public bool RandomRt { get; }
        public bool DecoyPeptide => Identifications.First().IsDecoy;

        public void CalculateIntensityForThisFeature(bool integrate)
        {
            if (IsotopicEnvelopes.Any())
            {
                Apex = IsotopicEnvelopes.MaxBy(p => p.Intensity);

                if (integrate)
                {
                    Intensity = IsotopicEnvelopes.Sum(p => p.Intensity);
                }
                else
                {
                    Intensity = Apex.Intensity;
                }

                MassError = double.NaN;

                foreach (Identification id in Identifications)
                {
                    double massErrorForId = ((ClassExtensions.ToMass(Apex.IndexedPeak.Mz, Apex.ChargeState) - id.PeakfindingMass) / id.PeakfindingMass) * 1e6;

                    if (double.IsNaN(MassError) || Math.Abs(massErrorForId) < Math.Abs(MassError))
                    {
                        MassError = massErrorForId;
                    }
                }

                NumChargeStatesObserved = IsotopicEnvelopes.Select(p => p.ChargeState).Distinct().Count();
            }
            else
            {
                Intensity = 0;
                MassError = double.NaN;
                NumChargeStatesObserved = 0;
                Apex = null;
            }
        }

        /// <summary>
        /// Merges ChromatographicPeaks by combining Identifications and IsotopicEnvelopes,
        /// then recalculates feature intensity.
        /// </summary>
        /// <param name="otherFeature"> Peak to be merged in. This peak is not modified</param>
        public void MergeFeatureWith(ChromatographicPeak otherFeature, bool integrate)
        {
            if (otherFeature != this)
            {
                var thisFeaturesPeaks = new HashSet<IIndexedMzPeak>(IsotopicEnvelopes.Select(p => p.IndexedPeak));
                this.Identifications = this.Identifications
                    .Union(otherFeature.Identifications)
                    .Distinct()
                    .ToList();
                ResolveIdentifications();
                this.IsotopicEnvelopes.AddRange(otherFeature.IsotopicEnvelopes
                    .Where(p => !thisFeaturesPeaks.Contains(p.IndexedPeak)));
                this.CalculateIntensityForThisFeature(integrate);
            }
        }

        /// <summary>
        /// Sets two ChromatographicPeak properties: NumIdentificationsByBaseSeq and NumIdentificationsByFullSeq
        /// </summary>
        public void ResolveIdentifications()
        {
            this.NumIdentificationsByBaseSeq = Identifications.Select(v => v.BaseSequence).Distinct().Count();
            this.NumIdentificationsByFullSeq = Identifications.Select(v => v.ModifiedSequence).Distinct().Count();
        }
        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("File Name" + "\t");
                sb.Append("Base Sequence" + "\t");
                sb.Append("Full Sequence" + "\t");
                sb.Append("Protein Group" + "\t");
                sb.Append("Organism" + '\t');
                sb.Append("Peptide Monoisotopic Mass" + "\t");
                sb.Append("MS2 Retention Time" + "\t");
                sb.Append("Precursor Charge" + "\t");
                sb.Append("Theoretical MZ" + "\t");
                sb.Append("Peak intensity" + "\t");
                sb.Append("Peak RT Start" + "\t");
                sb.Append("Peak RT Apex" + "\t");
                sb.Append("Peak RT End" + "\t");
                sb.Append("Peak MZ" + "\t");
                sb.Append("Peak Charge" + "\t");
                sb.Append("Num Charge States Observed" + "\t");
                sb.Append("Peak Detection Type" + "\t");
                sb.Append("PIP Q-Value" + "\t");
                sb.Append("PIP PEP" + "\t");
                sb.Append("PSMs Mapped" + "\t");
                sb.Append("Base Sequences Mapped" + "\t");
                sb.Append("Full Sequences Mapped" + "\t");
                sb.Append("Peak Split Valley RT" + "\t");
                sb.Append("Peak Apex Mass Error (ppm)" + "\t");
                sb.Append("Decoy Peptide" + "\t");
                sb.Append("Random RT");
                return sb.ToString();
            }
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(SpectraFileInfo.FilenameWithoutExtension + "\t");
            sb.Append(string.Join("|", Identifications.Select(p => p.BaseSequence).Distinct()) + '\t');
            sb.Append(string.Join("|", Identifications.Select(p => p.ModifiedSequence).Distinct()) + '\t');

            //The semi-colon here splitting the protein groups requires some explanation
            //During protein parsimony, you can get situations where all peptides are shared between two or more proteins. In other words, there is no unique peptide that could resolve the parsimony.
            //In this case you would see something like P00001 | P00002.

            //That’s the easy part and you already understand that.

            //    Now imagine another scenario where you have some other peptides(that are not in either P00001 or P00002) that give you a second group, like the one above.Let’s call it P00003 | P00004.
            // Everything is still fine her.

            //    Now you have two protein groups each with two proteins. 

            //    Here is where the semi - colon comes in.
            //Imagine you now find a new peptide(totally different from any of the peptides used to create the two original protein groups) that is shared across all four proteins.The original peptides
            //require that two different protein groups exist, but this new peptide could come from either or both.We don’t know. So, the quantification of that peptide must be allowed to be
            //either/ both groups. For this peptide, the protein accession in the output will be P00001 | P00002; P00003 | P00004.

            //    You could see an output that looks like P0000A; P0000B.Here there is only one protein in each protein group(as decided by parsimony).And you have a peptide that is shared.This would
            // not ever be reported as P0000A | P0000B because each protein has a unique peptide that confirms its existence.

            var t = Identifications.SelectMany(p => p.ProteinGroups.Select(v => v.ProteinGroupName)).Distinct().OrderBy(p => p);
            if (t.Any())
            {
                sb.Append(string.Join(";", t) + '\t');
                sb.Append(string.Join(";", Identifications.SelectMany(id => id.ProteinGroups).Select(p => p.Organism).Distinct()) + '\t');
            }
            else
            {
                sb.Append("" + '\t');
                sb.Append("" + '\t');
            }

            sb.Append("" + Identifications.First().MonoisotopicMass + '\t');
            if (!IsMbrPeak)
            {
                sb.Append("" + Identifications.First().Ms2RetentionTimeInMinutes + '\t');
            }
            else
            {
                sb.Append("" + '\t');
            }

            sb.Append("" + Identifications.First().PrecursorChargeState + '\t');
            sb.Append("" + ClassExtensions.ToMz(Identifications.First().MonoisotopicMass, Identifications.First().PrecursorChargeState) + '\t');
            sb.Append("" + Intensity + "\t");

            if (Apex != null)
            {
                sb.Append("" + IsotopicEnvelopes.Min(p => p.IndexedPeak.RetentionTime) + "\t");
                sb.Append("" + Apex.IndexedPeak.RetentionTime + "\t");
                sb.Append("" + IsotopicEnvelopes.Max(p => p.IndexedPeak.RetentionTime) + "\t");

                sb.Append("" + Apex.IndexedPeak.Mz + "\t");
                sb.Append("" + Apex.ChargeState + "\t");
            }
            else
            {
                sb.Append("" + "-" + "\t");
                sb.Append("" + "-" + "\t");
                sb.Append("" + "-" + "\t");

                sb.Append("" + "-" + "\t");
                sb.Append("" + "-" + "\t");
            }

            sb.Append("" + NumChargeStatesObserved + "\t");

            // temporary way to distinguish between MBR, MBR_IsoTrack, IsoTrack_Ambiguous and MSMS peaks
            if (IsMbrPeak && DetectionType == DetectionType.IsoTrack_MBR)
            {
                sb.Append("" + "MBR_IsoTrack" + "\t");
            }

            else if (IsMbrPeak && DetectionType == DetectionType.IsoTrack_Ambiguous)
            {
                sb.Append("" + "IsoTrack_Ambiguous" + "\t");
            }

            else if (IsMbrPeak)
            {
                sb.Append("" + "MBR" + "\t");
            }
            else
            {
                sb.Append("" + "MSMS" + "\t");
            }

            sb.Append("" + (IsMbrPeak ? MbrQValue.ToString() : "") + "\t");
            sb.Append("" + (IsMbrPeak ? MbrPep.ToString() : "") + "\t");

            sb.Append("" + Identifications.Count + "\t");
            sb.Append("" + NumIdentificationsByBaseSeq + "\t");
            sb.Append("" + NumIdentificationsByFullSeq + "\t");
            sb.Append("" + SplitRT + "\t");
            sb.Append("" + MassError);
            sb.Append("\t" + DecoyPeptide);
            sb.Append("\t" + RandomRt);

            return sb.ToString();
        }
    }
}