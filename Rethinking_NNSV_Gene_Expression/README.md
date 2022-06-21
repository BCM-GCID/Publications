This consists of the scripts used in the manuscript : Rethinking NNSV Gene Expression: modeling RSV and VSV transcription with ejective polymerase collisions and biased diffusion.


**Summary**

Infections by non-segmented negative-strand RNA viruses (NNSV) are widely thought to entail gradient gene expression from the well-established existence of a single promoter at the 3’ end of the viral genome and the assumption of constant transcriptional attenuation between genes. But multiple recent studies show viral mRNA levels in respiratory syncytial virus (RSV), a major human pathogen and member of NNSV, infections that are inconsistent with a simple gradient. Here we integrate known and newly predicted phenomena into a biophysically reasonable model of NNSV transcription. Our model succeeds in capturing published observations of RSV and vesicular stomatitis virus (VSV) mRNA levels. We therefore propose a novel understanding of NNSV transcription based on the possibility of ejective polymerase-polymerase collisions and, in the case of RSV, biased polymerase diffusion.      


**Running the Scripts**

**Input**:
- Values for the below parameters are required to run the script:
 1.	Events (# of events in 1 run of model).
 2.	Dscan (=diffusion speed of non-transcribing pols).
 3.	k_transc (=5' speed of transcribing pols).
 4.	Dbias (=multiplicative factor biasing 5' diffusion of non-transcribing pols).
 5.	pol_footprint (= pol footprint size in nt).
 6.	transc_prob (=array containing transcription initiation probabilities for gene start (GS) signals).
 7.	transc_term_prob (=array containing transcription termination probabilities for gene end (GE) signals)
 8.	pol_pos (=array defining position of each modeled pol with number of elements equal to maximum number of pols bound to genome).
 9.	pol_state (=array defining state [0 for non-transcribing, 1 for transcribing] of each modeled pol with number of elements equal to maximum number of pols bound to genome.

**Output**:
- Running the script generates the below outputs.
pol_pos and pol_state array values with each event of the simulation; # of pols bound to genome with each event of the simulation; transcription initiation @ gene x with event #; transcription termination at gene x with event #; pol ejection(s) from position x with event #.

- The following values are output at the end of each simulation: 
  1.	transc_events (=array containing total # of transcription initiation events at each GS signal). total_transc_events (=sum of transcription initiation events at all GS signals). transc_term_events (=array containing total # of transcription termination events at each GE signal).
  2.	Diff (=array of differences between # of transcription initiation events and # of termination events for each gene).
  3.	transc_RT_events (=array containing total # of transcription readthrough events at each GE signal).
  4.	Number of pols ejected.
  5.	a list of values corresponding to the mRNA level for gene x divided by the sum of mRNA levels for all genes (=a gene expression pattern).
  6.	error (=root-mean-square deviation [RMSD] for modeled gene expression pattern vs. experimentally observed gene expression pattern).  
