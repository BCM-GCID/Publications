This consists of all the scripts, used to assemble the  *Cryptosporidium parvum* genome and compare to older assemblies. All scripts and outputs have also been submitted to GigaScience as part of the manuscript submission.

**Publication Link**

bioRixv : https://www.biorxiv.org/content/10.1101/2021.07.07.451495v2

**Long-read Sequencing :**
  - The DNA from *Cryptosporidium parvum* oocyst was run on the PromethION :
    - Flow Cell : FLO-PRO002
    - Library was prepared using the kit : SQK-LSK110
    - Basecalling was done using the on-board installed Guppy software version : 4.0.11+f1071ce

**Short-read Sequencing :**
  - The DNA from the same extraction used for long-read sequencing was used to generate short-read data. The data was used for polishing of the assembly.
    - Library preparation : Kapa Hyper protocol was used without PCR amplification
    - Library name to track in the HGSC LIMS : ILWGS_GCIDWP_GCIDP3_crypto_202001116_273718_1
    - Instrument : NovaSeq 6000
    - Instrument id : A00733
    - Flow cell id : H5TTKDSX2
    - Lane barcode : H5TTKDSX2-4-IDUDI0001

All various sub-folders in this repo are listed below:

- assembly : Different assemblies performed on long reads with pass and fail ONT reads
- busco  : Buso analysis for Canu assembly
- clustal_alignment  : ClustalW alignments for the last result section
- comparison  : Genome alignments using Mummer
- Cryptosporidium_assemblies  : Holds the differnet versions of our assembly from raw Canu to renamed v2 that was submitted to NCBI
- genome_est  : Estimation of genome size using genomescope +jellyfish based on short reads
- polish  : Error correction of Canu assmebly
- read_len  : Summary of readlenght to produce Fig2



The assembled genome has been deposited at NCBI under
 - BioProject : PRJNA744539 ,  https://www.ncbi.nlm.nih.gov/bioproject/744539
 - BioSample : SAMN20223754 , https://www.ncbi.nlm.nih.gov/biosample/SAMN20223754/
 - The raw data from both the short-read and long-read sequencing are also being uploaded under the same BioProject.
