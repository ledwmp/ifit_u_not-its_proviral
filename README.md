# Data analysis pipeline for IFIT2 eCLIP and ribosome profiling

[CLIP](https://github.com/mehlelab/ifit_u_not-its_proviral/tree/master/clip) and [ribosome profiling](https://github.com/mehlelab/ifit_u_not-its_proviral/tree/master/ribosome-profiling) scripts to re-perform analysis as performed in the paper on [raw sequencing data](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP261790&o=acc_s%3Aa)

## IFIT2 eCLIP
Follows eCLIP pipeline from Yeo lab pretty closely, with much of the code translated to python. Detailed pipeline [here](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/clip.md)

## Ribosome profiling
Map to human rRNA contig (NT_167214.1) using hisat2
```bash
hisat2 -p <> -x <path to rRNA index> -U <index_number_trimmed.fastq> -S <index_aligned.sam>
```
Filter out rRNA
```bash
samtools fastq -f 4 <index_aligned.sam> > <index_unaligend.fastq>
```
