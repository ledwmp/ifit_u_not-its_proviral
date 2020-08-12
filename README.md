# Data analysis pipeline for IFIT2 eCLIP and ribosome profiling

[CLIP](https://github.com/mehlelab/ifit_u_not-its_proviral/tree/master/clip) and [ribosome profiling](https://github.com/mehlelab/ifit_u_not-its_proviral/tree/master/ribosome-profiling) scripts to re-perform analysis as performed in the paper on [raw sequencing data](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP261790&o=acc_s%3Aa)

## IFIT2 eCLIP
Trim multiplexed sequencing file with bbudk
```bash
bbduk.sh in=<seq_file.fastq> out=<seq_file_trimmed.fastq> qtrim=r literal=CTGTAGGCACCATCAATCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG ktrim=r k=21 mink=10 hdist=1
```
