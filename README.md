# Data analysis pipeline for IFIT2 eCLIP and ribosome profiling

[CLIP](https://github.com/mehlelab/ifit_u_not-its_proviral/tree/master/clip) and [ribosome profiling](https://github.com/mehlelab/ifit_u_not-its_proviral/tree/master/ribosome-profiling) scripts to re-perform analysis as performed in the paper on [raw sequencing data](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP261790&o=acc_s%3Aa)

## IFIT2 eCLIP
Trim multiplexed fastq with bbduk
```bash
bbduk.sh in=<seq_file.fastq> out=<seq_file_trimmed.fastq> qtrim=r literal=CTGTAGGCACCATCAATCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG ktrim=r k=21 mink=10 hdist=1
```
De-multiplex and slice UMIs with [debarcode.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/debarcode.py)
```bash
python debarcody.py <seq_file_trimmed.fastq> <barcode_file.txt>
```
where barcode file takes form:
```
barcode_in_primer_orientation\tindex_number\tbarcode_in_seq_orientation
```
Trim a second time with bbduk (if adapters were stacked on top of eachother)
```bash
bbduk.sh in=<index_number.fastq> out=<index_number_trimmed.fastq> qtrim=r literal=CTGTAGGCACCATCAATCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG ktrim=r k=11 mink=10 hdist=1 minlen=18
```
Map to GRCh38/WSN hisat2 index
```bash
hisast2 -p <> -x <index_location> -U <index_number_trimmed.fastq> -S <index.sam>
```
Filter out unmapped reads
```bash
samtools view -F 4 -h <index.sam> > <index_mapped.sam>
```
PCR de-duplicate with UMIs and read-start position using [PCR_dup_remover.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/PCR_dup_remover.py)
```bash
python pcr_dup_remover.py log.log <index_mapped.sam>
```
Call peaks using clipper with a modified gene annotation set including WSN
```bash
clipper -b <index_dedup.bam> -s GRCh38_flu -o <index.bed>
```
Concatenate adjacent and overlapping peaks with [merge_adjacent_peaks.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/merge_adjacent_peaks.py)
```bash
python merge_adjacent_peaks.py <index.bed>
```
Convert CLIP bam and size-matched input bam to bed files using bedtools 
```bash
bedtools bamtobed -split -i <who_gives_a_bam.bam> > <my_favorite_bam.bed>
```
Perform input normalization similar to original eCLIP paper but with [input_norm.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/input_norm.py)
```bash
python input_norm.py <CLIP.bed> <SMInput.bed> <index_concat.bed>
```
Collapse overlapping transcripts with [mergetranscriptpeaks.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/mergetranscriptpeaks.py)
```bash
python mergetranscriptpeaks.py <input_norm.bed>
```

## Ribosome profiling
Map to human rRNA contig (NT_167214.1) using hisat2
```bash
hisat2 -p <> -x <path to rRNA index> -U <index_number_trimmed.fastq> -S <index_aligned.sam>
```
Filter out rRNA
```bash
samtools fastq -f 4 <index_aligned.sam> > <index_unaligend.fastq>
```
