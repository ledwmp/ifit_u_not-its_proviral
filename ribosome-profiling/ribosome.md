## General pre-processing
Trim multiplexed fastq with bbduk
```bash
bbduk.sh in=<seq_file.fastq> out=<seq_file_trimmed.fastq> qtrim=r literal=CTGTAGGCACCATCAATCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG ktrim=r k=21 mink=10 hdist=1 minlen=31
```
De-multiplex and slice UMIs with [debarcode.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/ribosome-profiling/debarcode.py)
```bash
python debarcody.py <seq_file_trimmed.fastq> <barcode_file.txt>
```
where barcode file takes form:
```
barcode_in_primer_orientation\tindex_number\tbarcode_in_seq_orientation
```
Trim a second time with bbduk (if adapters were stacked on top of eachother)
```bash
bbduk.sh in=<index_number.fastq> out=<index_number_trimmed.fastq> qtrim=r literal=CTGTAGGCACCATCAATCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG ktrim=r k=11 mink=10 hdist=1 minlen=20 maxlen=50
```
Map to human rRNA contig (NT_167214.1) using hisat2
```bash
hisat2 -p <> -x <path to rRNA index> -U <index_number_trimmed.fastq> -S <index_aligned.sam>
```
Filter out rRNA
```bash
samtools fastq -f 4 <index_aligned.sam> > <index_unaligned.fastq>
```
Map to GRCh38/WSN hisat2 index with splice sites derived using HISAT2 and GRCh38/WSN gtf
```bash
hisast2 -p <> --known-splicesite-infile <splice_sites.txt> -x <index_location> -U <index_unaligned.fastq> -S <index.sam>
```
## Prepare for counting
Filter out unmapped reads
```bash
samtools view -F 4 -h <index.sam> > <index_mapped.sam>
```
PCR de-duplicate with UMIs and read-start position using [PCR_dup_remover.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/ribosome-profiling/PCR_dup_remover.py)
```bash
python pcr_dup_remover.py log.log <index_mapped.sam>
```
QC RPF libraries by phasing using plastid
```bash
phase_by_size <GRCh38_start_rois.txt <out_name> --count_files <dedup.bam> --fiveprime_variable --offset <psite_offsets.txt> --codon_buffer 5 --min_length 20 --max_length 40
```
## Counting and DE analysis 
Count RPFs and RNAs using featureCounts
```bash
featureCounts -T <> -s 1 -t exon -g gene_id -a <GRCh38_WSN.gtf> -o <dedup_counts.txt> <dedup.bam>
```
Perform DE using edgeR and 
