## General pre-processing
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
## Preparation for peak-calling
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
## Input normalization and replicate analysis
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
Run irreproducible discovery rate analysis on replicates using IDR with ranking performed on -log10(pvalue) from clipper
```bash
idr.sh -s rep1_input_norm.bed rep2_input_norm.bed --plot --rank 5 --input-file-type bed -o <idr.bed>
```
## Peak characterization
Intersect bed with [gene_biotype.bed](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/bed/Homo_sapiens.GRCh38.95_gene_biotype.bed) derived from GRCh38 gtf 
```bash
bedtools intersect -a <idr.bed> -b <gene_biotype.bed> -split -s -wo <idr_biotype.bed>
```
Count and display overlapping meta-features with [count_CLIP_metabiotype.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/count_CLIP_metabiotype_new.py)
```bash
python count_CLIP_metabiotype.py intersect.bed out_file
```
Filter protein-coding genes with bedtools intersect and [protein_coding.bed](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/bed/Homo_sapiens.GRCh38.95_merge_5_3_CDS_collapse.bed)
```bash
bedtools intersect -a <idr.bed> -b <protein_coding.bed> -wa -s -split > <idr_protein_coding.bed>
```
Calculate bam coverage over CLIP peaks with bedtools coverage
```bash
bedtools coverage -a <idr_protein_coding.bed> -b <index_dedup.bam> -d -split -s > <index_dedup.depth>
```
Filter peaks on 4-fold change (CLIP/input) and relative coverage with [plot_coverage_CLIPvexpression.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/plot_coverage_CLIPvexpression.py)
```bash
python plot_coverage_CLIPvexpression.py CLIP1.depth CLIP2.depth SMin1.depth SMin2.depth
```
Make heatmaps of peak coverage with [make_heatmap_peakcoverage.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/make_heatmap_peakcoverage.py)
```bash
python make_heatmap_peakcoverage.py CLIP1.depth CLIP2.depth SMin1.depth SMin2.depth CLIP3.depth CLIP4.depth out_file
```
Intersect and separate peaks with [5'UTR,3'UTR, and CDS](https://github.com/mehlelab/ifit_u_not-its_proviral/tree/master/clip/bed) derived from GRCh38 gtf
```bash
bedtools intersect -a <idr_protein_coding.bed> -b <5UTR.bed> <CDS.bed> <3UTR.bed> -s -wo > <idr_proteing_coding_intersect.bed> 
```
[Split](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/split_CLIP_meta.py) and [plot](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/count_CLIP_meta.py) the intersect
```bash
python split_CLIP_meta.py <idr_protein_coding_intersect.bed>
python count_CLIP_meta.py <idr_protein_coding_intersect.bed> <out_file>
```
Compare FPKM of CLIPed and un-CLIPped genes with [fpkm_distribution.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/fpkm_distribution.py)
```bash
python fpkm_distribution.py clip_list.txt out_file
```
Derive fastas from bed intervals with bedtools
```bash
bedtools getfasta -fi <genome_file.fa> -bed <idr_protein_coding.bed> -s  > <idr_protein_coding.fasta>
```
And calculate GC-content of bed intervals with [fasta_peak_to_gc_kde.py](https://github.com/mehlelab/ifit_u_not-its_proviral/blob/master/clip/fasta_peak_to_gc_kde.py)
```bash
fasta_peak_to_gc_kde.py CLIP_peaks.fasta sub_region.fasta out_file
```
