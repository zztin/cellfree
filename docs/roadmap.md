### Overview:
Cell free DNA are DNA outside of cells in a bio-liquid. The mapping location and length of cell-free DNA are meaningful,
as well as the genomic variants they contains. This package focus on sequencing data acquired by methods designed for 
cell-free DNA sequencing. There are multiple popular sequencing methods, each may require different preprocessing. 
- ssDNA preparation
- dsDNA preparation
- Nanopore Cyclomics preparation
- Nanopore short read preparation

The aim of this python package is to provide friendly tool for users to analyses cell-free DNA sequencing data with 
ease. We aim to update to include most analyses done for cell-free DNA with the focus on fragmentomics 
features. 

Input file format is [bam files](https://samtools.github.io/hts-specs/SAMv1.pdf). 
Use [samtools](http://www.htslib.org) to sort and create index.

(Later on we can discuss if we want to harness information from raw fastq files and add relevant functionalities.)

Functionalities of this package is listed below:


### Cell free DNA length
```bash
cellfree frag_length --tsv <out_per_read.tsv> --hist <out_per_length.tsv> --png <dist.png> <bamfile_list>
```
### Cell free DNA end motif
```bash
cellfree end_motif --tsv <out_per_read.tsv> --hist <out_per_motif.tsv> --png <dist.png> <bamfile_list>
```
### Cell free DNA overlaying with other (epi)genomic tracks
```bash
cellfree overlay --track_bed <track1.bed,track2.bed> --tsv <out_per_read.tsv> --hist <out_per_track.tsv> --png <dist.png> <bamfile_list>
```
alternative file format: track3.bigwig

### Cell free DNA methylation
This function is only applicable to cfDNA with methylation measured. Either with bi-sulfite sequencing or with 
nanopore sequencing. Determine the ratio of methylation of the cell-free DNA and surrounding region (per MB bin.)
```bash
cellfree call_methyl --tsv <out_per_read.tsv> --hist <out_per_bin.tsv> <bamfile_list>
```

### Cell free DNA assoicated nucleosome type
```bash
cellfree infer_nucleosome --tsv <out_per_read.tsv> --hist <dist.tsv> <bamfile_list>
```

### cell free DNA feature (per molecule)
```bash
cellfree features --feature_list <list_of_subcommand> --tsv <out_per_read.tsv> <bamfile_list>
```

### cell free DNA feature plotting
```bash
cellfree plot --feature_list <list_of_subcommand> --input <tsv or bamfile_list> --png <path_to_figure_folder>
```

### cell free DNA genomic features
This component integrates other available tools since it is universal not restricted to cell-free DNA analysis. 
#### CNV
#### SNV
#### SV


