# Each bam file record is a read unit (can contain fastq read1 + fastq read2??).
# Each read is a bam file record read. Multiple reads mapping to the same genomics location can be deduplicated and become a molecule.
# A molecule can contain multiple reads.
# A molecule that mapped to different locations are molecule that contains Structure Variants (SV). They contains reads with same read name but different mapping location.
# Read = Fragment in SingleCellMultiOmics



