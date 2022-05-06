# cellfree

`cellfree`is a analysis suite for cell-free DNA sequencing data.

See docs for more information.

## Feature table format:

| Column Name | Type | Description |  
| ----------- | ---- | ----------- |
| chrom       | str  | Chromosomeor scaffold name  (e.g. chr3, chrY, chr2_random, scaffold10671) |
| chromStart  | int  | Start coordinate on the chromosome or scaffold for the sequence considered (the first base on the chromosome is numbered 0) |
| chromEnd    | int  | End coordinate on the chromosome or scaffold for the sequence considered. **Inclusive (incompatible wit bed format)** |
| name        | str  | name of the read if available. This value is unique per bam file but is not unique across bam files.  |
| score       | float| confidence of the sequence quality and mapping quality   |
| strand      | str  | Which strand the original molecul is sequenced if available. ("+" for forward strand, "-" for reverse strand, "." for unknown.|
| rlen        | int  | Length of the molecule in basepairs (bp).  |
| end5upstream| str  | k-mer sequence compose of (A,T,G,C) of the l bases upstream of 5' end of the molecule. Length of the string is consistent based on user defined value. Default k=4 (4-mer). |
| end5        | str  | k-mer sequence compose of (A,T,G,C) of the k bases of 5' end termini of the molecule. Length of the string is consistent based on user defined value. Default k=4 (4-mer). |
| end3        | str  | k-mer sequence compose of (A,T,G,C) of the k bases of 3' end termini of the molecule. Length of the string is consistent based on user defined value. Default k=4 (4-mer). |
|end3downstream|str  | k-mer sequence compose of (A,T,G,C) of the k bases downstream of 3' end termini of the molecule. Length of the string is consistent based on user defined value. Default k=4 (4-mer). |
| extra_features | unknown  | There can be 0-n columns of extra features upon development of the project. |
| label    | str     | User defined label for each molecules in the supplied bam file.  |
