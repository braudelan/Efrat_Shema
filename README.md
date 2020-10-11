### Compare the number of recognition sites in two different sets of sequences.

The sets are given as *.bed* files. Sequences are extended by a constant number of base pairs on each side. This is done by using coordinates of the relevant sequences to extract the extended sequences from a reference genome.

Reference genome is stored in fasta files (one for each chromosome).

The number of recognition sites (including it's complementary sequence) and the number of total base pairs (bp) are counted for each set. The total number of recognition sites is normalized to the total number of bp and given as a percentage.

The repository also contains a short script to download the reference genome from the NCBI nucleotide database.
Mitochondrial chromosome and Ribosomal RNA are not included in the analysis.

*methylation.py* contains the main functions to extract extended sequences and count.
 

