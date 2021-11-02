# nucleotide_differences

A quick and dirty little script for pulling total nucleotide differences from a MSA in Fasta format.
The first portion of the script (outputting the one-line fasta file) was shared with me by someone who had cobbled it together from Google searches. I amended this to fit in with the bits that I wrote (line 39 on). 

The script takes an MSA as input and outputs a comparison log in text format with total missingness and gaps plus the known number of substitutions for each pairwise comparison in the fasta file. Gaps and missingness are removed from the tally of total substitutions. It also outputs a pairwise matrix containing the total number of substitutions.

The script should be run like this: \
python3 nucleotide_differences.py input.fasta


Two files will be output: \
compare.log \
substitution_matrix.csv
