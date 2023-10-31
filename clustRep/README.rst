clustRep: Create condensed representations of datasets of related genomes. Clustering can be performed within groups, such as viral subtypes.
===============================================================================================================================================

Requirements
^^^^^^^^^^^^^

mmseqs must be installed for the script to run. Please follow instrcutions at: https://github.com/soedinglab/MMseqs2

About
^^^^^^
This script performs clustering of input sequences with mmseqs linclust, parses the output and selects representative sequences. Additional metadata can be provided to restrict clustering to specific groups - only seqeuences with the same group ID (such as clade or subtype) will be clustered together. If no groups are provided, all sequences will be evaluated in a single group.

The script enables some simple filtering such as via the length or completeness of the genome sequence.

The script outputs a fasta file with representative sequences for each group.

Usage
^^^^^^^^^^^^^^

usage: clustrep.py [-h] -i I -o O [--meta META] [--threads THREADS] [--tmp TMP] [--keep_tmp] [--debug] [-c C] [--kmer_per_seq_scale KMER_PER_SEQ_SCALE] [--minlen MINLEN] [--maxlen MAXLEN]
                   [--nts NTS] [--perc_ambiguity PERC_AMBIGUITY]

Help Page

options:
  -h, --help            show this help message and exit
  -i I                  File containing all the sequnces to be clustered in FASTA format.
  -o O                  Output FASTA file containing the representative sequences.
  --meta META           File containing group information for each sequence in the input file. The file is expected to have 2 columns: 1. sequence name as it appears in the input file; 2.
                        group specification. Sequences will be clustered by the 2nd column and subsetting will be performed for each group. Sequences to be preserved will be added in
                        addition to the subsetted tree. If not enabled, any two sequences can be clustered together.
  --threads THREADS     Number of threads to use [1]. Corresponds to the `--threads` parameter in mmseqs2 easy-linclust.
  --tmp TMP             Directory in which temporary data will be stored [./tmp]
  --keep_tmp            If enabled - will keep all temporary data
  --debug               If enabled - will print all commands to stdout
  -c C                  List matches above this fraction of aligned (covered) residues (see --cov-mode) [0.800]. Corresponds to the `-c` parameter in mmseqs2 easey-linclust.
  --kmer_per_seq_scale KMER_PER_SEQ_SCALE
                        Scale the number of kmers per sequence by this factor [0.5]. Corresponds to the `--kmer-per-seq-scale` parameter in mmseqs2 easy-linclust.
  --minlen MINLEN       Minimum genome length to permit in the analysis
  --maxlen MAXLEN       Maximum genome length to permit in the analysis
  --nts NTS             String of nucleotides to consider when filtering sequences. If a sequnece contains a character that is not listed - it will be skipped. Both aupper and lower case
                        variants should be provided. For example: 'ATGCatgc' will only consider sequences that contain only A,T,G,C,a,t,g,c nucleotides. If not provided - all sequences
                        will be considered.
  --perc_ambiguity PERC_AMBIGUITY
                        Maximum fraction of ambiguous nucleotides (N) to permit in the analysis. Sequences with more than this fraction of Ns will be skipped [0.0].
