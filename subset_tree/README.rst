subset_tree: Select a subset of the phylogenetic tree maximizing diversity representation
===============================================================================================================================================

Requirements
^^^^^^^^^^^^^
pip install -r requirements.txt

About
^^^^^^
This script does not build phylogenetic trees. It is up to the user to build a phylogenetic tree and convert it into the NWK format for the script to run properly.

This script iterates over pairs of leaves in a tree, joining the two closest leaves on each iteration. The script terminates when it reaches N sequences as requested by the user.

The script outputs a subset of the tree in nwk format, a csv file describing which sequences are groupped together. An illustration of the tree can be requested to be dreawn and saved to disk as well. If FASTA is provided with sequences corresponding to the input tree, a subset of the FASTA will also be created using representative sequences

Usage
^^^^^^^^
usage: subset_tree.py [-h] --tree TREE [--fasta FASTA] --meta META --output OUTPUT [-N NUM_SEQS] [--draw] [--seed SEED] [--debug]

Help Page

options:
  -h, --help            show this help message and exit
  --tree TREE           File containing a tree to be subsetted in Newick format
  --fasta FASTA         File containing a fasta file with sequences corresponding to the tree. If provided, an output fasta file will be generated with sequences corresponding to the once
                        found in the subsetted tree.
  --meta META           CSV file containing metadata corresponding to the sequences in the tree. Metadata is expected to have 3 columns: 1. sequence name as it appears in the tree; 2.
                        group specification; 3. Whether the sequence has to be preserved in the output (1/0). The data will be groupped by the 2nd column and subsetting will be performed
                        for each group. Sequences to be preserved will be added in addition to the subsetted tree
  --output OUTPUT       Output directory in which all output and temporary data will be stored
  -N NUM_SEQS, --num_seqs NUM_SEQS
                        Number of genomes to preserve in each group
  --draw                If enabled - will draw the resulting phylogenetic tree
  --seed SEED           If set - will be used a seed for all pseudo-random choices
  --debug               If enabled - will print debug information
