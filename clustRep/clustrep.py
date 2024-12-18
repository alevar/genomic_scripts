#!/usr/bin/env python3

# ./represent_seqs.py 

import os
import sys
import random
import shutil
import argparse
import subprocess

# compute percent non-ACGTacgt characters in a sequence    
def count_ambiguity(seq):
    return sum([1 for x in seq if not x in "ACGTacgt"])

def run(args):
    if not shutil.which("mmseqs") is not None:
        print("mmseqs not found. Please download and install it prior to running this script.\n Additional information can be found at https:://github.com/soedinglab/MMseqs2")
        exit(1)

    # iterate over the input file and identify suitable sequences
    tmpdir = args.tmp if args.tmp[-1]=="/" else args.tmp+"/"
    tmpdir = os.path.abspath(tmpdir)+"/"
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    # make sure the script has write permissions to the tmp dir
    if not os.access(tmpdir,os.W_OK):
        print("Cannot write to the tmp directory: "+tmpdir)
        exit(1)


    # filter sequences before clustering
    filtered_fasta_fname = tmpdir+"filtered.fasta"

    nts = set(args.nts) if not args.nts is None else None # transform into set to speedup the search        

    with open(filtered_fasta_fname,"w+") as outFP:
        cur_seqid = ""
        cur_seq = ""
        with open(args.i,"r") as inFP:
            for line in inFP:
                if line[0]==">":
                    if len(cur_seq)>0 &\
                       (len(cur_seq)>=args.minlen) and \
                       (len(cur_seq)<=args.maxlen) and \
                       (count_ambiguity(cur_seq)/len(cur_seq)<=args.perc_ambiguity) and \
                       (nts is None or all([x in nts for x in cur_seq])):
                        outFP.write(">"+cur_seqid+"\n")
                        outFP.write(cur_seq+"\n")

                    cur_seqid = line[1:].strip()
                    cur_seq = ""
                else:
                    cur_seq+=line.strip()

        if len(cur_seq)>0 &\
            (len(cur_seq)>=args.minlen) and \
            (len(cur_seq)<=args.maxlen) and \
            (count_ambiguity(cur_seq)/len(cur_seq)<=args.perc_ambiguity) and \
            (nts is None or all([x in nts for x in cur_seq])):
            outFP.write(">"+cur_seqid+"\n")
            outFP.write(cur_seq+"\n")

    # if metafile is provided - we will process sequences by the group column
    # otherwise - we will process all sequences together
    groups_dir = tmpdir+"groups/"
    if not os.path.exists(groups_dir):
        os.mkdir(groups_dir)

    groups = dict()
    if not args.meta is None:
        with open(args.meta,"r") as inFP:
            for line in inFP:
                lcs = line.strip().split("\t")
                assert len(lcs)==2,"meta file must have 2 columns"
                groups[lcs[0]] = lcs[1]

    # iterate over the filtered set and divide sequences into groups
    group_files = dict()
    if len(groups)==0:
        group_files[0] = filtered_fasta_fname
    else:
        cur_seqid = ""
        with open(filtered_fasta_fname,"r") as inFP:
            for line in inFP:
                if line[0]==">":
                    cur_seqid = line[1:].strip()

                if cur_seqid in groups:
                    group = groups[cur_seqid]
                    with open(groups_dir+group+".fasta","a+") as outFP:
                        outFP.write(line)

                    group_files[group] = groups_dir+group+".fasta"
    
    for grp,grp_file in group_files.items():
        print(grp)
        linclust_cmd = ["mmseqs","easy-linclust",
                        grp_file,
                        groups_dir+str(grp)+".clu",
                        groups_dir+str(grp)+".clu.tmp",
                        "--cov-mode","1",
                        "-c",str(args.c),
                        "--min-seq-id",str(args.c),
                        "--kmer-per-seq-scale",str(args.kmer_per_seq_scale),
                        "--threads",str(args.threads)]
        if args.debug:
            print(" ".join(linclust_cmd))
        subprocess.call(linclust_cmd)

    # now we need to combine all representative sequences together
    total_reps = 0
    with open(args.o,"w+") as outFP:
        for grp,file in group_files.items():
            with open(groups_dir+str(grp)+".clu_rep_seq.fasta","r") as inFP:
                for line in inFP:
                    if line[0]==">":
                        total_reps+=1
                    outFP.write(line)

    print("total number of representative sequences: "+str(total_reps))

    # cleanup tmp files unless requested to keep
    if not args.keep_tmp:
        shutil.rmtree(tmpdir)


# from: https://stackoverflow.com/questions/55324449/how-to-specify-a-minimum-or-maximum-float-value-with-argparse
def float_range(mini,maxi):
    """Return function handle of an argument type function for 
       ArgumentParser checking a float range: mini <= arg <= maxi
         mini - minimum acceptable argument
         maxi - maximum acceptable argument"""

    # Define the function with default arguments
    def float_range_checker(arg):
        """New Type function for argparse - a float within predefined range."""

        try:
            f = float(arg)
        except ValueError:    
            raise argparse.ArgumentTypeError("must be a floating point number")
        if f < mini or f > maxi:
            raise argparse.ArgumentTypeError("must be in range [" + str(mini) + " .. " + str(maxi)+"]")
        return f

    # Return function handle to checking function
    return float_range_checker

def main(args):
    global all_stages
    parser = argparse.ArgumentParser(description='''Help Page''')

    ################################
    ######## GENERAL ###############
    ################################
    parser.add_argument('-i',
                        required=True,
                        type=str,
                        help="File containing all the sequnces to be clustered in FASTA format.")
    parser.add_argument('-o',
                        required=True,
                        type=str,
                        help="Output FASTA file containing the representative sequences.")
    parser.add_argument('--meta',
                        required=False,
                        type=str,
                        help="File containing group information for each sequence in the input file. The file is expected to have 2 columns: 1. sequence name as it appears in the input file; 2. group specification. Sequences will be clustered by the 2nd column and subsetting will be performed for each group. Sequences to be preserved will be added in addition to the subsetted tree. If not enabled, any two sequences can be clustered together.")
    parser.add_argument("--threads",
                         required=False,
                         type=int,
                         default=1,
                         help="Number of threads to use [1]. Corresponds to the `--threads` parameter in mmseqs2 easy-linclust.")
    parser.add_argument("--tmp",
                        required=False,
                        type=str,
                        default="./tmp",
                        help="Directory in which temporary data will be stored [./tmp]")
    parser.add_argument("--keep_tmp",
                        required=False,
                        action="store_true",
                        help="If enabled - will keep all temporary data")
    parser.add_argument("--debug",
                        required=False,
                        action="store_true",
                        help="If enabled - will print all commands to stdout")
    
    ################################
    ########## MMSEQS ##############
    ################################
    parser.add_argument("-c",
                        required=False,
                        type=float_range(0.0,1.0),
                        default=0.9,
                        help="List matches above this fraction of aligned (covered) residues (see --cov-mode) [0.800]. Corresponds to the `-c` parameter in mmseqs2 easey-linclust.")
    parser.add_argument("--kmer_per_seq_scale",
                        required=False,
                        type=float_range(0.0,1.0),
                        default=0.3,
                        help="Scale the number of kmers per sequence by this factor [0.5]. Corresponds to the `--kmer-per-seq-scale` parameter in mmseqs2 easy-linclust.")

    ################################
    ########## FILTERS #############
    ################################
    parser.add_argument("--minlen",
                        required=False,
                        type=int,
                        default = 0,
                        help="Minimum genome length to permit in the analysis")
    parser.add_argument("--maxlen",
                        required=False,
                        type=int,
                        default = sys.maxsize,
                        help="Maximum genome length to permit in the analysis")
    parser.add_argument("--nts",
                        required=False,
                        type=str,
                        help="String of nucleotides to consider when filtering sequences. If a sequnece contains a character that is not listed - it will be skipped. Both aupper and lower case variants should be provided. For example: 'ATGCatgc' will only consider sequences that contain only A,T,G,C,a,t,g,c nucleotides. If not provided - all sequences will be considered.")
    parser.add_argument("--perc_ambiguity",
                        required=False,
                        type=float_range(0.0,1.0),
                        default=0.0,
                        help="Maximum fraction of ambiguous nucleotides (N) to permit in the analysis. Sequences with more than this fraction of Ns will be skipped [0.0].")

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])


# tests:
# 1. what happens if two sequences are readlly different (nothing in common)?
