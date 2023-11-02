#!/usr/bin/env python3

# ./blast_ris.py 

import os
import sys
import random
import shutil
import argparse
import subprocess

def extract_donor_acceptor(fname):
    # parses a GTF file to extract information about all donor and acceptor sites
    # return two dictionaries, one for donors and one for acceptors
    # each dictionary maps position to the gene name

    donors = {}
    acceptors = {}



    return donors,acceptors

def process(name,m1,m2,args):
    # take two mappings of the same read and find the breakpoint between them


    # when approximate site has been found - scan a window around it to see if known donor/acceptor

    # if no donor/acceptor - report site which maximizes combined alignment score (or minimizes combined edit distance)
    return

def run(args):
    # iterate over r1 and r2 files in pairs (assumes the same order)

    name = ""
    r1_seq = ""
    r2_seq = ""
    for line1,line2 in zip(open(args.i1,"r"),open(args.i2,"r")):
        if line1[0]==">" and line2[0]==">":
            if not name=="":
                # process the previous pair
                process(name,r1_seq,r2_seq,args)
            
            assert line1[1:].strip()==line2[1:].strip(), "Sequence IDs do not match: "+line1[1:].strip()+" vs "+line2[1:].strip()
            name = line1[1:].strip()
            r1_seq = ""
            r2_seq = ""
        else:
            r1_seq+=line1.strip()
            r2_seq+=line2.strip()


def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

    ################################
    ######## GENERAL ###############
    ################################
    parser.add_argument('-i1',
                        required=True,
                        type=str,
                        help="File containing mapping of reads to genome #1.")
    parser.add_argument('-i2',
                        required=True,
                        type=str,
                        help="File containing mapping of reads to genome #2.")
    parser.add_argument('-o',
                        required=True,
                        type=str,
                        help="Output file.")
    
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])