#!/usr/bin/env python3

# ./blast_ris.py 

import re
import os
import sys
import copy
import random
import shutil
import argparse
import subprocess

def extract_donor_acceptor(fname):
    # parses a GTF file to extract information about all donor and acceptor sites
    # return two dictionaries, one for donors and one for acceptors
    # each dictionary maps position to the gene name

    # run gffread to standardize
    cmd = ["gffread","-T","-F","-o","tmp.gtf",fname]

    donors = {} # position to gene gene_id
    acceptors = {} # position to gene gene_id

    # iterate over lines of gtf to extract the donor acceptor sites
    # donor is defined as the last base of the previous exon
    # acceptor is defined as the first base of the current exon
    with open("tmp.gtf","r") as inFP:
        chain = []
        gid = ""
        for line in inFP:
            if line[0]=="#":
                continue
            lcs = line.strip().split("\t")
            
            if lcs[2] == "transcript":
                for i in range(1,len(chain)):
                    donors.setdefault(lcs[0],{})[chain[i][0]] = (gid,lcs[6])
                    acceptors.setdefault(lcs[0],{})[chain[i-1][0]] = (gid,lcs[6])
                
                # new transcript
                chain = []
                gid = lcs[8].split("gene_id ")[1].split(";")[0].strip("\"")

            if lcs[2] == "exon":
                chain.append((int(lcs[3]),int(lcs[4]))) # position: start,end,strand
    
    if os.path.exists("tmp.gtf"):
        os.remove("tmp.gtf")

    return donors,acceptors

def parse_btop(input_string):
    result = []
    pattern = re.compile(r'(\d+)|([A-Za-z-]{2})')
    matches = re.finditer(pattern, input_string)

    for match in matches:
        if match.group(1):
            result.append(int(match.group(1)))
        else:
            result.append(match.group(2))

    return result

def btop_to_list(btop,length,start_pos):
    binread = [0] * length
    index = start_pos-1

    for b in btop:
        if isinstance(b, int):
            for i in range(index,index+b,1):
                binread[i] += 1
            index += b
        elif isinstance(b, str):
            if b[0]=="-": # insertion
                # if insertion - decrement the score of the next element
                # this is equivalent to saying, by that position the score on the reference would have decreased by this much due to insertion penalty
                binread[index] -= 1
            elif b[1]=="-": # deletion
                binread[index] -= 1
                index += 1
            else: # mismatch
                binread[index] = 0
                index += 1

    return binread

def find_breakpoint(list1, list2):
    max_sum = float('-inf')
    max_indices = []
    for i in range(len(list1)):
        current_sum = sum(list1[:i]) + sum(list2[i:])
        if current_sum > max_sum:
            max_sum = current_sum
            max_indices = [i]
        elif current_sum == max_sum:
            max_indices.append(i)
    return sum(max_indices) // len(max_indices) if max_indices else -1

def get_sites(sites,btop,start_pos_genome,start_pos_read): # start position is on the template this time
    res = []
    
    index_read = start_pos_read-1
    index_genome = start_pos_genome-1

    for b in btop:
        if isinstance(b, int):
            for i in range(0,b,1):
                if index_genome+i in sites:
                    res.append((index_read+i,sites[index_genome+i]))
            index_genome += b
            index_read += b
        elif isinstance(b, str):
            if b[0]=="-": # insertion in read
                index_genome += 1
            elif b[1]=="-": # deletion from read
                index_read += 1
            else: # mismatch - treat as a regular match here
                if index_genome in sites:
                    res.append((index_read,sites[index_genome]))
                index_genome += 1
                index_read += 1

    return res

def match_donor_acceptor(donors,acceptors):
    # finds pairs of donors, acceptors that are adjacent to each other
    return [(x, y) for x in donors for y in acceptors if y[0] - x[0] == 1]

def find_breakpoint_with_sj(list1, list2, sj):
    tmp1 = copy.deepcopy(list1)
    tmp2 = copy.deepcopy(list2)
    tmp1[sj[0][0]]+=5
    tmp2[sj[1][0]]+=5
    breakpoint = find_breakpoint(tmp1,tmp2)
    # only add if breakpoint is actually at the junction site
    if breakpoint == sj[1][0]:
        res = tmp1[:breakpoint]+tmp2[breakpoint:]
        score = sum(res)
        return [score,res,breakpoint,sj]
    return None

def process(name,m1,m2,donors1,acceptors1,donors2,acceptors2,args):
    # take two mappings of the same read and find the breakpoint between them

    qseqid1,qlen1,sseqid1,qstart1,qend1,sstart1,send1,btop1 = m1.split("\t")
    qseqid2,qlen2,sseqid2,qstart2,qend2,sstart2,send2,btop2 = m2.split("\t")

    assert qlen2 == qlen1, "Read lengths do not match."
    assert qseqid1 == qseqid2, "Read names do not match."

    btopl1 = parse_btop(btop1)
    binread1 = btop_to_list(btopl1,qlen1,qstart1)
    btopl2 = parse_btop(btop2)
    binread2 = btop_to_list(btopl2,qlen2,qstart2)

    cur_donors1 = get_sites(donors1,btopl1,sstart1,qstart1)
    cur_acceptors1 = get_sites(acceptors1,btopl1,sstart1,qstart1)
    cur_donors2 = get_sites(donors2,btopl2,sstart2,qstart2)
    cur_acceptors2 = get_sites(acceptors2,btopl2,sstart2,qstart2)

    # given lists of donors and acceptors, find suitable pairs
    # good pair is defined as a combination of g1 donor and g2 acceptor or vice versa 
    # such that they are adjacent on the read (no nuceotides between them)
    d1a2 = match_donor_acceptor(donors1,acceptors2)
    d2a1 = match_donor_acceptor(donors2,acceptors1)

    # now for each combination of donor acceptors
    # we can modify the scores
    # find best brakepoint
    # find combination of donor/acceptor which yields optimal brakepoint
    results = []
    for sj in d1a2+d2a1:
        res = find_breakpoint_with_sj(binread1,binread2,sj)
        if res is not None:
            results.append(res)

    if len(results)==0:
        # add the breakpoint without junction
        breakpoint = find_breakpoint(binread1,binread2)
        res = binread1[:breakpoint]+binread2[breakpoint:]
        score = sum(res)
        return [score,res,breakpoint,"-"]
    else:
        # lastly, select the highest-scoring result
        best_breakpoint = max(results, key=lambda x: x[0])
        return best_breakpoint

def run(args):
    donors1,acceptors1 = extract_donor_acceptor(args.a1)
    donors2,acceptors2 = extract_donor_acceptor(args.a2)

    # what if instead of finding integration sites, we simply do the following:
    # 1. find all reads which map to both genomes
    # 2. search for the junction

    return


def main(args):
    # sample blast command to generate desired input
    # blastn -db macaque_blast_nt_db -query r1.fasta -out r1.host.blastn.6 -outfmt "6 qseqid qlen sseqid qstart qend sstart send btop"
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
    parser.add_argument('-a1',
                        required=True,
                        type=str,
                        help="GTF file containing gene annotations for genome #1.")
    parser.add_argument('-a2',
                        required=True,
                        type=str,
                        help="GTF file containing gene annotations for genome #2.")
    parser.add_argument('-o',
                        required=True,
                        type=str,
                        help="Output file.")
    
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])