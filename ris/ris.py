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
from intervaltree import Interval, IntervalTree

class Read:
    def __init__(self):
        self.qseqid = None
        self.qlen = None
        self.sseqid = None
        self.qstart = None
        self.qend = None
        self.sstart = None
        self.send = None
        self.btop = None

        self.btopl = None
        self.binread = None

    def from_line(self,line):
        self.qseqid,self.qlen,self.sseqid,self.qstart,self.qend,self.sstart,self.send,self.btop = line.strip().split("\t")
        self.qlen = int(self.qlen)
        self.qstart = int(self.qstart)
        self.qend = int(self.qend)
        self.sstart = int(self.sstart)
        self.send = int(self.send)

        self.btopl = self.parse_btop(self.btop)
        self.binread = self.btop_to_list(self.btopl,self.qlen,self.qstart)

    def parse_btop(btop_str):
        result = []
        pattern = re.compile(r'(\d+)|([A-Za-z-]{2})')
        matches = re.finditer(pattern, btop_str)

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
    


class Breakpoint:
    def __init__(self):
        self.read1 = Read()
        self.read2 = Read()

        self.binread = []
        self.breakpoint = None
        self.score = None
        self.sj1 = None
        self.sj2 = None

    def add_read1(self,line):
        assert qlen2 == qlen1, "Read lengths do not match."
        assert qseqid1 == qseqid2, "Read names do not match."
        self.read1.from_line(line)
    def add_read2(self,line):
        assert qlen2 == qlen1, "Read lengths do not match."
        assert qseqid1 == qseqid2, "Read names do not match."
        self.read2.from_line(line)

    def _find_breakpoint(self):
        max_sum = float('-inf')
        max_indices = []
        for i in range(len(list1)):
            current_sum = sum(list1[:i]) + sum(list2[i:])
            if current_sum > max_sum:
                max_sum = current_sum
                max_indices = [i]
            elif current_sum == max_sum:
                max_indices.append(i)

        breakpoint = sum(max_indices) // len(max_indices) if max_indices else -1
        binread = list1[:breakpoint]+list2[breakpoint:]
        return [breakpoint,binread,round(max_sum/len(binread),2)]

    def find_breakpoint(self):
        fw_pos,fw_binread,fw_score = _find_breakpoint(list1, list2)
        rv_pos,rv_binread,rv_score = _find_breakpoint(list2, list1)
        if fw_score > rv_score:
            return [fw_pos,fw_binread,fw_score]
        else:
            return [rv_pos,rv_binread,rv_score]
        
    def find_breakpoint_with_sj(self, site1, site2):
        tmp1 = copy.deepcopy(list1)
        tmp2 = copy.deepcopy(list2)
        if site1 and site2: # prioritize cases with matching known junctions
            tmp1[site1[0]]+=5
            tmp2[site2[0]]+=5
        else: # however if only odnor or acceptor is known - it is possible the junction is cryptic. give it a slight bump
            if site1:
                tmp1[site1[0]]+=3
            if site2:
                tmp2[site2[0]]+=3
        breakpoint,binread,score = find_breakpoint(tmp1,tmp2)
        # only add if breakpoint is actually at the junction site
        if breakpoint == site1[0]:
            return [score,binread,breakpoint,(site1,site2)]
        return None

    def read2genome(btop,start_pos_read,start_pos_genome,pos):
        # given a btop string and the start position of the read and the genome
        # return the position on the genome corresponding to the position on the read
        index_read = start_pos_read-1
        index_genome = start_pos_genome-1

        for b in btop:
            if index_read == pos:
                return index_genome
            if isinstance(b, int):
                for i in range(0,b,1):
                    index_genome += 1
                    index_read += 1
            elif isinstance(b, str):
                if b[0]=="-": # insertion in read
                    index_genome += 1
                elif b[1]=="-": # deletion from read
                    index_read += 1
                else: # mismatch - treat as a regular match here
                    index_genome += 1
                    index_read += 1

        if pos < start_pos_read:
            # if insertion is upstream of the read start position
            # just return the start position - the difference between the two
            return start_pos_genome - (start_pos_read - pos)
        else:
            # Otherwise if the position is further down the read than the btop string
            # just add the difference
            return index_genome + (pos - index_read)

def extract_genes(fname):
    # parses a GTF file to extract information about all genes
    # returns a dictionary mapping position to the gene name

    # run gffread to standardize
    cmd = ["gffread","-T","-F","-o","tmp.gtf",fname]
    subprocess.call(cmd)

    genes = {} # gene name to start-end coordinate
    gene_tree = IntervalTree()

    # iterate over lines of gtf to extract the genes
    with open("tmp.gtf","r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            lcs = line.strip().split("\t")
            
            if lcs[2] == "transcript":
                gid = lcs[8].split("gene_id ")[1].split(";")[0].strip("\"")
                genes.setdefault(gid,[int(lcs[3]),int(lcs[4])])
                # update start to minimum of current and recorded
                genes[gid][0] = min(genes[gid][0],int(lcs[3]))
                # update end to maximum of current and recorded
                genes[gid][1] = max(genes[gid][1],int(lcs[4]))

    # move to interval tree
    for gid,coords in genes.items():
        gene_tree.addi(coords[0], coords[1], gid)

    if os.path.exists("tmp.gtf"):
        os.remove("tmp.gtf")

    return genes

def extract_donor_acceptor(fname):
    # parses a GTF file to extract information about all donor and acceptor sites
    # return two dictionaries, one for donors and one for acceptors
    # each dictionary maps position to the gene name

    # run gffread to standardize
    cmd = ["gffread","-T","-F","-o","tmp.gtf",fname]
    subprocess.call(cmd)

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
                    donors.setdefault(lcs[0],{})[chain[i-1][1]] = (gid,lcs[6])
                    acceptors.setdefault(lcs[0],{})[chain[i][0]] = (gid,lcs[6])
                
                # new transcript
                chain = []
                gid = lcs[8].split("gene_id ")[1].split(";")[0].strip("\"")

            if lcs[2] == "exon":
                chain.append((int(lcs[3]),int(lcs[4]))) # position: start,end,strand
    
    if os.path.exists("tmp.gtf"):
        os.remove("tmp.gtf")

    return donors,acceptors

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
    res = []
    for x in donors:
        for y in acceptors:
            if y[0] - x[0] == 1:
                res.append((x,y))
            else:
                res.append((x,None))
                res.append((None,y))
    return res

def process(name,m1,m2,donors1,acceptors1,donors2,acceptors2,args):
    # take two mappings of the same read and find the breakpoint between them

    binread = Binread()
    binread.add_read1(m1)
    binread.add_read2(m2)

    cur_donors1 = get_sites(donors1[sseqid1],btopl1,sstart1,qstart1)
    cur_acceptors1 = get_sites(acceptors1[sseqid1],btopl1,sstart1,qstart1)
    cur_donors2 = get_sites(donors2[sseqid2],btopl2,sstart2,qstart2)
    cur_acceptors2 = get_sites(acceptors2[sseqid2],btopl2,sstart2,qstart2)

    # given lists of donors and acceptors, find suitable pairs
    # good pair is defined as a combination of g1 donor and g2 acceptor or vice versa 
    # such that they are adjacent on the read (no nuceotides between them)
    d1a2 = match_donor_acceptor(cur_donors1,cur_acceptors2)
    d2a1 = match_donor_acceptor(cur_donors2,cur_acceptors1)

    # now for each combination of donor acceptors
    # we can modify the scores
    # find best brakepoint
    # find combination of donor/acceptor which yields optimal brakepoint
    results = []
    for sj in d1a2:
        res = find_breakpoint_with_sj(binread1,binread2,sj[0],sj[1])
        if res is not None:
            results.append(res)
    for sj in d2a1:
        res = find_breakpoint_with_sj(binread1,binread2,sj[1],sj[0])
        if res is not None:
            results.append(res)

    best_breakpoint = []
    if len(results)==0:
        # add the breakpoint without junction
        best_breakpoint = find_breakpoint(binread1,binread2)
        best_breakpoint.append((None,None)) # set gene to None
    else:
        # lastly, select the highest-scoring result
        best_breakpoint = max(results, key=lambda x: x[0])
    
    # get breakpoint on the genomes
    g1breakpoint = read2genome(btopl1,qstart1,sstart1,best_breakpoint[2])
    g2breakpoint = read2genome(btopl2,qstart2,sstart2,best_breakpoint[2])
    best_breakpoint[2] = (best_breakpoint[2],g1breakpoint,g2breakpoint)
    return best_breakpoint

def next_read_group(fname1,fname2):
    # iterate over lines of two sorted files
    # group lines by read name (1st column)
    # yield two lists of lines, one from each file, that have the same read name
    with open(fname1, 'r') as inFP1, open(fname2, 'r') as inFP2:
        iter1, iter2 = iter(inFP1), iter(inFP2)
        line1, line2 = next(iter1, None), next(iter2, None)
        while line1 == "\n":
            line1 = next(iter1, None)
        while line2 == "\n":
            line2 = next(iter2, None)

        while line1 is not None or line2 is not None:
            current_read_name = None
            lines1, lines2 = [], []

            if line1 is not None and (line2 is None or line1.split("\t")[0] <= line2.split("\t")[0]):
                current_read_name = line1.split("\t")[0]
                while line1 is not None and line1.split("\t")[0] == current_read_name:
                    lines1.append(line1.strip())
                    line1 = next(iter1, None)

            if line2 is not None and (line1 is None or line2.split("\t")[0] <= line1.split("\t")[0]):
                current_read_name = line2.split("\t")[0] if current_read_name is None else current_read_name
                while line2 is not None and line2.split("\t")[0] == current_read_name:
                    lines2.append(line2.strip())
                    line2 = next(iter2, None)

            yield current_read_name, lines1, lines2

def ris(args,outFP):
    genes1 = extract_genes(args.a1)
    genes2 = extract_genes(args.a2)
    donors1,acceptors1 = extract_donor_acceptor(args.a1)
    donors2,acceptors2 = extract_donor_acceptor(args.a2)

    # sort inputs by read name
    cmd = ["sort","-k1,1",args.i1,"-o","tmp1"]
    subprocess.call(cmd)
    cmd = ["sort","-k1,1",args.i2,"-o","tmp2"]
    subprocess.call(cmd)

    # iterate over read groups
    for read_name, lines1, lines2 in next_read_group("tmp1","tmp2"):
        # create all unique combinations of two lists
        breakpoints = []
        for line1,line2 in [(x,y) for x in lines1 for y in lines2]:
            res = process(read_name,line1,line2,donors1,acceptors1,donors2,acceptors2,args)
            breakpoints.append(res)

        if len(breakpoints)==0:
            continue

        score,res,breakpoint,gene = max(breakpoints, key=lambda x: x[0])
        gene1 = gene[0][1][0] if gene[0] is not None else "-"
        gene2 = gene[1][1][0] if gene[1] is not None else "-"
        outFP.write("\t".join([read_name,
                                str(breakpoint[1]),
                                str(breakpoint[2]),
                                str(breakpoint[0]),
                                str(score),
                                gene1,
                                gene2,
                                ",".join([str(x) for x in res])])+"\n")

    return

def run(args):
    outFP = open(args.o,"w+")
    try:
        outFP.write("read_name\tgenome1_breakpoint\tgenome2_breakpoint\tread_breakpoint\tscore\tgene\tbinread\n")
        ris(args,outFP)
    except Exception as e:
        print(e)
        outFP.close()
        sys.exit(1)
    outFP.close()

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