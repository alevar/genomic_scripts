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

        self.donors = []
        self.acceptors = []

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

    def get_sites(self,sites): # start position is on the template this time
        res = []
        
        index_read = self.qstart-1
        index_genome = self.sstart-1

        for b in self.btopl:
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
    
    def load_donors(self,donors):
        self.donors = self.get_sites(donors[self.sseqid])
    def load_acceptors(self,acceptors):
        self.donors = self.get_sites(acceptors[self.sseqid])
    
    def read2genome(self,pos):
        # given a btop string and the start position of the read and the genome
        # return the position on the genome corresponding to the position on the read
        index_read = self.qstart-1
        index_genome = self.sstart-1

        for b in self.btopl:
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

        if pos < self.qstart:
            # if insertion is upstream of the read start position
            # just return the start position - the difference between the two
            return self.sstart - (self.qstart - pos)
        else:
            # Otherwise if the position is further down the read than the btop string
            # just add the difference
            return index_genome + (pos - index_read)

class Binread:
    def __init__(self):
        self.read1 = Read()
        self.read2 = Read()

        self.binread = []
        self.breakpoint = None
        self.score = None
        self.sj1 = None
        self.sj2 = None

    def __str__(self):
        
        return "\t".join([self.read1.qseqid,
                          str(self.breakpoint),
                          str(breakpoint[2]),
                          str(breakpoint[0]),
                          str(score),
                          gene1,
                          gene2,
                          ",".join([str(x) for x in res])])+"\n")
               

    def add_read1(self,line):
        self.read1.from_line(line)
        assert self.read2.qlen is None or self.read1.qlen == self.read2.qlen, "Read lengths do not match."
        assert self.read2.qseqid is None or self.read1.qseqid == self.read2.qseqid, "Read names do not match."
    def add_read2(self,line):
        self.read2.from_line(line)
        assert self.read1.qlen is None or self.read1.qlen == self.read2.qlen, "Read lengths do not match."
        assert self.read1.qseqid is None or self.read1.qseqid == self.read2.qseqid, "Read names do not match."

    def _find_breakpoint(list1, list2):
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
        fw_pos,fw_binread,fw_score = self._find_breakpoint(self.read1.binread, self.read2.binread)
        rv_pos,rv_binread,rv_score = self._find_breakpoint(self.read2.binread, self.read1.binread)
        if fw_score > rv_score:
            return [fw_pos,fw_binread,fw_score]
        else:
            return [rv_pos,rv_binread,rv_score]

    def break_at(self,pos):
        # returns the result if read is broken at specific position
        # can be used to score breaks at junction sites by explicitely splitting at a given position
        # and incrementing the returned score by the wight of the junction
        
        # evaluate both sides and return the highest score
        fw_binread = self.read1.binread[:pos]+self.read2.binread[pos:]
        fw_score = sum(fw_binread)/len(fw_binread)
        rv_binread = self.read2.binread[:pos]+self.read1.binread[pos:]
        rv_score = sum(rv_binread)/len(rv_binread)
        if fw_score > rv_score:
            return [pos,fw_binread,fw_score]
        else:
            return [pos,rv_binread,rv_score]
        
    def add_sj(self,site1,site2):
        self.sj1 = site1
        self.sj2 = site2

    def find_best_breakpoint(self):
        results = []

        # first process without junctions searching for the optimal breakpoint from the alignment alone
        bp = self.find_breakpoint()
        bp.append((None,None)) # set gene to None
        results.append(bp)

        # now process with junctions
        if self.sj1 is not None and self.sj2 is not None and self.sj2[0]-self.sj1[0] == 1:
            bp = self.break_at(self.sj1[0])
            bp[2] += 10/len(bp[1]) # add weight of junction
            bp.append((self.sj1,self.sj2))
            results.append(bp)
        if self.sj1 is not None:
            bp = self.break_at(self.sj1[0])
            bp[2] += 3/len(bp[1]) # add weight of junction
            bp.append((self.sj1,None))
            results.append(bp)
        if self.sj2 is not None:
            bp = self.break_at(self.sj2[0])
            bp[2] += 3/len(bp[1]) # add weight of junction
            bp.append((None,self.sj2))
            results.append(bp)
        
        # select the highest-scoring result
        self.breakpoint,self.binread,self.score = max(results, key=lambda x: x[2])

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

def process(name,m1,m2,donors1,acceptors1,donors2,acceptors2,args):
    # take two mappings of the same read and find the breakpoint between them

    binread = Binread()
    binread.add_read1(m1)
    binread.add_read2(m2)

    # add donor/acceptor sites
    binread.read1.load_donors(donors1)
    binread.read1.load_acceptors(acceptors1)
    binread.read2.load_donors(donors1)
    binread.read2.load_acceptors(acceptors1)

    binread.find_best_breakpoint()

    if binread.breakpoint is not None:
        return binread
    else:
        return None
    

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
            if res is not None:
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