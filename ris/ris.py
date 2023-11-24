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
        self.weights = []

    def to_interval(self):
        return Interval(min(self.sstart,self.send),max(self.sstart,self.send),self)

    def from_line(self,line):
        self.qseqid,self.qlen,self.sseqid,self.qstart,self.qend,self.sstart,self.send,self.btop = line.strip().split("\t")
        self.qlen = int(self.qlen)
        self.qstart = int(self.qstart)
        self.qend = int(self.qend)
        self.sstart = int(self.sstart)
        self.send = int(self.send)

        self.parse_btop()
        self.btop_to_list()

    def is_reversed(self):
        return self.sstart > self.send
    
    def reverse(self):
        self.sstart,self.send = self.send,self.sstart
        self.qstart,self.qend = self.qlen-self.qend+1,self.qlen-self.qstart+1
        self.btopl = self.btopl[::-1]
        self.binread = self.binread[::-1]
        self.donors = [(self.qlen-x[0],x[1]) for x in self.donors]
        self.acceptors = [(self.qlen-x[0],x[1]) for x in self.acceptors]
        self.weights = [(self.qlen-x[0],x[1]) for x in self.weights]

    def parse_btop(self):
        self.btopl = []
        pattern = re.compile(r'(\d+)|([A-Za-z-]{2})')
        matches = re.finditer(pattern, self.btop)

        for match in matches:
            if match.group(1):
                self.btopl.append(int(match.group(1)))
            else:
                self.btopl.append(match.group(2))

    def btop_to_list(self):
        self.binread = [0] * self.qlen
        index = self.qstart-1

        for b in self.btopl:
            if isinstance(b, int):
                for i in range(index,index+b,1):
                    self.binread[i] += 1
                index += b
            elif isinstance(b, str):
                if b[0]=="-": # insertion
                    # if insertion - decrement the score of the next element
                    # this is equivalent to saying, by that position the score on the reference would have decreased by this much due to insertion penalty
                    self.binread[index] -= 1
                elif b[1]=="-": # deletion
                    self.binread[index] -= 1
                    index += 1
                else: # mismatch
                    self.binread[index] = 0
                    index += 1

    def get_sites(self,sites): # start position is on the template this time
        res = []

        # handle the case when the insertion is upstream of the read start position
        for i in range(self.qstart-2,-1,-1):
            genome_pos = self.sstart - (self.qstart - i)
            if genome_pos in sites:
                res.append((i,sites[genome_pos]))
        
        index_read = self.qstart-1
        index_genome = self.sstart

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

        # handle the case when the insertion is downstream of the read end position
        for i in range(index_read,self.qlen,1):
            genome_pos = index_genome + (i - index_read)
            if genome_pos in sites:
                res.append((i,sites[genome_pos]))

        return res
    
    def load_donors(self,donors):
        if self.sseqid in donors:
            self.donors = self.get_sites(donors[self.sseqid])
    def load_acceptors(self,acceptors):
        if self.sseqid in acceptors:
            self.acceptors = self.get_sites(acceptors[self.sseqid])

    def load_weights(self,weights):
        # iterates over the positions on the read and extracts weights that belong to the read

        # handle the case when the insertion is upstream of the read start position
        for i in range(self.qstart-1,-1,-1):
            genome_pos = self.sstart - (self.qstart - i)
            weight = weights.get((self.sseqid,genome_pos),None)
            if weight is not None:
                self.weights.append((i,weight))

        index_read = self.qstart-1
        index_genome = self.sstart

        for b in self.btopl:
            if isinstance(b, int):
                for i in range(0,b,1):
                    weight = weights.get((self.sseqid,index_genome+i),None)
                    if weight is not None:
                        self.weights.append((index_read+i,weight))
                index_genome += b
                index_read += b

            elif isinstance(b, str):
                if b[0]=="-": # insertion in read
                    weight = weights.get((self.sseqid,index_genome),None)
                    if weight is not None:
                        self.weights.append((index_read,weight))
                    index_genome += 1

                elif b[1]=="-": # deletion from read
                    index_read += 1
                else: # mismatch - treat as a regular match here
                    weight = weights.get((self.sseqid,index_genome),None)
                    if weight is not None:
                        self.weights.append((index_read,weight))
                    index_genome += 1
                    index_read += 1

        # handle the case when the insertion is downstream of the read end position
        for i in range(index_read,self.qlen,1):
            genome_pos = index_genome + (i - index_read)
            weight = weights.get((self.sseqid,genome_pos),None)
            if weight is not None:
                self.weights.append((i,weight))

        return
    
    def read2genome(self,pos):
        # given a btop string and the start position of the read and the genome
        # return the position on the genome corresponding to the position on the read
        index_read = self.qstart-1
        index_genome = self.sstart

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

        self.sj1 = None
        self.sj2 = None

        self.binread = []
        self.breakpoint = None
        self.score = None
        self.genes = None

        self.orientation = None # value 1 means that in the RIS read1 is on the left, read2 is on the right. Value 2 means the opposite.

    def __str__(self):
        genome1_read_breakpoint = self.breakpoint if self.orientation==1 else self.breakpoint+1
        genome2_read_breakpoint = self.breakpoint+1 if self.orientation==1 else self.breakpoint

        genome1_pos = self.read1.read2genome(genome1_read_breakpoint)
        genome2_pos = self.read2.read2genome(genome2_read_breakpoint)
        gene1 = self.sj1[1][0] if self.sj1 is not None else "-"
        gene2 = self.sj2[1][0] if self.sj2 is not None else "-"
        
        return "\t".join([self.read1.qseqid,
                          str(genome1_read_breakpoint),
                          str(genome2_read_breakpoint),
                          self.read1.sseqid,
                          self.read2.sseqid,
                          str(genome1_pos),
                          str(genome2_pos),
                          str(self.score),
                          gene1,
                          gene2,
                          ";".join(list(self.genes[0])),
                          ";".join(list(self.genes[1])),
                          "".join([str(x) for x in self.binread])])
               

    def add_read1(self,line):
        self.read1.from_line(line)
        assert self.read2.qlen is None or self.read1.qlen == self.read2.qlen, "Read lengths do not match."
        assert self.read2.qseqid is None or self.read1.qseqid == self.read2.qseqid, "Read names do not match."
    def add_read2(self,line):
        self.read2.from_line(line)
        assert self.read1.qlen is None or self.read1.qlen == self.read2.qlen, "Read lengths do not match."
        assert self.read1.qseqid is None or self.read1.qseqid == self.read2.qseqid, "Read names do not match."

    def is_reversed(self):
        assert self.read1.is_reversed() == self.read2.is_reversed(), "Segments of the read are not oriented the same way."
        return self.read1.is_reversed() and self.read2.is_reversed()
    
    def reverse(self):
        self.read1.reverse()
        self.read2.reverse()

    @staticmethod
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
            self.orientation = 1
            self.breakpoint,self.binread,self.score,self.genes = fw_pos,fw_binread,fw_score,(None,None)
        else:
            self.orientation = 2
            self.breakpoint,self.binread,self.score,self.genes = rv_pos,rv_binread,rv_score,(None,None)

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
            return [pos,fw_binread,fw_score,1]
        else:
            return [pos,rv_binread,rv_score,2]
        
    def add_sj(self,site1,site2):
        self.sj1 = site1
        self.sj2 = site2

    def find_spliced_breakpoint(self):
        results = []

        # now process with junctions
        if self.sj1 is not None and self.sj2 is not None and self.sj2[0]-self.sj1[0] == 1:
            bp = self.break_at(self.sj1[0])
            bp[2] += 10/len(bp[1]) # add weight of junction
            bp.append((self.sj1,self.sj2))
            results.append(bp)
        elif self.sj1 is not None and self.sj2 is not None and self.sj1[0]-self.sj2[0] == 1:
            bp = self.break_at(self.sj2[0])
            bp[2] += 10/len(bp[1]) # add weight of junction
            bp.append((self.sj2,self.sj1))
            results.append(bp)
        elif self.sj1 is not None:
            bp = self.break_at(self.sj1[0])
            bp[2] += 3/len(bp[1]) # add weight of junction
            bp.append((self.sj1,None))
            results.append(bp)
        elif self.sj2 is not None:
            bp = self.break_at(self.sj2[0])
            bp[2] += 3/len(bp[1]) # add weight of junction
            bp.append((None,self.sj2))
            results.append(bp)
        else:
            print("unknown clause")
            sys.exit(1)
        
        if len(results)==0:
            self.breakpoint,self.binread,self.score,self.genes,self.orientation = None,None,None,None,None
        else:
            self.breakpoint,self.binread,self.score,self.genes,self.orientation = max(results, key=lambda x: x[2])

def extract_genes(fname,exonic=False):
    # parses a GTF file to extract information about all genes
    # returns a dictionary mapping position to the gene name

    # run gffread to standardize
    cmd = ["gffread","-T","-F","-o","tmp.gtf",fname]
    subprocess.call(cmd)

    genes = {} # gene name to start-end coordinate
    gene_trees = {}

    # iterate over lines of gtf to extract the genes
    with open("tmp.gtf","r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            lcs = line.strip().split("\t")
            if exonic:
                if lcs[2] == "exon":
                    gid = lcs[8].split("gene_id ")[1].split(";")[0].strip("\"")

                    gene_trees.setdefault((lcs[0],lcs[6]),IntervalTree()).addi(int(lcs[3]),int(lcs[4])+1, gid) # +1 because intervaltree is exclusive on the right
            else:
                if lcs[2] == "transcript":
                    gid = lcs[8].split("gene_id ")[1].split(";")[0].strip("\"")
                    genes.setdefault(gid,[lcs[0],lcs[6],int(lcs[3]),int(lcs[4])+1])
                    # update start to minimum of current and recorded
                    genes[gid][2] = min(genes[gid][2],int(lcs[3]))
                    # update end to maximum of current and recorded
                    genes[gid][3] = max(genes[gid][3],int(lcs[4])+1)

    if not exonic:
        # move to interval tree
        for gid,coords in genes.items():
            gene_trees.setdefault((coords[0],coords[1]),IntervalTree()).addi(coords[2], coords[3], gid)

    if os.path.exists("tmp.gtf"):
        os.remove("tmp.gtf")

    return gene_trees

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
                    if chain[i][0] - chain[i-1][1] < 1:
                        print("Intron too short. Discarding. ",chain[i-1],chain[i])
                        continue

                    # consider strand
                    if lcs[6] == "+":
                        donors.setdefault(lcs[0],{})[chain[i-1][1]] = (gid,lcs[6])
                        acceptors.setdefault(lcs[0],{})[chain[i][0]] = (gid,lcs[6])
                    else:
                        donors.setdefault(lcs[0],{})[chain[i][0]] = (gid,lcs[6])
                        acceptors.setdefault(lcs[0],{})[chain[i-1][1]] = (gid,lcs[6])
                        
                # new transcript
                chain = []
                gid = lcs[8].split("gene_id ")[1].split(";")[0].strip("\"")

            if lcs[2] == "exon":
                chain.append((int(lcs[3]),int(lcs[4]))) # position: start,end,strand
    
    if os.path.exists("tmp.gtf"):
        os.remove("tmp.gtf")

    return donors,acceptors

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

def process(m1,m2,donors1,acceptors1,donors2,acceptors2,args,pass1_bps=None):
    # take two mappings of the same read and find the breakpoint between them

    binread = Binread()
    binread.add_read1(m1)
    binread.add_read2(m2)

    # check if read needs to be reversed
    if binread.is_reversed():
        binread.reverse()

    # add weights from the first pass if available
    weight_pairs = [None]
    if not pass1_bps is None:
        binread.read1.load_weights(pass1_bps)
        binread.read2.load_weights(pass1_bps)
        weight_pairs = [(x,y) for x in binread.read1.weights for y in binread.read2.weights]
        
    # add donor/acceptor sites
    binread.read1.load_donors(donors1)
    binread.read1.load_acceptors(acceptors1)
    binread.read2.load_donors(donors2)
    binread.read2.load_acceptors(acceptors2)

    # match junctions:
    d1a2 = match_donor_acceptor(binread.read1.donors,binread.read2.acceptors)
    d2a1 = match_donor_acceptor(binread.read2.donors,binread.read1.acceptors)

    # now for each combination of donor acceptors
    # we can modify the scores
    # find best brakepoint
    # find combination of donor/acceptor which yields optimal brakepoint
    results = []

    for weight_pair in weight_pairs:
        base_binread = copy.deepcopy(binread) # create a master copy of the binread
        if not weight_pair is None:
            base_binread.read1.binread[weight_pair[0][0]] = weight_pair[0][1]
            base_binread.read1.binread[weight_pair[1][0]] = weight_pair[1][1]

        # find unspliced breakpoint first
        tmp = copy.deepcopy(base_binread)
        tmp.find_breakpoint()
        results.append(tmp)

        for sj in d1a2:
            tmp = copy.deepcopy(base_binread)
            tmp.add_sj(sj[0],sj[1])
            tmp.find_spliced_breakpoint()
            if not tmp.breakpoint is None:
                results.append(tmp)

        for sj in d2a1:
            tmp = copy.deepcopy(base_binread)
            tmp.add_sj(sj[1],sj[0])
            tmp.find_spliced_breakpoint()
            if not tmp.breakpoint is None:
                results.append(tmp)

    if len(results)==0:
        return None
    else:
        return max(results, key=lambda x: x.score)
    

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

def ris(args,in1_fname,in2_fname,out_fname,gene_trees1,gene_trees2,donors1,acceptors1,donors2,acceptors2,pass1_bps=None):
    outFP = open(out_fname,"w+")
    outFP.write("read_name\t" +
                "genome1_read_breakpoint\t" +
                "genome2_read_breakpoint\t" +
                "genome1_seqid\t" +
                "genome2_seqid\t" +
                "genome1_breakpoint\t" +
                "genome2_breakpoint\t" +
                "score\t" +
                "junction1\t" +
                "junction2\t" +
                "gene1\t" +
                "gene2\t" +
                "binread\n")

    # iterate over read groups
    for read_name, lines1, lines2 in next_read_group(in1_fname,in2_fname):
        if len(lines1)==0 or len(lines2)==0:
            continue
        
        # create all unique combinations of two lists
        breakpoints = []
        for line1,line2 in [(x,y) for x in lines1 for y in lines2]:
            res = process(line1,line2,donors1,acceptors1,donors2,acceptors2,args,pass1_bps)
            if res is not None:
                it1 = res.read1.to_interval()
                it2 = res.read2.to_interval()
                found_genes1 = set()
                found_genes2 = set()
                for seqid,genes in gene_trees1.items():
                    inter = genes.overlap(it1)
                    if len(inter)>0:
                        found_genes1.update([x[2] for x in inter])
                for seqid,genes in gene_trees2.items():
                    inter = genes.overlap(it2)
                    if len(inter)>0:
                        found_genes2.update([x[2] for x in inter])
                
                if len(found_genes1)==0:
                    found_genes1.add("-")
                if len(found_genes2)==0:
                    found_genes2.add("-")
                
                res.genes = (found_genes1,found_genes2)
                breakpoints.append(res)

        if len(breakpoints)==0:
            continue

        breakpoint = max(breakpoints, key=lambda x: x.score)
        outFP.write(str(breakpoint)+"\n")

    outFP.close()
    return

def ris_pass2(args,in_fname,in1_fname,in2_fname,out_fname,gene_trees1,gene_trees2,donors1,acceptors1,donors2,acceptors2):
    # read in the output of the first pass and collect all breakpoints
    # extra weight is assigned if the breakpoint matches donor/acceptor or both
    pass1_bps = {} # (seqid,bp):weighted count
    with open(in_fname,"r") as inFP:
        next(inFP) # skip header
        for line in inFP:
            lcs = line.strip().split("\t")
            seqid1 = lcs[3]
            seqid2 = lcs[4]
            bp1 = int(lcs[5])
            bp2 = int(lcs[6])
            sj1 = lcs[8]!="-"
            sj2 = lcs[9]!="-"
            weight1 = 1 + (1 if sj1 else 0) + (1 if sj2 else 0)
            weight2 = 1 + (1 if sj2 else 0) + (1 if sj1 else 0)
            pass1_bps.setdefault((seqid1,bp1),0)
            pass1_bps.setdefault((seqid2,bp2),0)
            pass1_bps[(seqid1,bp1)] += weight1
            pass1_bps[(seqid2,bp2)] += weight2

    # rerun ris detection, this time adding weights according to the previous round
    ris(args,in1_fname,in2_fname,out_fname,gene_trees1,gene_trees2,donors1,acceptors1,donors2,acceptors2,pass1_bps)

def group_breakpoints(args,in_fname,outFP):
    # read in the input, count the number of times each breakpoint is observed, output the result

    groups = {} # (seqid1,seqid2,bp1,bp2):(count,sj1,sj2,gene1,gene2)


    with open(in_fname,"r") as inFP:
        next(inFP)
        for line in inFP:
            lcs = line.strip().split("\t")
            seqid1 = lcs[3]
            seqid2 = lcs[4]
            bp1 = int(lcs[5])
            bp2 = int(lcs[6])
            sj1 = lcs[8]
            sj2 = lcs[9]
            gene1 = lcs[10]
            gene2 = lcs[11]

            groups.setdefault((seqid1,seqid2,bp1,bp2),[0,sj1,sj2,gene1,gene2])
            group_sj1 = groups[(seqid1,seqid2,bp1,bp2)][1]
            group_sj2 = groups[(seqid1,seqid2,bp1,bp2)][2]
            group_gene1 = groups[(seqid1,seqid2,bp1,bp2)][3]
            group_gene2 = groups[(seqid1,seqid2,bp1,bp2)][4]
            assert group_sj1 == sj1, "sj1 does not match"
            assert group_sj2 == sj2, "sj2 does not match"
            assert group_gene1 == gene1, "gene1 does not match"
            assert group_gene2 == gene2, "gene2 does not match"

            groups[(seqid1,seqid2,bp1,bp2)][0] += 1

    for k,v in groups.items():
        outFP.write(k[0]+"\t"+
                    k[1]+"\t"+
                    str(k[2])+"\t"+
                    str(k[3])+"\t"+
                    str(v[0])+"\t"+
                    v[1]+"\t"+
                    v[2]+"\t"+
                    v[3]+"\t"+
                    v[4]+"\n")


    return

def run(args):
    gene_trees1 = extract_genes(args.a1,args.exonic)
    gene_trees2 = extract_genes(args.a2,args.exonic)
    donors1,acceptors1 = extract_donor_acceptor(args.a1)
    donors2,acceptors2 = extract_donor_acceptor(args.a2)

    # sort inputs by read name
    cmd = ["sort","-k1,1",args.i1,"-o","tmp1"]
    subprocess.call(cmd)
    cmd = ["sort","-k1,1",args.i2,"-o","tmp2"]
    subprocess.call(cmd)

    group_in_fname = None
    try:
        
        ris(args,"tmp1","tmp2",args.o,gene_trees1,gene_trees2,donors1,acceptors1,donors2,acceptors2)

        group_in_fname = args.o

        if args.two_pass:
            
            ris_pass2(args,args.o,"tmp1","tmp2",args.o+".corrected",gene_trees1,gene_trees2,donors1,acceptors1,donors2,acceptors2)

            group_in_fname = args.o+".corrected"

        if args.group:
            group_outFP = open(args.o+".grouped","w+")
            group_outFP.write("genome1_seqid\t" +
                                "genome2_seqid\t" +
                                "genome1_breakpoint\t" +
                                "genome2_breakpoint\t" +
                                "count\t" +
                                "junction1\t" +
                                "junction2\t" +
                                "gene1\t" +
                                "gene2\n")
            group_breakpoints(args,group_in_fname,group_outFP)
            group_outFP.close()

    except Exception as e:
        print(e)
        sys.exit(1)

    if os.path.exists("tmp1"):
        os.remove("tmp1")
    if os.path.exists("tmp2"):
        os.remove("tmp2")

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
    parser.add_argument('-exonic',
                        required=False,
                        action='store_true',
                        help="Only consider exonic regions when extracting gene coordinates. Mapped coordinates of the read are matched against the gene coordinate set and if any position overlaps - the match will be reported. By default, exonic and intronic positions are considered. If this flag is enabled and read is contained entirely within an intron - it will not be reported as a match for the gene.")
    parser.add_argument('-two_pass',
                        required=False,
                        action='store_true',
                        help="Run two pass approach. First pass will find all possible breakpoints. Second pass will try to match breakpoints that are close to each other.")
    parser.add_argument('-group',
                        required=False,
                        action='store_true',
                        help="If enabled, will output a file with breakpoints groupped by position.")
    parser.add_argument('-o',
                        required=True,
                        type=str,
                        help="Output file.")
    
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])