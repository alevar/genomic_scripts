#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Author: Ales Varabyou
"""

# this script performs pairwise comparison between multiple annotations and computes and displays several key statistics

from upsetplot import generate_counts
from upsetplot import plot
import pandas as pd
import numpy as np
import subprocess
import itertools
import argparse
import shutil
import sys
import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sbn

gff3cols=["seqid","source","type","start","end","score","strand","phase","attributes"]

def get_combs(names):
    all_combinations = [list(itertools.combinations(names,x)) for x in range(1,len(names)+1)]
    res = list()
    for comb_set in all_combinations:
        for comb in comb_set:
            res.append([list(set(comb)),list(set(names)-set(comb))])

    return res

def get_mat(names):
    tmp = list(itertools.product([False, True], repeat=len(names)))
    index = pd.MultiIndex.from_tuples(tmp,names=names)
    series = pd.Series(np.zeros(len(tmp)),index=index)
    return series

def get_not_none(d):
    res = set()
    for k,v in d.items():
        if not v is None:
            res.add(k)

    return res

def get_set(tm,c,c_rev): # get transcripts in c but not in union of c_rev
    res = set(tm[c[0]][list(tm[c[0]])[0]])
    for name in c: # from requested union
        if name == c[0]: # skip itself
            continue

        res = res.intersection(get_not_none(tm[c[0]][name]))

    for name in c_rev: # remove unwanted
        res = res - get_not_none(tm[c[0]][name])

    return res

def read_setup(fname):
    res = list()
    assert os.path.exists(fname),"setup file not found: "+fname
    with open(fname,"r") as inFP:
        for line in inFP:
            label,path = line.strip().split(",")
            assert os.path.exists(path),"file not found: "+path
            res.append(tuple([label,path]))
    assert len(res)>=2,"setup file must have at least two entries"
    return res

def gffcmp_multi(args):
    gtf_list = read_setup(args.setup)
    input_label = args.input
    input_fname = None
    for x in gtf_list:
        if x[0]==input_label:
            input_fname=x[1]
            break
    assert input_fname is not None,"requested input label not found in the setup file"

    gtf_pairs = list(itertools.permutations(gtf_list,2))
    for gp in gtf_pairs:
        gffcmp_cmd = [args.gffcompare,"--no-merge",
                      "-o",args.output+gp[0][0]+"_"+gp[1][0],
                      "-r",gp[0][1],
                      gp[1][1]]
        subprocess.call(gffcmp_cmd)

        fp_dir = "/".join(gp[1][1].split("/")[:-1])+"/"
        fp_name_base = gp[1][1].split("/")[-1]
        out_base = args.output.split("/")[-1]
        assert os.path.exists(fp_dir+out_base+gp[0][0]+"_"+gp[1][0]+"."+fp_name_base+".tmap"),"tmap/refmap do not exist"
        shutil.move(fp_dir+out_base+gp[0][0]+"_"+gp[1][0]+"."+fp_name_base+".tmap",args.output+gp[0][0]+"_"+gp[1][0]+".tmap")
        shutil.move(fp_dir+out_base+gp[0][0]+"_"+gp[1][0]+"."+fp_name_base+".refmap",args.output+gp[0][0]+"_"+gp[1][0]+".refmap")

    # cycle through annotations and load all transcripts in
    tx_map = dict()
    for gp in gtf_pairs:
        gtf_fname = gp[0][1]
        rn = gp[0][0]
        tn = gp[1][0]
        tx_map.setdefault(rn,dict())
        tx_map[rn].setdefault(tn,dict())
        with open(gtf_fname,"r") as inFP:
            for line in inFP:
                cols = line.strip().split("\t")
                attrs = cols[8]
                tmp = attrs.split("transcript_id \"",1)
                if len(tmp)==1:
                    continue
                tid = tmp[-1].split("\"",1)[0]
                tx_map[rn][tn][tid]=None

    # get complete map between all transcripts in comparison
    for gp in gtf_pairs:
        rn = gp[0][0]
        tn = gp[1][0]
        refmap_fname = args.output+rn+"_"+tn+".refmap"
        with open(refmap_fname,"r") as inFP:
            for line in inFP:
                cols = line.strip().split("\t")
                code = cols[2]
                if not code == "=":
                    continue
                rgid = cols[0]
                rtid = cols[1]
                ttids = [x.split("|")[1] for x in cols[3].split(",")]
                assert rtid in tx_map[rn][tn],"tid not preloaded: "+rtid+" for "+rn+"/"+tn
                assert tx_map[rn][tn][rtid]==None,"duplicates in: "+rtid
                tx_map[rn][tn][rtid] = ttids

    combs = get_combs(list(tx_map))
    mat = get_mat(list(tx_map))
    index_names = mat.index.names

    for comb,comb_rev in combs:
        res = get_set(tx_map,comb,comb_rev)
        idx_loc = tuple([True if x in comb else False for x in index_names])
        mat.loc[idx_loc] = len(res)

    plt.close('all')
    plt.clf()
    plot(mat)
    plt.savefig(args.output+".png")

    ref_labels = set(tx_map)-set([input_label])
    with open(args.output+".gtf","w+") as outFP:
        with open(input_fname,"r") as inFP:
            for line in inFP:
                line = line.rstrip()
                if line[0]=="#":
                    outFP.write(line+"\n")
                else:
                    cols = line.split("\t")
                    if cols[2]=="transcript":
                        attrs = cols[8]
                        assert attrs[-1]==";","invalid format"
                        tmp = attrs.split("transcript_id \"",1)
                        assert len(tmp)>1,"invalid format?"
                        tid = tmp[-1].split("\"",1)[0]
                        for rl in ref_labels:
                            rtid = tx_map[input_label][rl][tid]
                            if not rtid is None:
                                rtid = ",".join(rtid)
                            attrs+=" "+rl+" \""+str(rtid)+"\";"
                        cols[8] = attrs
                        outFP.write("\t".join(cols)+"\n")
                    else:
                        outFP.write(line+"\n")


def main(args):
    parser=argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("-i",
                       "--input",
                       required=True,
                       help="Name of the file to be annotated as provided in the first column of the '-s/--setup' argument.")
    parser.add_argument("-s",
                        "--setup",
                        required=True,
                        help="Setup file with a table of all file names to be used in CSV format. The first column is used to specify the custom label for each source. The second column is used to provide full file path for each source")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        help="basename for the outputs")
    parser.add_argument("--gffcompare",
                        required=False,
                        type=str,
                        default="gffcompare",
                        help="path to the gffcompare executable")
    parser.set_defaults(func=gffcmp_multi)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])
