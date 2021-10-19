#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Author: Ales Varabyou
"""

# this script simply pools the latest human annotation sources - unzips and standardizes them into GTF-formatted UCSC nomenclature files

from contextlib import closing
import pandas as pd
import subprocess
import argparse
import requests
import zipfile
import tarfile
import shutil
import urllib
import gzip
import sys
import os

def get_nomenclature(ann_fname,chr_map,out_fname):
    nomenclature = None

    seqids = set()
    with open(ann_fname) as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            seqid = line.split("\t",1)[0]
            seqids.add(seqid)

    wrong_seqids = set()
    sub_map = dict()
    for seqid in seqids:
        found_unique = False
        for col in list(chr_map.columns):
            tmp = set(chr_map[chr_map[col]==seqid]["ucsc"])
            if len(tmp)==1:
                sub_map[seqid]=list(tmp)[0]
                found_unique=True
                break
        if not found_unique:
            print("multiple choices for sequence: "+seqid+" in source: "+ann_fname)
            wrong_seqids.add(seqid)

    with open(out_fname,"w+") as outFP:
        for k,v in sub_map.items():
            outFP.write(k+" "+v+"\n")

    return wrong_seqids

def remove_gtf_seqids(fname,out_fname,wrong_seqids): # removes records with the specified IDs from the GTF/GFF file inplace
    with open(out_fname,"w+") as outFP:
        with open(fname,"r") as inFP:
            for line in inFP:
                seqid = line.split("\t")[0]
                if seqid in wrong_seqids:
                    continue
                outFP.write(line)

def fetch_annotations(args):
    
    # load mapping file
    chrMapCols=["name","role","molecule","type","genbank","rel","refseq","unit","seqLen","ucsc"]
    chrMap=pd.read_csv(args.map,sep="\t",names=chrMapCols,comment="#")

    # standardize
    splits = args.input.split(".")
    assert len(splits)>=2,"no extension found in the filename: "+args.input
    extension = splits[-1].lower()
    assert extension in ["gtf","gff","gff3","gff2"],"unknown file extension detected from url: "+item["url"]


    # detect which nomenclature is used
    tmp_map_fname = args.output+".map.tsv"
    wrong_seqids = get_nomenclature(args.input,chrMap,tmp_map_fname)
    remove_gtf_seqids(args.input,args.output+".tmp",wrong_seqids)


    gffread_cmd = [args.gffread,"-T",
                   "-m",tmp_map_fname,
                   "-o",args.output,
                   args.output+".tmp"]
    ret = subprocess.call(gffread_cmd)
    assert ret==0,"non-0 return code from gffread"

    if not args.keep_tmp:
        os.remove(tmp_map_fname)
        os.remove(args.output+".tmp")

    return

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("-i",
                        "--input",
                        required=True,
                        help="input annotation")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        help="Output annotation")
    parser.add_argument("-m",
                        "--map",
                        required=True,
                        help="Mapfile")
    parser.add_argument("--gffread",
                        required=False,
                        type=str,
                        default="gffread",
                        help="Path to the gffread executable")
    parser.add_argument("--keep-tmp",
                        action="store_true",
                        help="do not remove temporary data")
    parser.set_defaults(func=fetch_annotations)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
