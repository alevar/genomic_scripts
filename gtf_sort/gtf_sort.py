#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Author: Ales Varabyou
"""

# this script sorts transcripts in a GTF file by chromosome/strand/start/end
# preserving the hierarchy of the transcripts
# the script will only retain features: "transcript","exon","CDS"

import pandas as pd
import numpy as np
import argparse
import csv
import sys
import os

def gtf_sort(args):
    gff3cols=["seqid","source","type","start","end","score","strand","phase","attributes"]
    df=pd.read_csv(args.input,sep="\t",names=gff3cols,comment="#")

    df = df[df["type"].isin(["transcript","exon","CDS"])].reset_index(drop=True)

    # check if GFF or GTF
    gff3 = False
    if "Parent=" in df.iloc[0].attributes:
        gff3=True

    assert not gff3,"Error: detected GFF3 format. This script only supports GTF. You may convert between formats using gffread"

    df["parent"]=df["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    df["type"]=pd.Categorical(df["type"],categories=["transcript","exon","CDS"],ordered=True)
    df.sort_values(by=["seqid","strand","parent","type","start","end"],ascending=True,inplace=True)
    df[gff3cols].to_csv(args.output,sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)

def main(args):
    parser=argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("-i",
                       "--input",
                       required=True,
                       help="annotation in a GFF/GTF format. The script will only retain features of types: 1. transcript; 2. exon; 3. CDS. Everything else will be discarded")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        help="output file for the intronic annotation")
    parser.set_defaults(func=gtf_sort)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])