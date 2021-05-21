#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Author: Ales Varabyou
"""

# this script simply pools the latest human annotation sources - unzips and standardizes them into GTF-formatted UCSC nomenclature files

import urllib.request as request
from contextlib import closing
import pandas as pd
import subprocess
import argparse
import requests
import zipfile
import tarfile
import shutil
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

    return seqid

def remove_gtf_seqids(fname,seqids): # removes records with the specified IDs from the GTF/GFF file inplace
    with open(fname+"tmp","w+") as outFP:
        with open(fname,"r") as inFP:
            for line in inFP.readlines():
                seqid = line.split("\t")[0]
                if seqid in seqids:
                    continue
                outFP.write(line)
    shutil.move(fname+"tmp",fname)

def fetch_annotations(args):
    output_dir = args.output.rstrip("/")+"/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # fetch mapping file for nomenclature conversions
    assembly_stats_url = 'https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_assembly_report.txt'
    req = requests.get(assembly_stats_url,allow_redirects=False)
    if req.status_code==200:
        with open(output_dir+"mapfile.raw.txt","wb") as outFP:
            outFP.write(req.content)
    else:
        assert False,"Invalid mapfile URL: "+assembly_stats_url

    # load mapping file
    chrMapCols=["name","role","molecule","type","genbank","rel","refseq","unit","seqLen","ucsc"]
    chrMap=pd.read_csv(output_dir+"mapfile.raw.txt",sep="\t",names=chrMapCols,comment="#")

    # read the setup
    urls = dict()
    with open(args.setup,"r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            cols = line.rstrip("\n").split(",")
            assert len(cols)==2,"incorrect setup file format"
            filename = output_dir+cols[1].split("/")[-1]
            label = cols[0]
            url = cols[1]
            urls[cols[0]]={"label":label,"filename":filename,"url":url}


    # fetch data
    for label,item in urls.items():
        if os.path.exists(item["url"]):
            shutil.copy(item["url"],item["filename"])
        elif item["url"].startswith("http"):
            req = requests.get(item["url"],allow_redirects=False)
            if req.status_code == 200:
                with open(item["filename"],"wb") as outFP:
                    outFP.write(req.content)
                assert os.path.exists(item["filename"])
            else:
                assert False,'Invalid URL in the setup file: '+item["url"]
        elif item["url"].startswith("ftp"):
            with closing(request.urlopen(item["url"])) as inFP:
                with open(item["filename"],"wb") as outFP:
                    shutil.copyfileobj(inFP,outFP)
        else:
            assert False,"unrecognized url format: "+item["url"]


    # standardize
    for label,item in urls.items():
        splits = item["filename"].split(".")
        assert len(splits)>=2,"no extension found in the filename: "+item["filename"]
        extension = splits[-1].lower()
        if extension in ["gtf","gff","gff3","gff2"]:
            continue
        else:
            if len(item["filename"].split("/")[-1])>7 and item["filename"][-7:].lower()==".tar.gz":
                # run tar -xvzf
                with tarfile.open(item["filename"], "r:gz") as inFP:
                    namelist = inFP.getnames()
                    found=False
                    for n in namelist:
                        if n==item["filename"].split("/")[-1].rstrip(".tar.gz"):
                            inFP.extract(item["filename"].split("/")[-1].rstrip(".tar.gz"),output_dir)
                            found=True
                    assert found,"did not locate valid file in the archive: "+item["filename"]

                os.remove(item["filename"])
                item["filename"]=item["filename"].rstrip(".tar.gz")

            elif len(item["filename"].split("/")[-1])>4 and item["filename"][-4:].lower()==".tar":
                with tarfile.open(item["filename"], "r") as inFP:
                    namelist = inFP.getnames()
                    found=False
                    for n in namelist:
                        if n==item["filename"].split("/")[-1].rstrip(".tar"):
                            inFP.extract(item["filename"].split("/")[-1].rstrip(".tar"),output_dir)
                            found=True
                    assert found,"did not locate valid file in the archive: "+item["filename"]

                os.remove(item["filename"])
                item["filename"]=item["filename"].rstrip(".tar")

            elif len(item["filename"].split("/")[-1])>3 and item["filename"][-3:].lower()==".gz":
                # run gunzip
                with gzip.open(item["filename"],'rb') as inFP:
                    with open(item["filename"].rstrip(".gz"),'wb') as outFP:
                        shutil.copyfileobj(inFP, outFP)

                os.remove(item["filename"])
                item["filename"]=item["filename"].rstrip(".gz")

            elif len(item["filename"].split("/")[-1])>4 and item["filename"][-4:].lower()==".zip":
                # run unzip
                with zipfile.ZipFile(item["filename"],'r') as inFP:
                    namelist = inFP.namelist()
                    found = False
                    for n in namelist:
                        if n == item["filename"].split("/")[-1].rstrip(".zip"):
                            inFP.extract(item["filename"].split("/")[-1].rstrip(".zip"),output_dir)
                            found=True
                    assert found,"did not locate valid file in the archive: "+item["filename"]

                os.remove(item["filename"])
                item["filename"]=item["filename"].rstrip(".zip")

            else:
                assert False,"unknown file extension detected from url: "+item["url"]

    # verify all files were extracted correctly and standardize them
    for label,item in urls.items():
        assert os.path.exists(item["filename"]),"extracted path does not exist: "+item["filename"]

        # detect which nomenclature is used
        tmp_map_fname = output_dir+"tmp.tsv"
        wrong_seqids = get_nomenclature(item["filename"],chrMap,tmp_map_fname)
        if len(wrong_seqids)>0: # wrong IDs found - remove them from the annotation
            remove_gtf_seqids(item["filename"],wrong_seqids)


        gffread_cmd = [args.gffread,"-T",
                       "-m",tmp_map_fname,
                       "-o",output_dir+label+".gtf",
                       item["filename"]]
        ret = subprocess.call(gffread_cmd)
        assert ret==0,"non-0 return code from gffread"

        os.remove(output_dir+"tmp.tsv")
        os.remove(item["filename"])
        item["filename"] = output_dir+label+".gtf"

    os.remove(output_dir+"mapfile.raw.txt")

    return

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("-s",
                        "--setup",
                        required=True,
                        help="File with a list of urls and labels in a CSV format. The first column of each column is the label to be used throughout analysis and the second column is the url. Lines that start with # are interpreted as comments and skipped")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        help="Output directory")
    parser.add_argument("--gffread",
                        required=False,
                        type=str,
                        default="gffread",
                        help="Path to the gffread executable")
    parser.set_defaults(func=fetch_annotations)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
