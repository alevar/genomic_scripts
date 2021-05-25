#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Author: Ales Varabyou
"""

# this script performs pairwise comparison between multiple annotations and computes and displays several key statistics

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
import matplotlib.pyplot as plt

mpl.use('Agg')

gff3cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

def get_combs(names):
    all_combinations = [list(itertools.combinations(names, x)) for x in range(1, len(names) + 1)]
    res = list()
    for comb_set in all_combinations:
        for comb in comb_set:
            res.append([list(set(comb)), list(set(names) - set(comb))])

    return res


def get_mat(names,keep_all_false=True):
    tmp = list(itertools.product([False, True], repeat=len(names)))
    index = pd.MultiIndex.from_tuples(tmp, names=names)
    series = pd.Series(np.zeros(len(tmp)), index=index)
    if not keep_all_false:
        series = series[~(series==False)] # remove the first row with all Falses
    return series


def get_not_none(d):
    res = set()
    for k, v in d.items():
        if v is not None:
            res.add(k)

    return res


def get_set(tm, c, c_rev):  # get transcripts in c but not in union of c_rev
    res = set(tm[c[0]][list(tm[c[0]])[0]])
    for name in c:  # from requested union
        if name == c[0]:  # skip itself
            continue

        res = res.intersection(get_not_none(tm[c[0]][name]))

    for name in c_rev:  # remove unwanted
        res = res - get_not_none(tm[c[0]][name])

    return res

def get_set_def(m, c, c_rev):  # get entities in c but not in union of c_rev
    res = set(m[c[0]])
    for name in c:  # from requested union
        res = res.intersection(m[name])

    for name in c_rev:  # remove unwanted
        res = res - m[name]

    return res

def get_introns_exons(fname,coding_only=False): # if set to True "coding_only" will use features with type=="CDS" (default is type=="exon")
    introns = set()
    exons = set()

    with open(fname,"r") as inFP:
        cur_tid = ""
        cur_intron = ""
        for line in inFP:
            line = line.rstrip("\n")
            cols = line.split("\t")

            tmp = cols[8].split("transcript_id \"", 1)
            if len(tmp) == 1:
                continue
            tid = tmp[-1].split("\"", 1)[0]

            if cols[2]=="transcript":
                cur_tid = tid
                cur_intron = ""

            elif cols[2]=="exon" and not coding_only:
                assert tid==cur_tid,"incorrect hierarchy: lower-level tid does not match upper-level tid (current exon tid != last transcript tid)"
                exons.add(cols[0]+cols[6]+cols[3]+"-"+cols[4])
                if len(cur_intron)==0:
                    cur_intron+=cols[0]+cols[6]+cols[4]+"-"
                else:
                    cur_intron+=cols[3]
                    introns.add(cur_intron)
                    cur_intron = ""

            elif cols[2]=="CDS" and coding_only:
                assert tid==cur_tid,"incorrect hierarchy: lower-level tid does not match upper-level tid (current exon tid != last transcript tid)"
                exons.add(cols[0]+cols[6]+cols[3]+"-"+cols[4])
                if len(cur_intron)==0:
                    cur_intron+=cols[0]+cols[6]+cols[4]+"-"
                else:
                    cur_intron+=cols[3]
                    introns.add(cur_intron)
                    cur_intron = ""

            else:
                continue

    return introns,exons

def read_setup(fname):
    res = list()
    assert os.path.exists(fname), "setup file not found: " + fname
    with open(fname, "r") as inFP:
        for line in inFP:
            label, path = line.strip().split(",")
            assert os.path.exists(path), "file not found: " + path
            res.append(tuple([label, path]))
    assert len(res) >= 2, "setup file must have at least two entries"
    return res


def gffcmp_multi(args):
    output_dir = "/".join(args.output.split("/")[:-1])+"/"
    assert os.path.exists(output_dir),"output path is incorrect"

    fig_dir = "/".join(args.output.split("/")[:-1])+"/figs/"
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_out = fig_dir+args.output.split("/")[-1]

    raw_gtf_list = read_setup(args.setup)

    # run gffread on all inputs to standardize them
    gtf_list = list()
    labels = []
    for label,path in raw_gtf_list:
        assert os.path.exists(path),"input file does not exist: "+path
        gffread_cmd = [args.gffread,"-T",
                       "-o",output_dir+label+".gtf",
                       path]
        subprocess.call(gffread_cmd)
        gtf_list.append(tuple([label,output_dir+label+".gtf"]))
        labels.append(label)

    gtf_pairs = list(itertools.permutations(gtf_list, 2))
    for gp in gtf_pairs:
        gffcmp_cmd = [args.gffcompare, "--no-merge",
                      "-o", args.output + gp[0][0] + "_" + gp[1][0],
                      "-r", gp[0][1],
                      gp[1][1]]
        subprocess.call(gffcmp_cmd)

        fp_dir = "/".join(gp[1][1].split("/")[:-1]) + "/"
        fp_name_base = gp[1][1].split("/")[-1]
        out_base = args.output.split("/")[-1]
        assert os.path.exists(
            fp_dir + out_base + gp[0][0] + "_" + gp[1][0] + "." + fp_name_base + ".tmap"), "tmap/refmap do not exist"
        shutil.move(fp_dir + out_base + gp[0][0] + "_" + gp[1][0] + "." + fp_name_base + ".tmap",
                    args.output + gp[0][0] + "_" + gp[1][0] + ".tmap")
        shutil.move(fp_dir + out_base + gp[0][0] + "_" + gp[1][0] + "." + fp_name_base + ".refmap",
                    args.output + gp[0][0] + "_" + gp[1][0] + ".refmap")

    # compute exon overlaps and plot upset
    introns = dict()
    exons = dict()
    for label,path in gtf_list:
        introns[label],exons[label] = get_introns_exons(path)
    combs = get_combs(labels)
    mat_introns = get_mat(labels,False)
    mat_exons = get_mat(labels,False)
    index_names = mat_introns.index.names
    for comb,comb_rev in combs:
        res_introns = get_set_def(introns,comb,comb_rev)
        res_exons = get_set_def(exons,comb,comb_rev)
        idx_loc = tuple([True if x in comb else False for x in index_names])
        mat_introns.loc[idx_loc] = len(res_introns)
        mat_exons.loc[idx_loc] = len(res_exons)
    plt.close('all')
    plt.clf()
    plot(mat_exons,show_percentages=True,show_counts=False,sort_by="degree")
    plt.suptitle('Exon Set Comparison Between Sources')
    plt.savefig(fig_out + ".exons.eps")
    mat_exons.to_csv(fig_out+".exons.csv")

    # compute intron overlaps and plot upset
    plt.close('all')
    plt.clf()
    plot(mat_introns,show_percentages=True,show_counts=False,sort_by="degree")
    plt.suptitle('Intron Set Comparison Between Sources')
    plt.savefig(fig_out + ".introns.eps")
    mat_introns.to_csv(fig_out+".introns.csv")

    # cycle through annotations and load all transcripts in
    tx_map = dict()
    for gp in gtf_pairs:
        gtf_fname = gp[0][1]
        rn = gp[0][0]
        tn = gp[1][0]
        tx_map.setdefault(rn, dict())
        tx_map[rn].setdefault(tn, dict())
        with open(gtf_fname, "r") as inFP:
            for line in inFP:
                cols = line.strip().split("\t")
                attrs = cols[8]
                tmp = attrs.split("transcript_id \"", 1)
                if len(tmp) == 1:
                    continue
                tid = tmp[-1].split("\"", 1)[0]
                tx_map[rn][tn][tid] = None

    # get complete map between all transcripts in comparison
    for gp in gtf_pairs:
        rn = gp[0][0]
        tn = gp[1][0]
        refmap_fname = args.output + rn + "_" + tn + ".refmap"
        with open(refmap_fname, "r") as inFP:
            for line in inFP:
                cols = line.strip().split("\t")
                code = cols[2]
                if not code == "=":
                    continue
                rtid = cols[1]
                ttids = [x.split("|")[1] for x in cols[3].split(",")]
                assert rtid in tx_map[rn][tn], "tid not preloaded: " + rtid + " for " + rn + "/" + tn
                assert tx_map[rn][tn][rtid] is None, "duplicates in: " + rtid
                tx_map[rn][tn][rtid] = ttids

    combs = get_combs(list(tx_map))
    mat = get_mat(list(tx_map),False)
    index_names = mat.index.names

    for comb, comb_rev in combs:
        res = get_set(tx_map, comb, comb_rev)
        idx_loc = tuple([True if x in comb else False for x in index_names])
        mat.loc[idx_loc] = len(res)

    plt.close('all')
    plt.clf()
    plot(mat,show_percentages=True,show_counts=False,sort_by="degree")
    plt.suptitle('Transcript Set Comparison Between Sources')
    plt.savefig(fig_out + ".tx.eps")
    mat.to_csv(fig_out+".tx.csv")

    

    if args.input is not None: # annotate
        input_label = args.input
        input_fname = None
        for x in gtf_list:
            if x[0] == input_label:
                input_fname = x[1]
                break
        assert input_fname is not None, "requested input label not found in the setup file"

        ref_labels = list(set(tx_map) - {input_label})
        ref_labels = sorted(ref_labels)
        out_gtf_fp = open(args.output + ".gtf", "w+")
        out_map_fp = open(args.output + ".map.tsv", "w+")
        map_line = "\t".join([input_label] + ref_labels) + "\n"
        out_map_fp.write(map_line)
        with open(input_fname, "r") as inFP:
            for line in inFP:
                line = line.rstrip()
                if line[0] == "#":
                    out_gtf_fp.write(line + "\n")
                else:
                    cols = line.split("\t")
                    if cols[2] == "transcript":
                        attrs = cols[8]
                        assert attrs[-1] == ";", "invalid format"
                        tmp = attrs.split("transcript_id \"", 1)
                        assert len(tmp) > 1, "invalid format?"
                        tid = tmp[-1].split("\"", 1)[0]
                        map_line = tid + "\t"
                        for rl in ref_labels:
                            rtid = tx_map[input_label][rl][tid]
                            if rtid is not None:
                                rtid = ",".join(rtid)
                                map_line += rtid + "\t"
                            else:
                                map_line += "-\t"
                            attrs += " " + rl + " \"" + str(rtid) + "\";"
                        out_map_fp.write(map_line.rstrip("\t") + "\n")
                        cols[8] = attrs
                        out_gtf_fp.write("\t".join(cols) + "\n")
                    else:
                        out_gtf_fp.write(line + "\n")

        out_gtf_fp.close()
        out_map_fp.close()

    # cleanup

    for label,path in gtf_list: # cleanup gffread conversions
        os.remove(path)

    if not args.keep_tmp:
        for gp in gtf_pairs:
            os.remove(args.output + gp[0][0] + "_" + gp[1][0] + ".tmap")
            os.remove(args.output + gp[0][0] + "_" + gp[1][0] + ".refmap")
            os.remove(args.output + gp[0][0] + "_" + gp[1][0] + ".annotated.gtf")
            os.remove(args.output + gp[0][0] + "_" + gp[1][0] + ".stats")
            os.remove(args.output + gp[0][0] + "_" + gp[1][0] + ".tracking")
            os.remove(args.output + gp[0][0] + "_" + gp[1][0] + ".loci")


def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("-i",
                        "--input",
                        required=False,
                        help="Name of the file to be annotated as provided in the first column of the '-s/--setup' "
                             "argument.")
    parser.add_argument("-s",
                        "--setup",
                        required=True,
                        help="Setup file with a table of all file names to be used in CSV format. The first column is "
                             "used to specify the custom label for each source. The second column is used to provide "
                             "full file path for each source")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        help="basename for the outputs")
    parser.add_argument("--gffcompare",
                        required=False,
                        type=str,
                        default="gffcompare",
                        help="path to the gffcompare executable")
    parser.add_argument("--gffread",
                        required=False,
                        type=str,
                        default="gffread",
                        help="path to the gffread executable")
    parser.add_argument("--keep-tmp",
                        action="store_true",
                        help="keep all final temporary outputs (individual gffcompare results)")
    parser.set_defaults(func=gffcmp_multi)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
