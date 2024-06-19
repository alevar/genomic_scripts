#!/usr/bin/env python3

import os
import subprocess
import sys
import shutil
import argparse
import datetime
import pandas as pd

def gtf_or_gff(file_path):
    """
    Checks whether a file is in GTF or GFF format.

    Args:
    file_path (str): Path to the file to check.

    Returns:
    str: 'GTF', 'GFF', or 'Unknown' based on the file format.
    """
    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Skip comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                
                # Split the line into fields
                lcs = line.strip().split('\t')
                
                # Check the number of fields
                if len(lcs) != 9:
                    return None
                
                # Check for GTF-specific attributes (field 9 contains key-value pairs separated by semicolons)
                if 'gene_id \"' in lcs[8] and 'transcript_id \"' in lcs[8]:
                    return 'gtf'
                
                # Check for GFF-specific format (attribute field contains key-value pairs separated by semicolons, but without "gene_id" or "transcript_id")
                if 'ID=' in lcs[8] or 'Parent=' in lcs[8]:
                    return 'gff'
                
            return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def get_comments(number: str) -> str:
    assert len(number.split(".")) == 3, "Invalid release number"
    comment_str = "##NAME: CHESS\n" + \
                  "##VERSION: "+number+"\n" + \
                  "##DESCRIPTION: Comprehensive Human Expressed SequenceS\n" + \
                  "##DATE: "+datetime.datetime.now().strftime("%Y-%m-%d")+"\n" + \
                  "##CONTACT: ales[dot]varabyou[at]jhu[dot][edu],mpertea[at]jhu[dot][edu]\n"
    return comment_str

def extract_attributes(attribute_str:str,gff=False)->dict:
    """
    This function extracts attributes from an attribute string. Assumes that the only information passed is the 9th column of the GTF file;
    
    Parameters:
    attribute_str (str): The attribute string to extract attributes from.
    gff (bool, optional): A flag indicating whether the attribute string is from a GFF file. Defaults to False.

    Returns:
    dict: A dictionary of attributes extracted from the attribute string.
    """
    attrs = attribute_str.rstrip().rstrip(";").split(";")
    attrs = [x.strip() for x in attrs]
    attrs = [x.strip("\"") for x in attrs]
    attrs_dict = dict()
    sep = " \""
    if gff:
        sep = "="
    for at in attrs:
        try:
            k,v = at.split(sep)
            attrs_dict.setdefault(k,v)
        except:
            continue
        
    return attrs_dict

def to_attribute_string(attrs: dict, for_gff: bool, feature_type: str) -> str:
    return "; ".join([f'{key}="{value}"' for key, value in attrs.items()])

def gtf2gff(gtf_fname: str, gff_fname: str):
    with open(gff_fname, "w+") as outFP:
        outFP.write("##gff-version 3\n")
        outFP.write("#!gff-spec-version 1.21\n")
        with open(gtf_fname, "r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if lcs[0] == "#":
                    outFP.write(line)
                if not len(lcs) == 9:
                    continue
                attrs = extract_attributes(lcs[8])
                res_line = "\t".join(lcs[:-1]) + "\t" + to_attribute_string(attrs, True, lcs[2])
                outFP.write(res_line + "\n")

def get_attribute(gtf_fname:str,attrs,cols:list=None,feature:str="transcript",gff:bool=False,as_dict=False) -> pd.DataFrame:
    """
    This function extracts attributes from a GTF file and loads them into a dictionary grouped by transcript ID.

    Parameters:
    gtf_fname (str): The name of the GTF file to load.
    attrs (list): The attributes to extract.
    cols (list, optional): The columns to include in the output DataFrame. If None, all columns will be included. Defaults to None.
    feature (str, optional): The feature to filter lines by. Only lines where the third column equals this value will be processed. Defaults to "transcript".
    gff (bool, optional): A flag indicating whether the file is in GFF format. Defaults to False.

    Returns:
    pd.DataFrame: A Pandas DataFrame containing the extracted attributes.
    """
    clean_attrs = []
    if type(attrs)==list:
        clean_attrs = attrs
    else: # is string
        clean_attrs = [attrs]
        
    tid2name = dict()
    with open(gtf_fname,"r") as inFP:
        for line in inFP:
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            if lcs[2]==feature:
                tid = ""
                if not gff:
                    tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                else:
                    tid = lcs[8].split("ID=",1)[1].split(";",1)[0]
                
                tid2name[tid]=[]
                if cols is not None:
                    for c in cols:
                        tid2name[tid].append(lcs[c])
                cur_attrs = extract_attributes(lcs[8],gff)
                for a in clean_attrs:
                    gn = cur_attrs.get(a,"-")
                    tid2name[tid].append(gn)
                    
    if as_dict:
        return tid2name
    
    res = pd.DataFrame.from_dict(tid2name,orient="index").reset_index()
    tmp_cols = ["tid"]
    if cols is not None:
        for c in cols:
            tmp_cols.append(c)
    for a in clean_attrs:
        tmp_cols.append(a)
    res.columns = tmp_cols
    return res

def load_gene_descriptions(fname: str) -> dict:
    # get  gene names from the refseq gtf
    gene_desc = get_attribute(fname,["Name","description"],cols=None,feature="gene",gff=True)
    gene_desc = gene_desc[gene_desc["description"]!="-"][["Name","description"]].reset_index(drop=True)
    gene_desc.drop_duplicates(inplace=True)
    gene_desc_dict = gene_desc.set_index('Name').to_dict()['description']
    return gene_desc_dict

def run(args):
    # setup output directory
    outdir = os.path.abspath(args.outdir)+"/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # check that the input is in the GTF format
    assert os.path.exists(args.gtf_file), "Input file not found"
    assert gtf_or_gff(args.gtf_file) == 'gtf', "Input file is not in GTF format"

    # load gene descriptions from the reference file
    assert os.path.exists(args.gene_description_reference), "Gene description reference file not found"
    gene_desc_dict = load_gene_descriptions(args.gene_description_reference)

    # assert that contigs lengths file exists
    assert os.path.exists(args.contig_lengths), "Contig lengths file not found"

    # load a list of genes
    genes = dict()

    with open(args.gtf_file, "r") as inFP:
        for line in inFP:
            if line[0] == "#":
                continue

            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            attrs = extract_attributes(lcs[8])
            gid = attrs["gene_id"]
            tid = attrs["transcript_id"]

            if lcs[2] == "transcript":
                genes.setdefault(gid, [dict({"seqid": None,
                                             "start": sys.maxsize,
                                             "end": 0,
                                             "strand": None,
                                             "name": set(),
                                             "type": set()}), dict()])
                genes[gid][0]["start"] = min(genes[gid][0]["start"], int(lcs[3]))
                genes[gid][0]["end"] = max(genes[gid][0]["end"], int(lcs[4]))

                if genes[gid][0]["seqid"] is not None:
                    assert genes[gid][0]["seqid"] == lcs[0], "wrong seqid: " + gid
                else:
                    genes[gid][0]["seqid"] = lcs[0]

                assert lcs[6] in ["+", "-"], "wrong strand: " + tid + " : " + lcs[6]
                if genes[gid][0]["strand"] is not None:
                    if genes[gid][0]["strand"] != lcs[6]:
                        genes[gid][0]["strand"] = "."
                else:
                    genes[gid][0]["strand"] = lcs[6]

                if "gene_name" in attrs:
                    genes[gid][0]["name"].add(attrs["gene_name"])

                if "gene_type" in attrs:
                    genes[gid][0]["type"].add(attrs["gene_type"])

                genes[gid][1].setdefault(tid, "")

            assert tid in genes[gid][1], "wrong tid (unsorted?): " + tid
            genes[gid][1][tid] += line

    # write out a version of the GTF file with gene information
    with_genes_gtf_fname = outdir + "1.with_genes.gtf"
    with open(with_genes_gtf_fname, "w+") as outFP:
        outFP.write(get_comments(args.release_number))
        with open(args.gtf_file, "r") as inFP:
            for line in inFP:
                if line[0] == "#":
                    outFP.write(line)
                else:
                    break

        for gid, gv in genes.items():
            gline = gv[0]["seqid"] + "\t" + \
                    "CHESS" + "\t" + \
                    "gene" + "\t" + \
                    str(gv[0]["start"]) + "\t" + \
                    str(gv[0]["end"]) + "\t" + \
                    "." + "\t" + \
                    gv[0]["strand"] + "\t" + \
                    "." + "\t" + \
                    "gene_id \"" + gid + "\";"

            if len(gv[0]["type"]) > 0:
                gline += " gene_type \"" + ", ".join(list(gv[0]["type"])) + "\";"

            if len(gv[0]["name"]) > 0:
                gline += " gene_name \"" + ", ".join(list(gv[0]["name"])) + "\";"
                desc = list()
                for n in gv[0]["name"]:
                    if n in gene_desc_dict:
                        desc.append([n, gene_desc_dict[n]])

                if len(desc) == 1:
                    gline += " description \"" + desc[0][1] + "\";"
                if len(desc) > 1:
                    print(gv[0]["name"])
                    gline += " description \"" + ", ".join([x[0] + " : " + x[1] for x in desc]) + "\";"

            outFP.write(gline + "\n")
            for tid, tv in gv[1].items():
                outFP.write(tv)

    # convert the GTF file to GFF
    with_genes_gff_fname = outdir + "1.with_genes.gff"
    gtf2gff(with_genes_gtf_fname, with_genes_gff_fname)

    os.environ['LD_LIBRARY_PATH'] = "/ccb/sw/lib/:LD_LIBRARY_PATH"

    genePred_fname = outdir + "2.with_genes.genePred"
    cmd = ["gff3ToGenePred", with_genes_gff_fname, genePred_fname]
    subprocess.call(cmd)

    bedPlus_fname = outdir + "2.with_genes.bedPlus"
    cmd = ["genePredToBigGenePred", genePred_fname, bedPlus_fname]
    subprocess.call(cmd)

    sorted_bedPlus_fname = outdir + "2.with_genes.srt.bedPlus"
    cmd = ["bedSort", bedPlus_fname, sorted_bedPlus_fname]
    subprocess.call(cmd)

    os.rename(sorted_bedPlus_fname, bedPlus_fname)

    bb_fname = outdir + "2.with_genes.bb"
    cmd = ["bedToBigBed", "-type=bed12+8", "-tab", "-as="+args.bed_as, bedPlus_fname,
           args.contig_lengths, bb_fname, "-extraIndex=name"]
    subprocess.call(cmd)

    base_name = args.prefix+args.release_number+"."+args.suffix
    rdir = os.path.join(outdir, base_name+"/")
    if not os.path.exists(rdir):
        os.makedirs(rdir)

    shutil.copy(args.gtf_file, os.path.join(rdir, base_name+".gtf"))
    os.rename(with_genes_gff_fname, os.path.join(rdir, base_name+".gff"))
    os.rename(bb_fname, os.path.join(rdir, base_name+".bb"))

    with open(os.path.join(rdir, base_name+".primary.gtf"), "w+") as outFP:
        with open(os.path.join(rdir, base_name+".gtf"), "r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if lcs[0] == "#":
                    outFP.write(line)
                if not len(lcs) == 9:
                    continue

                if "_alt" in lcs[0]:
                    continue

                outFP.write(line)

    with open(os.path.join(rdir, base_name+".primary.gff"), "w+") as outFP:
        with open(os.path.join(rdir, base_name+".gff"), "r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if lcs[0] == "#":
                    outFP.write(line)
                if not len(lcs) == 9:
                    continue

                if "_alt" in lcs[0]:
                    continue

                outFP.write(line)

    cmd = ["gffread", "-g", args.genome, "-y",
           os.path.join(rdir, base_name+".protein.fa"), os.path.join(rdir, base_name+".primary.gtf")]
    subprocess.call(cmd)

def main(args):
    parser = argparse.ArgumentParser(description="Process GTF file to generate various derivative files and create a new release of the CHESS annotation.")
    parser.add_argument("--gtf_file", help="Input GTF file")
    parser.add_argument("--outdir", help="Base directory for output files")
    parser.add_argument("--release_number", help="Release number")
    parser.add_argument("--gene_description_reference", help="GFF3 file containing gene records which can be used to provide gene descriptions for CHESS genes")
    parser.add_argument("--contig_lengths", help="File containing the lengths of the contigs in the reference genome")
    parser.add_argument("--prefix", help="Prefix for the output files")
    parser.add_argument("--suffix", help="Suffix for the output files")
    parser.add_argument("--genome", help="Genome file")
    parser.add_argument("--bed_as", help="Non-standard bed format for running bedToBigBed (consult bedToBigBed documentation for more information)")
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)
    
    # TODO:
    # 1. Compute summaries (# tx, cds, gid, etc)
if __name__=="__main__":
    main(sys.argv[1:])
