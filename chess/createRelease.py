#!/usr/bin/env python3

import os
import subprocess
import sys
import shutil
import argparse
import datetime
import pandas as pd

from definitions import *

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

def to_attribute_string(attrs:dict,gff:bool=False,feature_type:str=None)->str:
    """
    This function converts a dictionary of attributes to an attribute string. Guarantees order of essential attributes.

    Parameters:
    attrs (dict): A dictionary of attributes.
    gff (bool, optional): A flag indicating whether to convert to GFF format. Defaults to False.
    feature_type (str, optional): The feature type of the GFF file. Defaults to None.

    Returns:
    str: An attribute string.
    """
    order = ["ID","Parent","transcript_id","gene_id","gene_name","gene_type","db_xref","description","max_TPM","sample_count","assembly_id","tag"]
    res = ""
    sep = " "
    quote = "\""
    end = "; "
    if gff:
        assert feature_type in ["gene","transcript","exon","CDS"],"wrong type: "+str(feature_type)
        sep = "="
        quote = ""
        end = ";"
        
    for k in order:
        if k in attrs:
            if gff:
                assert ";" not in attrs[k],"invalid character in attribute: "+attrs[k]
            
            if gff and feature_type=="gene" and k=="transcript_id":
                continue
            elif gff and feature_type=="gene" and k=="gene_id":
                res+="ID="+quote+attrs[k]+quote+end
            elif gff and feature_type=="transcript" and k=="transcript_id":
                res+="ID="+quote+attrs[k]+quote+end
            elif gff and feature_type=="transcript" and k=="gene_id":
                res+="Parent="+quote+attrs[k]+quote+end
            elif gff and feature_type in ["exon","CDS"] and k=="transcript_id":
                res+="Parent="+quote+attrs[k]+quote+end
            elif gff and feature_type in ["exon","CDS"] and k=="gene_id":
                continue
            else:        
                res+=k+sep+quote+attrs[k]+quote+end
    
    # add any other attributes in sorted order
    for k in sorted(list(attrs)):
        if k not in order:
            if gff:
                assert ";" not in attrs[k],"invalid character in attribute: "+attrs[k]
            res+=k+sep+quote+attrs[k]+quote+end
    
    if not gff:
        res = res.rstrip()
    if gff:
        res = res.rstrip(";")
    return res

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

def deduplicate_gtf(in_gtf_fname, out_gtf_fname, borf_fname):
    # removes duplicate transcripts (full-matching exon chain) and keeps a single copy
    # copy to keep is chosen if:
    # 1. CDS matches MANE
    # 2. CDS is closest to MANE
    # 3. CDS has optimal median of medians
    # 4. random choice

    # run borf on chess
    cmd = [borf_fname,"--use_geneid","--score","-g",in_gtf_fname,"-o",out_gtf_fname+".borf"]
    subprocess.call(cmd)
    
    c3 = get_chains(in_gtf_fname,"exon",True)
    c3["echain"] = c3["seqid"]+c3["strand"]+c3.apply(lambda row: ",".join(str(x[0])+"-"+str(x[1]) for x in row["chain"]),axis=1)
    c3["ichain"] = c3["seqid"]+c3["strand"]+c3.apply(lambda row: ",".join(str(x[0])+"-"+str(x[1]) for x in chain_inv(row["chain"])),axis=1)
    c3["ichain"] = np.where(c3["ichain"].str.endswith(("-","+")),"-",c3["ichain"])
    c3 = c3[["tid","has_cds","echain","ichain"]]

    c3 = c3.merge(get_attribute(in_gtf_fname,["gene_id"]),on="tid",how="left")
    c3 = c3.merge(get_attribute(out_gtf_fname+".borf.gtf",["mom"]),on="tid",how="left")

    # attach cds chain to it as well
    c3_cds = get_chains(in_gtf_fname,"CDS",True)
    c3_cds["cchain"] = c3_cds["seqid"]+c3_cds["strand"]+c3_cds.apply(lambda row: ",".join(str(x[0])+"-"+str(x[1]) for x in row["chain"]),axis=1)
    c3_cds = c3_cds[["tid","cchain"]]
    c3_cds["cchain"] = np.where(c3_cds["cchain"].str.endswith(("-","+")),"-",c3_cds["cchain"])

    # merge
    c3 = c3.merge(c3_cds,on="tid")
    
    # compute similarity to mane and decide which ones to keep
    # for every gene, we need to get MANE transcript and CDS chain
    # load transcripts with source and gene_id
    mane_cmp_df = get_attribute(in_gtf_fname,["gene_id"],cols=[1])
    mane_cmp_df.rename({1:"source"},axis=1,inplace=True)

    # split off mane and join on the dataframe to mark mane_tid for each transcript
    mane_cmp_df = mane_cmp_df.merge(mane_cmp_df[mane_cmp_df["source"]=="MANE"][["tid","gene_id"]].rename({"tid":"tid_mane"},axis=1),on="gene_id",how="left")

    mane_cmp_df = mane_cmp_df[~(mane_cmp_df["tid_mane"].isna())].reset_index(drop=True)
    
    # load chains for all transcripts
    cds_chains = get_chains(in_gtf_fname,"CDS",True)
    
    # merge onto the dataframe mane_cmp
    mane_cmp_df = mane_cmp_df.merge(cds_chains,on="tid",how="left")
    mane_cmp_df = mane_cmp_df.merge(cds_chains[["tid","chain"]].rename({"chain":"chain_mane"},axis=1),left_on="tid_mane",right_on="tid",how="left")
    mane_cmp_df.drop("tid_y",axis=1,inplace=True)
    
    mane_cmp_df[["mod_chain","c1len","c2len","match_start","match_end","num_bp_extra","num_bp_missing","num_bp_inframe","num_bp_match","num_bp_outframe","lpi","ilpi","mlpi"]] = mane_cmp_df.apply(lambda row: compare_and_extract(row["chain"],row["chain_mane"],row["strand"]),axis=1)
    mane_cmp_df.rename({"tid_x":"tid"},axis=1,inplace=True)
    c3 = c3.merge(mane_cmp_df[["tid","tid_mane","ilpi"]],on="tid",how="left")
    
    # iterate over duplicates and decide what to do based on the comparisons
    # keep group information (all tids so we can add tag with duplicates)

    keep_tids = {}

    for (gid,echain), grp in c3.groupby(by=['gene_id','echain']):
        if len(set(grp["tid"].tolist()))==1:
            keep_tids[grp["tid"].tolist()[0]] = []
        
        # find one with highest ilpi to mane
        max_ilpi = 0
        max_ilpi_tid = None

        max_mom = 0
        max_mom_tid = None
        for idx, row in grp.iterrows():
            ilpi = 0
            try:
                ilpi = float(row["ilpi"])
            except:
                ilpi = 0
            if ilpi>max_ilpi:
                max_ilpi = ilpi
                max_ilpi_tid = row["tid"]

            mom = 0
            try:
                mom = float(row["mom"])
            except:
                mom = 0
            if mom>max_mom:
                max_mom = mom
                max_mom_tid = row["tid"]

        if max_ilpi > 0: # found match
            keep_tids[max_ilpi_tid] = []
            for idx, row in grp.iterrows():
                keep_tids[max_ilpi_tid].append(row["tid"])
            # continue to next group
            continue

        if max_mom > 0: # found match
            keep_tids[max_mom_tid] = []
            for idx, row in grp.iterrows():
                keep_tids[max_mom_tid].append(row["tid"])
            # continue to next group
            continue

        # for everything else - pick at random
        rnd_tid = random.choice(grp["tid"].tolist())
        keep_tids[rnd_tid] = []
        for idx, row in grp.iterrows():
            keep_tids[rnd_tid].append(row["tid"])
            
        with open(out_gtf_fname,"w+") as outFP:
            with open(in_gtf_fname,"r") as inFP:
                for line in inFP:
                    lcs = line.split("\t")
                    if lcs[0]=="#":
                        outFP.write(line)
                    if not len(lcs)==9:
                        outFP.write(line)
                        
                    tid = lcs[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                    if tid in keep_tids:
                        attrs = extract_attributes(lcs[8])
                        if len(keep_tids[tid])>1:
                            attrs["duplicates"] = ",".join(keep_tids[tid])
                        res_line = "\t".join(lcs[:-1]) + "\t" + to_attribute_string(attrs, True, lcs[2])
                        outFP.write(res_line + "\n")
                        
def build_with_genes_gtf(in_gtf_fname, out_gtf_fname, gene_desc_dict):
    # load a list of genes
    genes = dict()

    with open(in_gtf_fname, "r") as inFP:
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
    with open(out_gtf_fname, "w+") as outFP:
        outFP.write(get_comments(args.release_number))
        with open(in_gtf_fname, "r") as inFP:
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
                
def build_conversions(gtf_file,out_base_fname,contig_lengths,gene_description_reference,bed_as):
    # check that the input is in the GTF format
    assert os.path.exists(gtf_file), "Input file not found"
    assert gtf_or_gff(gtf_file) == 'gtf', "Input file is not in GTF format"

    # load gene descriptions from the reference file
    assert os.path.exists(gene_description_reference), "Gene description reference file not found"
    gene_desc_dict = load_gene_descriptions(gene_description_reference)

    # assert that contigs lengths file exists
    assert os.path.exists(contig_lengths), "Contig lengths file not found"

    with_genes_gtf_fname = out_base_fname + "1.with_genes.gtf"
    build_with_genes_gtf(gtf_file, with_genes_gtf_fname, gene_desc_dict)

    # convert the GTF file to GFF
    with_genes_gff_fname = out_base_fname + "1.with_genes.gff"
    gtf2gff(with_genes_gtf_fname, with_genes_gff_fname)

    os.environ['LD_LIBRARY_PATH'] = "/ccb/sw/lib/:LD_LIBRARY_PATH"

    genePred_fname = out_base_fname + "2.with_genes.genePred"
    cmd = ["gff3ToGenePred", with_genes_gff_fname, genePred_fname]
    subprocess.call(cmd)

    bedPlus_fname = out_base_fname + "2.with_genes.bedPlus"
    cmd = ["genePredToBigGenePred", genePred_fname, bedPlus_fname]
    subprocess.call(cmd)

    sorted_bedPlus_fname = out_base_fname + "2.with_genes.srt.bedPlus"
    cmd = ["bedSort", bedPlus_fname, sorted_bedPlus_fname]
    subprocess.call(cmd)

    os.rename(sorted_bedPlus_fname, bedPlus_fname)

    bb_fname = out_base_fname + "2.with_genes.bb"
    cmd = ["bedToBigBed", "-type=bed12+8", "-tab", "-as="+bed_as, bedPlus_fname,
           contig_lengths, bb_fname, "-extraIndex=name"]
    subprocess.call(cmd)

    shutil.copy(gtf_file, out_base_fname+".gtf")
    os.rename(with_genes_gff_fname, out_base_fname+".gff")
    os.rename(bb_fname, out_base_fname+".bb")
                
def run(args):
    # setup output directory
    outdir = os.path.abspath(args.outdir)+"/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    rdir = os.path.join(args.outdir,args.prefix+args.release_number+"."+args.suffix)
    if not os.path.exists(rdir):
        os.makedirs(rdir)
    
    base_name = rdir+"/"+args.prefix+args.release_number+"."+args.suffix

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
    parser.add_argument("--borf", help="Path to the borf executable")
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)
    
    # TODO:
    # 1. Compute summaries (# tx, cds, gid, etc)
if __name__=="__main__":
    main(sys.argv[1:])
