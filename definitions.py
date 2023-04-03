# contains reusable funcitons used throughout the experiments


# 1s cds chains of two transcripts and returns a set of statistics:
# 1. number match
#    - inframe
#    - outframe
# 2. number mismatch
# 3. match start
# 4. match stop

import os
import random
import pyBigWig
import upsetplot
import subprocess
import numpy as np
import pandas as pd

gff3cols=["seqid","source","type","start","end","score","strand","phase","attributes"]

def load_fasta_dict(fa_fname,rev=False):
    res = dict()
    with open(fa_fname,"r") as inFP:
        cur_nm = None
        
        for line in inFP:
            if line[0]==">":
                cur_nm = line.strip()[1:].split()[0]
                assert cur_nm not in res,"duplicate record name: "+nm
                res[cur_nm] = ""
            else:
                assert cur_nm is not None,"empty record name"
                res[cur_nm]+=line.strip()
                
    if rev:
        im = dict()
        for k, v in res.items():
            im[v] = im.get(v, []) + [k]
        
        res = im
    return res

def get_gff3cols():
    return gff3cols

def test_defs(): # simple test to check whether the file has been loaded
    print("test passed")

def intersect(s1,s2):
    res = [0,-1,0]
    tis = max(s1[0],s2[0])
    tie = min(s1[1],s2[1])
    if(tis<=tie):
        res[0] = tis
        res[1] = tie
        return (tie-tis)+1,res
    return 0,res

def split(s1,s2):
    left  = [0,-1,-1]
    right = [0,-1,-1]
    
    il,inter = intersect(s1,s2)
    if il>0:
        if inter[0]>s1[0]:
            left[0] = s1[0]
            left[1] = inter[0]-1
            left[2] = s1[2]
        if inter[1]<s1[1]:
            right[0] = inter[1]+1
            right[1] = s1[1]
            right[2] = s1[2]
    else:
        if s1[0]<s2[0]:
            left = s1
        else:
            right = s1
        
    return left,inter,right

def slen(s):
    return (s[1]-s[0])+1

def clen(chain):
    res = 0
    for c in chain:
        res+=slen(c)
    return res

def compare(i1,i2):
    intervals = []
    for i in i1:
        intervals.append([i[0],i[1]])
        intervals[-1].append(-1)
    for i in i2:
        intervals.append([i[0],i[1]])
        intervals[-1].append(1)
    intervals.sort()
    
    if len(i1)==0 and len(i2)==0:
        return []
    
    stack = []
    stack.append(intervals[0])
    for i in intervals[1:]:
        
        left,inter,right = split(stack[-1],i)
        if slen(right)>0:
            assert slen(inter)==slen(i) # must be intirely contained within
        else:
            tmp,inter2,right = split(i,stack[-1])
            if(slen(tmp)>0):
                t2 = stack[-1]
                stack[-1]=tmp
                stack.append(t2)
                
            else:
                assert slen(tmp)<=0,str(tmp)+","+str(inter2)+","+str(right)
            assert inter==inter2
        
        stack.pop()
        
        if slen(left)>0:
            stack.append(left)
        if slen(inter)>0:
            inter[2] = 0
            stack.append(inter)
        if slen(right)>0:
            stack.append(right)
        
    return stack

def subset_gtf(in_gtf_fname,out_gtf_fname,tids):
    with open(out_gtf_fname,"w+") as outFP:
        with open(in_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue

                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid in tids:
                    outFP.write(line)
                    
def load_gtf(gtf_fname):
    assert os.path.exists(gtf_fname),"input file does not exist: "+gtf_fname
    res = pd.read_csv(gtf_fname,sep="\t",names=gff3cols,comment="#")
    return res

def extract_gtf_attribute(gtf_df,k,column="attributes"):
    assert column in set(gtf_df.columns), "please speify the correct column to extract attributes from"
    return gtf_df[column].str.split(k+" \"",expand=True)[1].str.split("\"",expand=True)[0]

def load_map(gtf_fname, qname, tname, pass_codes):
    assert os.path.exists(gtf_fname),"input file does not exist: "+gtf_fname
    
    res = dict()
    
    with open(gtf_fname,"r") as inFP:
        for line in inFP.readlines():
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            if not lcs[2]=="transcript":
                continue

            code = lcs[8].split("class_code \"",1)[1].split("\"",1)[0]
            if code not in pass_codes:
                continue

            qtid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
            ttid = lcs[8].split("cmp_ref \"", 1)[1].split("\"", 1)[0]

            assert qtid not in res,"duplicate qtid: "+qtid

            res[qtid] = [code,ttid,lcs[0],lcs[6]]
    res = pd.DataFrame.from_dict(res,orient="index").reset_index()
    res.columns = [qname,"code",tname,"seqid","strand"]
    return res

def get_attribute(gtf_fname,attrs):
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

            if lcs[2]=="transcript":
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                tid2name[tid]=[]
                for a in clean_attrs:
                    gn = "-"
                    if a in lcs[8]:
                        gn = lcs[8].split(a+" \"", 1)[1].split("\"", 1)[0]
                    tid2name[tid].append(gn)
    res = pd.DataFrame.from_dict(tid2name,orient="index").reset_index()
    tmp_cols = ["tid"]
    for a in clean_attrs:
        tmp_cols.append(a)
    res.columns = tmp_cols
    return res

def get_chains(gtf_fname,feature_type,coords,phase=False):
    res = dict()
    with open(gtf_fname,"r") as inFP:
        for line in inFP:
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            if lcs[2]=="transcript":
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]

                if(coords):
                    res.setdefault(tid,[0,lcs[0],lcs[6],lcs[0]+":"+str(lcs[3])+"-"+str(lcs[4]),[]])
                else:
                    res.setdefault(tid,[0,[]])

            if not lcs[2]==feature_type:
                continue
                
            tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]

            res[tid][-1].append([int(lcs[3]),int(lcs[4])])
            if phase:
                res[tid][-1][-1].append(lcs[7])
            res[tid][0]=1
            
    for k,v in res.items():
        res[k] = v[:-1]
        res[k].append(sorted(v[-1]))
    res = pd.DataFrame.from_dict(res,orient="index").reset_index()
    if(coords):
        res.columns = ["tid","has_cds","seqid","strand","coords","chain"]
    else:
        res.columns = ["tid","has_cds","chain"]
    return res

# runs compare() funciton and labels all matches as in and out of frame accordingly
def compare_label_frame(chain1,chain2,strand):
    if chain2 is np.nan or len(chain2)==0:
        [[x[0],x[1],-1] for x in chain1]
    if chain1 is np.nan or len(chain1)==0:
        [[x[0],x[1],1] for x in chain2]

    mod_chain = compare(chain1,chain2)

    if strand=="-":
        mod_chain.reverse()

    t_frame = 0
    q_frame = 0
    
    for mc in mod_chain:
        if (mc[2] == -1): # extra positions in the query
            q_frame += slen(mc)
        elif (mc[2] == 1): # template positions missing from the query
            t_frame += slen(mc)
        elif (mc[2] == 0): # matching positions between query and template
            if (q_frame % 3 == t_frame % 3):
                mc[2] = 100 # inframe
            else:
                mc[2] = -100 # outframe
        else:
            print("wrong code")
            return

    return mod_chain

def compare_and_extract(chain1,chain2,strand):
    if chain2 is np.nan or len(chain2)==0:
        return pd.Series([[[x[0],x[1],-1] for x in chain1],-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
    if chain1 is np.nan or len(chain1)==0:
        return pd.Series([[[x[0],x[1],1] for x in chain2],-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])

    # 1. compute the total number of matching positions between query and template
    # 2. compute the number of matching positions in frame between query and template
    mod_chain = compare(chain1,chain2)
    
    c1len = clen(chain1)
    c2len = clen(chain2)

    if strand=="-":
        mod_chain.reverse()

    num_bp_extra = 0
    num_bp_missing = 0
    num_bp_inframe = 0
    num_bp_match = 0
    num_bp_outframe = 0

    t_frame = 0
    q_frame = 0
    
    for mc in mod_chain:
        if (mc[2] == -1): # extra positions in the query
            num_bp_extra += slen(mc)
            q_frame += slen(mc)
        elif (mc[2] == 1): # template positions missing from the query
            num_bp_missing += slen(mc)
            t_frame += slen(mc)
        elif (mc[2] == 0): # matching positions between query and template
            num_bp_match += slen(mc)
            if (q_frame % 3 == t_frame % 3):
                num_bp_inframe += slen(mc) # TODO: shouldn't this be stranded?
            else:
                num_bp_outframe += slen(mc)
        else:
            print("wrong code")
            return

    # compute lpd, ilpd, mlpd, etc
    lpd = int((100.0 * (float(c1len) / float(c2len))))
    ilpd = int((100.0 * (float(num_bp_inframe) / float(c2len))))
    mlpd = int((100.0 * (float(num_bp_match) / float(c2len))))

    match_start = chain1[0][0]==chain2[0][0] if strand=='+' else chain1[-1][1]==chain2[-1][1]
    match_end = chain1[-1][1]==chain2[-1][1] if strand=='+' else chain1[0][0]==chain2[0][0]

    return pd.Series([mod_chain,c1len,c2len,match_start,match_end,num_bp_extra,num_bp_missing,num_bp_inframe,num_bp_match,num_bp_outframe,lpd,ilpd,mlpd])

def load_tid2aa(fname):
    tid2aa = dict()

    with open(fname,"r") as inFP:
        cur_tid = ""
        cur_aa = ""
        for line in inFP:
            if line[0]==">":

                if not len(cur_tid)==0:
                    tid2aa[cur_tid]=cur_aa

                cur_tid = line[1:].rstrip()
                cur_aa = ""
            else:
                cur_aa += line.rstrip()

        if not len(cur_tid)==0:
            tid2aa[cur_tid]=cur_aa

    res = pd.DataFrame.from_dict(tid2aa,orient="index").reset_index()
    res.columns = ["tid","aa"]
    return res

def merge(segs):
    segs.sort()
    res = [[segs[0][0],segs[0][1]]]
    for s in segs:
        prev = res[-1]
        if s[0] <= prev[1]:
            prev[1] = max(prev[1], s[1])
        else:
            res.append([s[0],s[1]])
            
    return res

def load_segments(fname,feature_type,strandless):
    
    res = dict({"+":dict(),
                "-":dict()})
    if strandless:
        res = dict()
    
    with open(fname,"r") as inFP:
        for line in inFP:
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            if not lcs[2]==feature_type:
                continue
                
            if strandless:
                res.setdefault(lcs[0],set())
                res[lcs[0]].add((int(lcs[3]),int(lcs[4])))
            else:
                res[lcs[6]].setdefault(lcs[0],set())
                res[lcs[6]][lcs[0]].add((int(lcs[3]),int(lcs[4])))
            
    for k,v in res.items():
        if strandless:
            res[k] = merge(list(v))
        else:
            for k2,v2 in v.items():
                res[k][k2] = merge(list(v2))
            
    return res

def extract_from_comp(segs): # separated "shared,left,right" into separate objects
    left=[]
    shared=[]
    right=[]
    for s in segs:
        if s[2]==-1:
            left.append(s[:2])
        if s[2]==1:
            right.append(s[:2])
        if s[2]==0:
            shared.append(s[:2])
            
    return left,shared,right

# extract three-way comparisons
def upset_data(segs):
    res = {111:dict(),
           100:dict(),
             1:dict(),
            10:dict(),
           101:dict(),
           110:dict(),
            11:dict()}
    
    for seqid in segs[0]:
        for k,v in res.items():
            v.setdefault(seqid,list())
            
        s0 = segs[0][seqid]
        s1 = segs[1][seqid]
        s2 = segs[2][seqid]
        
        s0s1 = compare(s0,s1)
        s2s1 = compare(s2,s1)
        s0s2 = compare(s0,s2)
        
        s0_not_s1,s0_and_s1,s1_not_s0 = extract_from_comp(s0s1)
        s2_not_s1,s2_and_s1,s1_not_s2 = extract_from_comp(s2s1)
        s0_not_s2,s0_and_s2,s2_not_s0 = extract_from_comp(s0s2)
        
        # shared between all three
        tmp = compare(s0_and_s1,s2)
        x,s0_and_s1_and_s2,y = extract_from_comp(tmp)
        # in s0 only
        tmp = compare(s0_not_s1,s2)
        s0_not_s1_not_s2,x,y = extract_from_comp(tmp)
        # in s2 only
        tmp = compare(s2_not_s1,s0)
        s2_not_s1_not_s0,x,y = extract_from_comp(tmp)
        # in s1 only
        tmp = compare(s1_not_s0,s2)
        s1_not_s0_not_s2,x,y = extract_from_comp(tmp)
        # s0 and s2
        tmp = compare(s0_and_s2,s1)
        s0_and_s2_not_s1,x,y = extract_from_comp(tmp)
        # s0 and s1
        tmp = compare(s0_and_s1,s2)
        s0_and_s1_not_s2,x,y = extract_from_comp(tmp)
        # s1 and s2
        tmp = compare(s2_and_s1,s0)
        s2_and_s1_not_s0,x,y = extract_from_comp(tmp)

        res[111][seqid] = s0_and_s1_and_s2
        res[100][seqid] = s0_not_s1_not_s2
        res[1  ][seqid] = s2_not_s1_not_s0
        res[10 ][seqid] = s1_not_s0_not_s2
        res[101][seqid] = s0_and_s2_not_s1
        res[110][seqid] = s0_and_s1_not_s2
        res[11 ][seqid] = s2_and_s1_not_s0
        
    return res

# extract scores
def get_scores(scores_fname,ud,num_random=None):
    res = dict()
    min_score = 0
    max_score = 0
    
    bw = pyBigWig.open(scores_fname)   
    
    for k,v in ud.items():
        res.setdefault(k,[])
        for seqid,v2 in v.items():
            if not seqid in bw.chroms():
                continue
            for r in v2:
                res[k].extend(bw.values(seqid,r[0],r[1]+1))
        if len(res[k])>0:
            min_score = min(min_score,min(res[k]))
            max_score = max(max_score,max(res[k]))
        
    if num_random is not None:
        for k,v in res.items():
            if len(v)>num_random:
                res[k] = random.sample(v,num_random)
                
    return res,min_score,max_score

# get position count matrix
def pos_mat(ud,labels,scores=None):
    res =  upsetplot.from_memberships( [[],
                                        [labels[0]],
                                        [labels[1]],
                                        [labels[1], labels[0]],
                                        [labels[2]],
                                        [labels[2], labels[0]],
                                        [labels[2], labels[1]],
                                        [labels[2], labels[1], labels[0]],],data=[0,
                                                                                sum([clen(v) for k,v in ud[100].items()]),
                                                                                sum([clen(v) for k,v in ud[10].items()]),
                                                                                sum([clen(v) for k,v in ud[110].items()]),
                                                                                sum([clen(v) for k,v in ud[1].items()]),
                                                                                sum([clen(v) for k,v in ud[101].items()]),
                                                                                sum([clen(v) for k,v in ud[11].items()]),
                                                                                sum([clen(v) for k,v in ud[111].items()])])
    
    res_scores = [(sum([clen(v) for k,v in ud[100].items()]),[]),
                (sum([clen(v) for k,v in ud[10].items()]),[]),
                (sum([clen(v) for k,v in ud[110].items()]),[]),
                (sum([clen(v) for k,v in ud[1].items()]),[]),
                (sum([clen(v) for k,v in ud[101].items()]),[]),
                (sum([clen(v) for k,v in ud[11].items()]),[]),
                (sum([clen(v) for k,v in ud[111].items()]),[])]
    if scores is not None:
        res_scores = [(sum([clen(v) for k,v in ud[100].items()]),scores[100]),
                    (sum([clen(v) for k,v in ud[10].items()]),scores[10]),
                    (sum([clen(v) for k,v in ud[110].items()]),scores[110]),
                    (sum([clen(v) for k,v in ud[1].items()]),scores[1]),
                    (sum([clen(v) for k,v in ud[101].items()]),scores[101]),
                    (sum([clen(v) for k,v in ud[11].items()]),scores[11]),
                    (sum([clen(v) for k,v in ud[111].items()]),scores[111])]
    
    return [res,res_scores]

# get_mean_score
def mean_score(score_fname):
    means = []
    bw = pyBigWig.open(score_fname)
    for seqid,l in bw.chroms().items():
        means.append(bw.stats(seqid,0,l-1))
    return np.mean(means)

# extract sashimi and gtf for a specified set of transcripts based on several annotations
def extract_sashimi(sbin,cmp_gtf_fname,ref_gtf_fname,q_gtf_fname,out_base_fname,cmp_tid,qtid,title_str):
    out_gtf_fname = out_base_fname+".gtf"
    out_svg_fname = out_base_fname+".svg"
    
    with open(out_gtf_fname,"w+") as outFP:
        # first write out MANE
        with open(cmp_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid == cmp_tid:
                    outFP.write(line)
        # next write out ORFanage
        with open(q_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid == qtid:
                    lcs[8] = "transcript_id \"ORFanage:"+tid+"\""
                    line = "\t".join(lcs)+"\n"
                    outFP.write(line)
        # lastly write out regular RefSeq
        with open(ref_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid == qtid:
                    outFP.write(line)
    
    sashimi_cmd = [sbin,
                   "--compare",cmp_tid,
                   "--title",title_str,
                   "--gtf",out_gtf_fname,
                   "-o",out_svg_fname]
    print(" ".join(sashimi_cmd))
    subprocess.call(sashimi_cmd)
    
def get_poly_gids(gtf_fname):
    df = get_chains(gtf_fname,"CDS",True)
    df = df[df["has_cds"]==1].reset_index(drop=True)
    df["seqid"]=df["coords"].str.split(":",n=1,expand=True)[0]
    df["start"] = df["chain"].apply(lambda row: row[0][0])
    df["end"] = df["chain"].apply(lambda row: row[-1][1])
    # add gene ids
    gid=pd.read_csv(gtf_fname,sep="\t",names=gff3cols,comment="#")
    gid=gid[gid["type"]=="transcript"].reset_index(drop=True)
    gid["tid"]=gid["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    gid["gid"]=gid["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    gid = gid[["gid","tid"]]

    df = df.merge(gid,on="tid",how="left",indicator=False)

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    df.sort_values(by=["seqid","strand","start","end"],ascending=True,inplace=True)

    df = df.groupby(by=["seqid","strand","gid"]).agg({"start":min,"end":max}).reset_index()
    df.sort_values(by=["seqid","strand","start","end"],ascending=True,inplace=True)
    df["nc"]=df.seqid.shift(-1)
    df["nt"]=df.strand.shift(-1)
    df["ns"]=df.start.shift(-1)
    df["nid"]=df.gid.shift(-1)
    df.fillna(0,inplace=True)
    df["od"] = np.where((df["seqid"]==df["nc"]) & 
                               (df["strand"]==df["nt"]) & 
                               (df["end"]>df["ns"]),1,0)
    pids = set(df[df["od"]==1]["gid"]).union(set(df[df["od"]==1]["nid"]))
    
    return pids

def clean(in_gtf_fname,out_gtf_fname,gffread_bin,fa_fname,remove_single_exon,keep_seqids = None):
    # deal with stops
    tmp1_fname = out_gtf_fname+".tmp1.gtf"
    cmd = [gffread_bin,
           "-g",fa_fname,
           "--adj-stop","-V","-T","-F","-J",
           "-o",tmp1_fname,
           in_gtf_fname]

    print(" ".join(cmd))
    subprocess.call(cmd)
    
    seids = set()
    if remove_single_exon:
        df=pd.read_csv(tmp1_fname,sep="\t",names=gff3cols,comment="#")
        df=df[df["type"]=="exon"].reset_index(drop=True)
        df["tid"]=df["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
        df["gid"]=df["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
        gdf = df[["tid","seqid"]].groupby(by="tid").count().reset_index()
        gdf.columns = ["tid","num_exons"]
        seids = set(gdf[gdf["num_exons"]==1]["tid"])
        seids = set(df[df["tid"].isin(seids)]["gid"])
        print("number of single exon genes discarded: "+str(len(seids)))
        
        
    # deal with poly
    pids = get_poly_gids(tmp1_fname)
    
    # deal with exceptions
    df=pd.read_csv(tmp1_fname,sep="\t",names=gff3cols,comment="#")
    df=df[df["type"]=="transcript"].reset_index(drop=True)
    df["tid"]=df["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    df["gid"]=df["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]

    # retain chromosomes
    cids = set()
    if not keep_seqids is None:
        cids = df[~(df["seqid"].isin(keep_seqids))]["gid"]
        print("number of transcripts not on selected seqids: "+str(len(cids)))
    
    try:
        df["attr"]=df["attributes"].str.split("exception \"",expand=True)[1].str.split("\"",expand=True)[0]
    except:
        df["attr"]=None
    df["seleno"] = df["attributes"].str.lower().str.contains("selen")

    sids = set(df[df["seleno"]]["gid"])
    print("number of seleno: "+str(len(sids)))

    print("exceptions: "+", ".join(list(set(df[~(df["attr"].isna())]["attr"].tolist()))))
    eids = set(df[~(df["attr"].isna())]["gid"])
    print("number of exceptions: "+str(len(eids)))
    
    dirty_gids = seids.union(pids).union(sids).union(eids).union(cids)
    print("number of genes to discard: "+str(len(dirty_gids)))
    
    with open(out_gtf_fname,"w+") as outFP:
        with open(tmp1_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                gid = lcs[8].split("gene_id \"",1)[1].split("\"",1)[0]
                if not gid in dirty_gids:
                    outFP.write(line)
                    
def get_coding_gids(gtf_fname):
    df=pd.read_csv(gtf_fname,sep="\t",names=gff3cols,comment="#")
    df = df[df["type"]=="CDS"].reset_index(drop=True)
    df["gid"] = df["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    return set(df["gid"])

def subset_gtf_by_seqid(in_gtf_fname,out_gtf_fname,seqids):
    with open(out_gtf_fname,"w+") as outFP:
        with open(in_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.rstrip().split("\t")
                if not seqids == False and lcs[0] not in seqids:
                    continue
                    
                outFP.write(line)
    

def subset_gtf(in_gtf_fname,out_gtf_fname,gids,tids):
    writing_tid = ""
    with open(out_gtf_fname,"w+") as outFP:
        with open(in_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.rstrip().split("\t")
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if lcs[2]=="transcript":
                    gid = lcs[8].split("gene_id \"", 1)[1].split("\"", 1)[0]
                    if not gids == False and gid in gids:
                        outFP.write(line)
                        writing_tid = tid
                        continue
                
                if not tids == False and tid in tids:
                    outFP.write(line)
                    continue
                    
                # handle non transcript ffeatures without gene_id for whihc a transcript was found based on gene-id
                if writing_tid == tid:
                    outFP.write(line)
                    continue
                
                
def get_coding_stats(gtf_fname):
    # average number of isoforms per gene

    cgids = get_coding_gids(gtf_fname)
    t2g = get_attribute(gtf_fname,"gene_id")
    t2g.columns = ["tid","gid"]
    chains = get_chains(gtf_fname,"CDS",False)
    chains = chains[chains["has_cds"]==1].reset_index(drop=True)
    chains["chain_str"] = chains.apply(lambda row: ",".join(str(x[0])+"-"+str(x[1]) for x in row["chain"]),axis=1)

    t2g_c = t2g[t2g["gid"].isin(cgids)].reset_index(drop=True)
    t2g_c_g = t2g_c.groupby(by="gid").count().reset_index()
    print("mean number of transcripts per coding gene: "+str(t2g_c_g["tid"].mean()))

    chains = chains.merge(t2g_c,on="tid",how="inner")
    chains = chains[["chain_str","gid"]]
    chains.drop_duplicates(inplace=True)
    chains_g = chains.groupby(by="gid").count().reset_index()
    print("mean number of unique proteins per coding gene: "+str(chains_g["chain_str"].mean()))