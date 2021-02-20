#!/usr/bin/env python3

# ./subset_tree.py --tree data/nextstrain_ncov_global_timetree_170221.nwk --meta data/nextstrain_meta.csv -N 5 --output ./data/subtree.nwk --draw

import os
import sys
import random
import argparse
import numpy as np
import pandas as pd

from ete3 import Tree, TreeStyle,PhyloTree,faces, AttrFace, NodeStyle

def run(args):
    # load metadata
    meta_df = pd.read_csv(args.meta)
    meta_df.columns = ["seqname","grp","preserve"]
    preserved_genomes = set(meta_df[meta_df["preserve"]==1]["seqname"])

    # load tree
    t = Tree(args.tree,format=1)

    good_seqids = set(meta_df["seqname"]).intersection(set(t.get_leaf_names()))
    t.prune(good_seqids)
    meta_df = meta_df[meta_df["seqname"].isin(good_seqids)].reset_index(drop=True)


    # get groups of genomes by clade
    all_kept_genomes = set()

    clade_groups = meta_df[["seqname","grp"]].groupby(by="grp").apply(lambda grp: list(grp["seqname"]))
    for clade,strains in clade_groups.items():
        # find first genome of the clade in the tree
        strain_node = t.get_leaves_by_name(strains[0])
        assert len(strain_node)==1,"more than one node"
        strain_node = strain_node[0]

        # find the common ancestor for a set of genomes from the same clade
        clade_t = t.copy()
        clade_t.prune(strains)

        clade_mat,seq_order = clade_t.cophenetic_matrix()
        clade_mat = np.array([np.array(xi) for xi in clade_mat])

        np.fill_diagonal(clade_mat,1) # replace diagonal with 1 to search for argmin

        keep_genomes = seq_order
        srtarting_num_genomes = len(seq_order)

        while len(seq_order)>args.num_seqs:
            # get two two points with the closest distance between them
            min_dist = np.min(clade_mat)
            idxs = np.where(clade_mat==min_dist)

            # order 2d indices to remove duplicates (from lower and upper triangular forms)
            idxs = set(tuple(sorted(list(x))) for x in idxs)

            for i in idxs: # for  each occurrence of the minimum value
                # get a random one of the two selected genome to remove
                if args.seed is not None:
                    random.seed(args.seed)
                rm_seqidx = random.choice(i)

                # remove the genome from the matrix and the seq_order
                del keep_genomes[rm_seqidx]
                clade_mat = np.delete(clade_mat,rm_seqidx,0) # remove row
                clade_mat = np.delete(clade_mat,rm_seqidx,1) # remove column

        print("retained "+str(len(keep_genomes))+" out of "+str(srtarting_num_genomes)+" for group: "+str(clade))
        assert len(all_kept_genomes.intersection(set(keep_genomes)))==0,"duplicate genomes"
        all_kept_genomes = all_kept_genomes.union(keep_genomes)

    for seq in preserved_genomes:
        all_kept_genomes.add(seq)

    # final pruning to select the tree with only selected genomes
    subt = t.copy()
    subt.prune(list(all_kept_genomes))

    print("The final subsetted tree contains "+str(len(subt.get_leaf_names()))+" out of "+str(len(t.get_leaf_names())))

    # write to disk
    subt.sort_descendants()
    subt.write(features=[],outfile=args.output,format=1)
    if args.draw is True:
        os.environ['QT_QPA_PLATFORM']='offscreen'
        from matplotlib import cm
        from matplotlib.colors import to_hex

        viridis = cm.get_cmap('viridis',12)

        clades = dict([(i,x) for x,i in enumerate(sorted(list(set(meta_df["grp"]))))]) # assign IDs for clades
        meta_df["grpid"] = meta_df["grp"]
        meta_df = meta_df.replace({"grpid":clades})

        # normalize cluids and generate palette
        clade_colors = dict()
        for clu,cluid in clades.items():
            min_cluid = min(clades.values())
            max_cluid = max(clades.values())
            norm_cluid = (((0.9999-0.0001)*(cluid-min_cluid))/(max_cluid-min_cluid))+0.0001
            norm_color = to_hex(viridis(norm_cluid))
            clade_colors[cluid] = norm_color

        ts = TreeStyle()

        seen_clades = set()
        for node in subt.get_leaves():
            seqid = node.name
            # for everything else we just need to set the color of the cluster
            if node.is_leaf():
                nstyle = NodeStyle()
                cluid = meta_df[meta_df["seqname"]==seqid]["grpid"].iloc[0]
                c = clade_colors[cluid]
                nstyle["bgcolor"] = c
                node.set_style(nstyle)

                if not cluid in seen_clades: # add to legend if haven't seen the clade before
                    clu = meta_df[meta_df["grpid"]==cluid]["grp"].iloc[0]
                    ts.legend.add_face(faces.RectFace(height=10,width=10,fgcolor=c,bgcolor=c),column=0)
                    ts.legend.add_face(faces.TextFace(clu),column=1)
                    seen_clades.add(cluid)

        ts.show_leaf_name = False
        ts.allow_face_overlap=True

        ts.mode = "c"
        ts.root_opening_factor = 0.25
        subt.render(args.output+".png",tree_style=ts,w=600,h=1000)

    return

def main(args):
    global all_stages
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('--tree',
                        required=True,
                        type=str,
                        help="File containing a tree to be subsetted in Newick format")
    parser.add_argument('--meta',
                        required=True,
                        type=str,
                        help="CSV file containing metadata corresponding to the sequences in the tree. Metadata is expected to have 3 columns: 1. sequence name as it appears in the tree; 2. group specification; 3. Whether the sequence has to be preserved in the output (1/0). The data will be groupped by the 2nd column and subsetting will be performed for each group. Sequences to be preserved will be added in addition to the subsetted tree")
    parser.add_argument('--output',
                        required=True,
                        type=str,
                        help="Output directory in which all output and temporary data will be stored")
    parser.add_argument("-N",
                        "--num_seqs",
                        required=False,
                        type=int,
                        default=10,
                        help="Number of genomes to preserve in each group")
    parser.add_argument("--draw",
                        required=False,
                        action="store_true",
                        help="If enabled - will draw the resulting phylogenetic tree")
    parser.add_argument("--seed",
                        required=False,
                        type=int,
                        help="If set - will be used a seed for all pseudo-random choices")

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])


