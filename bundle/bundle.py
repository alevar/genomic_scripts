#!/usr/bin/env python3

import os
import sys
import copy
import argparse
import numpy as np


class Node():
    def __init__(self, index, weight):
        self.index = index
        self.weight = weight
        self.children = list()
        self.parents = list()
        self.top = True

    def add_child(self, child_idx):
        self.children.append(child_idx)

    def add_parent(self, parent_idx):
        self.parents.append(parent_idx)

    def set_top(self, top):
        self.top = top

    def is_top(self):
        return self.top == True

    def get_weight(self):
        return self.weight

    def get_idx(self):
        return self.index

    def get_parent(self):
        assert len(self.parents) == 1, "wrong parents"
        return self.parents[0]

    def get_children(self):
        return self.children

    def is_orphan(self):  # checks if the transcript is orphan (no neighbors/parents/children)
        if len(self.parents) == 0 and len(self.children) == 0 and self.top:
            return True
        return False

    def is_bottom(self):
        if len(self.children) == 0:
            return True
        return False

    def remove_child(self, child_idx):
        self.children.remove(child_idx)  # remove by value

    def remove_parent(self, parent_idx):
        self.parents.remove(parent_idx)

    def __str__(self):
        return "node: " + str(self.index) + \
               " top: " + str(self.top) + \
               " weight: " + str(self.weight) + \
               " children: " + ",".join([str(x) for x in self.children]) + \
               " parents: " + ",".join([str(x) for x in self.parents])


class Tree():
    def __init__(self):
        self.nodes = list()

    def add_node(self, node):
        self.nodes.append(node)

    def common_top_node(self, idx1, idx2):
        idx1_trace = list()
        idx2_trace = list()

        top_node = self.nodes[idx1]
        if not top_node.is_top():
            while True:
                parent_idx = top_node.get_parent()
                top_node = self.nodes[parent_idx]
                idx1_trace.append(parent_idx)
                if top_node.is_top():
                    break

        top_node = self.nodes[idx2]
        if not top_node.is_top():
            while True:
                parent_idx = top_node.get_parent()
                top_node = self.nodes[parent_idx]
                idx2_trace.append(parent_idx)
                if top_node.is_top():
                    break

        commons = set(idx1_trace).intersection(set(idx2_trace))
        if len(commons) == 0:
            return None
        else:
            common_top = list()
            for n in commons:
                if self.nodes[n].is_top():
                    common_top.append(n)
            assert len(common_top) == 1, "wrong number of top nodes: " + str(len(common_top))
            return common_top[0]

    def connect_nodes(self, idx1, idx2, weight):
        self.nodes.append(Node(len(self.nodes), weight))
        self.nodes[-1].set_top(True)

        common_node_idx = self.common_top_node(idx1, idx2)
        if common_node_idx is not None:
            common_node = self.nodes[common_node_idx]
            self.nodes[-1].add_child(common_node.get_idx())
            self.nodes[common_node.get_idx()].add_parent(self.nodes[-1].get_idx())
            self.nodes[common_node.get_idx()].set_top(False)
        else:
            if self.nodes[idx1].is_top():
                self.nodes[-1].add_child(idx1)
                self.nodes[idx1].add_parent(self.nodes[-1].get_idx())
                self.nodes[idx1].set_top(False)
            else:
                # has a parent - need to find for linking
                top_node = self.nodes[idx1]
                while True:
                    top_node = self.nodes[top_node.get_parent()]
                    if top_node.is_top():
                        break

                self.nodes[-1].add_child(top_node.get_idx())
                self.nodes[top_node.get_idx()].add_parent(self.nodes[-1].get_idx())
                self.nodes[top_node.get_idx()].set_top(False)

            if self.nodes[idx2].is_top():
                self.nodes[-1].add_child(idx2)
                self.nodes[idx2].add_parent(self.nodes[-1].get_idx())
                self.nodes[idx2].set_top(False)
            else:
                # has a parent - need to find for linking
                top_node = self.nodes[idx2]
                while True:
                    top_node = self.nodes[top_node.get_parent()]
                    if top_node.is_top():
                        break

                self.nodes[-1].add_child(top_node.get_idx())
                self.nodes[top_node.get_idx()].add_parent(self.nodes[-1].get_idx())
                self.nodes[top_node.get_idx()].set_top(False)

        # what if nodes are already joined? we need to introduce an intermediate node

    def get_top_idxs(self):
        idxs = list()
        for n in self.nodes:
            if n.is_top():
                idxs.append(n.get_idx())
        return idxs

    def is_orphan(self, idx):  # checks if the node is an orphan - transcript without a cluster/neighbors
        return self.nodes[idx].is_orphan()

    def is_bottom(self, idx):  # checks if the node has no children
        return self.nodes[idx].is_bottom()

    def top_down_score(self, idx):
        res = []
        for c in self.nodes[idx].get_children():
            res += self.top_down_score(c)
        if self.is_orphan(idx):
            res += [self.nodes[idx].get_weight()]
        elif self.is_bottom(
                idx):  # if is bottom and not orphan - we should not include the default value in the calculations
            pass
        else:
            res += [self.nodes[idx].get_weight()]
        return res

    def get_score(
            self):  # needs to combine all weights not just the top... (top-down recursive add) to compute an average of averages (mean within each top-down path)
        clu_scores = list()
        for i in self.get_top_idxs():
            tmp_scores = self.top_down_score(i)
            clu_scores.append(np.mean(tmp_scores))
        return np.mean(clu_scores)

    def get_num_clus(self):  # count the number of clusters
        return

    def get_bottom_idxs(self, idx):  # return a list of all lowest-level/bottom ids
        res = []
        for c in self.nodes[idx].get_children():
            res += self.get_bottom_idxs(c)
        if self.is_bottom(
                idx):  # if is bottom and not orphan - we should not include the default value in the calculations
            res += [idx]
        return res

    def get_clus(self):  # return a dictionary of clusters and transcripts
        clus = dict()
        for i in self.get_top_idxs():
            clus[i] = self.get_bottom_idxs(i)
        return clus


class Transcript():
    def __init__(self, seqid, strand, start, end, attrs):
        start = int(start)
        end = int(end)
        assert start < end, "Error (Coordinates): start>=end: "+attrs
        assert strand == "+" or strand == "-", "Error (Coordinates): unknown strand"
        self.seqid = seqid
        self.strand = strand
        self.start = int(start)
        self.end = int(end)
        self.attrs = dict()
        for pair in attrs.rstrip("\n").rstrip(";").split("\";"):
            k, v = pair.split(" \"")
            v = v.rstrip("\"")
            if k == "transcript_id":
                self.tid = v
            self.attrs[k] = v
        assert "transcript_id" in self.attrs, "Error (Attributes): no 'transcript_id' key present in the attributes"

        self.exons = list()
        self.introns = list()
        self.donors = list()
        self.acceptors = list()
        self.num_exons = 0
        self.num_introns = 0
        self.elen = 0  # effective length

    def get_tid(self):
        return self.tid

    def __str__(self):
        res = self.seqid + ":" + str(self.start) + "-" + str(self.end) + "\t" + self.strand + "\t"
        for k, v in self.attrs.items():
            res += k + " \"" + v + "\";"

        res += "num_exons \"" + str(self.num_exons) + "\";"
        res += "num_introns \"" + str(self.num_introns) + "\";"
        res += "elen \"" + str(self.elen) + "\";"
        return res

    def __repr__(self):
        res = self.seqid + ":" + str(self.start) + "-" + str(self.end) + "\t" + self.strand + "\t"
        for k, v in self.attrs.items():
            res += k + " \"" + v + "\";"

        res += "num_exons \"" + str(self.num_exons) + "\";"
        res += "num_introns \"" + str(self.num_introns) + "\";"
        res += "elen \"" + str(self.elen) + "\";"
        return res

    def add_exon(self, seqid, strand, start, end, attrs):
        start = int(start)
        end = int(end)
        assert seqid == self.seqid, "Error (Coordinates): incorrect seqid"
        assert strand == self.strand, "Error (Coordinates): incorrect strand"
        assert start >= self.start, "Error (Coordinates): start<tx_start"
        assert end <= self.end, "Error (Coordinates): end>tx_end: "+str(self.tid)
        for pair in attrs.rstrip("\n").rstrip(";").split("\";"):
            k, v = pair.split(" \"")
            v = v.rstrip("\"")
            if k == "transcript_id":
                assert v == self.tid, "Error (Attributes): transcript ids don't match"

        if self.num_exons > 0:
            assert start > self.exons[-1][0], "Error (Coordinates): exons not sorted"
            if self.strand == "+":
                self.donors.append(self.exons[-1][1])
                self.acceptors.append(start)
            elif self.strand == "-":
                self.donors.append(start)
                self.acceptors.append(self.exons[-1][1])
            else:
                raise Exception("Error (Coordinates): unknown strand")
        self.exons.append((start, end))

        exon_len = (end - start) + 1

        self.num_exons += 1
        self.num_introns = self.num_exons - 1
        self.elen += exon_len

    def get_seqid(self):
        return self.seqid

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_strand(self):
        return self.strand

    def get_acceptors(self):
        return self.acceptors

    def get_donors(self):
        return self.donors

    def get_exons(self):
        return self.exons

    def get_introns(self):
        return self.introns

    def normalize_bp(self, max_val):  # normalize the number of bps of overlap to a custom range where max_val==100%
        return 1

    # https://codereview.stackexchange.com/questions/178427/given-2-disjoint-sets-of-intervals-find-the-intersections
    @staticmethod
    def single_intersection(interval1, interval2):
        start = max(interval1[0], interval2[0])
        end = min(interval1[1], interval2[1])
        if start <= end:
            return (start, end), (end - start) + 1
        return None, 0

    @staticmethod
    def intersection(intervals1, intervals2):
        lengths = 0

        start = 0
        for interval1 in intervals1:
            found_yet = False
            for j in range(start, len(intervals2)):
                intersection, matches = Transcript.single_intersection(interval1, intervals2[j])
                if intersection:
                    lengths += matches
                    found_yet = True
                    start += 1
                elif found_yet:
                    break

        return lengths

    @staticmethod
    def union(intervals):
        res = []
        length = 0
        for s, e in sorted(intervals):
            if res and res[-1][1] > s - 1:
                res[-1][1] = max(res[-1][1], e)
            else:
                if len(res) > 0:
                    length += (res[-1][1] - res[-1][0]) + 1
                res.append([s, e])
        if len(res) > 0:
            length += (res[-1][1] - res[-1][0]) + 1
        return res, length

    def get_distance(self, tx):
        num_match = Transcript.intersection(self.exons, tx.get_exons())
        _, num_union = Transcript.union(self.exons + tx.get_exons())
        num_mismatch = num_union - num_match
        assert num_mismatch >= 0, "Error (Distance): wrong overlaps computed"

        # bases
        norm_perc_match = (num_match / num_union) * 20
        norm_perc_mismatch = (num_mismatch / num_union) * 20

        # # acceptors
        # num_intersection_acceptors = len(set(self.acceptors).intersection(set(tx.get_acceptors())))
        # num_union_acceptors = len(set(self.acceptors).union(set(tx.get_acceptors())))
        # norm_shared_acceptors=0
        # norm_nonshared_acceptors=0
        # if not num_union_acceptors == 0:
        #     norm_shared_acceptors = (num_intersection_acceptors / num_union_acceptors) * 4
        #     norm_nonshared_acceptors = ((num_union_acceptors-num_intersection_acceptors) / num_union_acceptors) * 4
        #
        # # donors
        # num_intersection_donors = len(set(self.donors).intersection(set(tx.get_donors())))
        # num_union_donors = len(set(self.donors).union(set(tx.get_donors())))
        # norm_shared_donors=0
        # norm_nonshared_donors=0
        # if not num_union_donors == 0:
        #     norm_shared_donors = (num_intersection_donors / num_union_donors) * 4
        #     norm_nonshared_donors = ((num_union_donors-num_intersection_donors) / num_union_donors) * 4
        #
        # # introns
        # num_intersection_introns = len(set(self.introns).intersection(set(tx.get_introns())))
        # num_union_introns = len(set(self.introns).union(set(tx.get_introns())))
        # norm_shared_introns=0
        # norm_nonshared_introns=0
        # if not num_union_introns == 0:
        #     norm_shared_introns = (num_intersection_introns / num_union_introns) * 6
        #     norm_nonshared_introns = ((num_union_introns-num_intersection_introns) / num_union_introns) * 6
        #
        # # exons
        # num_intersection_exons = len(set(self.exons).intersection(set(tx.get_exons())))
        # num_union_exons = len(set(self.exons).union(set(tx.get_exons())))
        # norm_shared_exons=0
        # norm_nonshared_exons=0
        # if not num_union_exons == 0:
        #     norm_shared_exons = (num_intersection_exons / num_union_exons) * 8
        #     norm_nonshared_exons = ((num_union_exons-num_intersection_exons) / num_union_exons) * 8
        #
        # return (norm_perc_match - norm_perc_mismatch) + norm_shared_acceptors + norm_shared_donors + norm_shared_exons + norm_shared_introns

        num_shared_acceptors = len(set(self.acceptors).intersection(set(tx.get_acceptors()))) * 4
        num_shared_donors = len(set(self.donors).intersection(set(tx.get_donors()))) * 4
        num_shared_exons = len(set(self.exons).intersection(set(tx.get_exons()))) * 10
        num_shared_introns = len(set(self.introns).intersection(set(tx.get_introns()))) * 7
        return (norm_perc_match - norm_perc_mismatch) + num_shared_acceptors + num_shared_donors + num_shared_exons + num_shared_introns

        # instead of using num_shared what if we also used normalization intersection/union?

        # need to investigate weird cases where non-overlapping things are groupped together


class Bundle():
    def __init__(self):
        self.txs = list()
        self.size = 0
        self.seqid = ""
        self.strand = ""
        self.start = np.inf  # minimum start coordinate of contained transcripts
        self.end = 0  # maximum end coordinate of contained transcripts

        # scores
        self.bp = 5
        self.exon = 10
        self.intron = 8
        self.donor = 4
        self.acceptor = 4

        # distances
        self.distmat = list()

    def add_tx(self, tx):
        assert self.can_add(tx), "Error (Bundle Coordinates): can not add tx to the bundle"
        self.txs.append(tx)
        self.seqid = tx.get_seqid()
        self.strand = tx.get_strand()
        self.start = min(self.start, tx.get_start())
        self.end = max(self.end, tx.get_end())
        self.size += 1

        # add empty entry to the distmat
        self.distmat.append([0 if x < self.size - 1 else 20 for x in range(self.size)])

    def empty(self):
        return self.size == 0

    def can_add(self, tx):
        if self.empty():
            return True
        else:
            if self.seqid == tx.get_seqid() and self.strand == tx.get_strand():
                if self.overlap(self.start, self.end, tx.get_start(), tx.get_end()) > 0:
                    return True
                else:
                    return False
            else:
                return False

    def overlap(self, start_x, end_x, start_y, end_y):
        return max(0, min(end_x, end_y) - max(start_x, start_y))

    def __str__(self):
        res = "BUNDLE\n"
        for tx in self.txs:
            res += str(tx) + "\n"
        return res

    def set_params(self, bp, exon, intron, donor, acceptor):
        self.bp = bp
        self.exon = exon
        self.intron = intron
        self.donor = donor
        self.acceptor = acceptor

    def dist_mat(self):
        # all vs all comparison to compute distances between transcripts in the bundle

        # upper triangular w/o diagonal
        for col in range(self.size - 1):
            for row in range(col + 1, self.size, 1):
                dist = self.txs[col].get_distance(self.txs[row])
                print(self.txs[col].tid,self.txs[row].tid,dist)
                self.distmat[row][col] = dist
        return

    def cluster(self):
        dist_np = self.to_np()
        if dist_np.shape[0]==1:
            clus = dict()
            clus[0] = [0]
            clus_tids = self.idxs2txs(clus)
            return clus_tids

        scores = []
        pairs = []
        for col in range(dist_np.shape[0] - 1):
            for row in range(col + 1, dist_np.shape[0], 1):
                scores.append(dist_np[col, row])
                pairs.append((col, row))

        zsp = zip(scores, pairs)
        szsp = sorted(zsp)[::-1]
        tuples = zip(*szsp)
        pair_scores, pairs = [list(x) for x in tuples]

        tree = Tree()
        for i in range(dist_np.shape[0]):
            tree.add_node(Node(i, 0))

        best_score = 0
        for i in range(len(pairs)):
            tree_copy = copy.deepcopy(tree)
            highest_pair_idx = pairs[i]
            highest_pair_score = pair_scores[i]
            # if highest_pair_score<=0:
            #     break
            tree_copy.connect_nodes(highest_pair_idx[0], highest_pair_idx[1], highest_pair_score)
            new_score = tree_copy.get_score()
            # tree = copy.deepcopy(tree_copy)
            # best_score = new_score
            if new_score >= best_score:
                # can re-write the tree
                tree = copy.deepcopy(tree_copy)
                best_score = new_score
        clus = tree.get_clus()
        clus_tids = self.idxs2txs(clus)
        return clus_tids

    def idxs2txs(self, clus):
        res = dict()
        for clu, txs in clus.items():
            tx_names = [self.txs[tx].get_tid() for tx in txs]
            res[clu] = tx_names
        return res

    def to_np(self):
        res = [[0 for i in range(self.size)] for j in range(self.size)]
        for col in range(self.size):
            for row in range(col, self.size, 1):
                res[row][col] = self.distmat[row][col]
                res[col][row] = self.distmat[row][col]
        return np.array([np.array(x) for x in res]).T

    def write(self, clus_tids, outFP):
        outFP.write(str(len(clus_tids))+"\t"+str(self.seqid)+":"+str(self.start)+"-"+str(self.end))
        for clu,txs in clus_tids.items():
            outFP.write("\t" + str(clu)+":"+",".join([str(tx) for tx in txs]))
        outFP.write("\n")


def run(args):
    with open(args.output,"w+") as outFP:
        with open(args.gtf, "r") as inFP:
            tx = None
            bundle = Bundle()
            for line in inFP.readlines():
                line = line.strip()
                lineCols = line.split("\t")
                if lineCols[2] == "transcript":
                    tx = Transcript(lineCols[0], lineCols[6], lineCols[3], lineCols[4], lineCols[8])
                    if bundle.can_add(tx):
                        bundle.add_tx(tx)
                    else:
                        bundle.dist_mat()
                        print(repr(bundle.to_np()))
                        clus_tids = bundle.cluster()
                        bundle.write(clus_tids,outFP)
                        bundle = Bundle()
                        bundle.add_tx(tx)
                if lineCols[2] == "exon":
                    tx.add_exon(lineCols[0], lineCols[6], lineCols[3], lineCols[4], lineCols[8])
            # process last bundle
            bundle.dist_mat()
            print(repr(bundle.to_np()))
            clus_tids = bundle.cluster()
            bundle.write(clus_tids,outFP)


def main(argv):
    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument('--gtf',
                        required=True,
                        type=str,
                        help="GTF file to analyze")
    parser.add_argument('--output',
                        required=True,
                        type=str,
                        help="output file name")
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
