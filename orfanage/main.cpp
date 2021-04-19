//
// Created by Ales Varabyou, Beril Erdogdu, Natalia Rincon on 03/30/21.
//

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>

#include <gclib/gff.h>
#include <fstream>
#include <gff_utils.h>
#include "arg_parse.h"

//#include <libBigWig/bigWig.h>
//#include <libBigWig/bwValues.c>
//#include <libBigWig/io.c>
//#include <libBigWig/bwWrite.c>
//#include <libBigWig/bwRead.c>

#define DBUF_LEN=1024;

//typedef std::vector<std::pair<uint,uint>> CHAIN_TYPE;
typedef std::vector<std::tuple<uint,uint,uint>> CDS_CHAIN_TYPE; // start,end,phase

struct Globals{
    std::ofstream out_stats_fp;
    std::ofstream out_cds_stats_fp;
    std::ofstream out_gtf_fp;
    std::ofstream out_gtf_perfect_fp;
    std::ofstream out_gtf_imperfect_fp;
    std::ofstream out_gtf_nonoverlap_fp;

    std::map<std::string,uint> fai; // seqid to length
    std::pair<std::map<std::string,uint>::iterator,bool> fai_it;
} globals;

std::string chain2str(CDS_CHAIN_TYPE& chain){
    std::string res = "";
    for(auto& c : chain){
        res+=std::to_string(std::get<0>(c))+"-"+std::to_string(std::get<1>(c))+",";
    }
    if(!res.empty()){
        res.pop_back();
    }
    if(res.empty()){
        res = "-";
    }
    return res;
}

int chain_len(CDS_CHAIN_TYPE& chain){
    int len=0;
    for(auto& c : chain) {
        len += (std::get<1>(c) + 1) - std::get<0>(c);
    }
    return len;
}

struct Mods{
    CDS_CHAIN_TYPE orig_chain;
    CDS_CHAIN_TYPE new_chain;
    bool missing_start = false;
    bool missing_end = false;

    CDS_CHAIN_TYPE missing;
    CDS_CHAIN_TYPE extra;

    int num_bp_extra = 0;
    int num_bp_missing = 0;
    int num_bp_match = 0;
    int num_bp_outframe = 0;
    int num_bp_inframe = 0;

    std::string orig_cds_tid = ""; // tid of the original CDS

    std::string cds_nt = ""; // nucleotide sequence
    std::string cds_aa = ""; // translated amino-acid sequence
    bool adjusted = false; // whether the chain has been adjusted to the first stop codon

    int get_score(){
        return missing_start*10+missing_end*10+num_bp_extra+num_bp_missing+num_bp_outframe;
    }
    std::string get_desc(){
        return "CDS_tid \""+this->orig_cds_tid+"\"; "+
               "missing_start \""+std::to_string(missing_start)+"\"; "+
               "missing_end \""+std::to_string(missing_end)+"\"; "+
               "num_bp_outframe \""+std::to_string(num_bp_outframe)+"\"; "+
               "num_bp_inframe \""+std::to_string(num_bp_inframe)+"\"; "+
               "num_bp_extra \""+std::to_string(num_bp_extra)+"\"; "+
               "num_bp_missing \""+std::to_string(num_bp_missing)+"\"; "+
               "num_bp_match \""+std::to_string(num_bp_match)+"\";";
    }
};

const Mods empty_mod; // placeholder for whenever a comparison is needed to an empty Mods object

// TODO: need start codon in addition to stop codon

//void run_annotate_with_tracks(std::string& seqid,char strand,std::tuple<uint,uint,uint>& coords,bigWigFile_t *bw_files[7],Globals& globals){
//
//    float transcript_sum = 0.0;
//    float transcript_power_sum = 0.0;
//    float transcript_all_power_sum = 0.0;
//    uint64_t transcript_power_count = 0;
//
////    // check whether chromosome exists in tracks
////    const auto chrom_sizes_it = params.chrom_sizes.find(seqid);
////    if (chrom_sizes_it == params.chrom_sizes.end())
////    {
////        t.phylo_score = NAN;
////        t.phylo_power = NAN;
////        // only report a missing sequence once
////        if (missing_sequences.find(t.chr) == missing_sequences.end())
////        {
////            missing_sequences.insert(t.chr);
////            printf(OUT_DEL OUT_INFO "Sequence %s from the GFF file does not occur in the tracks. Skipping ...\n" OUT_RESET, t.chr.c_str());
////        }
////    }
////    else
//    {
//        globals.fai_it.first = globals.fai.find(seqid);
//        uint64_t chr_len = 0;
//        if(globals.fai_it.first!=globals.fai.end()){
//            chr_len = globals.fai_it.first->second;
//        }
//
//        uint8_t wig_phase;
//        if (strand == '+')
//            wig_phase = 0 + (c.phase + c.begin - 1) % 3;
//        else
//            wig_phase = 3 + (chr_len - c.end - 1 + c.phase + 1) % 3;
//
//        float cds_sum = 0.0;
//
//        float cds_power_sum = 0.0;
//        float cds_all_power_sum = 0.0;
//        uint64_t cds_power_count = 0;
//
//        // cds_sum and cds_power_sum are for the weighted computation (and only consider positions where a score is available)
//        count_weighted_scores(cds_sum, cds_power_sum, cds_all_power_sum, cds_power_count, params.bw_files[wig_phase], params.bw_files[6], seqid.c_str(), c.begin - 1, c.end);
//        c.phylo_score = cds_sum / cds_power_sum;
//        // here we consider all positions, even if the confidence is below the threshold (when the tracks were computed)
//        c.phylo_power = (cds_power_count == 0) ? 0 : cds_all_power_sum / cds_power_count;
//
//        transcript_sum += cds_sum;
//        transcript_power_sum += cds_power_sum;
//        transcript_all_power_sum += cds_all_power_sum;
//        transcript_power_count += cds_power_count;
//    }
//}

void evaluate_fasta(Mods& res,const char* subseq,int bundle_start,char bundle_strand,uint start_offset,uint end_offset){
    int stop_coord = 0;
    std::string stop_nts = "";
    if(bundle_strand=='+'){
        for(auto& c : res.new_chain){
            uint cs = std::get<0>(c)-bundle_start; // start coordinate of the chain with respect to the bundle start
            uint ce = std::get<1>(c)-bundle_start; // end coordinate of the chain with respect to the bundle start
            for(uint i=cs;i<=ce;i++){
                res.cds_nt+=subseq[i];
            }
        }
    }
    else{ // strand=='-'
        for(int ci=res.new_chain.size()-1;ci>=0;ci--){
            uint cs = start_offset+(std::get<0>(res.new_chain[ci])-bundle_start);
            uint ce = start_offset+(std::get<1>(res.new_chain[ci])-bundle_start);
            for (uint i=ce;i>=cs;i--) {
                res.cds_nt+=ntComplement(subseq[i]);
            }
        }
    }
    int cds_len = 0;
    char* cdsaa = translateDNA(res.cds_nt.c_str(),cds_len,res.cds_nt.size());
    res.cds_aa = cdsaa;

    return;
    // TODO: also for each chain - output a separate stats file describing
    //     1. fraction of transcripts in the bundle it fits
    //     2. list of transcripts it fits

    // TODO: remove duplicate transcripts?
}

void adjust_stop(Mods& m,int bundle_start,char bundle_strand){ // adjust the chain coordinates to stop at the first stop codon
    if(m.cds_aa.empty()){
        std::cerr<<"amino chain is empty"<<std::endl;
        exit(2);
    }
    size_t stop_aa_pos = m.cds_aa.find('.');
    if(stop_aa_pos!=std::string::npos){ // stop codon was found
        size_t stop_nt_pos = stop_aa_pos*3;
        size_t left_to_stop = stop_nt_pos; // number of positions left until the stop codon
        // find CDS piece in which the first nt of the stop codon was found
        if(bundle_strand=='+'){
            for(int i=0;i<m.new_chain.size();i++){
                uint cs = std::get<0>(m.new_chain[i])-bundle_start; // start coordinate of the chain with respect to the bundle start
                uint ce = std::get<1>(m.new_chain[i])-bundle_start; // end coordinate of the chain with respect to the bundle start
                size_t clen = (ce+1)-cs;
                if(left_to_stop<clen){ // found the cds segment with the stop codon
                    std::get<1>(m.new_chain[i]) = (std::get<0>(m.new_chain[i])+left_to_stop)-1;
                    if(std::get<1>(m.new_chain[i])<std::get<0>(m.new_chain[i])){ // this case happens when the stop codon starts at the start of the 'i' exon
                        i--;
                    }
                    int chain_len = m.new_chain.size();
                    if(i<chain_len-1){
                        m.new_chain.erase(m.new_chain.begin()+i+1,m.new_chain.end());
                    }
                    m.cds_nt.erase(m.cds_nt.begin()+stop_nt_pos,m.cds_nt.end());
                    m.cds_aa.erase(m.cds_aa.begin()+stop_aa_pos,m.cds_aa.end());
                    m.adjusted = true;
                    break;
                }
                left_to_stop-=clen;
            }
        }
        else{ // strand=='-
            for(int i=m.new_chain.size()-1;i>=0;i--){
                uint cs = std::get<0>(m.new_chain[i])-bundle_start; // start coordinate of the chain with respect to the bundle start
                uint ce = std::get<1>(m.new_chain[i])-bundle_start; // end coordinate of the chain with respect to the bundle start
                size_t clen = (ce+1)-cs;
                if(left_to_stop<clen){ // found the cds segment with the stop codon
                    std::get<0>(m.new_chain[i]) = (std::get<1>(m.new_chain[i])-left_to_stop)+1;
                    if(std::get<1>(m.new_chain[i])<std::get<0>(m.new_chain[i])){ // this case happens when the stop codon starts at the start of the 'i' exon
                        i++;
                    }
                    if(i>0) {
                        if(i==0){
                            m.new_chain.erase(m.new_chain.begin());
                        }
                        else if(i==m.new_chain.size()){
                            m.new_chain.clear();
                        }
                        else{
                            m.new_chain.erase(m.new_chain.begin(),m.new_chain.begin()+i);
                        }
                    }
                    m.cds_nt.erase(m.cds_nt.begin()+stop_nt_pos,m.cds_nt.end());
                    m.cds_aa.erase(m.cds_aa.begin()+stop_aa_pos,m.cds_aa.end());
                    m.adjusted = true;
                    break;
                }
                left_to_stop-=clen;
            }
        }
    }
}

void adjust_start(Mods& m,int bundle_start,char bundle_strand){ // adjust the chain coordinates to start at the first start codon

}

void compare_aa(Mods& ref_res, Mods& new_res){
    // check if starts with start codon - if not - check the first one

    // if the original chain contains an in-frame stop codon at the same position as the fitted chain - do not re-adjust

    // check if ends in stop codon - if not - check the first occurrence
    // TODO: this may require examining the 3 bases after the CDS frame - since stop codons are typically not included in the CDS
    return;
}

struct TX{
    TX(GffObj* tx,uint idx):
                            cds_start(tx->hasCDS() ? tx->CDstart : 0),cds_end(tx->hasCDS() ? tx->CDend : 0){
        if(tx->hasCDS()){
            switch(tx->CDphase){
                case '0':
                    this->cds_phase=0;
                    break;
                case '1':
                    this->cds_phase=1;
                    break;
                case '2':
                    this->cds_phase=2;
                    break;
                case '3':
                    this->cds_phase=3;
                    break;
                case '.':
                    this->cds_phase=0;
                    break;
                default:
                    std::cerr<<"unknown CDS phase for the initial segment"<<std::endl;
                    exit(-1);
            }
        }

        this->id = idx;
        this->tid = tx->getID();
        this->seqid = tx->getGSeqName();
        this->strand = tx->strand;
        this->source = tx->getTrackName();
        for(int i=0;i<tx->exons.Count();i++){
            this->exons.push_back(std::make_tuple(tx->exons.Get(i)->start,tx->exons.Get(i)->end,0));
        }
        if(tx->hasCDS()){
            this->is_coding = true;
        }
        // store attributes
        if (tx->attrs!=NULL) {
            bool trId=false;
            //bool gId=false;
            for (int i=0;i<tx->attrs->Count();i++) {
                const char* attrname=tx->names->attrs.getName(tx->attrs->Get(i)->attr_id);
                const char* attrval=tx->attrs->Get(i)->attr_val;
                if (attrval==NULL || attrval[0]=='\0') continue;
                if (strcmp(attrname, "transcriptID")==0) {
                    if (trId) continue;
                    trId=true;
                }
                if (strcmp(attrname, "transcript_id")==0 && !trId) {
                    attrname="transcriptID";
                    trId=true;
                }
                if (Gstrcmp(attrname, "geneID")==0 || strcmp(attrname, "gene_id")==0){
                    continue;
                }
                this->ait = this->attrs.insert(std::make_pair(attrname,tx->attrs->Get(i)->attr_val));
                if(!this->ait.second){
                    std::cerr<<"Attribute: "<<attrname<<" already exists"<<std::endl;
                    exit(-1);
                }
            }
        }
    }
    ~TX()=default;

    std::string get_attributes(){
        std::string res = "";
        for(auto& a : this->attrs){
            res+=a.first+" \""+a.second+"\";";
        }
        return res;
    }

    void assign_phase(Mods& m) {
        int cdsacc=0;
        if (this->strand=='-') { //reverse strand
            for(int i=m.new_chain.size()-1;i>=0;i--){
                std::get<2>(m.new_chain[i])=(3-cdsacc%3)%3;
                cdsacc+=std::get<1>(m.new_chain[i])-std::get<0>(m.new_chain[i])+1;
            }
        }
        else { //forward strand
            for(auto& cds : m.new_chain){
                std::get<2>(cds)=(3-cdsacc%3)%3;
                cdsacc+=std::get<1>(cds)-std::get<0>(cds)+1;
            }
        }
    }

    int get_pppscore(std::tuple<uint,uint,uint>& coord){

    }

    int exon_count() const{
        return this->exons.size();
    }

    bool operator< (const TX& tx) const{
        if(this->get_seqid()!=tx.get_seqid()){
            return this->get_seqid()<tx.get_seqid();
        }
        else if(this->get_strand()!=tx.get_strand()){
            return this->get_strand()<tx.get_strand();
        }
        else if(this->get_start()!=tx.get_start()){
            return this->get_start()<tx.get_start();
        }
        else if(this->get_end()!=tx.get_end()){
            return this->get_end()<tx.get_end();
        }
        else{ // doesn't matter - they definitely overlap
            return true;
        }
    }
    bool operator> (const TX& tx) const{
        if(this->get_seqid()!=tx.get_seqid()){
            return this->get_seqid()>tx.get_seqid();
        }
        else if(this->get_strand()!=tx.get_strand()){
            return this->get_strand()>tx.get_strand();
        }
        else if(this->get_start()!=tx.get_start()){
            return this->get_start()>tx.get_start();
        }
        else if(this->get_end()!=tx.get_end()){
            return this->get_end()>tx.get_end();
        }
        else{ // doesn't matter - they definitely overlap
            return false;
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const TX& t)
    {
        os << t.get_tid() << "\t" <<std::endl;
        return os;
    }


    bool has_cds(){
        return this->is_coding;
    }

    void build_cds_chain(CDS_CHAIN_TYPE& chain){
        if(!this->is_coding){
            return;
        }
        if(this->cds_phase!=0){
            if(this->strand=='+'){
                this->cds_start+=cds_phase;
            }
            else{
                this->cds_end-=cds_phase;
            }
        }
        for(auto& e : this->exons){

            uint es=std::get<0>(e);
            uint ee=std::get<1>(e);
            if(this->cds_end<es || this->cds_start>ee){
                continue;
            }
            if(this->cds_start>=es && this->cds_start<=ee){
                es = cds_start;
            }
            if(this->cds_end>=es && this->cds_end<=ee){
                ee=this->cds_end;
            }

            chain.push_back(std::make_tuple(es,ee,0));
        }
    }

    int get_end() const{
        return std::get<1>(this->exons.back());
    }
    int get_start() const{
        return std::get<0>(this->exons.front());
    }
    std::string get_seqid() const{
        return this->seqid;
    }
    char get_strand() const{
        return this->strand;
    }
    std::string get_tid() const{
        return this->tid;
    }
    int get_id(){
        return this->id;
    }

    int single_intersection(std::tuple<uint,uint,uint>& i1,std::tuple<uint,uint,uint>& i2,std::tuple<uint,uint,uint>& res){
        uint start = std::max(std::get<0>(i1),std::get<0>(i2));
        uint end = std::min(std::get<1>(i1),std::get<1>(i2));
        if(start<=end){
            std::get<0>(res) = start;
            std::get<1>(res) = end;
            return (end-start)+1;
        }
        std::get<0>(res)=NULL;
        std::get<1>(res)=NULL;
        return 0;
    }

    int intersection(CDS_CHAIN_TYPE& chain1,CDS_CHAIN_TYPE& chain2,CDS_CHAIN_TYPE& res){
        bool found_yet = false;
        int start = 0;
        std::tuple<uint,uint,uint> inter;
        int full_inter_len = 0;
        int inter_len = 0;
        for(auto& i1 : chain1){
            found_yet = false;
            for(int j=start;j<chain2.size();j++){
                inter_len = single_intersection(i1,chain2[j],inter);
                if(inter_len>0){
                    res.push_back(inter);
                    full_inter_len+=inter_len;
                    found_yet = true;
                    if(std::get<1>(inter)==std::get<1>(chain2[j])){
                        start++;
                    }
                }
                else{
                    if(found_yet){ // found the last piece in the intersection for the current i1
                        break;
                    }
                }
            }
        }
        return full_inter_len;
    }

    int cut(uint start,uint end,CDS_CHAIN_TYPE& res){ // cuts it's own exon chain - result contains everything between the start/end
        CDS_CHAIN_TYPE start_end;
        start_end.push_back(std::make_tuple(start,end,0));
        uint len = intersection(this->exons,start_end,res);
        return len;
    }

    int compare(CDS_CHAIN_TYPE& chain1,CDS_CHAIN_TYPE& chain2,Mods& res){
        int c1_i=0,c2_i=0;
        std::tuple<uint,uint,uint> cur_c1 = chain1[c1_i];
        std::tuple<uint,uint,uint> cur_c2 = chain2[c2_i];

        std::vector<std::pair<int,int>> tmp_mods;
        std::vector<std::pair<int,int>> tmp_mods_clean;

        std::tuple<uint,uint,uint> inter;
        while(true){
            int inter_len = single_intersection(cur_c1, cur_c2, inter);

            if(inter_len<=0){
                if(std::get<1>(cur_c2)<std::get<0>(cur_c1)){
                    c2_i++;
                    int right = (std::get<1>(cur_c2)-std::get<0>(cur_c2))+1;
                    if(right!=0){
                        res.missing.push_back(std::make_tuple(std::get<0>(cur_c2),std::get<1>(cur_c2),0));
                    }
                    tmp_mods.push_back(std::make_pair(right,1));

                    if(c2_i == chain2.size()){ // done
                        int left = 0-(((int)std::get<1>(cur_c1)+1) - (int)std::get<0>(cur_c1));
                        if(left!=0){
                            res.extra.push_back(std::make_tuple(std::get<0>(cur_c1),std::get<1>(cur_c1),0));
                        }
                        for(auto cc= chain1.begin() + c1_i + 1; cc != chain1.end(); cc++){
                            left-=(((int)std::get<1>(*cc)+1)-(int)std::get<0>(*cc));
                            if(left!=0){
                                res.extra.push_back(std::make_tuple(std::get<0>(*cc),std::get<1>(*cc),0));
                            }
                        }
                        tmp_mods.push_back(std::make_pair(left,1));
                        break;
                    }
                    else{
                        cur_c2 = chain2[c2_i];
                        continue;
                    }
                }
                if(std::get<1>(cur_c1)<std::get<0>(cur_c2)){
                    c1_i++;
                    int left = 0-((std::get<1>(cur_c1)+1)-std::get<0>(cur_c1));
                    if(left!=0){
                        res.extra.push_back(std::make_tuple(std::get<0>(cur_c1),std::get<1>(cur_c1),0));
                    }
                    tmp_mods.push_back(std::make_pair(left,1));

                    if(c1_i == chain1.size()){ // done
                        int right = ((int)std::get<1>(cur_c2) - (int)std::get<0>(cur_c2)) + 1;
                        if(right!=0){
                            res.missing.push_back(std::make_tuple(std::get<0>(cur_c2),std::get<1>(cur_c2),0));
                        }
                        for(auto cf= chain2.begin() + c2_i + 1; cf != chain2.end(); cf++){
                            right+=(((int)std::get<1>(*cf)-(int)std::get<0>(*cf))+1);
                            if(right!=0){
                                res.missing.push_back(std::make_tuple(std::get<0>(*cf),std::get<1>(*cf),0));
                            }
                        }
                        tmp_mods.push_back(std::make_pair(right,1));
                        break;
                    }
                    else{
                        cur_c1 = chain1[c1_i];
                        continue;
                    }
                }
            }

            int left_start = std::min((int)std::get<0>(cur_c1),(int)std::get<0>(inter));
            int left_end = std::max((int)std::get<0>(cur_c1),(int)std::get<0>(inter));
            int left = 0-(left_end - left_start);
            if(left!=0){
                res.extra.push_back(std::make_tuple(left_start,left_end,0));
            }
            int right_start = std::min((int)std::get<0>(cur_c2), (int)std::get<0>(inter));
            int right_end = std::max((int)std::get<0>(cur_c1), (int)std::get<0>(inter));
            int right = right_end - right_start;
            if(right!=0){
                res.missing.push_back(std::make_tuple(right_start,right_end,0));
            }
            int inters = ((int)std::get<1>(inter)-(int)std::get<0>(inter))+1;

            if(left!=0 && right!=0){
                std::cerr<<"something wrong with operations"<<std::endl;
                exit(-1);
            }
            if(left!=0){
                tmp_mods.push_back(std::make_pair(left,1));
            }
            if(right!=0){
                tmp_mods.push_back(std::make_pair(right,1));
            }
            if(inters!=0){
                tmp_mods.push_back(std::make_pair(inters,0));
            }

            cur_c1 = std::make_tuple(std::get<1>(inter) + 1, std::get<1>(cur_c1),0);
            cur_c2 = std::make_tuple(std::get<1>(inter) + 1, std::get<1>(cur_c2),0);

            if(((int)std::get<1>(cur_c1) - (int)std::get<0>(cur_c1)) < 0){
                c1_i++;
                if(c1_i == chain1.size()){ // done
                    right = ((int)std::get<1>(cur_c2) - (int)std::get<0>(cur_c2)) + 1;
                    if(right!=0){
                        res.missing.push_back(std::make_tuple(std::get<0>(cur_c2),std::get<1>(cur_c2),0));
                    }
                    for(auto cf= chain2.begin() + c2_i + 1; cf != chain2.end(); cf++){
                        right+=(((int)std::get<1>(*cf)-(int)std::get<0>(*cf))+1);
                        if(right!=0){
                            res.missing.push_back(std::make_tuple(std::get<0>(*cf),std::get<1>(*cf),0));
                        }
                    }
                    tmp_mods.push_back(std::make_pair(right,1));
                    break;
                }
                cur_c1 = chain1[c1_i];
            }

            if(((int)std::get<1>(cur_c2) - (int)std::get<0>(cur_c2)) < 0){
                c2_i++;
                if(c2_i == chain2.size()){ // done
                    left = 0-(((int)std::get<1>(cur_c1)+1) - (int)std::get<0>(cur_c1));
                    if(left!=0){
                        res.extra.push_back(std::make_tuple(std::get<0>(cur_c1),std::get<1>(cur_c1),0));
                    }
                    for(auto cc= chain1.begin() + c1_i + 1; cc != chain1.end(); cc++){
                        left-=(((int)std::get<1>(*cc)+1)-(int)std::get<0>(*cc));
                        if(left!=0){
                            res.extra.push_back(std::make_tuple(std::get<0>(*cc),std::get<1>(*cc),0));
                        }
                    }
                    tmp_mods.push_back(std::make_pair(left,1));
                    break;
                }
                cur_c2 = chain2[c2_i];
            }
        }

        for(auto& m : tmp_mods){
            if(tmp_mods_clean.size()>0 && m.second==tmp_mods_clean.back().second){
                tmp_mods_clean.back().first+=m.first;
            }
            else{
                if(m.first==0){
                    continue;
                }
                tmp_mods_clean.push_back(m);
            }
        }

        // extract modifications
        extract_mods(tmp_mods_clean,this->strand=='-',res);

        return res.get_score();
    }

    void extract_mods(std::vector<std::pair<int,int>>& mod_chain,bool rev,Mods& res){
        if(rev){
            std::reverse(mod_chain.begin(),mod_chain.end());
        }

        res.missing_start = mod_chain.front().second==1;
        res.missing_end = mod_chain.back().second==1;

        int frame = 0;
        bool inframe = true;

        for(auto& m : mod_chain){
            if(m.second==1 && m.first<0){ // mod and extra
                res.num_bp_extra+=(0-m.first);
                frame+=std::abs(m.first)%3;
                inframe = frame%3==0;
            }
            if(m.second==1 && m.first>0){ // mod and missing
                res.num_bp_missing+=(m.first);
                frame+=std::abs(m.first)%3;
                inframe = frame%3==0;
            }
            if(m.second==0){
                res.num_bp_match+=(m.first);
                if(inframe){
                    res.num_bp_inframe+=m.first;
                }
                else{
                    res.num_bp_outframe+=m.first;
                }
            }
        }

        return;
    }

    Mods fit_new(Mods& orig_cds_mod,const char* bundle_seq, int bundle_start,uint start_offset,uint end_offset){ // computes the intersection of own exons and the chain
        this->mods.push_back(Mods());
        this->mods.back().orig_cds_tid = orig_cds_mod.orig_cds_tid;
        this->mods.back().orig_chain = orig_cds_mod.new_chain;
        uint cut_len = cut(std::get<0>(orig_cds_mod.new_chain.front()),std::get<1>(orig_cds_mod.new_chain.back()),this->mods.back().new_chain);
        if(cut_len<3){ // no overlap found - create a dummy chain which will never overlap anything. This way the Mods will still get populated with data and will make it possible to filter
            this->mods.back().new_chain = CDS_CHAIN_TYPE{std::make_tuple(0,0,0)};
        }
        else{
            if(bundle_seq!=NULL){ // check_ref flag toggled
                evaluate_fasta(this->mods.back(),bundle_seq,bundle_start,this->strand,start_offset,end_offset);
                adjust_stop(this->mods.back(),bundle_start,this->strand);
                compare_aa(orig_cds_mod,this->mods.back());
            }
        }

        if(this->mods.back().new_chain.empty()){ // no matching chain - likely due to stop-codon
            this->mods.back().new_chain = CDS_CHAIN_TYPE{std::make_tuple(0,0,0)};
        }

        // now we can take the cut piece and directly compare it to the desired CDS counting any changes
        int score = compare(orig_cds_mod.new_chain,this->mods.back().new_chain,this->mods.back());

        if(cut_len!=0){ // no overlap found - create a dummy chain which will never overlap anything. This way the Mods will still get populated with data and will make it possible to filter
            assign_phase(this->mods.back());
        }

        return this->mods.back();
    }

    bool overlap(uint s, uint e) {
        if (s>e){
            std::cerr<<s<<"\t"<<e<<std::endl;
            std::cerr<<"start>end"<<std::endl;
            exit(-1);
        }
        return (this->get_start()<=e && this->get_end()>=s);
    }

    std::string get_stats_str(){
        std::string res = "";

        // since we don't know which one is the correct CDS yet (multiple might fit)
        // we shall create duplicate of the transcript with different CDSs
        std::string out_tid_name = this->get_tid();
        uint chain_i=0;
        for(auto& m : this->mods){
            res += this->get_tid()+"\t"+
                   out_tid_name+":-:"+std::to_string(chain_i)+"\t"+
                   m.orig_cds_tid+"\t"+
                   this->get_seqid()+"\t"+
                   this->get_strand()+"\t"+
                   chain2str(m.orig_chain)+"\t"+
                   chain2str(m.new_chain)+"\t"+
                   std::to_string(m.get_score())+"\t"+
                   std::to_string(m.missing_start)+"\t"+
                   std::to_string(m.missing_end)+"\t"+
                   std::to_string(m.num_bp_outframe)+"\t"+
                   std::to_string(m.num_bp_inframe)+"\t"+
                   std::to_string(m.num_bp_extra)+"\t"+
                   std::to_string(m.num_bp_missing)+"\t"+
                   std::to_string(m.num_bp_match)+"\t"+
                   chain2str(m.missing)+"\t"+
                   chain2str(m.extra)+"\n";
            chain_i++;
        }
        res.pop_back(); // remove last new-line
        return res;
    }

    std::string get_gtf(std::string& tid,Mods& m){
        std::string gtf_str = "";
        // get transcript line
        gtf_str+=this->seqid+"\t"
                 +this->source+"\t"
                 +"transcript"+"\t"
                 +std::to_string(this->get_start())+"\t"
                 +std::to_string(this->get_end())+"\t"
                 +"."+"\t"
                 +this->strand+"\t"
                 +"."+"\t"
                 +"transcript_id \""+tid+"\"; "
                 +this->get_attributes();
        if(m.new_chain.empty()){
            gtf_str+="\n";
        }
        else{
            gtf_str+=m.get_desc()+"\n";
        }
        // get exon lines
        for(auto& e : this->exons){
            gtf_str+=this->seqid+"\t"
                     +this->source+"\t"
                     +"exon"+"\t"
                     +std::to_string(std::get<0>(e))+"\t"
                     +std::to_string(std::get<1>(e))+"\t"
                     +"."+"\t"
                     +this->strand+"\t"
                     +"."+"\t"
                     +"transcript_id \""+tid+"\"; "
                     +this->get_attributes()+"\n";
        }
        if(!m.new_chain.empty()){
            for(auto& c : m.new_chain){
                gtf_str+=this->seqid+"\t"
                         +this->source+"\t"
                         +"CDS"+"\t"
                         +std::to_string(std::get<0>(c))+"\t"
                         +std::to_string(std::get<1>(c))+"\t"
                         +"."+"\t"
                         +this->strand+"\t"
                         +std::to_string(std::get<2>(c))+"\t"
                         +"transcript_id \""+tid+"\"; "
                         +this->get_attributes()+"\n";
            }
        }
        gtf_str.pop_back();
        return gtf_str;
    }

    std::string get_all_gtf(){
        if(this->mods.empty()){ // no mods available - just output what's already there
            return this->get_gtf(this->tid, const_cast<Mods &>(empty_mod));
        }

        std::string gtf_str = "";
        uint chain_i=0;
        std::string cur_tid = "";
        for(auto& m : this->mods){
            cur_tid = this->tid+":-:"+std::to_string(chain_i);
            gtf_str+=get_gtf(cur_tid,m)+"\n";
            chain_i++;
        }

        gtf_str.pop_back(); // remove trailing new-line

        return gtf_str;
    }

private:
    int id = -1;
    std::string tid = "";
    std::string seqid = "";
    char strand = 0;
    CDS_CHAIN_TYPE exons;
    uint cds_start;
    uint cds_end;
    char cds_phase;
    bool is_coding = false;
    std::string source = "";
    std::string mods_gtf_str = "";

    std::vector<Mods> mods;

    std::map<std::string,std::string> attrs;
    std::pair<std::map<std::string,std::string>::iterator,bool> ait;
};

class Bundle{
public:
    Bundle()=default;
    ~Bundle(){
//        delete this->seqid_seq;
    }

    void set_ref(const std::string& ref_fa_fname){
        this->ref_fa_fname = ref_fa_fname;
        this->check_ref = true;
        gfasta.init(ref_fa_fname.c_str());
    }

    bool can_add(TX& t){
        if(this->size==0){
            return true;
        }
        else{
            if(std::strcmp(this->seqid.c_str(),t.get_seqid().c_str())==0 && this->strand==t.get_strand()){
                if(t.overlap(this->start,this->end)){
                    return true;
                }
                else{
                    return false;
                }
            }
            return false;
        }
    }

    bool add_tx(TX& t){
        if(!this->can_add(t)){
            return false;
        }
        this->txs.push_back(t);
        this->seqid = t.get_seqid();
        this->strand = t.get_strand();
        this->start = std::min(this->start,t.get_start());
        this->end = std::max(this->end,t.get_end());
        this->size++;
    }

    void clear(){
        this->txs.clear();
        this->size = 0;
        this->seqid = "";
        this->strand = 0;
        this->start = MAX_INT;
        this->end = 0;
    }

    void print(){
        std::cout<<"process"<<std::endl;
        // construct a searchable of all CDSs in the bundle
        for(auto& tx : this->txs){
            std::cout<<tx.get_tid()<<std::endl;
            if(tx.has_cds()){
                CDS_CHAIN_TYPE cur_cds_chain;
                tx.build_cds_chain(cur_cds_chain);
                std::pair<std::string,CDS_CHAIN_TYPE> cds_chain = {"",cur_cds_chain};
                for(auto& c : cds_chain.second){
                    std::cout<<std::get<0>(c)<<"-"<<std::get<1>(c)<<" ; ";
                }
                std::cout<<std::endl;
            }
        }
        std::cout<<"------"<<std::endl;
    }

    std::string get_nt(CDS_CHAIN_TYPE& chain, const char* subseq,uint start_offset,uint end_offset){
        std::string cds_nt = "";
        if(this->strand=='+'){
            for(auto& c : chain){
                uint cs = std::get<0>(c)-this->start; // start coordinate of the chain with respect to the bundle start
                uint ce = std::get<1>(c)-this->start; // end coordinate of the chain with respect to the bundle start
                for(uint i=cs;i<=ce;i++){
                    cds_nt+=subseq[i];
                }
            }
        }
        else{ // strand=='-'
            for(int ci=chain.size()-1;ci>=0;ci--){
                uint cs = start_offset+(std::get<0>(chain[ci])-this->start); // -3 to adjust for the stop codon on -; start coordinate of the chain with respect to the bundle start
                uint ce = start_offset+(std::get<1>(chain[ci])-this->start); // -3 to adjust for the stop codon on -; end coordinate of the chain with respect to the bundle start
                for (uint i=ce;i>=cs;i--) {
                    cds_nt+=ntComplement(subseq[i]);
                }
            }
        }
        return cds_nt;
    }

    std::string get_aa(CDS_CHAIN_TYPE& chain, const char* subseq,uint start_offset,uint end_offset){
        std::string cds_nt = get_nt(chain,subseq,start_offset,end_offset);
        int cds_len = 0;
        std::string cds_aa = translateDNA(cds_nt.c_str(),cds_len,cds_nt.size());
        return cds_aa;
    }

    int process(Globals& globals){
        const char* bundle_seq = NULL;
        uint start_offset = 0;
        uint end_offset = 0;
        if(this->check_ref){
            this->seqid_seq=fastaSeqGet(gfasta, this->seqid.c_str());
            if(this->seqid_seq==NULL){
                std::cerr<<"Error: no genomic sequence available (check -r option!)."<<std::endl;
            }
            int bundle_start = this->start;
            start_offset = std::min(3,this->start-1); // if start < 3 - will compensate; -1 is to compensate for the 1-based used by gfaseqget
            end_offset = std::min(3,this->seqid_seq->getseqlen()-this->end);
//            if(start_offset<3 || end_offset<3){
//                std::cout<<"found"<<std::endl;
//            }
            if(this->strand=='-'){
                bundle_start-=start_offset;
            }
            if(start_offset>3 || end_offset>3){
                std::cerr<<"shouldn't happen"<<std::endl;
                exit(-1);
            }
            int bundle_len = this->end-this->start;
            if(this->strand=='+'){
                bundle_len+=end_offset;
            }
            else{
                bundle_len+=start_offset;
            }
            if(bundle_start+bundle_len>seqid_seq->getseqlen() || bundle_start<0){
                std::cerr<<"bundle+stop codon extends past the end of the reference sequence"<<std::endl;
                exit(2);
            }
            bundle_seq=this->seqid_seq->subseq(bundle_start,bundle_len); // +3 is to account for the last stop codon if the end of the last CDS corresponds to the end of the last exon
        }
        // built a dict of all ORFs for searching
        std::vector<Mods> cds_chains; // outer pair int is the TID of the transcript
        std::set<CDS_CHAIN_TYPE> cur_cds_chains; // for duplicate removal
        std::pair<std::set<CDS_CHAIN_TYPE>::iterator,bool> cit;
        for(auto& tx : this->txs){
            if(tx.has_cds()){
                if(std::strcmp(tx.get_tid().c_str(),"ENST00000379101.8")==0){
                    std::cout<<"found"<<std::endl;
                }
                CDS_CHAIN_TYPE cur_cds_chain;
                tx.build_cds_chain(cur_cds_chain);
                cit = cur_cds_chains.insert(cur_cds_chain);
                if(cit.second){ // successfully inserted
                    Mods m;
                    m.orig_chain = cur_cds_chain;
                    m.new_chain = cur_cds_chain;
                    m.orig_cds_tid = tx.get_tid();
                    if(bundle_seq!=NULL){
                        evaluate_fasta(m,bundle_seq,this->start,this->strand,start_offset,end_offset);
                    }
                    cds_chains.push_back(m);
                }
            }
        }

        // iterate over each transcript-ORF pair to gauge compatibility - assign compatibility scores
        for(auto& tx : this->txs){
            if(std::strcmp(tx.get_tid().c_str(),"ENST00000379101.8")==0){
                std::cout<<"found"<<std::endl;
            }
            std::string cur_tid = tx.get_tid();
            if(cds_chains.empty()){
                globals.out_gtf_fp<<tx.get_gtf(cur_tid,const_cast<Mods &>(empty_mod))<<std::endl;
                globals.out_gtf_nonoverlap_fp<<tx.get_gtf(cur_tid, const_cast<Mods &>(empty_mod))<<std::endl;

                globals.out_stats_fp<<tx.get_tid()<<"\t"
                         <<tx.get_tid()<<"\t"
                         <<"-\t"
                         <<tx.get_seqid()<<"\t"
                         <<tx.get_strand()<<"\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-\t"
                         <<"-"<<std::endl;
            }
            else{
                // TODO: when a premature stop-codon is found - needs to be stored as a separate stat

                // find the best chain fit
                Mods longest_perfect_mod;
                int longest_perfect_mod_len = 0;
                for(auto& chain : cds_chains){
                    Mods new_mod = tx.fit_new(chain,bundle_seq,this->start,start_offset,end_offset);
                    int new_mod_len = chain_len(new_mod.new_chain);
                    if(new_mod.get_score()==0 && new_mod_len>longest_perfect_mod_len){
                        longest_perfect_mod = new_mod;
                        longest_perfect_mod_len = new_mod_len;
                    }
                }
                globals.out_gtf_fp<<tx.get_all_gtf()<<std::endl;
                globals.out_stats_fp<<tx.get_stats_str()<<std::endl;
                if(longest_perfect_mod_len>0){ // found perfect fit
                    globals.out_gtf_perfect_fp<<tx.get_gtf(cur_tid,longest_perfect_mod)<<std::endl;
                }
                else{ // no perfect fits found
                    globals.out_gtf_imperfect_fp<<tx.get_gtf(cur_tid, const_cast<Mods &>(empty_mod))<<std::endl;
                }
            }
        }

        return 0;
    }
private:
    std::vector<TX> txs;
    int size = 0;
    std::string seqid = "";
    char strand = 0;
    int start = MAX_INT;
    int end = 0;

    bool check_ref = false;
    std::string ref_fa_fname = "";
    GFaSeqGet* seqid_seq = NULL;
};

// uses gffReader to read-in all transcripts
// and then sorts them using custom rules
struct Transcriptome{
public:
    Transcriptome(const std::string& gtf_fname){
        FILE* gff_file = fopen(gtf_fname.c_str(), "r");
        if (gff_file == nullptr)
        {
            std::cerr << "@ERROR::Couldn't open the GTF: " << gtf_fname << std::endl;
            exit(1);
        }
        GffReader gffReader(gff_file,false,false);
        gffReader.readAll();

        GffObj *pGffObj;
        for(int i=0;i<gffReader.gflst.Count();++i) {
            pGffObj = gffReader.gflst.Get(i);
            tx_vec.push_back(TX(pGffObj,i));
        }
    }
    ~Transcriptome()=default;

    void add(const std::string& gtf_fname){
        FILE* gff_file = fopen(gtf_fname.c_str(), "r");
        if (gff_file == nullptr)
        {
            std::cerr << "@ERROR::Couldn't open the GTF: " << gtf_fname << std::endl;
            exit(1);
        }
        GffReader gffReader(gff_file,false,false);
        gffReader.readAll();

        GffObj *pGffObj;
        for(int i=0;i<gffReader.gflst.Count();++i) {
            pGffObj = gffReader.gflst.Get(i);
            tx_vec.push_back(TX(pGffObj,i));
        }
    }

    void sort(){
        std::sort(this->tx_vec.begin(),this->tx_vec.end());
    }

    typedef std::vector<TX>::iterator it;
    typedef std::vector<TX>::const_iterator cit;
    it begin() {return this->tx_vec.begin();}
    cit cbegin() const { return this->tx_vec.cbegin();}
    it end() {return this->tx_vec.end();}
    cit cend() const { return this->tx_vec.cend();}
private:
    std::vector<TX> tx_vec;
};

int run(const std::vector<std::string>& known_gtf_fnames, const std::string& novel_gtf_fname,const std::string& out_fname,const std::string& ref_fa_fname,const std::string& ppp_track_fname){
    globals.out_stats_fp.open(out_fname+".stats");
    globals.out_cds_stats_fp.open(out_fname+".cds");
    globals.out_gtf_fp.open(out_fname+".all.gtf");
    globals.out_gtf_perfect_fp.open(out_fname+".perfect.gtf");
    globals.out_gtf_imperfect_fp.open(out_fname+".imperfect.gtf");
    globals.out_gtf_nonoverlap_fp.open(out_fname+".nonoverlapping.gtf");

    // read from both GTF streams simultaneously
    // only consider CDS from the reference and exons from the novel
    // form bundles (all overlapping

    // load the reference GFF
    Transcriptome transcriptome(novel_gtf_fname); // TODO: if transcript has a cds - just report it unless better CDS available
    for(auto& k : known_gtf_fnames){
        transcriptome.add(k);
    }
    transcriptome.sort();

    Bundle bundle;

    bool check_ref = !ref_fa_fname.empty();
    if(check_ref){
        bundle.set_ref(ref_fa_fname);
    }

    globals.out_stats_fp<<"tid\t"
               "new_tid\t"
               "cds_tid\t"
               "seqid\t"
               "strand\t"
               "original_cds\t"
               "fitted_cds\t"
               "score\t"
               "missing_start\t"
               "missing_end\t"
               "num_bp_outframe\t"
               "num_bp_inframe\t"
               "num_bp_extra\t"
               "num_bp_missing\t"
               "num_bp_match\t"
               "missing_frags\t"
               "extra_frags"<<std::endl;

    globals.out_cds_stats_fp<<"cds_tid\t"
                      "chain\t"
                      "overlapping_txs\t"
                      "perfect_txs\t"
                      "num_overlapping_txs\t"
                      "num_perfect_fit_txs\t"
                      "start\t"
                      "stop\t"
                      "premature_stop\t"<<std::endl;

    for(auto& t : transcriptome){
        if(!bundle.can_add(t)){
            bundle.process(globals);
            bundle.clear();
            bundle.add_tx(t);
        }
        else{
            bundle.add_tx(t);
        }
    }
    // handle the final bundle
    bundle.process(globals);

    // questions:
    // 1. any transcripts contain multiple CDSs? which one to chose?
    // 2. what to do with minor modifications? (no frame shifts)

    globals.out_stats_fp.close();
    globals.out_cds_stats_fp.close();
    globals.out_gtf_fp.close();
    globals.out_gtf_perfect_fp.close();
    globals.out_gtf_imperfect_fp.close();
    globals.out_gtf_nonoverlap_fp.close();
    return 0;
}

enum Opt {CDS       = 'c',
          INPUT     = 'i',
          OUTPUT    = 'o',
          REFERENCE = 'r',
          CLEANREF  = 'c',
          FILTER    = 'f',
          PPPTRACK  = 'p'};

int main(int argc, char** argv) {

    ArgParse args("orfanage");
    args.add_multi_string(Opt::CDS,"cds","","Comma-separated list of GTF filenames with known transcripts and CDS annotations",true);
    args.add_string(Opt::INPUT, "input", "", "Input GTF with transcripts to which CDSs are to be ported", true);
    args.add_string(Opt::OUTPUT, "output", "", "Output name", true);
    args.add_string(Opt::REFERENCE,"reference","","Reference fasta",false);
    args.add_flag(Opt::CLEANREF,"cleanref","Remove transcripts which contain mistakes in pre-annotated CDS",false);
    args.add_flag(Opt::FILTER,"filter","Select the best fitting CDS where possible",false);
    args.add_string(Opt::PPPTRACK,"ppptrack","","PhyloCSF track",false);

    if(argc <= 1 || strcmp(argv[1], "--help") == 0){
        std::cerr << args.get_help() << std::endl;
        exit(1);
    }

    args.parse_args(argc, argv);

    // first create the execution string
    std::string cl = "orfanage ";
    for (int i = 0; i < argc; i++) {
        if (i == 0) {
            cl += argv[i];
        } else {
            cl += " ";
            cl += argv[i];
        }
    }

    // check that all files exist
    std::ifstream if_ss;
    if_ss.open(args.get_string(Opt::INPUT));
    if (!if_ss){
        std::cerr << "Input file does not exist! "<<args.get_string(Opt::INPUT)<<std::endl;
        exit(2);
    }
    if_ss.close();

    const std::vector<std::string> knowns = args.get_multi_string(Opt::CDS);
    for(auto& k : knowns){
        if_ss.open(k);
        if(!if_ss){
            std::cerr<<"Input file does not exist! "<<k<<std::endl;
            exit(2);
        }
        if_ss.close();
    }

    if(args.is_set(Opt::REFERENCE)){ // much from gpertea bamcons.cpp
        std::ifstream rf_ss;
        rf_ss.open(args.get_string(Opt::REFERENCE));
        if (!rf_ss){
            std::cerr << "Reference FASTA file does not exist!\n";
            exit(2);
        }
        rf_ss.close();

        // get potential fasta index file name
        std::string fa_idx_fname = args.get_string(Opt::REFERENCE)+".fai";
        GFastaIndex faIdx(args.get_string(Opt::REFERENCE).c_str(), fa_idx_fname.c_str());
        if (!faIdx.hasIndex()){
            std::cerr<<"No fasta index found for "<<fa_idx_fname<<". Building now"<<std::endl;
            faIdx.buildIndex();
            if (faIdx.getCount() == 0){
                std::cerr<<"Error: no fasta records found!"<<std::endl;
                exit(2);
            }
            FILE* fcreate = fopen(fa_idx_fname.c_str(), "w");
            if (fcreate == NULL){
                std::cerr<<"Error: cannot create fasta index: "<<fa_idx_fname<<std::endl;
                exit(2);
            }
            if (faIdx.storeIndex(fcreate) < faIdx.getCount()){
                std::cerr<<"Warning: error writing the index file!"<<std::endl;
                exit(2);
            }
            std::cerr<<"FASTA index rebuilt."<<std::endl;
        }
    }

//    // check phylocsf track if set
//    if(args.is_set(Opt::PPPTRACK)){
//        std::ifstream ppp_ss;
//        ppp_ss.open(args.get_string(Opt::PPPTRACK));
//        if (!ppp_ss){
//            std::cerr << "PhyloCSF track file does not exist! "<<args.get_string(Opt::PPPTRACK)<<std::endl;
//            exit(2);
//        }
//        ppp_ss.close();


//        // TESTS
//        bigWigFile_t *bw_files[7] = { NULL };
//        std::string bw_path = args.get_string(Opt::PPPTRACK);
//        const size_t bw_path_suffix_pos = bw_path.find("+1");
//        if (bw_path_suffix_pos == std::string::npos)
//        {
//            std::cerr<<"Could not find '+1' in tracks file name. Expecting a name like 'PhyloCSF+1.bw'."<<std::endl;
//            exit(3);
//        }
//        for (uint16_t i = 0; i < 7; ++i)
//        {
//            std::string suffix;
//            if (i == 6)
//                suffix = "power";
//            else
//                suffix = ((i < 3) ? "+" : "-") + std::to_string((i % 3) + 1);
//            bw_path.replace(bw_path_suffix_pos, 2, suffix); // NOTE: length of "+1" is 2
//
//            bw_files[i] = bwOpen(const_cast<char * >(bw_path.c_str()), NULL, "r");
//            if (!bw_files[i])
//            {
//                // check whether the user has used an unindexed wig file, then print a useful hint
//                if (access(bw_path.c_str(), F_OK) == 0 &&
//                    bw_path.size() >= 4 && bw_path.compare(bw_path.size() - 4, 4, ".wig") == 0)
//                {
//                    std::cerr<<"An error occurred while opening the PhyloCSF file '%s'."<<bw_path.c_str()<<std::endl;
//                    std::cerr<<"It seems you provided a *.wig file. You need to simply index them first with wigToBigWig and then use the *.bw files."<<std::endl;
//                    exit(3);
//                }
//                else
//                {
//                    std::cerr<<"Could not find PhyloCSF track file '%s'."<<bw_path.c_str()<<std::endl;
//                    exit(3);
//                }
//            }
//        }

//        chromList_t *chr_list = bwReadChromList(bw_files[0]);
//        for (int64_t i = 0; i < chr_list->nKeys; ++i)
//        {
//            std::map<std::string, uint64_t> chrom_sizes;
//            chrom_sizes.insert(std::make_pair(chr_list->chrom[i], chr_list->len[i]));
//        }
//    }

    // run
    run(knowns,args.get_string(Opt::INPUT),args.get_string(Opt::OUTPUT),args.get_string(Opt::REFERENCE),args.get_string(Opt::PPPTRACK));

    return 0;
}