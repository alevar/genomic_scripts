//
// Created by Ales Varabyou, Beril Erdogdu, Natalia Rincon on 03/30/21.
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <map>

#include <gclib/gff.h>
#include <fstream>
#include <set>
#include "arg_parse.h"

typedef std::vector<std::array<int,3>> CDS_CHAIN_TYPE; // start,end,phase

struct Globals{
    std::ofstream out_reader_gtf_fp;
    int po_threshold;
} globals;

void assign_phase(CDS_CHAIN_TYPE& chain,char strand,int start_phase) {
    int cdsacc=start_phase;
    if (strand=='-') { //reverse strand
        for(int i=chain.size()-1;i>=0;i--){
            std::get<2>(chain[i])=(3-cdsacc%3)%3;
            cdsacc+=std::get<1>(chain[i])-std::get<0>(chain[i])+1;
        }
    }
    else { //forward strand
        for(auto& cds : chain){
            std::get<2>(cds)=(3-cdsacc%3)%3;
            cdsacc+=std::get<1>(cds)-std::get<0>(cds)+1;
        }
    }
}

std::string chain2str(CDS_CHAIN_TYPE& chain){
    std::string res;
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

int single_intersection(std::array<int,3>& i1,std::array<int,3>& i2,std::array<int,3>& res){
    int start = std::max(std::get<0>(i1),std::get<0>(i2));
    int end = std::min(std::get<1>(i1),std::get<1>(i2));
    if(start<=end){
        std::get<0>(res) = start;
        std::get<1>(res) = end;
        return (end-start)+1;
    }
    std::get<0>(res)=0;
    std::get<1>(res)=0;
    return 0;
}

int intersection(CDS_CHAIN_TYPE& chain1,CDS_CHAIN_TYPE& chain2,CDS_CHAIN_TYPE& res){
    bool found_yet = false;
    int start = 0;
    std::array<int,3> inter{};
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

class TX{
public:
    TX() = default;
    TX(GffObj* tx,int idx,bool is_templ,std::string tag=""){
        if(tx->hasCDS()){
            cds_start = (int)tx->CDstart;
            cds_end = (int)tx->CDend;
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

        this->tag = tag;
        this->is_templ=is_templ;
        this->id = idx;
        this->tid = tx->getID();
        this->seqid = tx->getGSeqName();
        this->strand = tx->strand;
        this->source = tx->getTrackName();
        for(int i=0;i<tx->exons.Count();i++){
            this->exons.push_back(std::array<int,3>{(int)tx->exons.Get(i)->start,(int)tx->exons.Get(i)->end,0});
        }
        if(tx->hasCDS()){
            this->is_coding = true;
        }

        // store attributes
        const char* gid = tx->getGeneID();
        if(!(gid==NULL || gid[0]=='\0')){
            ait = this->attrs.insert(std::make_pair("gene_id",gid));
            if(!ait.second){
                std::cerr<<"Gene id already assigned"<<std::endl;
                exit(-1);
            }
        }
        else{
            ait = this->attrs.insert(std::make_pair("gene_id",tid));
        }
        if (tx->attrs!=NULL) {
            bool trId=false;
            //bool gId=false;
            for (int i=0;i<tx->attrs->Count();i++) {
                const char* attrname=GffObj::names->attrs.getName(tx->attrs->Get(i)->attr_id);
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
//                if (Gstrcmp(attrname, "geneID")==0 || strcmp(attrname, "gene_id")==0){
//                    continue;
//                }
                ait = this->attrs.insert(std::make_pair(attrname,tx->attrs->Get(i)->attr_val));
                if(!ait.second){
                    std::cerr<<"Attribute: "<<attrname<<" already exists"<<std::endl;
                    exit(-1);
                }
            }
        }
    }
    ~TX()=default;

    std::string get_igv_str(bool include_strand=true){
        if(include_strand){
            return this->seqid+this->strand+":"+std::to_string(this->get_start())+"-"+std::to_string(this->get_end());
        }
        else{
            return this->seqid+":"+std::to_string(this->get_start())+"-"+std::to_string(this->get_end());
        }
    }

    void set_cds_start(int cs){
        this->cds_start = cs;
    }
    void set_cds_end(int ce){
        this->cds_end = ce;
    }

    int get_cds_start(){
        return this->cds_start;
    }
    int get_cds_end(){
        return this->cds_end;
    }

    bool is_tag(std::string t){
        return this->tag==t;
    }

    void add_attribute(std::string k, std::string v){
        this->ait = this->attrs.insert(std::make_pair(k,v));
        if(!this->ait.second){
            std::cerr<<"attribute key already exists"<<std::endl;
            exit(-1);
        }
    }

    void update_attribute(std::string k,std::string v, char delim=',',bool replace = false){
        this->ait = this->attrs.insert(std::make_pair(k,v));
        if(!this->ait.second && replace){
            this->attrs.erase(this->ait.first);
            this->ait = this->attrs.insert(std::make_pair(k,v));
            return;
        }
        if(!this->ait.second && !replace){ // if already exists - add to the value list with the delimiter
            this->ait.first->second = this->ait.first->second+delim+v;
            return;
        }
    }

    std::string get_attributes(bool skip_gene_id = false){
        std::string res;
        for(auto& a : this->attrs){
            if(skip_gene_id && a.first=="gene_id"){
                continue;
            }
            res+=a.first+" \""+a.second+"\"; ";
        }
        res.pop_back();
        return res;
    }

    std::string get_attribute(const std::string key){ // returns NULL if key not found
        std::string val;
        this->ait.first = this->attrs.find(key);
        if(this->ait.first!=this->attrs.end()){
            return this->ait.first->second;
        }
        return "";
    }

    int exon_count() const{
        return this->exons.size();
    }

    int intron_count() const{
        return this->exons.size()-1;
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
            return false;
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

    bool intron_eq(const TX& tx) const{ // always returns false for single-exon transcripts
        if(this->get_seqid()!=tx.get_seqid()){
            return false;
        }
        else if(this->get_strand()!=tx.get_strand()){
            return false;
        }
        else if(this->intron_count()<=0 || tx.intron_count()<=0){ // always returns false for single-exon transcripts
            return false;
        }
        else if(this->intron_count()!=tx.intron_count()){
            return false;
        }

        // compare introns
        for(int i=0;i<this->intron_count();i++){
            if(std::get<1>(this->exons[i])!=std::get<1>(tx.exons[i]) ||
                    std::get<0>(this->exons[i+1])!=std::get<0>(tx.exons[i+1])){
                return false;
            }
        }
        return true;
    }

    int num_introns_shared (const TX& tx) const{ // returns the number of matching introns
        if(this->get_seqid()!=tx.get_seqid()){
            return 0;
        }
        else if(this->get_strand()!=tx.get_strand()){
            return 0;
        }
        else if(this->intron_count()<=0 || tx.intron_count()<=0){ // always returns false for single-exon transcripts
            return 0;
        }
        else{
            int nit=0;
            std::set<std::pair<int,int>> its;
            std::set<std::pair<int,int>>::iterator its_it;
            for (int i = 0; i < this->intron_count(); i++) {
                its.insert(std::make_pair(std::get<1>(this->exons[i]),std::get<0>(this->exons[i+1])));
            }
            for (int i = 0; i < tx.intron_count(); i++) {
                its_it = its.find(std::make_pair(std::get<1>(tx.exons[i]),std::get<0>(tx.exons[i+1])));
                if(its_it!=its.end()){
                    nit+=1;
                }
            }
            return nit;
        }
    }

    int num_overlapping_exons(const TX& tx) const{ // returns the number of exons in the current transcript which overlap the other transcript;
        int count = 0;
        std::array<int,3> i2{tx.get_start(),tx.get_end(),0};
        std::array<int,3> res{0,0,0};
        for(std::array<int,3> e : this->exons){
            if(single_intersection(e,i2,res)!=0){
                count+=1;
            }
        }
        return count;
    }

    int poverlap(TX& tx,bool consider_strand, bool request_contained_in_exon){ // computes percent overlap between two transcripts; contained_in_exon if set to true will only compute the number of matching bases if the "this->" transcript is entirely contained within the argument
        int nmbp = bpoverlap(tx,consider_strand,request_contained_in_exon);
        if(nmbp == 0){
            return 0;
        }
        else{
            float pop = (static_cast<float>(nmbp)*100.0)/static_cast<float>(this->length());
            return static_cast<int>(pop);
        }
    }

    bool contains(TX& tx,bool consider_strand) const{
        if(consider_strand && this->strand != tx.strand){
            return false;
        }
        if(this->get_seqid()!=tx.get_seqid()){
            return false;
        }
        return this->get_start()<=tx.get_start() && this->get_end()>=tx.get_end();
    }

    int bpoverlap(TX& tx,bool consider_strand, bool request_contained_in_exon){
        if(this->get_seqid()!=tx.get_seqid()){
            return 0;
        }
        else if(consider_strand && this->get_strand()!=tx.get_strand()){
            return 0;
        }
        else{ // count number of matching bases
            if(request_contained_in_exon && !is_contained_in_exon(tx)){
                return 0;
            }
            int nmbp = this->num_matching_bp(tx.exons);
            return nmbp;
        }
    }

    bool is_contained_in_exon(TX& tx) const{
        for(auto& e : tx.exons){
            if(std::get<0>(e)<=this->get_start() && std::get<1>(e)>=this->get_end()){
                return true;
            }
        }
        return false;
    }

    int length(){
        int rl = 0;
        for(auto& e : this->exons){
            rl+=(std::get<1>(e)-std::get<0>(e))+1;
        }
        return rl;
    }

    friend std::ostream& operator<<(std::ostream& os, const TX& t){
        os << t.get_tid() << "\t" <<std::endl;
        return os;
    }

    CDS_CHAIN_TYPE get_exons(){
        return this->exons;
    }

    CDS_CHAIN_TYPE get_exons(int start_idx,int end_idx){
        CDS_CHAIN_TYPE::const_iterator first = this->exons.begin() + start_idx;
        CDS_CHAIN_TYPE::const_iterator last = this->exons.begin() + end_idx;
        CDS_CHAIN_TYPE subvec(first, last);
        return subvec;
    }

    std::array<int,3> get_exon(int pos){
        return this->exons[pos];
    }

    std::pair<int,int> get_intron(int first_exon_pos){
        if(first_exon_pos<0 || first_exon_pos>=this->exon_count()-1){
            std::cerr<<"requesting intron at non-existing position"<<std::endl;
        }
        return std::make_pair(std::get<1>(this->exons[first_exon_pos]),std::get<0>(this->exons[first_exon_pos+1]));
    }

    bool has_cds() const{
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
            this->cds_phase = 0;
        }
        int cp = this->cds_phase;
        for(auto& e : this->exons){
            int es=std::get<0>(e);
            int ee=std::get<1>(e);
            if(this->cds_end<es || this->cds_start>ee){
                continue;
            }
            if(this->cds_start>=es && this->cds_start<=ee){
                es = cds_start;
            }
            if(this->cds_end>=es && this->cds_end<=ee){
                ee=this->cds_end;
            }
            chain.push_back(std::array<int,3>{es,ee,cp});
        }
        assign_phase(chain,this->strand,0);
    }

    int get_end() const{
        return std::get<1>(this->exons.back());
    }
    int get_start() const{
        return std::get<0>(this->exons.front());
    }
    void set_end(int new_end){
        if(new_end<std::get<0>(this->exons[this->exons.size()-1])){ // if smaller than the start of last exon
            std::cerr<<"incompatible coordinate for the new end"<<std::endl;
            exit(-1);
        }
        std::get<1>(this->exons[this->exons.size()-1]) = new_end;
    }
    void set_start(int new_start){
        if(new_start>std::get<1>(this->exons[0])){ // if greater than the end of the first exon
            std::cerr<<"incompatible coordinate for the new start"<<std::endl;
            exit(-1);
        }
        std::get<0>(this->exons[0]) = new_start;
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
    int get_id() const{
        return this->id;
    }

    // cut strictly - means start and end coordinates must be in the original chain
    // this is different from cut, where start and end can be in reference intronic regions
    int strict_cut(CDS_CHAIN_TYPE& template_chain,CDS_CHAIN_TYPE& res){
        CDS_CHAIN_TYPE tmp_res;
        int tmp_len = intersection(this->exons,template_chain,tmp_res);
        if(tmp_len>0){
            CDS_CHAIN_TYPE start_end;
            int trs = std::get<0>(tmp_res.front());
            int tre = std::get<1>(tmp_res.back());
            start_end.push_back(std::array<int,3>({trs,tre,0}));
            int len = intersection(this->exons,start_end,res);
            return len;
        }
        else{
            res = tmp_res;
            return tmp_len;
        }
    }

    // in this mode only the intersection is computed
    // the longest set of consecutive intervals in which each base belongs to both template and query is reported
    // instead of returning everything between matching start/end coordinates
    int very_strict_cut(CDS_CHAIN_TYPE& template_chain,CDS_CHAIN_TYPE& res){
        std::vector<CDS_CHAIN_TYPE> tmp_res;
        bool found_yet = false;
        int start = 0;
        std::array<int,3> inter{};
        int inter_len = 0;
        bool continue_prev = false;
        for(auto& i1 : this->exons){
            found_yet = false;
            for(int j=start;j<template_chain.size();j++){
                inter_len = single_intersection(i1,template_chain[j],inter);
                if(inter_len>0){
                    // check whether intersection follows the same chain as before - if not a new intersection chain is initiated
                    int inter_e = std::get<1>(inter);
                    int i1_e = std::get<1>(i1);
                    int t_e = std::get<1>(template_chain[j]);
                    bool ends_match = inter_e==i1_e && i1_e==t_e;

                    int inter_s = std::get<0>(inter);
                    int i1_s = std::get<0>(i1);
                    int t_s = std::get<0>(template_chain[j]);
                    bool starts_match = inter_s==i1_s && i1_s==t_s;
                    if(!continue_prev || !starts_match) { // if previous ends didn't match or current starts don't match - create new intersection chain
                        tmp_res.push_back(CDS_CHAIN_TYPE{});
                    }
                    tmp_res.back().push_back(inter);
                    continue_prev = ends_match;

                    found_yet = true;
                    if(std::get<1>(inter)==std::get<1>(template_chain[j])){
                        start++;
                    }
                }
                else{
                    if(found_yet){ // found the last piece in the intersection for the current i1
                        break;
                    }
                    else{
                        continue_prev = false;
                    }
                }
            }
            if(!found_yet && continue_prev){
                continue_prev = false;
            }
        }

        // find the longest chain
        int max_len = 0;
        for(auto& t : tmp_res){
            int cur_len = chain_len(t);
            if(cur_len>max_len){
                res = t;
                max_len = cur_len;
            }
        }
        return max_len;
    }

    // same as very_strict_cut but only computes the number of overlapping bases - not just for the longest chain
    int num_matching_bp(CDS_CHAIN_TYPE& template_chain){
        int res = 0;

        bool found_yet = false;
        int start = 0;
        std::array<int,3> inter{};
        int inter_len = 0;
        bool continue_prev = false;
        for(auto& i1 : this->exons){
            found_yet = false;
            for(int j=start;j<template_chain.size();j++){
                inter_len = single_intersection(i1,template_chain[j],inter);
                if(inter_len>0){
                    // check whether intersection follows the same chain as before - if not a new intersection chain is initiated
                    int inter_e = std::get<1>(inter);
                    int i1_e = std::get<1>(i1);
                    int t_e = std::get<1>(template_chain[j]);
                    bool ends_match = inter_e==i1_e && i1_e==t_e;

                    int inter_s = std::get<0>(inter);
                    int i1_s = std::get<0>(i1);
                    int t_s = std::get<0>(template_chain[j]);
                    bool starts_match = inter_s==i1_s && i1_s==t_s;

                    res+=((std::get<1>(inter) + 1) - std::get<0>(inter));
                    continue_prev = ends_match;

                    found_yet = true;
                    if(std::get<1>(inter)==std::get<1>(template_chain[j])){
                        start++;
                    }
                }
                else{
                    if(found_yet){ // found the last piece in the intersection for the current i1
                        break;
                    }
                    else{
                        continue_prev = false;
                    }
                }
            }
            if(!found_yet && continue_prev){
                continue_prev = false;
            }
        }
        return res;
    }

    int cut(int start,int end,CDS_CHAIN_TYPE& res){ // cuts it's own exon chain - result contains everything between the start/end
        CDS_CHAIN_TYPE start_end;
        start_end.push_back(std::array<int,3>{start,end,0});
        int len = intersection(this->exons,start_end,res);
        return len;
    }

    bool overlap(int s, int e) const {
        if (s>e){
            std::cerr<<s<<"\t"<<e<<std::endl;
            std::cerr<<"start>end"<<std::endl;
            exit(-1);
        }
        return (this->get_start()<=e && this->get_end()>=s);
    }

    bool overlap(TX& tx) const{
        int s = tx.get_start();
        int e = tx.get_end();
        return this->overlap(s,e);
    }

    std::string get_gtf(){
        std::string gtf_str;
        // get transcript line
        gtf_str+=this->seqid+"\t"
                 +this->source+"\t"
                 +"transcript"+"\t"
                 +std::to_string(this->get_start())+"\t"
                 +std::to_string(this->get_end())+"\t"
                 +"."+"\t"
                 +this->strand+"\t"
                 +"."+"\t"
                 +"transcript_id \""+tid+"\"; "+"gene_id \""+this->get_attribute("gene_id")+"\"; "
                 +this->get_attributes(true);
        gtf_str.pop_back();
        gtf_str+="\n";
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
                     +"transcript_id \""+tid+"\";\n";
        }
        if(this->has_cds()){
            CDS_CHAIN_TYPE cds_chain;
            this->build_cds_chain(cds_chain);
            for(auto& c : cds_chain){
                gtf_str+=this->seqid+"\t"
                         +this->source+"\t"
                         +"CDS"+"\t"
                         +std::to_string(std::get<0>(c))+"\t"
                         +std::to_string(std::get<1>(c))+"\t"
                         +"."+"\t"
                         +this->strand+"\t"
                         +std::to_string(std::get<2>(c))+"\t"
                         +"transcript_id \""+tid+"\";\n";
            }
        }
        gtf_str.pop_back();
        return gtf_str;
    }

    bool is_template() const{
        return this->is_templ;
    }

    std::string get_source(){
        return this->source;
    }

private:
    bool is_templ = false;
    std::string tag; // used to set specific tags to distinguish transcripts by type
    int id = -1;
    std::string tid;
    std::string seqid;
    char strand = '.';
    CDS_CHAIN_TYPE exons;
    int cds_start = 0;
    int cds_end = 0;
    int cds_phase;
    bool is_coding = false;
    std::string source;

    std::map<std::string,std::string> attrs;
    std::pair<std::map<std::string,std::string>::iterator,bool> ait;
};

class Bundle{
public:
    Bundle()=default;
    ~Bundle()= default;

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
        if(t.is_template()){
            this->contains_ref = true;
        }
        this->txs.push_back(t);
        this->seqid = t.get_seqid();
        this->strand = t.get_strand();
        this->start = std::min(this->start,t.get_start());
        this->end = std::max(this->end,t.get_end());

        // assign bundle-specific IDs to the transcripts within
        this->t2id.insert(std::make_pair(t.get_id(),this->size));

        this->size++;
        return true;
    }

    int get_bundle_tid(TX& t){
        this->t2id_it.first = this->t2id.find(t.get_id());
        if(this->t2id_it.first != this->t2id.end()){
            return this->t2id_it.first->second;
        }
        std::cerr<<"ID not found in the map"<<std::endl;
        exit(-1);
    }

    bool has_ref() const{
        return this->contains_ref;
    }

    void clear(){
        this->txs.clear();
        this->size = 0;
        this->seqid = "";
        this->strand = 0;
        this->start = MAX_INT;
        this->end = 0;
        this->t2id.erase(this->t2id.begin(),this->t2id.end());
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

    // TODO: why did gene KRTAP6-2 chr21:30,530,659-30,621,497 (one of many examples) not get included???!!! It's in both gencode and refseq (not in mane) and also in chess 2
    //     could it be due to gffcompare? also why are these not showing on the main IGV?

    void build_compat_mat(std::vector<std::vector<std::pair<int,int>>>& compat){
        compat = std::vector<std::vector<std::pair<int,int>>>(this->size,std::vector<std::pair<int,int>>(this->size));

        std::string tx1_gid,tx2_gid;
        int bid1,bid2;
        for(auto& tx1 : this->txs){
            bid1 = this->get_bundle_tid(tx1);
            tx1_gid = tx1.get_attribute("gene_id");
            for(auto& tx2 : this->txs) {
                bid2 = this->get_bundle_tid(tx2);
                int status = 0;
                int nmbp = evaluate(tx1,tx2,status); // evaluate returns the number of bases of overlap
                compat[bid1][bid2] = std::make_pair(status,nmbp);
            }
        }
    }

    void group_identical_genes(std::vector<std::vector<std::pair<int,int>>>& compat,std::map<std::string,std::set<std::string>>& onetomany){ // find identical genes
        TX tx1,tx2;
        std::string tx1_gid,tx2_gid;
        int idx1=0,idx2=0;

        std::pair<std::map<std::string,std::set<std::string>>::iterator,bool> otm_it;
        idx1 = 0;
        for(auto& txv : compat) {
            tx1 = this->txs[idx1];
            if (tx1.is_template()) { // found template transcript - add to the map appropriately
                // check if already in the map
                tx1_gid = tx1.get_attribute("gene_id");
                otm_it = onetomany.insert(std::make_pair(tx1_gid,std::set<std::string>{tx1_gid}));

                idx2 = 0;
                for(auto& val: txv) {
                    tx2 = this->txs[idx2];
                    if (tx2.is_template()) {
                        tx2_gid = tx2.get_attribute("gene_id");
                        if (val.first == EVAL_STATUS::st_full_intron_match ||
                            (val.first == EVAL_STATUS::st_se_overlap &&
                             ((val.second * 100) / tx1.length() >= globals.po_threshold ||
                              (val.second * 100) / tx2.length() >= globals.po_threshold))) { // TODO: move globals to the private member
                            otm_it.first->second.insert(tx2_gid);
                        }
                    }
                    idx2++;
                }
            }
            idx1++;
        }
    }

    int merge(int* parent, int x){
        if (parent[x] == x)
            return x;
        return merge(parent, parent[x]);
    }
    void connectedComponents(int n, std::vector<std::vector<int> >& edges,std::map<int,std::vector<int>>& res)
    {
        if(edges.size()==0){
            return;
        }
        int parent[n];
        for (int i = 0; i < n; i++) {
            parent[i] = i;
        }
        std::vector<bool> e_indices(n,false);
        for (auto x : edges) {
            parent[merge(parent, x[0])] = merge(parent, x[1]);
            e_indices[x[0]]=true;
            e_indices[x[1]]=true;
        }
        for (int i = 0; i < n; i++) {
            if(!e_indices[i]){
                parent[i]=-1;
                continue;
            }
            parent[i] = merge(parent, parent[i]);
        }
        for (int i = 0; i < n; i++) {
            if(!e_indices[i]){
                continue;
            }
            res[parent[i]].push_back(i);
        }
        return;
    }

    int generate_refmap(std::vector<std::vector<std::pair<int,int>>>& compat,std::vector<int>& res){
        std::vector<std::vector<int>> edges;

        TX tx1,tx2;
        int idx1=0,idx2=0;

        idx1 = 0;
        for(auto& txv : compat) {
            tx1 = this->txs[idx1];
            if (tx1.is_template()) {
                idx2 = 0;
                for(auto& val: txv) {
                    if(idx2<idx1){
                        idx2++;
                        continue;
                    }
                    tx2 = this->txs[idx2];
                    if (tx2.is_template()) {
                        if (val.first!=EVAL_STATUS::st_no_overlap) {
                            edges.push_back(std::vector<int>{idx1,idx2});
                        }
                    }
                    idx2++;
                }
            }
            idx1++;
        }

        std::map<int, std::vector<int>> ccs; //conected components;
        connectedComponents(compat.size(),edges,ccs);

        res = std::vector<int>(compat.size(),-1);
        int kv_idx=0;
        for(auto& kv : ccs){
            for(auto& v : kv.second){
                res[v]=kv_idx;
            }
            kv_idx+=1;
        }

        return 0;
    }

    int process(Globals& globals_local){ // returns true if new gene_id has been assigned

        std::vector<std::vector<std::pair<int,int>>> compat(this->size,std::vector<std::pair<int,int>>(this->size));
        build_compat_mat(compat);

        std::vector<int> gene_map; // each position is the index of transcript in this->txs and the value is the gene assignment
        int res = generate_refmap(compat,gene_map);

        // evaluate queries
        int qidx=0;
        for(auto& txv : compat){
            TX& q = this->txs[qidx];
            if(!q.is_template()){
                int tidx=0;
                // collect IDs of all overlapping genes
                std::set<int> ids;
                for(auto& val: txv){
                    TX& t = this->txs[tidx];
                    if(t.is_template() && val.first != EVAL_STATUS::st_no_overlap && !(q.exon_count()>1 && val.first==EVAL_STATUS::st_overlap_first_last_exon)){
                        ids.insert(gene_map[tidx]);
                    }
//                    if(t.is_template() && val.first != EVAL_STATUS::st_no_overlap){
//                        ids.insert(gene_map[tidx]);
//                    }
                    tidx+=1;
                }
                // evaluate collected IDs and assign genes or remove readthroughs

                if(ids.size()==1){ // no readthroughs here
                    // figure out if there is a previous transcript id assignment
                    q.add_attribute("reader_status","pass");
                }
                else if(ids.size()==0){ // novel gene
                    q.add_attribute("reader_status","unknown");
                }
                else{
                    std::cerr<<"skippping readthrough: "<<q.get_tid()<<"\t"<<q.get_igv_str(true)<<std::endl;
                    q.add_attribute("reader_status","fail");
                }
            }
            qidx+=1;
        }

        // lastly we can write to the output
        for(auto& t : this->txs){
            if(!t.is_template()){
                globals_local.out_reader_gtf_fp << t.get_gtf() << std::endl;
            }
        }

        return 1;
    }

    // TODO: another status - start/end coordinates overlap some reference transcript by poverlap
    //      chr7:122319330-122886500 - handle alt splicing of known readthroughs
    // what we can do is once we have a map of genes we can further refine it into a map of intervals
    // any query that is within the interval can then be included
    // this can be done by merging genes if a readthrough is identified in the reference

    // compare to transcripts to determine the type of overlap they share)
    // return percent overlap
    int evaluate(TX& tx1, TX& tx2, int& status){
        int nmbp = tx1.bpoverlap(tx2,true,false);
        bool containment = tx1.contains(tx2,true) || tx2.contains(tx1,true);
        if(nmbp==0){
            status = containment ? EVAL_STATUS::st_other_overlap : EVAL_STATUS::st_no_overlap;
            return 0;
        }

        if(tx1.exon_count()==1 && tx2.exon_count()==1){
            status = EVAL_STATUS::st_se_overlap;
        }
        else if(tx1.intron_eq(tx2)){
            status = EVAL_STATUS::st_full_intron_match;
        }
        else if(tx1.num_introns_shared(tx2)>0){
            status = EVAL_STATUS::st_partial_intron_match;
        }
        else{
            // check if only overlap between last/first exons
            if(tx1.num_overlapping_exons(tx2)==1 || tx2.num_overlapping_exons(tx1)==1){
                status = EVAL_STATUS::st_overlap_first_last_exon;
            }
            else {
                status = EVAL_STATUS::st_overlap;
            }
        }

        return nmbp;
    }

    // we need some quantifiable measure here
    int get_overlapping_txs(TX& q,std::vector<TX>& ops,bool template_only){ // returns the count of overlapping transcripts
        for(auto& t : this->txs){
            if(!template_only || t.is_template()){
                if(t.overlap(q)){

                }
            }
        }
        return ops.size();
    }

    struct intron_tx_cmp {
        bool operator()(TX a, TX b) {
            return a.intron_eq(b);
        }
    };
private:

    std::vector<TX> txs;
    std::map<int,int> t2id; // transcript ID to bundle index lookup
    std::pair<std::map<int,int>::iterator,bool> t2id_it;
    int size = 0;
    std::string seqid;
    char strand = 0;
    int start = MAX_INT;
    int end = 0;

    bool contains_ref = false; // one or more of the transcripts in the bundle are in the reference

    enum EVAL_STATUS {
        st_full_intron_match = 1,
        st_partial_intron_match = 2,
        st_overlap = 3,
        st_se_overlap = 4,
        st_other_overlap = 5, // could be contained in an intron or contains another transcript in it's intron
        st_overlap_first_last_exon = 6, // the only overlap found between transcripts is with one or both in the last/first exon
        st_no_overlap = 0
    };
};

// TODO: - exonic overlaps - only consider readthrough if >threshold

// uses gffReader to read-in all transcripts
// and then sorts them using custom rules
struct Transcriptome{
public:
    Transcriptome(const std::string& gtf_fname,bool is_templ,std::string tag=""){
        FILE* gff_file = fopen(gtf_fname.c_str(), "r");
        if (gff_file == nullptr)
        {
            std::cerr << "@ERROR::Couldn't open the GTF: " << gtf_fname << std::endl;
            exit(1);
        }
        GffReader gffReader(gff_file,false,false);
        gffReader.readAll(true);


        for(int i=0;i<gffReader.gflst.Count();++i) {
            GffObj *pGffObj = gffReader.gflst.Get(i);
            TX tmp(pGffObj,tx_vec.size(),is_templ,tag);
            tx_vec.push_back(tmp);
        }
    }
    ~Transcriptome()=default;

    void add(const std::string& gtf_fname,bool is_templ,std::string tag=""){
        FILE* gff_file = fopen(gtf_fname.c_str(), "r");
        if (gff_file == nullptr)
        {
            std::cerr << "@ERROR::Couldn't open the GTF: " << gtf_fname << std::endl;
            exit(1);
        }
        GffReader gffReader(gff_file,false,false);
        gffReader.readAll(true);

        for(int i=0;i<gffReader.gflst.Count();++i) {
            GffObj *pGffObj = gffReader.gflst.Get(i);
            TX tmp(pGffObj,tx_vec.size(),is_templ,tag);
            tx_vec.push_back(tmp);
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

// this function takes a gene or transcript id in CHESS format (PREFIX.GENEID.TRANSSCRIPTID) and stores geneid and tid in the res pair as ints
int parse_chs_id(std::string& id,std::pair<int,int>& res){
    int gstart_pos = id.find('.');
    if(gstart_pos==std::string::npos){
        std::cerr<<"gene id: "<<id<<" does not appear to be in the correct format. cannot extract numerical number"<<std::endl;
        exit(-1);
    }
    gstart_pos+=1;
    int gend_pos = id.find('.',gstart_pos);
    int gid_len = (gend_pos==std::string::npos) ? id.size()-gstart_pos : gend_pos-gstart_pos;

    std::string gene_id_number = id.substr(gstart_pos,gid_len);
    try {
        res.first = std::stoi(gene_id_number);
    }
    catch (...) {
        std::cerr << "gene id: " << id
                  << " does not appear to be in the correct format. cannot extract numerical number" << std::endl;
        exit(-1);
    }

    if(gend_pos!=std::string::npos){ // check whether gene or transcript
        int tstart_pos = gend_pos+1;
        int tid_len = id.size()-tstart_pos;
        std::string t_id_number = id.substr(tstart_pos,tid_len);
        try {
            res.second = std::stoi(t_id_number);
        }
        catch (...) {
            std::cerr << "transcript id: " << id
                      << " does not appear to be in the correct format. cannot extract numerical number" << std::endl;
            exit(-1);
        }
    }

    return 1;
}

int run(const std::vector<std::string>& known_gtf_fnames, const std::string& novel_gtf_fname,const std::string& out_fname){
    globals.out_reader_gtf_fp.open(out_fname);

    // read from all the reference streams as well as the query stream
    // form bundles (all overlapping transcripts)
    // process each bundle to identify best suitable 3' and 5' ends

    // load the reference GFF
    Transcriptome transcriptome(novel_gtf_fname,false);
    for(auto& k : known_gtf_fnames){
        transcriptome.add(k,true);
    } // TODO: need to verify that the duplicate template transcripts are not just appended to the vector but rather the tag is appended to the structure
    //         or at least they should be bundled correctly - but better yet to have them be added to the map with appropriate modifiers, etc
    transcriptome.sort();

    Bundle bundle;
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

    return 0;
}

enum Opt {  REFERENCE   = 'r',
            INPUT       = 'i',
            OUTPUT      = 'o',
            PO          = 'p'};

int main(int argc, char** argv) {

    ArgParse args("reader");
    args.add_multi_string(Opt::REFERENCE,"r","","Comma-separated list of GTF filenames with known transcripts",true);
    args.add_string(Opt::INPUT, "input", "", "Input GTF with transcripts to be corrected", true);
    args.add_string(Opt::OUTPUT, "output", "", "Output name", true);
    args.add_int(Opt::PO,"po",100,"Percent overlap to consider for marking single-exon transcripts",false);
    if(argc <= 1 || strcmp(argv[1], "--help") == 0){
        std::cerr << args.get_help() << std::endl;
        exit(1);
    }

    args.parse_args(argc, argv);

    // first create the execution string
    std::string cl = "reader ";
    for (int i = 0; i < argc; i++) {
        if (i == 0) {
            cl += argv[i];
        } else {
            cl += " ";
            cl += argv[i];
        }
    }

    globals.po_threshold = args.get_int(Opt::PO);

    // check that all files exist
    std::ifstream if_ss;
    if_ss.open(args.get_string(Opt::INPUT));
    if (!if_ss){
        std::cerr << "Input file does not exist! "<<args.get_string(Opt::INPUT)<<std::endl;
        exit(2);
    }
    if_ss.close();

    const std::vector<std::string> knowns = args.get_multi_string(Opt::REFERENCE);
    for(auto& k : knowns){
        if_ss.open(k);
        if(!if_ss){
            std::cerr<<"Reference file does not exist! "<<k<<std::endl;
            exit(2);
        }
        if_ss.close();
    }

    // run
    run(knowns,args.get_string(Opt::INPUT),args.get_string(Opt::OUTPUT));

    return 0;
}