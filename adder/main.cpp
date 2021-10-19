//
// Created by Ales Varabyou, Beril Erdogdu, Natalia Rincon on 03/30/21.
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <gclib/gff.h>
#include <fstream>
#include <gff_utils.h>
#include <set>
#include "arg_parse.h"

typedef std::vector<std::array<int,3>> CDS_CHAIN_TYPE; // start,end,phase

struct Globals{
    std::ofstream out_adder_gtf_fp;
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
    std::array<int,3> inter;
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
    TX(GffObj* tx,int idx,bool is_templ){
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
        std::pair<std::map<std::string,std::string>::iterator,bool> ait;
        if(is_templ){
            const char* gid = tx->getGeneID();
            if(!(gid==NULL || gid[0]=='\0')){
                ait = this->attrs.insert(std::make_pair("gene_id",gid));
                if(!ait.second){
                    std::cerr<<"Gene id already assigned"<<std::endl;
                    exit(-1);
                }
            }
            else{
                std::cerr<<"gene_id field missing"<<std::endl;
                exit(-1);
            }
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

    void add_attribute(std::string k, std::string v){
        this->ait = this->attrs.insert(std::make_pair(k,v));
        if(!this->ait.second){
            std::cerr<<"attribute key already exists"<<std::endl;
            exit(-1);
        }
    }

    void update_attribute(std::string k,std::string v, char delim=','){
        this->ait = this->attrs.insert(std::make_pair(k,v));
        if(!this->ait.second){ // if already exists - add to the value list with the delimiter
            this->ait.first->second = this->ait.first->second+delim+v;
        }
    }

    std::string get_attributes(){
        std::string res = "";
        for(auto& a : this->attrs){
            res+=a.first+" \""+a.second+"\";";
        }
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

    int poverlap(TX& tx,bool consider_strand = true, bool exon_by_exon = true){ // computes percent overlap between two transcripts
        int po = 0;

        if(this->get_seqid()!=tx.get_seqid()){
            return 0;
        }
        else if(consider_strand && this->get_strand()!=tx.get_strand()){
            return 0;
        }
        else{ // count number of matching bases
            int nmbp = this->num_matching_bp(tx.exons); // TODO: after fixed - need to fix adder and other methods that use this method....
            float pop = (static_cast<float>(nmbp)*100.0)/static_cast<float>(this->length());
            return static_cast<int>(pop);
        }
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

    std::array<int,3> get_exon(int pos){
        return this->exons[pos];
    }

    std::pair<int,int> get_intron(int first_exon_pos){
        if(first_exon_pos<0 || first_exon_pos>=this->exon_count()-1){
            std::cerr<<"requesting intron at non-existing position"<<std::endl;
        }
        return std::make_pair(std::get<1>(this->exons[first_exon_pos]),std::get<0>(this->exons[first_exon_pos+1]));
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
    int get_id(){
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
                    continue_prev = ends_match ? true : false;

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
                    continue_prev = ends_match ? true : false;

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

    bool overlap(int s, int e) {
        if (s>e){
            std::cerr<<s<<"\t"<<e<<std::endl;
            std::cerr<<"start>end"<<std::endl;
            exit(-1);
        }
        return (this->get_start()<=e && this->get_end()>=s);
    }

    std::string get_gtf(){
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

    bool is_template(){
        return this->is_templ;
    }

    std::string get_source(){
        return this->source;
    }

private:
    bool is_templ = false;
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
    ~Bundle(){}

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
        this->size++;
        return true;
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

    int process(Globals& globals){
        // search for intron-compatible transcripts and label accordingly
        std::string t_tid,t_gid,t_name;
        for(auto& q : this->txs){
            if(!q.is_template()){
                for(auto& t : this->txs){
                    if(t.is_template()){
                        if(q.intron_eq(t) || (q.exon_count()==1 && t.exon_count()==1 && (q.poverlap(t)>=globals.po_threshold || t.poverlap(q)>=globals.po_threshold))){ // either intron_chain_match or entirely contained or containment of single-exon transcripts in reference or query
                            t_tid = t.get_tid();
                            t_gid = t.get_attribute("gene_id");
                            t_name = t.get_attribute("adder_name");
                            if(t_tid.empty() || t_gid.empty() || t_name.empty()){
                                std::cerr<<"Required attributes not found for the matching template"<<std::endl;
                                exit(-1);
                            }

                            q.update_attribute(t_name+"_tid",t_tid);
                            q.update_attribute(t_name+"_gid",t_gid);
                        }
                    }
                }
                globals.out_adder_gtf_fp<<q.get_gtf()<<std::endl;
            }
        }
        return 0;
    }

    // we need som e quantifiable measure here
    int get_overlapping_txs(TX& q,std::vector<TX>& ops,bool template_only){ // returns the count of overlapping transcripts
        for(auto& t : this->txs){
            if(!template_only || t.is_template()){
                if(t.overlap(q.get_start(),q.get_end())){

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
    int size = 0;
    std::string seqid;
    char strand = 0;
    int start = MAX_INT;
    int end = 0;

    bool contains_ref = false; // one or more of the transcripts in the bundle are in the reference
};

// uses gffReader to read-in all transcripts
// and then sorts them using custom rules
struct Transcriptome{
public:
    Transcriptome(const std::string& gtf_fname,std::string name = "", bool is_templ=false){
        this->add(gtf_fname,name,is_templ);
    }
    ~Transcriptome()=default;

    void add(const std::string& gtf_fname,std::string name = "", bool is_templ=false){
        FILE* gff_file = fopen(gtf_fname.c_str(), "r");
        if (gff_file == nullptr){
            std::cerr << "@ERROR::Couldn't open the GTF: " << gtf_fname << std::endl;
            exit(1);
        }
        GffReader gffReader(gff_file,false,false);
        gffReader.readAll(true);

        for(int i=0;i<gffReader.gflst.Count();++i) {
            GffObj *pGffObj = gffReader.gflst.Get(i);
            TX tmp(pGffObj,tx_vec.size(),is_templ);
            if(!name.empty()){
                tmp.add_attribute("adder_name",name);
            }
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

int run(std::vector<std::pair<std::string,std::string>>& known_gtf_fnames, const std::string& novel_gtf_fname,const std::string& out_fname){
    globals.out_adder_gtf_fp.open(out_fname);

    // load the reference GFF
    Transcriptome transcriptome(novel_gtf_fname,"",false);
    for(auto& k : known_gtf_fnames){
        transcriptome.add(k.second,k.first,true);
    }
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

// TODO: add ability to compare annotations and add 

int main(int argc, char** argv) {

    ArgParse args("adder");
    args.add_multi_string(Opt::REFERENCE,"r","","Comma-separated list of GTF filenames with known transcripts. "
                                                "Each filename must be preceeded by a keyname (with ':' separator) to be used in the attributes of matching transcripts. "
                                                                   "Eg. MANE_ID:/path/to/mane.gtf",true);
    args.add_string(Opt::INPUT, "input", "", "Input GTF with transcripts to be corrected", true);
    args.add_string(Opt::OUTPUT, "output", "", "Output name", true);
    args.add_int(Opt::PO,"po",100,"Percent overlap to consider for marking single-exon transcripts",false);
    if(argc <= 1 || strcmp(argv[1], "--help") == 0){
        std::cerr << args.get_help() << std::endl;
        exit(1);
    }

    args.parse_args(argc, argv);

    // first create the execution string
    std::string cl = "adder ";
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

    const std::vector<std::string> knowns_tmp = args.get_multi_string(Opt::REFERENCE);
    std::vector<std::pair<std::string,std::string>> knowns;
    std::string key,val;
    int col_pos = 0;
    for(auto& k : knowns_tmp){
        col_pos = k.find(":");
        if(col_pos<0){
            std::cerr<<"Reference arguments provided in the incorrect format. Each filename must be preceeded by a keyname "
                       "(with ':' separator) to be used in the attributes of matching transcripts. "
                       "Eg. MANE_ID:/path/to/mane.gtf"<<std::endl;
            exit(-1);
        }
        key = k.substr(0,col_pos);
        val = k.substr(col_pos+1,k.size());
        if_ss.open(val);
        if(!if_ss){
            std::cerr<<"Reference file does not exist! "<<k<<std::endl;
            exit(2);
        }
        knowns.push_back(std::make_pair(key,val));
        if_ss.close();
    }

    // run
    run(knowns,args.get_string(Opt::INPUT),args.get_string(Opt::OUTPUT));

    return 0;
}