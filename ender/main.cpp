//
// Created by Ales Varabyou, Beril Erdogdu, Natalia Rincon on 03/30/21.
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>

#include <gclib/gff.h>
#include <fstream>
#include <gff_utils.h>
#include "arg_parse.h"

typedef std::vector<std::array<int,3>> CDS_CHAIN_TYPE; // start,end,phase

struct Globals{
    bool extra = false;
    bool equal = false;
    bool use_cds = false;
    bool adjust_single = false;
    int extra_thresh = MAX_INT;

    std::ofstream out_ender_gtf_fp;
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

    std::string get_attributes(){
        std::string res = "";
        for(auto& a : this->attrs){
            res+=a.first+" \""+a.second+"\";";
        }
        return res;
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
        for(int i=0;i<=this->intron_count();i++){ // -1 for indexing
            if(std::get<1>(this->exons[i])!=std::get<1>(tx.exons[i]) ||
                    std::get<0>(this->exons[i+1])!=std::get<0>(tx.exons[i+1])){
                return false;
            }
        }
        return true;
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
    int set_end(int new_end){
        if(new_end<std::get<0>(this->exons[this->exons.size()-1])){ // if smaller than the start of last exon
            std::cerr<<"incompatible coordinate for the new end"<<std::endl;
            exit(-1);
        }
        std::get<1>(this->exons[this->exons.size()-1]) = new_end;
    }
    int set_start(int new_start){
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
        std::array<int,3> inter;
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
                 +this->get_attributes()+
                 +"ender_adjustment \""+std::to_string(this->adjustment)+"\"; ";
        if(this->is_adjusted){
            if(this->adjusted_start){
                gtf_str+="ender_start_source \""+this->adjusted_start_source+"\"; ";
            }
            if(this->adjusted_end){
                gtf_str+="ender_end_source \""+this->adjusted_end_source+"\"; ";
            }
        }
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
                     +"transcript_id \""+tid+"\"; "
                     +this->get_attributes()+"\n";
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
                         +"transcript_id \""+tid+"\"; "
                         +this->get_attributes()+"\n";
            }
        }
        gtf_str.pop_back();
        return gtf_str;
    }

    bool is_template(){
        return this->is_templ;
    }

    void set_is_adjusted(){
        this->is_adjusted = true;
    }
    void set_adjustment(int val){
        this->adjustment = val;
    }
    void set_adjustment_start_source(std::string val){
        this->adjusted_start = true;
        this->adjusted_start_source = val;
    }
    void set_adjustment_end_source(std::string val){
        this->adjusted_end = true;
        this->adjusted_end_source = val;
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

    bool is_adjusted = false;
    int adjustment = -1;
    bool adjusted_start = false;
    bool adjusted_end = false;
    std::string adjusted_start_source = "";
    std::string adjusted_end_source = "";

    std::map<std::string,std::string> attrs;
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

    bool has_ref(){
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
        // do any transcripts match the reference?
        TX ref_tx = TX();
        TX nov_tx = TX();

        std::array<std::vector<TX>,5> matches;
        for (auto& nov_tx : this->txs) { // TODO: need to handle single-exon novel transcripts separately
            matches.fill(std::vector<TX>());
//            if(std::strcmp(nov_tx.get_tid().c_str(),"ALL_15973932")==0){
//                std::cout<<"found"<<std::endl;
//            }

            if(!nov_tx.is_template()) {
                if(!globals.adjust_single && nov_tx.exon_count()==1){ // if not explicitely requested - skip single-exon query transcripts
                    continue;
                }
                for (auto& ref_tx : this->txs) {
                    if(ref_tx.is_template()){
                        // put everything into the ELSE container for a fall-back (can examine the same transcript for both ends)
                        matches[EQST::ELSE].push_back(ref_tx);
                        int equality_status = evaluate_chains(ref_tx,nov_tx);
                        switch (equality_status) {
                            case EQST::EQ: // found full intron-chain match
                                matches[equality_status].push_back(ref_tx);
                                break;
                            case EQST::EQ1m1: // found matching first and last introns
                                matches[equality_status].push_back(ref_tx);
                                break;
                            case EQST::EQ1: // found matching by first intron
                                matches[equality_status].push_back(ref_tx);
                                break;
                            case EQST::EQm1: // found matching by last intron
                                matches[equality_status].push_back(ref_tx);
                                break;
//                            case EQST::ELSE: // add everything else to a separate container - to be evaluated based on closest end in case of no intron matches
//                                matches[equality_status].push_back(ref_tx);
                            default:
                                break;
                        }
                    }
                }

                // after the transcripts have been sorted - pick the one with the closest end without violating CDS
                int res = correct_ends(matches,nov_tx);
                if(res==-1){
                    std::cerr<<"could not perform correction of the transcript: "<<nov_tx.get_tid()<<std::endl;
                }
                globals.out_ender_gtf_fp<<nov_tx.get_gtf()<<std::endl;
            }
        }

        return 0;
    }

    bool can_trim(TX& tx,int start,int end,bool check_thresh){ // checks CDS (and potentially other features to assert the trimming wont disrupt anything
        if(check_thresh){
            bool thresh_check = std::abs(start-std::get<0>(tx.get_exon(0)))<=globals.extra_thresh && std::abs(end>=std::get<1>(tx.get_exon(tx.exon_count()-1)))<=globals.extra_thresh;
            if(!thresh_check){
                return false;
            }
        }
        bool coord_check = start<=std::get<1>(tx.get_exon(0)) && end>=std::get<0>(tx.get_exon(tx.exon_count()-1));
        if(!globals.use_cds || !tx.has_cds()){
            return coord_check;
        }
        return coord_check && start<=tx.get_cds_start() && end>=tx.get_cds_end();
    }

    int _correct_start(std::vector<TX>& matches,TX& nov_tx,bool check_thresh){
        bool trimmed = false;

        // begin by sorting the elements by their proximity to the nov_tx start
        if(!matches.empty()){ // found perfect match - let's trim
            sort(matches.begin( ), matches.end( ), [nov_tx]( const TX& lhs, const TX& rhs){
                return std::abs(lhs.get_start()-nov_tx.get_start()) < std::abs(rhs.get_start()-nov_tx.get_start());
            });
            // now that they are sorted accordingly - we can iterate through to get the closest one and check for any assertions
            for(auto& m: matches){
                if(!trimmed && can_trim(nov_tx,m.get_start(),nov_tx.get_end(),check_thresh)){
                    nov_tx.set_start(m.get_start());
                    trimmed = true;
                    nov_tx.set_is_adjusted();
                    nov_tx.set_adjustment_start_source(m.get_tid());
                }
            }
        }
        return trimmed;
    }

    int _correct_end(std::vector<TX>& matches,TX& nov_tx,bool check_thresh){
        bool trimmed = false;

        // begin by sorting the elements by their proximity to the nov_tx start
        if(!matches.empty()){ // found perfect match - let's trim
            sort(matches.begin( ), matches.end( ), [nov_tx]( const TX& lhs, const TX& rhs){
                return std::abs(lhs.get_end()-nov_tx.get_end()) < std::abs(rhs.get_end()-nov_tx.get_end());
            });
            // now that they are sorted accordingly - we can iterate through to get the closest one and check for any assertions
            for(auto& m: matches){
                if(!trimmed && can_trim(nov_tx,nov_tx.get_start(),m.get_end(),check_thresh)){
                    nov_tx.set_end(m.get_end());
                    trimmed = true;
                    nov_tx.set_is_adjusted();
                    nov_tx.set_adjustment_end_source(m.get_tid());
                }
            }
        }
        return trimmed;
    }

    int correct_ends(std::array<std::vector<TX>,5>& matches,TX& nov_tx){

        bool trimmed_start = _correct_start(std::get<EQST::EQ>(matches),nov_tx,false);
        if(!trimmed_start){
            trimmed_start = _correct_start(std::get<EQST::EQ1m1>(matches),nov_tx,false);
        }
        if(!trimmed_start){
            trimmed_start = _correct_start(std::get<EQST::EQ1>(matches),nov_tx,false);
        }
        if(!trimmed_start){
            trimmed_start = _correct_start(std::get<EQST::ELSE>(matches),nov_tx,true);
        }


        // now the same for the end
        bool trimmed_end = _correct_end(std::get<EQST::EQ>(matches),nov_tx,false);
        if(!trimmed_end){
            trimmed_end = _correct_end(std::get<EQST::EQ1m1>(matches),nov_tx,false);
        }
        if(!trimmed_end){
            trimmed_end = _correct_end(std::get<EQST::EQm1>(matches),nov_tx,false);
        }
        if(!trimmed_end) {
            trimmed_end = _correct_end(std::get<EQST::ELSE>(matches), nov_tx,true);
        }


        int status = -1;
        if(trimmed_start && trimmed_end){
            status = EQST::EQ;
        }
        else if(trimmed_start){
            status = EQST::EQ1;
        }
        else if(trimmed_end){
            status = EQST::EQm1;
        }
        else{
            status = -1;
        }

        nov_tx.set_adjustment(status);

        return status;
    }

    // evaluate
    int evaluate_chains(TX& tx1,TX& tx2){
        if(tx1.intron_count()==0 || tx2.intron_count()==0){
            return EQST::ELSE;
        }

        if(tx1.intron_eq(tx2)){
            return EQST::EQ;
        }

        bool eq1 = tx1.get_intron(0)==tx2.get_intron(0);
        bool eqm1 = tx1.get_intron(tx1.intron_count()-1)==tx2.get_intron(tx2.intron_count()-1); // -1 for 0-based
        if(eq1 && eqm1){
            return EQST::EQ1m1;
        }
        else if(eq1){
            return EQST::EQ1;
        }
        else if(eqm1){
            return EQST::EQm1;
        }
        else{
            return EQST::ELSE; // everything else goes here
        }
    }
private:
    enum EQST {  EQ, // full intron chain equality
                 EQ1m1, // first and last intron equality
                 EQ1, // first intron equal
                 EQm1, // last intron equal
                 ELSE, // everything else goes here
                 };

    std::vector<TX> txs;
    int size = 0;
    std::string seqid = "";
    char strand = 0;
    int start = MAX_INT;
    int end = 0;

    bool contains_ref = false; // one or more of the transcripts in the bundle are in the reference
};

// uses gffReader to read-in all transcripts
// and then sorts them using custom rules
struct Transcriptome{
public:
    Transcriptome(const std::string& gtf_fname,bool is_templ){
        FILE* gff_file = fopen(gtf_fname.c_str(), "r");
        if (gff_file == nullptr)
        {
            std::cerr << "@ERROR::Couldn't open the GTF: " << gtf_fname << std::endl;
            exit(1);
        }
        GffReader gffReader(gff_file,false,false);
        gffReader.readAll();


        for(int i=0;i<gffReader.gflst.Count();++i) {
            GffObj *pGffObj = gffReader.gflst.Get(i);
            TX tmp(pGffObj,tx_vec.size(),is_templ);
            tx_vec.push_back(tmp);
        }
    }
    ~Transcriptome()=default;

    void add(const std::string& gtf_fname,bool is_templ){
        FILE* gff_file = fopen(gtf_fname.c_str(), "r");
        if (gff_file == nullptr)
        {
            std::cerr << "@ERROR::Couldn't open the GTF: " << gtf_fname << std::endl;
            exit(1);
        }
        GffReader gffReader(gff_file,false,false);
        gffReader.readAll();

        for(int i=0;i<gffReader.gflst.Count();++i) {
            GffObj *pGffObj = gffReader.gflst.Get(i);
            TX tmp(pGffObj,tx_vec.size(),is_templ);
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

int run(std::vector<std::string> known_gtf_fnames, std::string novel_gtf_fname,std::string out_fname){
    globals.out_ender_gtf_fp.open(out_fname+".ender.gtf");

    // read from all the reference streams as well as the query stream
    // form bundles (all overlapping transcripts)
    // process each bundle to identify best suitable 3' and 5' ends

    // load the reference GFF
    Transcriptome transcriptome(novel_gtf_fname,false);
    for(auto& k : known_gtf_fnames){
        transcriptome.add(k,true);
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
            CDS         = 'c',
            INPUT       = 'i',
            OUTPUT      = 'o',
            EXTRA       = 'e',
            EQUAL       = 'q',
            SINGLE      = 's',
            EXTRA_THRESH= 't'};

int main(int argc, char** argv) {

    ArgParse args("ender");
    args.add_multi_string(Opt::REFERENCE,"r","","Comma-separated list of GTF filenames with known transcripts",true);
    args.add_string(Opt::INPUT, "input", "", "Input GTF with transcripts to be corrected", true);
    args.add_string(Opt::OUTPUT, "output", "", "Output name", true);
    args.add_flag(Opt::CDS,"cds","use CDS features in the reference set",false);
    args.add_flag(Opt::EXTRA,"extra","perform adjustment based on additional data provided (chip-seq, etc)",false);
    args.add_flag(Opt::EQUAL, "equal","only adjust the coordinates of the novel transcripts which match intron chain of any reference transcripts",false);
    args.add_flag(Opt::SINGLE,"single","perform adjustment of the single-exon query transcripts based on the closest 3' and 5' ends from both single and multi-exon transcripts",false);
    args.add_int(Opt::EXTRA_THRESH,"thresh",MAX_INT,"the maximum distance between ends when performing correction of transcripts with no intron match to the reference",false);

    if(argc <= 1 || strcmp(argv[1], "--help") == 0){
        std::cerr << args.get_help() << std::endl;
        exit(1);
    }

    args.parse_args(argc, argv);

    // first create the execution string
    std::string cl = "ender ";
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

    const std::vector<std::string> knowns = args.get_multi_string(Opt::REFERENCE);
    for(auto& k : knowns){
        if_ss.open(k);
        if(!if_ss){
            std::cerr<<"Reference file does not exist! "<<k<<std::endl;
            exit(2);
        }
        if_ss.close();
    }

    globals.extra = args.get_flag(Opt::EXTRA);
    globals.equal = args.get_flag(Opt::EQUAL);
    globals.use_cds = args.get_flag(Opt::CDS);
    globals.adjust_single = args.get_flag(Opt::SINGLE);
    globals.extra_thresh = args.get_int(Opt::EXTRA_THRESH);

    // run
    run(knowns,args.get_string(Opt::INPUT),args.get_string(Opt::OUTPUT));

    return 0;
}