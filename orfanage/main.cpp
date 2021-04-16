//
// Created by Beril Erdogdu and Ales Varabyou on 03/30/21.
//

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include <gclib/gff.h>
#include <fstream>
#include <gff_utils.h>
#include "arg_parse.h"

#define DBUF_LEN=1024;

typedef std::vector<std::pair<uint,uint>> CHAIN_TYPE;
typedef std::vector<std::tuple<uint,uint,uint>> CDS_CHAIN_TYPE; // start,end,phase

struct Mods{
    CHAIN_TYPE new_chain;
    bool missing_start = false;
    bool missing_end = false;

    CHAIN_TYPE missing;
    CHAIN_TYPE extra;

    int num_bp_extra = 0;
    int num_bp_missing = 0;
    int num_bp_match = 0;
    int num_bp_outframe = 0;
    int num_bp_inframe = 0;

    std::string orig_cds_tid = ""; // tid of the original CDS

    std::string cds_nt = ""; // nucleotide sequence
    std::string cds_aa = ""; // translated amino-acid sequence
    char expected_stop_codon = '.'; // amino acid coded by the three bases after the end of the CDS
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

struct TX{
    TX(GffObj* tx,uint idx){
        this->id = idx;
        this->tid = tx->getID();
        this->seqid = tx->getGSeqName();
        this->strand = tx->strand;
        this->source = tx->getTrackName();
        for(int i=0;i<tx->exons.Count();i++){
            this->exons.push_back(std::make_pair(tx->exons.Get(i)->start,tx->exons.Get(i)->end));
        }
        if(tx->hasCDS()){
            this->is_coding = true;
            this->cds_start = tx->CDstart;
            this->cds_end = tx->CDend;
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

    void assign_phase() {
        int cdsacc=0;
        if (this->strand=='-') { //reverse strand
            for(int i=this->cdss.size()-1;i>=0;i--){
                std::get<2>(this->cdss[i])=(3-cdsacc%3)%3;
                cdsacc+=std::get<1>(this->cdss[i])-std::get<0>(this->cdss[i])+1;
            }
        }
        else { //forward strand
            for(auto& cds : this->cdss){
                std::get<2>(cds)=(3-cdsacc%3)%3;
                cdsacc+=std::get<1>(cds)-std::get<0>(cds)+1;
            }
        }
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

    void build_cds_chain(CHAIN_TYPE& chain){
        if(!this->is_coding){
            return;
        }
        bool stop=false; // signals to break iteration since the end was found
        for(auto& e : this->exons){
            if(this->cds_start > e.second){ // skip non-coding exons
                continue;
            }
            uint cur_exon_coding_start = e.first;
            uint cur_exon_coding_end = e.second;

            if(this->cds_start > e.first){
                cur_exon_coding_start = this->cds_start;
            }

            if(this->cds_end < e.second){
                cur_exon_coding_end = this->cds_end;
                stop=true;
            }

            chain.push_back(std::make_pair(cur_exon_coding_start,cur_exon_coding_end));
            if(stop){
                break;
            }
        }
    }

    int get_end() const{
        return this->exons.back().second;
    }
    int get_start() const{
        return this->exons.front().first;
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

    int single_intersection(std::pair<uint,uint>& i1,std::pair<uint,uint>& i2,std::pair<uint,uint>& res){
        uint start = std::max(i1.first,i2.first);
        uint end = std::min(i1.second,i2.second);
        if(start<=end){
            res.first = start;
            res.second = end;
            return (end-start)+1;
        }
        res.first=NULL;
        res.second=NULL;
        return 0;
    }

    int intersection(CHAIN_TYPE& chain1,CHAIN_TYPE& chain2,CHAIN_TYPE& res){
        bool found_yet = false;
        int start = 0;
        std::pair<uint,uint> inter;
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
                    if(inter.second==chain2[j].second){
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

    int cut(uint start,uint end,CHAIN_TYPE& res){ // cuts it's own exon chain - result contains everything between the start/end
        CHAIN_TYPE start_end;
        start_end.push_back(std::make_pair(start,end));
        uint len = intersection(this->exons,start_end,res);
        return len;
    }

    int recompare(CHAIN_TYPE& chain1,CHAIN_TYPE& chain2,Mods& res){ // clear the contents of the Mods and redo comparison
        res.missing.clear();
        res.extra.clear();
        res.num_bp_match=0;
        res.num_bp_missing=0;
        res.num_bp_extra=0;
        res.num_bp_inframe=0;
        res.num_bp_outframe=0;
        res.missing_start=false;
        res.missing_end=false;

        return this->compare(chain1,chain2,res);
    }

    int compare(CHAIN_TYPE& chain1,CHAIN_TYPE& chain2,Mods& res){
        int c1_i=0,c2_i=0;
        std::pair<uint,uint> cur_c1 = chain1[c1_i];
        std::pair<uint,uint> cur_c2 = chain2[c2_i];

        std::vector<std::pair<int,int>> mods;
        std::vector<std::pair<int,int>> mods_clean;

        std::pair<uint,uint> inter;
        while(true){
            int inter_len = single_intersection(cur_c1, cur_c2, inter);

            if(inter_len<=0){
                if(cur_c2.second<cur_c1.first){
                    c2_i++;
                    int right = (cur_c2.second-cur_c2.first)+1;
                    if(right!=0){
                        res.missing.push_back(std::make_pair(cur_c2.first,cur_c2.second));
                    }
                    mods.push_back(std::make_pair(right,1));

                    if(c2_i == chain2.size()){ // done
                        int left = 0-(((int)cur_c1.second+1) - (int)cur_c1.first);
                        if(left!=0){
                            res.extra.push_back(std::make_pair(cur_c1.first,cur_c1.second));
                        }
                        for(auto cc= chain1.begin() + c1_i + 1; cc != chain1.end(); cc++){
                            left-=(((int)cc->second+1)-(int)cc->first);
                            if(left!=0){
                                res.extra.push_back(std::make_pair(cc->first,cc->second));
                            }
                        }
                        mods.push_back(std::make_pair(left,1));
                        break;
                    }
                    else{
                        cur_c2 = chain2[c2_i];
                        continue;
                    }
                }
                if(cur_c1.second<cur_c2.first){
                    c1_i++;
                    int left = 0-((cur_c1.second+1)-cur_c1.first);
                    if(left!=0){
                        res.extra.push_back(std::make_pair(cur_c1.first,cur_c1.second));
                    }
                    mods.push_back(std::make_pair(left,1));

                    if(c1_i == chain1.size()){ // done
                        int right = ((int)cur_c2.second - (int)cur_c2.first) + 1;
                        if(right!=0){
                            res.missing.push_back(std::make_pair(cur_c2.first,cur_c2.second));
                        }
                        for(auto cf= chain2.begin() + c2_i + 1; cf != chain2.end(); cf++){
                            right+=(((int)cf->second-(int)cf->first)+1);
                            if(right!=0){
                                res.missing.push_back(std::make_pair(cf->first,cf->second));
                            }
                        }
                        mods.push_back(std::make_pair(right,1));
                        break;
                    }
                    else{
                        cur_c1 = chain1[c1_i];
                        continue;
                    }
                }
            }

            int left_start = std::min((int)cur_c1.first,(int)inter.first);
            int left_end = std::max((int)cur_c1.first,(int)inter.first);
            int left = 0-(left_end - left_start);
            if(left!=0){
                res.extra.push_back(std::make_pair(left_start,left_end));
            }
            int right_start = std::min((int)cur_c2.first, (int)inter.first);
            int right_end = std::max((int)cur_c1.first, (int)inter.first);
            int right = right_end - right_start;
            if(right!=0){
                res.missing.push_back(std::make_pair(right_start,right_end));
            }
            int inters = ((int)inter.second-(int)inter.first)+1;

            if(left!=0 && right!=0){
                std::cerr<<"something wrong with operations"<<std::endl;
                exit(-1);
            }
            if(left!=0){
                mods.push_back(std::make_pair(left,1));
            }
            if(right!=0){
                mods.push_back(std::make_pair(right,1));
            }
            if(inters!=0){
                mods.push_back(std::make_pair(inters,0));
            }

            cur_c1 = std::make_pair(inter.second + 1, cur_c1.second);
            cur_c2 = std::make_pair(inter.second + 1, cur_c2.second);

            if(((int)cur_c1.second - (int)cur_c1.first) < 0){
                c1_i++;
                if(c1_i == chain1.size()){ // done
                    right = ((int)cur_c2.second - (int)cur_c2.first) + 1;
                    if(right!=0){
                        res.missing.push_back(std::make_pair(cur_c2.first,cur_c2.second));
                    }
                    for(auto cf= chain2.begin() + c2_i + 1; cf != chain2.end(); cf++){
                        right+=(((int)cf->second-(int)cf->first)+1);
                        if(right!=0){
                            res.missing.push_back(std::make_pair(cf->first,cf->second));
                        }
                    }
                    mods.push_back(std::make_pair(right,1));
                    break;
                }
                cur_c1 = chain1[c1_i];
            }

            if(((int)cur_c2.second - (int)cur_c2.first) < 0){
                c2_i++;
                if(c2_i == chain2.size()){ // done
                    left = 0-(((int)cur_c1.second+1) - (int)cur_c1.first);
                    if(left!=0){
                        res.extra.push_back(std::make_pair(cur_c1.first,cur_c1.second));
                    }
                    for(auto cc= chain1.begin() + c1_i + 1; cc != chain1.end(); cc++){
                        left-=(((int)cc->second+1)-(int)cc->first);
                        if(left!=0){
                            res.extra.push_back(std::make_pair(cc->first,cc->second));
                        }
                    }
                    mods.push_back(std::make_pair(left,1));
                    break;
                }
                cur_c2 = chain2[c2_i];
            }
        }

        for(auto& m : mods){
            if(mods_clean.size()>0 && m.second==mods_clean.back().second){
                mods_clean.back().first+=m.first;
            }
            else{
                if(m.first==0){
                    continue;
                }
                mods_clean.push_back(m);
            }
        }

        // extract modifications
        extract_mods(mods_clean,this->strand=='-',res);

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

    int fit(CHAIN_TYPE& chain,Mods& res){ // computes the intersection of own exons and the chain
        uint cut_len = cut(chain.front().first,chain.back().second,res.new_chain);
        if(cut_len==0){ // no overlap found - create a dummy chain which will never overlap anything. This way the Mods will still get populated with data and will make it possible to filter
            res.new_chain = CHAIN_TYPE{std::make_pair(0,0)};
        }

        // now we can take the cut piece and directly compare it to the desired CDS counting any changes
        int score = compare(chain,res.new_chain,res);

        if(cut_len!=0){ // no overlap found - create a dummy chain which will never overlap anything. This way the Mods will still get populated with data and will make it possible to filter
            if(this->cdss.size()>0){
                std::cerr<<"cds is already filled"<<std::endl;
                exit(-1);
            }
            for(auto& cd : res.new_chain){
                this->cdss.push_back(std::make_tuple(cd.first,cd.second,0));
            }
            assign_phase();
            this->mods_gtf_str = res.get_desc();
            this->cds_start = res.new_chain.front().first;
            this->cds_end = res.new_chain.back().second;
        }

        return score;
    }

    void clear_fitted_cds(){
        this->cdss.clear();
        this->cds_start = 0;
        this->cds_end = 0;
        this->mods_gtf_str = "";
    }

    bool overlap(uint s, uint e) {
        if (s>e){
            std::cerr<<s<<"\t"<<e<<std::endl;
            std::cerr<<"start>end"<<std::endl;
            exit(-1);
        }
        return (this->get_start()<=e && this->get_end()>=s);
    }

    std::string get_gtf(std::string ci=""){
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
                +"transcript_id \""+this->tid+ci+"\"; "
                +this->get_attributes();
        if(this->cdss.size()>0){
            gtf_str.append(this->mods_gtf_str);
        }
        gtf_str.append("\n");
        // get exon lines
        for(auto& e : this->exons){
            gtf_str+=this->seqid+"\t"
                     +this->source+"\t"
                     +"exon"+"\t"
                     +std::to_string(e.first)+"\t"
                     +std::to_string(e.second)+"\t"
                     +"."+"\t"
                     +this->strand+"\t"
                     +"."+"\t"
                     +"transcript_id \""+this->tid+ci+"\"; "
                     +this->get_attributes()+"\n";
        }

        // if CDS is fitted - get CDS coordinates based on CDS start and end
        if(this->cdss.size()>0){
            for(auto& c : this->cdss){
                gtf_str+=this->seqid+"\t"
                         +this->source+"\t"
                         +"CDS"+"\t"
                         +std::to_string(std::get<0>(c))+"\t"
                         +std::to_string(std::get<1>(c))+"\t"
                         +"."+"\t"
                         +this->strand+"\t"
                         +std::to_string(std::get<2>(c))+"\t"
                         +"transcript_id \""+this->tid+ci+"\"; "
                         +this->get_attributes()+"\n";
            }
        }

        return gtf_str;
    }

private:
    int id = -1;
    std::string tid = "";
    std::string seqid = "";
    char strand = 0;
    CHAIN_TYPE exons;
    CDS_CHAIN_TYPE cdss;
    uint cds_start = 0;
    uint cds_end = 0;
    bool is_coding = false;
    std::string source = "";
    std::string mods_gtf_str = "";

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
                CHAIN_TYPE cur_cds_chain;
                tx.build_cds_chain(cur_cds_chain);
                std::pair<std::string,CHAIN_TYPE> cds_chain = {"",cur_cds_chain};
                for(auto& c : cds_chain.second){
                    std::cout<<c.first<<"-"<<c.second<<" ; ";
                }
                std::cout<<std::endl;
            }
        }
        std::cout<<"------"<<std::endl;
    }

    std::string chain2str(CHAIN_TYPE& chain){
        std::string res = "";
        for(auto& c : chain){
            res+=std::to_string(c.first)+"-"+std::to_string(c.second)+",";
        }
        if(!res.empty()){
            res.pop_back();
        }
        if(res.empty()){
            res = "-";
        }
        return res;
    }

    std::string get_nt(CHAIN_TYPE& chain, const char* subseq){
        std::string cds_nt = "";
        if(this->strand=='+'){
            for(auto& c : chain){
                uint cs = c.first-this->start; // start coordinate of the chain with respect to the bundle start
                uint ce = c.second-this->start; // end coordinate of the chain with respect to the bundle start
                for(uint i=cs;i<=ce;i++){
                    cds_nt+=subseq[i];
                }
            }
        }
        else{ // strand=='-'
            for(int ci=chain.size()-1;ci>=0;ci--){
                uint cs = 3+(chain[ci].first-this->start); // -3 to adjust for the stop codon on -; start coordinate of the chain with respect to the bundle start
                uint ce = 3+(chain[ci].second-this->start); // -3 to adjust for the stop codon on -; end coordinate of the chain with respect to the bundle start
                for (uint i=ce;i>=cs;i--) {
                    cds_nt+=ntComplement(subseq[i]);
                }
            }
        }
        return cds_nt;
    }

    std::string get_aa(CHAIN_TYPE& chain, const char* subseq){
        std::string cds_nt = get_nt(chain,subseq);
        int cds_len = 0;
        std::string cds_aa = translateDNA(cds_nt.c_str(),cds_len,cds_nt.size());
        return cds_aa;
    }

    void evaluate_fasta(Mods& res,const char* subseq){
        int stop_coord = 0;
        std::string stop_nts = "";
        if(this->strand=='+'){
            for(auto& c : res.new_chain){
                uint cs = c.first-this->start; // start coordinate of the chain with respect to the bundle start
                uint ce = c.second-this->start; // end coordinate of the chain with respect to the bundle start
                for(uint i=cs;i<=ce;i++){
                    res.cds_nt+=subseq[i];
                }
            }
            // now add the next three bases - that's where the stop codon is expected to be
            uint stop_codon_start = res.new_chain.back().second-this->start;
            uint stop_codon_end = 3+(res.new_chain.back().second-this->start);
            for(uint i=stop_codon_start+1;i<=stop_codon_end;i++){
                stop_nts+=subseq[i];
            }
        }
        else{ // strand=='-'
            for(int ci=res.new_chain.size()-1;ci>=0;ci--){
                uint cs = 3+(res.new_chain[ci].first-this->start); // -3 to adjust for the stop codon on -; start coordinate of the chain with respect to the bundle start
                uint ce = 3+(res.new_chain[ci].second-this->start); // -3 to adjust for the stop codon on -; end coordinate of the chain with respect to the bundle start
                for (uint i=ce;i>=cs;i--) {
                    res.cds_nt+=ntComplement(subseq[i]);
                }
            }
            // now add the next three bases - that's where the stop codon is expected to be
            uint stop_codon_start = 3+(res.new_chain.front().first-this->start);
            uint stop_codon_end = stop_codon_start-3;
            for(uint i=stop_codon_start-1;i>=stop_codon_end;i--){
                stop_nts+=ntComplement(subseq[i]);
            }
        }
        int cds_len = 0;
        char* cdsaa = translateDNA(res.cds_nt.c_str(),cds_len,res.cds_nt.size());
        res.cds_aa = cdsaa;

        int stop_len = 0;
        char* stopaa = translateDNA(stop_nts.c_str(),stop_len,3);
        res.expected_stop_codon = stopaa[0];

        return;
        // TODO: also for each chain - output a separate stats file describing
        //     1. fraction of transcripts in the bundle it fits
        //     2. list of transcripts it fits

        // TODO: remove duplicate transcripts?
    }

    void adjust_stop(Mods& m){ // adjust the chain coordinates to stop at the first stop codon
        if(m.cds_aa.empty()){
            std::cerr<<"amino chain is empty"<<std::endl;
            exit(2);
        }
        size_t stop_aa_pos = m.cds_aa.find('.');
        if(stop_aa_pos!=std::string::npos){ // stop codon was found
            size_t stop_nt_pos = stop_aa_pos*3;
            size_t left_to_stop = stop_nt_pos; // number of positions left until the stop codon
            // find CDS piece in which the first nt of the stop codon was found
            if(this->strand=='+'){
                for(int i=0;i<m.new_chain.size();i++){
                    uint cs = m.new_chain[i].first-this->start; // start coordinate of the chain with respect to the bundle start
                    uint ce = m.new_chain[i].second-this->start; // end coordinate of the chain with respect to the bundle start
                    size_t clen = (ce+1)-cs;
                    if(left_to_stop<clen){ // found the cds segment with the stop codon
                        m.new_chain[i].second = (m.new_chain[i].first+left_to_stop)-1;
                        if(i<m.new_chain.size()-1){
                            m.new_chain.erase(m.new_chain.begin()+i+1,m.new_chain.end());
                        }
                        m.cds_nt.erase(m.cds_nt.begin()+stop_nt_pos,m.cds_nt.end());
                        m.cds_aa.erase(m.cds_aa.begin()+stop_aa_pos,m.cds_aa.end());
                        m.expected_stop_codon = '.'; // since we found a stop codon and adjust with respect to it - we can reset here
                        m.adjusted = true;
                        break;
                    }
                    left_to_stop-=clen;
                }
            }
            else{ // strand=='-
                for(int i=m.new_chain.size()-1;i>=0;i--){
                    uint cs = m.new_chain[i].first-this->start; // start coordinate of the chain with respect to the bundle start
                    uint ce = m.new_chain[i].second-this->start; // end coordinate of the chain with respect to the bundle start
                    size_t clen = (ce+1)-cs;
                    if(left_to_stop<clen){ // found the cds segment with the stop codon
                        m.new_chain[i].first = m.new_chain[i].second-left_to_stop; // TODO: test
                        if(i>0) {
                            uint last_to_remove = m.new_chain.size()-i-1;
                            m.new_chain.erase(m.new_chain.begin(),m.new_chain.begin()+last_to_remove); // TODO: test
                        }
//                        m.cds_nt.erase(m.cds_nt.begin(),m.cds_nt.begin()+(m.cds_nt.size()-stop_nt_pos-1));
//                        m.cds_aa.erase(m.cds_aa.begin(),m.cds_aa.begin()+(m.cds_aa.size()-stop_aa_pos-1));
                        m.cds_nt.erase(m.cds_nt.begin()+stop_nt_pos,m.cds_nt.end());
                        m.cds_aa.erase(m.cds_aa.begin()+stop_aa_pos,m.cds_aa.end());
                        m.expected_stop_codon = '.'; // since we found a stop codon and adjust with respect to it - we can reset here
                        m.adjusted = true;
                        break;
                    }
                    left_to_stop-=clen;
                }
            }

            // TODO: test that the chain is broken correctly  when the stop codon is the first in the exon/last in the exon/etc
            //    also on the + and - strands

            // update info:
            // 1. length
            // 2. stop nt
            // 3. stop aa
            // 4. extra/missing
            // 5. nt sequence
            // 6. aa sequence

            // TODO: add length of the original and fitted ORFs to the stats

            // TODO: need a check if the CDS is no longer valid (0 aa) after stop codon detection (stop codon is first aa basically)
        }
    }

    void compare_aa(Mods& ref_res, Mods& new_res){
        // check if starts with start codon - if not - check the first one
        if(ref_res.expected_stop_codon == new_res.expected_stop_codon){

        }

        // if the original chain contains an in-frame stop codon at the same position as the fitted chain - do not re-adjust

        // check if ends in stop codon - if not - check the first occurrence
        // TODO: this may require examining the 3 bases after the CDS frame - since stop codons are typically not included in the CDS
        return;
    }

    int process(std::ofstream& out_al_fp,std::ofstream& out_cds_stats_fp,std::ofstream& out_gtf_fp){
        const char* bundle_seq = NULL;
        if(this->check_ref){
            this->seqid_seq=fastaSeqGet(gfasta, this->seqid.c_str());
            if(this->seqid_seq==NULL){
                std::cerr<<"Error: no genomic sequence available (check -r option!)."<<std::endl;
            }
            int bundle_start = this->start;
            if(this->strand=='-'){
                bundle_start-=3;
            }
            int bundle_len = ((this->end+1)-this->start)+3;
            if(bundle_start+bundle_len>seqid_seq->getseqlen() || bundle_start<0){
                std::cerr<<"bundle+stop codon extends past the end of the reference sequence"<<std::endl;
                exit(2);
            }
            bundle_seq=this->seqid_seq->subseq(bundle_start,bundle_len); // +3 is to account for the last stop codon if the end of the last CDS corresponds to the end of the last exon
        }
        // built a dict of all ORFs for searching
        std::vector<Mods> cds_chains; // outer pair int is the TID of the transcript
        std::set<CHAIN_TYPE> cur_cds_chains; // for duplicate removal
        std::pair<std::set<CHAIN_TYPE>::iterator,bool> cit;
        for(auto& tx : this->txs){
            if(tx.has_cds()){
                CHAIN_TYPE cur_cds_chain;
                tx.build_cds_chain(cur_cds_chain);
                cit = cur_cds_chains.insert(cur_cds_chain);
                if(cit.second){ // successfully inserted
                    Mods m;
                    m.new_chain = cur_cds_chain;
                    m.orig_cds_tid = tx.get_tid();
                    this->evaluate_fasta(m,bundle_seq);
                    cds_chains.push_back(m);
                }
            }
        }

        // iterate over each transcript-ORF pair to gauge compatibility - assign compatibility scores
        for(auto& tx : this->txs){
            if(cds_chains.empty()){
                out_gtf_fp<<tx.get_gtf();

                out_al_fp<<tx.get_tid()<<"\t"
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
            // TODO: when a premature stop-codon is found - needs to be stored as a separate stat
            int chain_i = 0;
            for(auto& chain : cds_chains){
                Mods mods_res;
                mods_res.orig_cds_tid = chain.orig_cds_tid;
                tx.fit(chain.new_chain,mods_res);
                if(this->check_ref){
                    this->evaluate_fasta(mods_res,bundle_seq);
                    this->adjust_stop(mods_res);
                    // redo comparison
                    tx.recompare(chain.new_chain,mods_res.new_chain,mods_res);
                    std::cout<<">"<<tx.get_tid()<<std::endl;
                    std::cout<<mods_res.cds_nt<<std::endl;
                    std::cout<<mods_res.cds_aa<<std::endl;
                    std::cout<<chain2str(mods_res.new_chain)<<std::endl;
                    this->compare_aa(chain,mods_res);
                }

                // since we don't know which one is the correct CDS yet (multiple might fit)
                // we shall create duplicate of the transcript with different CDSs
                std::string out_tid_name = tx.get_tid();
                if(cds_chains.size()>1 && !mods_res.new_chain.empty()){ // TODO: needs work - could lead to duplicates
                    out_gtf_fp<<tx.get_gtf(":-:"+std::to_string(chain_i));
                    out_tid_name+=":-:"+std::to_string(chain_i);
                    chain_i++;
                }
                else{
                    out_gtf_fp<<tx.get_gtf();
                }

                out_al_fp<<tx.get_tid()<<"\t"
                         <<out_tid_name<<"\t"
                         <<chain.orig_cds_tid<<"\t"
                         <<tx.get_seqid()<<"\t"
                         <<tx.get_strand()<<"\t"
                         <<chain2str(chain.new_chain)<<"\t"
                         <<chain2str(mods_res.new_chain)<<"\t"
                         <<mods_res.get_score()<<"\t"
                         <<mods_res.missing_start<<"\t"
                         <<mods_res.missing_end<<"\t"
                         <<mods_res.num_bp_outframe<<"\t"
                         <<mods_res.num_bp_inframe<<"\t"
                         <<mods_res.num_bp_extra<<"\t"
                         <<mods_res.num_bp_missing<<"\t"
                         <<mods_res.num_bp_match<<"\t"
                         <<chain2str(mods_res.missing)<<"\t"
                         <<chain2str(mods_res.extra)<<std::endl;
                tx.clear_fitted_cds(); // prepare for the new fitting since we wrote everything we needed
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

int run(const std::string& gtf_fname,const std::string& out_fname,const std::string& ref_fa_fname){
    std::ofstream out_stats_fp(out_fname+".stats");
    std::ofstream out_cds_stats_fp(out_fname+".cds");
    std::ofstream out_gtf_fp(out_fname+".gtf");

    // read from both GTF streams simultaneously
    // only consider CDS from the reference and exons from the novel
    // form bundles (all overlapping

    // load the reference GFF
    Transcriptome transcriptome(gtf_fname); // TODO: if transcript has a cds - just report it unless better CDS available
    transcriptome.sort();

    Bundle bundle;

    bool check_ref = !ref_fa_fname.empty();
    if(check_ref){
        bundle.set_ref(ref_fa_fname);
    }

    out_stats_fp<<"tid\t"
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

    out_cds_stats_fp<<"cds_tid\t"
                      "chain\t"
                      "overlapping_txs\t"
                      "perfect_txs\t"
                      "num_overlapping_txs\t"
                      "num_perfect_fit_txs\t"
                      "start\t"
                      "stop\t"
                      "premature_stop\t"<<std::endl;

    for(auto& t : transcriptome){ // TODO: now that we are adding transcripts to the custom transcriptome anyways - we can easily process separate GTFs (one input novel and a list of knowns)
        if(!bundle.can_add(t)){
            bundle.process(out_stats_fp,out_cds_stats_fp,out_gtf_fp);
            bundle.clear();
            bundle.add_tx(t);
        }
        else{
            bundle.add_tx(t);
        }
    }
    // handle the final bundle
    bundle.process(out_stats_fp,out_cds_stats_fp,out_gtf_fp);

    // questions:
    // 1. any transcripts contain multiple CDSs? which one to chose?
    // 2. what to do with minor modifications? (no frame shifts)

    out_stats_fp.close();
    out_cds_stats_fp.close();
    out_gtf_fp.close();
    return 0;
}

enum Opt {INPUT     = 'i',
          OUTPUT    = 'o',
          REFERENCE = 'r',
          CLEANREF  = 'c',
          FILTER    = 'f'};

int main(int argc, char** argv) {

    ArgParse args("orfanage");
    args.add_string(Opt::INPUT, "input", "", "Input GTF", true);
    args.add_string(Opt::OUTPUT, "output", "", "Output name", true);
    args.add_string(Opt::REFERENCE,"reference","","Reference fasta",false);
    args.add_flag(Opt::CLEANREF,"cleanref","Remove transcripts which contain mistakes in pre-annotated CDS",false);
    args.add_flag(Opt::FILTER,"filter","Select the best fitting CDS where possible",false);

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
        std::cerr << "Input file does not exist!\n";
        exit(2);
    }
    if_ss.close();

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

    // run
    run(args.get_string(Opt::INPUT),args.get_string(Opt::OUTPUT),args.get_string(Opt::REFERENCE));

    return 0;
}