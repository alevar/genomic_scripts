//
// Created by beril Erdogdu and Ales Varabyou on 03/30/21.
//

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include <gclib/gff.h>
#include <fstream>
#include "arg_parse.h"

typedef std::vector<std::pair<uint,uint>> CHAIN_TYPE;

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

    int get_score(){
        return missing_start*10+missing_end*10+num_bp_extra+num_bp_missing+num_bp_outframe;
    }
    std::string get_desc(){
        return "missing_start \""+std::to_string(missing_start)+"\"; "+
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
        for(int i=0;i<tx->exons.Count();i++){
            this->exons.push_back(std::make_pair(tx->exons.Get(i)->start,tx->exons.Get(i)->end));
        }
        if(tx->hasCDS()){
            this->is_coding = true;
            this->cds_start = tx->CDstart;
            this->cds_end = tx->CDend;
        }
    }
    ~TX()=default;

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

    int get_end(){
        return std::get<1>(this->exons[-1]);
    }
    int get_start(){
        return std::get<0>(this->exons[0]);
    }
    std::string get_seqid(){
        return this->seqid;
    }
    char get_strand(){
        return this->strand;
    }
    std::string get_tid(){
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

        return score;
    }

private:
    int id = -1;
    std::string tid = "";
    std::string seqid = "";
    char strand = 0;
    CHAIN_TYPE exons;
    uint cds_start = 0;
    uint cds_end = 0;
    bool is_coding = false;
};

class Bundle{
public:
    Bundle()=default;
    ~Bundle()=default;

    bool can_add(GffObj* tx){
        if(this->size==0){
            return true;
        }
        else{
            if(std::strcmp(this->seqid.c_str(),tx->getGSeqName())==0 && this->strand==tx->strand){
                if(tx->overlap(this->start,this->end)){
                    return true;
                }
                else{
                    return false;
                }
            }
            return false;
        }
    }

    bool add_tx(GffObj* tx,uint idx){
        if(!this->can_add(tx)){
            return false;
        }
        this->txs.push_back(TX(tx,idx));
        this->seqid = tx->getGSeqName();
        this->strand = tx->strand;
        this->start = std::min(this->start,tx->start);
        this->end = std::max(this->end,tx->end);
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
        return res;
    }

    int process(std::ofstream& out_al_fp){
        // built a dict of all ORFs for searching
        std::vector<std::pair<std::string,CHAIN_TYPE>> cds_chains; // outer pair int is the TID of the transcript
        std::set<CHAIN_TYPE> cur_cds_chains; // for duplicate removal
        std::pair<std::set<CHAIN_TYPE>::iterator,bool> cit;
        for(auto& tx : this->txs){
            if(tx.has_cds()){
                CHAIN_TYPE cur_cds_chain;
                tx.build_cds_chain(cur_cds_chain);
                cit = cur_cds_chains.insert(cur_cds_chain);
                if(cit.second){ // successfully inserted
                    cds_chains.push_back(std::make_pair(tx.get_tid(),cur_cds_chain));
                }
            }
        }

        // iterate over each transcript-ORF pair to gauge compatibility - assign compatibility scores
        for(auto& tx : this->txs){
            if(cds_chains.empty()){
                out_al_fp<<tx.get_tid()<<"\t"
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
            for(auto& chain : cds_chains){
                Mods mods_res;
                tx.fit(chain.second,mods_res);
                out_al_fp<<tx.get_tid()<<"\t"
                         <<chain.first<<"\t"
                         <<tx.get_seqid()<<"\t"
                         <<tx.get_strand()<<"\t"
                         <<chain2str(chain.second)<<"\t"
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
            }
        }

        return 0;
    }
private:
    std::vector<TX> txs;
    int size = 0;
    std::string seqid = "";
    char strand = 0;
    uint start = MAX_INT;
    uint end = 0;
};

int run(const std::string& gtf_fname,const std::string& out_fname){
    std::ofstream out_al_fp(out_fname);

    // read from both GTF streams simultaneously
    // only consider CDS from the reference and exons from the novel
    // form bundles (all overlapping

    // load the reference GFF
    FILE* gff_file = fopen(gtf_fname.c_str(), "r");
    if (gff_file == nullptr)
    {
        std::cerr << "@ERROR::Couldn't open the GTF: " << gtf_fname << std::endl;
        exit(1);
    }
    GffReader gffReader(gff_file,false,false);
    gffReader.readAll();

    Bundle bundle;

    GffObj *pGffObj;
    out_al_fp<<"tid\t"
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
    for(int i=0;i<gffReader.gflst.Count();++i){
        pGffObj = gffReader.gflst.Get(i);
        if(!bundle.can_add(pGffObj)){
            bundle.process(out_al_fp);
            bundle.clear();
            bundle.add_tx(pGffObj,i);
        }
        else{
            bundle.add_tx(pGffObj,i);
        }
    }
    // handle the final bundle
    bundle.process(out_al_fp);

    // questions:
    // 1. any transcripts contain multiple CDSs? which one to chose?
    // 2. what to do with minor modifications? (no frame shifts)

    out_al_fp.close();
    return 0;
}

enum Opt {INPUT     = 'i',
          OUTPUT    = 'o'};

int main(int argc, char** argv) {

    ArgParse args("orfanage");
    args.add_string(Opt::INPUT, "input", "", "input GTF", true);
    args.add_string(Opt::OUTPUT, "output", "", "output name", true);

    if (argc <= 1 || strcmp(argv[1], "--help") == 0) {
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

    run(args.get_string(Opt::INPUT),args.get_string(Opt::OUTPUT));

    return 0;
}