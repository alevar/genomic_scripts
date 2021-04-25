//
// Created by Ales Varabyou, Beril Erdogdu, Natalia Rincon on 03/30/21.
//

#include <cstdlib>
//#include <cstring>
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
    std::ofstream out_stats_fp;
    std::ofstream out_cds_stats_fp;
    std::ofstream out_gtf_fp;
    std::ofstream out_gtf_perfect_fp;
    std::ofstream out_gtf_imperfect_fp;
    std::ofstream out_gtf_nonoverlap_fp;
} globals;

std::map<std::string,char> codon_map = {{"AAA",'K'},{"AAC",'N'},{"AAG",'K'},{"AAR",'K'},{"AAT",'N'},{"AAY",'N'},{"ACA",'T'},{"ACB",'T'},{"ACC",'T'},{"ACD",'T'},{"ACG",'T'},{"ACH",'T'},{"ACK",'T'},{"ACM",'T'},{"ACN",'T'},{"ACR",'T'},{"ACS",'T'},{"ACT",'T'},{"ACV",'T'},{"ACW",'T'},{"ACY",'T'},{"AGA",'R'},{"AGC",'S'},{"AGG",'R'},{"AGR",'R'},{"AGT",'S'},{"AGY",'S'},{"ATA",'I'},{"ATC",'I'},{"ATG",'M'},{"ATH",'I'},{"ATM",'I'},{"ATT",'I'},{"ATW",'I'},{"ATY",'I'},{"CAA",'Q'},{"CAC",'H'},{"CAG",'Q'},{"CAR",'Q'},{"CAT",'H'},{"CAY",'H'},{"CCA",'P'},{"CCB",'P'},{"CCC",'P'},{"CCD",'P'},{"CCG",'P'},{"CCH",'P'},{"CCK",'P'},{"CCM",'P'},{"CCN",'P'},{"CCR",'P'},{"CCS",'P'},{"CCT",'P'},{"CCV",'P'},{"CCW",'P'},{"CCY",'P'},{"CGA",'R'},{"CGB",'R'},{"CGC",'R'},{"CGD",'R'},{"CGG",'R'},{"CGH",'R'},{"CGK",'R'},{"CGM",'R'},{"CGN",'R'},{"CGR",'R'},{"CGS",'R'},{"CGT",'R'},{"CGV",'R'},{"CGW",'R'},{"CGY",'R'},{"CTA",'L'},{"CTB",'L'},{"CTC",'L'},{"CTD",'L'},{"CTG",'L'},{"CTH",'L'},{"CTK",'L'},{"CTM",'L'},{"CTN",'L'},{"CTR",'L'},{"CTS",'L'},{"CTT",'L'},{"CTV",'L'},{"CTW",'L'},{"CTY",'L'},{"GAA",'E'},{"GAC",'D'},{"GAG",'E'},{"GAR",'E'},{"GAT",'D'},{"GAY",'D'},{"GCA",'A'},{"GCB",'A'},{"GCC",'A'},{"GCD",'A'},{"GCG",'A'},{"GCH",'A'},{"GCK",'A'},{"GCM",'A'},{"GCN",'A'},{"GCR",'A'},{"GCS",'A'},{"GCT",'A'},{"GCV",'A'},{"GCW",'A'},{"GCY",'A'},{"GGA",'G'},{"GGB",'G'},{"GGC",'G'},{"GGD",'G'},{"GGG",'G'},{"GGH",'G'},{"GGK",'G'},{"GGM",'G'},{"GGN",'G'},{"GGR",'G'},{"GGS",'G'},{"GGT",'G'},{"GGV",'G'},{"GGW",'G'},{"GGY",'G'},{"GTA",'V'},{"GTB",'V'},{"GTC",'V'},{"GTD",'V'},{"GTG",'V'},{"GTH",'V'},{"GTK",'V'},{"GTM",'V'},{"GTN",'V'},{"GTR",'V'},{"GTS",'V'},{"GTT",'V'},{"GTV",'V'},{"GTW",'V'},{"GTY",'V'},{"MGA",'R'},{"MGG",'R'},{"MGR",'R'},{"NNN",'X'},{"RAY",'B'},{"SAR",'Z'},{"TAA",'.'},{"TAC",'Y'},{"TAG",'.'},{"TAR",'.'},{"TAT",'Y'},{"TAY",'Y'},{"TCA",'S'},{"TCB",'S'},{"TCC",'S'},{"TCD",'S'},{"TCG",'S'},{"TCH",'S'},{"TCK",'S'},{"TCM",'S'},{"TCN",'S'},{"TCR",'S'},{"TCS",'S'},{"TCT",'S'},{"TCV",'S'},{"TCW",'S'},{"TCY",'S'},{"TGA",'.'},{"TGC",'C'},{"TGG",'W'},{"TGT",'C'},{"TGY",'C'},{"TRA",'.'},{"TTA",'L'},{"TTC",'F'},{"TTG",'L'},{"TTR",'L'},{"TTT",'F'},{"TTY",'F'},{"XXX",'X'},{"YTA",'L'},{"YTG",'L'},{"YTR",'L'}};

std::string translate(std::string& nts) {
    if(nts.size()%3){
        std::cerr<<"cannot translate sequence of length%3!=0"<<std::endl;
        exit(-1);
    }
    if(nts.size()==0){
        std::cerr<<"cannot translate empty sequence"<<std::endl;
        exit(-1);
    }

    std::string res;
    std::string cur_nts;
    for(int i=0;i+2<nts.size();i+=3) {
        for(auto& nt : nts.substr(i,3)){
            cur_nts+=toupper(nt);
        }
        res+=codon_map[cur_nts];
        cur_nts.clear();
    }
    return res;
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

    std::string orig_cds_tid; // tid of the original CDS

    std::string cds_nt; // nucleotide sequence
    std::string cds_aa; // translated amino-acid sequence
    char next_codon = '-';

    std::tuple<std::string,char,std::vector<int>> start_codon;

    bool adjusted = false;

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

bool check_valid_aa(Mods& m){
    if(m.cds_nt.size()%3!=0){
        return false;
    }
    if(m.cds_aa.front()!='M'){
        return false;
    }
    if(m.next_codon!='.'){
        return false;
    }
    size_t stop_aa_pos = m.cds_aa.find('.');
    if(stop_aa_pos!=std::string::npos){ // stop codon was found
        return false;
    }
    return true;
}

// return nts and the coordinates
std::pair<std::string,std::vector<int>> next_codon_nt_nonrevcomp(int last_position,Mods& m,CDS_CHAIN_TYPE& exons,const char* subseq,int bundle_start,int start_offset,int subseq_len,bool forward){
    std::string nts = "";
    std::vector<int> coords;
    int cp = last_position;
    if(forward){
//        int cs = std::get<0>(m.new_chain.back());
//        int ce = std::get<1>(m.new_chain.back());
        bool found_last_coding_exon_start = false;
        for(auto& e : exons){
            int es = std::get<0>(e);
            int ee = std::get<1>(e);
            if(cp>=es && cp<=ee){
                found_last_coding_exon_start = true;
                cp+=1; // step 1 base past the coding portion of the orf
            }

            if(found_last_coding_exon_start){
                cp=std::max(cp,es); // will either be the same or cp
                for(int i=cp;i<=ee;i++){
                    if(i-bundle_start>=subseq_len){ // exceeds sequence length
                        return std::make_pair(nts,coords);
                    }

                    coords.push_back(i);
                    nts+=subseq[i-bundle_start];
                    if(nts.size()==3){
                        break;
                    }
                }
                if(nts.size()==3){
                    break;
                }
            }
        }
        if(!found_last_coding_exon_start){
            std::cerr<<"last coding exon not found"<<std::endl;
            exit(-1);
        }
    }
    else{
//        int cs = std::get<0>(m.new_chain[0]);
//        int ce = std::get<1>(m.new_chain[0]);
        bool found_last_coding_exon_start = false;
        for(int eidx=exons.size()-1;eidx>=0;eidx--){
            int es = std::get<0>(exons[eidx]);
            int ee = std::get<1>(exons[eidx]);
            if(cp >= es && cp <= ee){
                found_last_coding_exon_start = true;
                cp-=1;
            }
            if(found_last_coding_exon_start){
                cp=std::min(cp, ee);
                for(int i=cp; i >= es; i--){
                    if(start_offset+(i-bundle_start)<0){
                        return std::make_pair(nts,coords);
                    }
                    coords.insert(coords.begin(),i);
                    nts.insert(nts.begin(),subseq[start_offset+(i-bundle_start)]);
                    if(nts.size()==3){
                        break;
                    }
                }
                if(nts.size()==3){
                    break;
                }
            }
        }
        if(!found_last_coding_exon_start){
            std::cerr<<"last coding exon not found"<<std::endl;
            exit(-1);
        }
    }
    return std::make_pair(nts,coords);
}

std::tuple<std::string,char,std::vector<int>> next_codon(Mods& m,CDS_CHAIN_TYPE& exons,const char* subseq,int bundle_start,char strand,int start_offset,int subseq_len){ // extracts the codon that follows the annotated orf - should be STOP ('.')
    int last_position = strand=='+' ? std::get<1>(m.new_chain.back()) : std::get<0>(m.new_chain[0]);
    std::pair<std::string,std::vector<int>> nts_ncs = next_codon_nt_nonrevcomp(last_position,m,exons,subseq,bundle_start,start_offset,subseq_len,strand=='+');
    std::string nc;
    if(nts_ncs.first.size()!=3){
        if(strand=='+'){
            return std::make_tuple(nts_ncs.first,'-',nts_ncs.second);
        }
        else{
            std::string nts_revcomp = "";
            for(int i=nts_ncs.first.size()-1;i>=0;i--){
                nts_revcomp+=ntComplement(nts_ncs.first[i]);
            }
            return std::make_tuple(nts_revcomp,'-',nts_ncs.second);
        }
    }

    if(strand=='+'){
        if(nts_ncs.first.size()!=3){
            std::cout<<"found in next"<<std::endl;
        }

        nc = translate(nts_ncs.first);
        return std::make_tuple(nts_ncs.first,nc[0],nts_ncs.second);
    }
    else{
        std::string nts_revcomp = "";
        for(int i=nts_ncs.first.size()-1;i>=0;i--){
            nts_revcomp+=ntComplement(nts_ncs.first[i]);
        }
        if(nts_revcomp.size()!=3){
            std::cout<<"found in next"<<std::endl;
        }
        nc = translate(nts_revcomp);
        return std::make_tuple(nts_revcomp,nc[0],nts_ncs.second);
    }
}

std::tuple<std::string,char,std::vector<int>> prev_codon(Mods& m,CDS_CHAIN_TYPE& exons,const char* subseq,int bundle_start,char strand,int start_offset,int subseq_len){ // extracts the codon that follows the annotated orf - should be STOP ('.')
    int last_position = strand=='-' ? std::get<1>(m.new_chain.back()) : std::get<0>(m.new_chain[0]);
    std::pair<std::string,std::vector<int>> nts_ncs = next_codon_nt_nonrevcomp(last_position,m,exons,subseq,bundle_start,start_offset,subseq_len,strand=='-');
    std::string pc;
    if(nts_ncs.first.size()!=3){
        if(strand=='+'){
            return std::make_tuple(nts_ncs.first,'-',nts_ncs.second);
        }
        else{
            std::string nts_revcomp = "";
            for(int i=nts_ncs.first.size()-1;i>0;i--){
                nts_revcomp+=ntComplement(nts_ncs.first[i]);
            }
            return std::make_tuple(nts_revcomp,'-',nts_ncs.second);
        }
    }

    if(strand=='+'){
        if(nts_ncs.first.size()!=3){
            std::cout<<"found in prev"<<std::endl;
        }
        pc = translate(nts_ncs.first);
        return std::make_tuple(nts_ncs.first,pc[0],nts_ncs.second);
    }
    else{
        std::string nts_revcomp = "";
        for(int i=nts_ncs.first.size()-1;i>=0;i--){
            nts_revcomp+=ntComplement(nts_ncs.first[i]);
        }
        if(nts_revcomp.size()!=3){
            std::cout<<"found in prev"<<std::endl;
        }
        pc = translate(nts_revcomp);
        return std::make_tuple(nts_revcomp,pc[0],nts_ncs.second);
    }
}

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

bool extend_chain(CDS_CHAIN_TYPE& exons, CDS_CHAIN_TYPE& cds,int num_bp,bool forward){
    int num_extended = 0;

    bool found_last_coding_exon_start = false;
    bool next_exon = false;
    if(forward){
        int cs = std::get<0>(cds.back());
        int ce = std::get<1>(cds.back());
        for(auto& e : exons){
            int es = std::get<0>(e);
            int ee = std::get<1>(e);
            if(cs>=es && ce<=ee){
                found_last_coding_exon_start = true;
                ce+=1; // step 1 base past the coding portion of the orf
            }

            if(found_last_coding_exon_start){
                ce=std::max(ce,es); // will either be the same or ce
                if(next_exon){
                    cds.push_back(std::array<int,3>({ce,ce-1,0}));
                }
                for(int i=ce;i<=ee;i++){
                    std::get<1>(cds.back())+=1;
                    num_extended+=1;
                    if(num_bp==num_extended){
                        break;
                    }
                }
                if(num_bp==num_extended){
                    break;
                }
                next_exon = true;
            }
        }
        if(!found_last_coding_exon_start){
            std::cerr<<"last coding exon not found"<<std::endl;
            exit(-1);
        }
    }
    else{
        int cs = std::get<0>(cds[0]);
        int ce = std::get<1>(cds[0]);
        for(int eidx=exons.size()-1;eidx>=0;eidx--){
            int es = std::get<0>(exons[eidx]);
            int ee = std::get<1>(exons[eidx]);
            if(cs >= es && ce <= ee){
                found_last_coding_exon_start = true;
                cs-=1;
            }
            if(found_last_coding_exon_start){
                cs=std::min(cs, ee);
                if(next_exon){
                    cds.insert(cds.begin(),std::array<int,3>({cs+1,cs,0}));
                }
                for(int i=cs; i >= es; i--){
                    std::get<0>(cds.front())-=1;
                    num_extended+=1;
                    if(num_bp==num_extended){
                        break;
                    }
                }
                if(num_bp==num_extended){
                    break;
                }
                next_exon=true;
            }
        }
        if(!found_last_coding_exon_start){
            std::cerr<<"last coding exon not found"<<std::endl;
            exit(-1);
        }
    }
    return num_extended==num_bp;
}

// searches downstream of the CDS for the next stop codon in the same frame
void extend_to_stop(Mods& m,CDS_CHAIN_TYPE& exons,const char* subseq,int bundle_start,char strand,int start_offset,int subseq_len){
    std::tuple<std::string,char,std::vector<int>> nt_na_nc;
    Mods tmp = m; // temp copy
    bool found_stop = false;
    bool res;
    while(true){
        nt_na_nc = next_codon(tmp,exons,subseq,bundle_start,strand,start_offset,subseq_len);
        if(std::get<1>(nt_na_nc)=='-'){ // no longer available
            break;
        }
        else if(std::get<1>(nt_na_nc)=='.'){
            found_stop = true;
            tmp.next_codon = '.';
            break;
        }
        else{ // valid aa found - store in tmp
            if(strand=='+'){
                res = extend_chain(exons,tmp.new_chain,3,true);
            }
            else{
                res = extend_chain(exons,tmp.new_chain,3,false);
            }
            if(!res){
                std::cerr<<"couldn't extend - something wrong"<<std::endl;
                exit(1);
            }
            tmp.cds_aa+=std::get<1>(nt_na_nc);
            tmp.cds_nt+=std::get<0>(nt_na_nc);
        }
    }
    if(found_stop){
        m = tmp;
    }
}

// searches upstream of the CDS for the next start codon in the same frame
void extend_to_start(Mods& m,Mods& t,CDS_CHAIN_TYPE& exons,const char* subseq,int bundle_start,char strand,int start_offset,int subseq_len){ // m - to be modified; t - original mod
    while(true){ // outer loop allows us to search for the farthest start codon
        std::tuple<std::string,char,std::vector<int>> nt_na_nc; // nt,aa,coords
        Mods tmp = m; // temp copy
        bool found_start = false;
        bool res;
        while(true){ // inner loop searches for the next available start codon and breaks the chain if the stop codon was reached to notify the outer loop
            nt_na_nc = prev_codon(tmp,exons,subseq,bundle_start,strand,start_offset,subseq_len);
            for(auto& n : std::get<0>(nt_na_nc)){ // even if nt.size()%3!=0 - should still be added to cds_nt unlike to cds_aa
                tmp.cds_nt.insert(tmp.cds_nt.begin(),n);
            }
            if(std::get<1>(nt_na_nc)=='-'){ // no longer available
                break;
            }
            if(std::get<1>(nt_na_nc)=='.'){ // stop codon found
                break;
            }
            if(std::get<2>(nt_na_nc)==std::get<2>(t.start_codon)){ // reached the same start codon as the template
                break;
            }
            else{ // valid aa found - store in tmp
                bool forward = strand=='+' ? false : true;
                res = extend_chain(exons,tmp.new_chain,3,forward);
                if(!res){
                    std::cerr<<"couldn't extend - something wrong"<<std::endl;
                    exit(1);
                }
                if(std::get<1>(nt_na_nc)=='M'){
                    found_start = true;
                    break;
                }
                tmp.cds_aa=std::get<1>(nt_na_nc);
            }
        }
        if(found_start){
            m = tmp;
        }
        else{ // conditions met - cannot search any further (in the current frame at least)
            break;
        }
    }
}

void evaluate_fasta(Mods& m, const char* subseq, int bundle_start, char strand, int start_offset, int end_offset){
    if(strand=='+'){
        for(auto& c : m.new_chain){
            int cs = std::get<0>(c)-bundle_start; // start coordinate of the chain with respect to the bundle start
            int ce = std::get<1>(c)-bundle_start; // end coordinate of the chain with respect to the bundle start
            for(int i=cs;i<=ce;i++){
                m.cds_nt+=subseq[i];
            }
        }
    }
    else{ // strand=='-'
        for(int ci= m.new_chain.size() - 1; ci >= 0; ci--){
            int cs = start_offset+(std::get<0>(m.new_chain[ci]) - bundle_start);
            int ce = start_offset+(std::get<1>(m.new_chain[ci]) - bundle_start);
            for (int i=ce;i>=cs;i--) {
                m.cds_nt+=ntComplement(subseq[i]);
            }
        }
    }
    int cds_len = 0;
    if(m.cds_nt.size()%3!=0){
        std::cout<<"found in eval"<<std::endl;
    }
    m.cds_aa = translate(m.cds_nt);

    return;
    // TODO: also for each chain - output a separate stats file describing
    //     1. fraction of transcripts in the bundle it fits
    //     2. list of transcripts it fits

    // TODO: remove duplicate transcripts?
}

int nt2chain_pos(CDS_CHAIN_TYPE& chain,int nt_pos,char strand){ // tells which coordinate on the chain corresponds to the coordinate on the nt
    int chain_pos = -1;
    int left_to_stop = nt_pos;
    bool found_pos = false;
    if(strand=='+'){
        for(int i=0;i<chain.size();i++){
            int cs = std::get<0>(chain[i]);
            int ce = std::get<1>(chain[i]);
            size_t clen = (ce+1)-cs;
            if(left_to_stop<=clen){ // found the cds segment with the stop codon
                chain_pos = cs+left_to_stop;
                found_pos = true;
                break;
            }
            left_to_stop-=clen;
        }
        if(!found_pos){ // return the last position
            chain_pos = std::get<1>(chain.back());
        }
    }
    else{
        for(int i=chain.size()-1;i>=0;i--){
            int cs = std::get<0>(chain[i]);
            int ce = std::get<1>(chain[i]);
            size_t clen = (ce+1)-cs;
            if(left_to_stop<=clen){ // found the cds segment with the stop codon
                chain_pos = (ce+1)-left_to_stop;
                found_pos = true;
                break;
            }
            left_to_stop-=clen;
        }
        if(!found_pos){ // return the last position
            chain_pos = std::get<0>(chain.front());
        }
    }
    if(chain_pos<0){
        std::cerr<<"unexpected chain_pos<0"<<std::endl;
        exit(-1);
    }
    return chain_pos;
}

// retrieves nt, aa and coordinates of a codon starting at position 'start_pos' on the m.cds_nt string
std::tuple<std::string,char,std::vector<int>> get_codon(Mods& m,int start_pos,char strand){
    if(start_pos>m.cds_nt.size()-3){
        std::cerr<<"requested coordinate is past the boundaries of the recorded sequence"<<std::endl;
        exit(-1);
    }
    std::string nts = "";
    std::vector<int> coords;

    for(int i=0;i<3;i++){
        nts+=m.cds_nt[i];
        int chain_pos = nt2chain_pos(m.new_chain,start_pos+i,strand);
        coords.push_back(chain_pos);
    }

    if(nts.size()!=3){
        std::cout<<"found in get"<<std::endl;
    }
    std::string aa = translate(nts);

    return std::make_tuple(nts,aa[0],coords);
}

void trim_chain_to_pos(Mods& m, int start_nt_pos, int end_nt_pos,char strand){ // from stop indicates whether we are trimming the end or the start of the sequence
    int num_bp = end_nt_pos-start_nt_pos;

    if(end_nt_pos>m.cds_nt.size()){
        std::cerr<<"nt shorter than expected"<<std::endl;
        exit(-1);
    }

    int tmp_chain_start = nt2chain_pos(m.new_chain,start_nt_pos,strand);
    int tmp_chain_end = nt2chain_pos(m.new_chain,end_nt_pos,strand);
    int start_chain_pos = strand=='+' ? tmp_chain_start : tmp_chain_end;
    int end_chain_pos = strand=='+' ? tmp_chain_end-1 : tmp_chain_start-1;

    if(start_chain_pos<std::get<0>(m.new_chain.front())){
        std::cerr<<"start pos<chain"<<std::endl;
        exit(-1);
    }
    if(end_chain_pos>std::get<1>(m.new_chain.back())){
        std::cerr<<"end pos>chain"<<std::endl;
        exit(-1);
    }

    CDS_CHAIN_TYPE start_end,res;
    start_end.push_back(std::array<int,3>{start_chain_pos,end_chain_pos,0});
    int len = intersection(m.new_chain,start_end,res);
    int new_start_phase = start_nt_pos%3;
    assign_phase(res,strand,new_start_phase);
    m.new_chain = res;
}

// only if the very last codon is stop - trim it from the nt, aa and chain so that it does not result in mismatching chains after adjust_stop on the query
int trim_stop(Mods& m,char strand){ // return adjusted end
    if(m.cds_aa.back()=='.'){
        size_t stop_nt_pos = m.cds_nt.size()-3;
        size_t stop_aa_pos = m.cds_aa.size()-1;
        trim_chain_to_pos(m,0,stop_nt_pos,strand);
        m.cds_nt.erase(m.cds_nt.begin()+stop_nt_pos,m.cds_nt.end());
        m.cds_aa.erase(m.cds_aa.begin()+stop_aa_pos,m.cds_aa.end());
        if(chain_len(m.new_chain)!=m.cds_nt.size()){
            std::cerr<<"lengths do not match"<<std::endl;
            exit(-1);
        }
        m.adjusted = true;
    }
    return std::get<1>(m.new_chain.back());
}

void adjust_stop(Mods& m,char strand){ // adjust the chain coordinates to stop at the first stop codon
    if(m.cds_aa.empty()){
        std::cerr<<"amino chain is empty"<<std::endl;
        exit(2);
    }
    size_t stop_aa_pos = m.cds_aa.find('.');
    if(stop_aa_pos!=std::string::npos){ // stop codon was found
        size_t stop_nt_pos = stop_aa_pos*3;
        // find CDS piece in which the first nt of the stop codon was found
        trim_chain_to_pos(m,0,stop_nt_pos,strand);
        m.cds_nt.erase(m.cds_nt.begin()+stop_nt_pos,m.cds_nt.end());
        m.cds_aa.erase(m.cds_aa.begin()+stop_aa_pos,m.cds_aa.end());
        m.adjusted = true;
    }
}

void adjust_start(Mods& m,char strand){ // adjust the chain coordinates to start at the first start codon
    if(m.cds_aa.empty()){
        std::cerr<<"amino chain is empty"<<std::endl;
        exit(2);
    }
    size_t start_aa_pos = m.cds_aa.find('M');
    if(start_aa_pos!=0){
        if(start_aa_pos!=std::string::npos){
            size_t start_nt_pos = (start_aa_pos*3);
            trim_chain_to_pos(m,start_nt_pos,m.cds_nt.size(),strand);
            m.cds_nt.erase(m.cds_nt.begin(),m.cds_nt.begin()+start_nt_pos);
            m.cds_aa.erase(m.cds_aa.begin(),m.cds_aa.begin()+start_aa_pos);
            m.adjusted = true;
        }
    }
}

void compare_aa(Mods& ref_res, Mods& new_res){
    // check if starts with start codon - if not - check the first one

    // if the original chain contains an in-frame stop codon at the same position as the fitted chain - do not re-adjust

    // check if ends in stop codon - if not - check the first occurrence
    // TODO: this may require examining the 3 bases after the CDS frame - since stop codons are typically not included in the CDS
    return;
}

int get_phase(int pos,char strand,CDS_CHAIN_TYPE& phased_chain){ // returns the phase of coordinate in the provided chain
    int cdsacc = strand=='+' ? std::get<2>(phased_chain.front()) : std::get<2>(phased_chain.back());
    int phase = strand=='+' ? std::get<2>(phased_chain.front()) : std::get<2>(phased_chain.back());
    bool found_pos = false;
    if(strand=='+'){
        for(int i=0;i<phased_chain.size();i++){
            int cs = std::get<0>(phased_chain[i]);
            int ce = std::get<1>(phased_chain[i]);
            if(pos>=cs && pos<=ce){ // found the cds segment with the stop codon
                for(int j=cs;j<=pos;j++){
                    phase=(3-cdsacc%3)%3;
                    cdsacc+=1;
                }
                found_pos = true;
                break;
            }
            cdsacc += (ce-cs)+1;
        }
    }
    else{
        for(int i=phased_chain.size()-1;i>=0;i--){
            int cs = std::get<0>(phased_chain[i]);
            int ce = std::get<1>(phased_chain[i]);
            if(pos>=cs && pos<=ce){ // found the cds segment with the stop codon
                for(int j=ce;j>=pos;j--) {
                    phase=(3-cdsacc%3)%3;
                    cdsacc+=1;
                }
                found_pos = true;
                break;
            }
            cdsacc+=(ce-cs)+1;
        }
    }

    if(!found_pos){
        std::cerr<<"position not found"<<std::endl;
        exit(-1);
    }

    return phase;
}

class TX{
public:
    TX() = default;
    TX(GffObj* tx,int idx){
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

    friend std::ostream& operator<<(std::ostream& os, const TX& t){
        os << t.get_tid() << "\t" <<std::endl;
        return os;
    }

    CDS_CHAIN_TYPE get_exons(){
        return this->exons;
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

    int trim_to_null_phase(Mods& m){
        int new_start_pos = std::get<0>(m.new_chain.front());
        int new_end_pos = std::get<1>(m.new_chain.back());
        if(this->strand=='+'){
            int start_phase = get_phase(std::get<0>(m.new_chain.front()),this->strand,m.orig_chain);
            int x = get_phase(std::get<1>(m.new_chain.back()),this->strand,m.orig_chain);
            int end_phase = (3-(x-1)%3)%3;

            new_start_pos+=start_phase;
            new_end_pos-=end_phase;
        }
        else{
            int start_phase = get_phase(std::get<1>(m.new_chain.back()),this->strand,m.orig_chain);
            int x = get_phase(std::get<0>(m.new_chain.front()),this->strand,m.orig_chain);
            int end_phase = (3-(x-1)%3)%3;

            new_start_pos+=end_phase;
            new_end_pos-=start_phase;
        }
        m.new_chain.clear();
        int cut_len = cut(new_start_pos,new_end_pos,m.new_chain); // cut again to the adjusted coordinates
        return cut_len;
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

    int compare(CDS_CHAIN_TYPE& chain1,CDS_CHAIN_TYPE& chain2,Mods& res){
        int c1_i=0,c2_i=0;
        std::array<int,3> cur_c1 = chain1[c1_i];
        std::array<int,3> cur_c2 = chain2[c2_i];

        std::vector<std::pair<int,int>> tmp_mods;
        std::vector<std::pair<int,int>> tmp_mods_clean;

        std::array<int,3> inter;
        while(true){
            int inter_len = single_intersection(cur_c1, cur_c2, inter);

            if(inter_len<=0){
                if(std::get<1>(cur_c2)<std::get<0>(cur_c1)){
                    c2_i++;
                    int right = (std::get<1>(cur_c2)-std::get<0>(cur_c2))+1;
                    if(right!=0){
                        res.missing.push_back(std::array<int,3>{std::get<0>(cur_c2),std::get<1>(cur_c2),0});
                    }
                    tmp_mods.push_back(std::make_pair(right,1));

                    if(c2_i == chain2.size()){ // done
                        int left = 0-(((int)std::get<1>(cur_c1)+1) - (int)std::get<0>(cur_c1));
                        if(left!=0){
                            res.extra.push_back(std::array<int,3>{std::get<0>(cur_c1),std::get<1>(cur_c1),0});
                        }
                        for(auto cc= chain1.begin() + c1_i + 1; cc != chain1.end(); cc++){
                            left-=(((int)std::get<1>(*cc)+1)-(int)std::get<0>(*cc));
                            if(left!=0){
                                res.extra.push_back(std::array<int,3>{std::get<0>(*cc),std::get<1>(*cc),0});
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
                        res.extra.push_back(std::array<int,3>{std::get<0>(cur_c1),std::get<1>(cur_c1),0});
                    }
                    tmp_mods.push_back(std::make_pair(left,1));

                    if(c1_i == chain1.size()){ // done
                        int right = ((int)std::get<1>(cur_c2) - (int)std::get<0>(cur_c2)) + 1;
                        if(right!=0){
                            res.missing.push_back(std::array<int,3>{std::get<0>(cur_c2),std::get<1>(cur_c2),0});
                        }
                        for(auto cf= chain2.begin() + c2_i + 1; cf != chain2.end(); cf++){
                            right+=(((int)std::get<1>(*cf)-(int)std::get<0>(*cf))+1);
                            if(right!=0){
                                res.missing.push_back(std::array<int,3>{std::get<0>(*cf),std::get<1>(*cf),0});
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
                res.extra.push_back(std::array<int,3>{left_start,left_end,0});
            }
            int right_start = std::min((int)std::get<0>(cur_c2), (int)std::get<0>(inter));
            int right_end = std::max((int)std::get<0>(cur_c1), (int)std::get<0>(inter));
            int right = right_end - right_start;
            if(right!=0){
                res.missing.push_back(std::array<int,3>{right_start,right_end,0});
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

            cur_c1 = std::array<int,3>{std::get<1>(inter) + 1, std::get<1>(cur_c1),0};
            cur_c2 = std::array<int,3>{std::get<1>(inter) + 1, std::get<1>(cur_c2),0};

            if(((int)std::get<1>(cur_c1) - (int)std::get<0>(cur_c1)) < 0){
                c1_i++;
                if(c1_i == chain1.size()){ // done
                    right = ((int)std::get<1>(cur_c2) - (int)std::get<0>(cur_c2)) + 1;
                    if(right!=0){
                        res.missing.push_back(std::array<int,3>{std::get<0>(cur_c2),std::get<1>(cur_c2),0});
                    }
                    for(auto cf= chain2.begin() + c2_i + 1; cf != chain2.end(); cf++){
                        right+=(((int)std::get<1>(*cf)-(int)std::get<0>(*cf))+1);
                        if(right!=0){
                            res.missing.push_back(std::array<int,3>{std::get<0>(*cf),std::get<1>(*cf),0});
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
                        res.extra.push_back(std::array<int,3>{std::get<0>(cur_c1),std::get<1>(cur_c1),0});
                    }
                    for(auto cc= chain1.begin() + c1_i + 1; cc != chain1.end(); cc++){
                        left-=(((int)std::get<1>(*cc)+1)-(int)std::get<0>(*cc));
                        if(left!=0){
                            res.extra.push_back(std::array<int,3>{std::get<0>(*cc),std::get<1>(*cc),0});
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

    Mods fit_new(Mods& orig_cds_mod,const char* bundle_seq, int bundle_start,int start_offset,int end_offset,int bundle_len){ // computes the intersection of own exons and the chain
//        if(std::strcmp(this->tid.c_str(),"rna-NM_130762.3")==0){
//            if(std::strcmp(orig_cds_mod.orig_cds_tid.c_str(),"rna-NM_130760.3")==0){
//                std::cout<<"found"<<std::endl;
//            }
//        }

        this->mods.push_back(Mods());
        this->mods.back().orig_cds_tid = orig_cds_mod.orig_cds_tid;
        this->mods.back().orig_chain = orig_cds_mod.new_chain;

        int cut_len = 0;
        if(bundle_seq!=NULL){
            cut_len = very_strict_cut(orig_cds_mod.new_chain,this->mods.back().new_chain);
        }
        else{
            cut_len = strict_cut(orig_cds_mod.new_chain,this->mods.back().new_chain);
        }
        if(cut_len<3){ // no overlap found - create a dummy chain which will never overlap anything. This way the Mods will still get populated with data and will make it possible to filter
            this->mods.back().new_chain = CDS_CHAIN_TYPE{std::array<int,3>{0,0,0}};
        }
        else{
            cut_len = trim_to_null_phase(this->mods.back());
            if(cut_len<3){ // no overlap found - create a dummy chain which will never overlap anything. This way the Mods will still get populated with data and will make it possible to filter
                this->mods.back().new_chain = CDS_CHAIN_TYPE{std::array<int,3>{0,0,0}};
            }
            else {
                if (bundle_seq != NULL) { // check_ref flag toggled
                    evaluate_fasta(this->mods.back(), bundle_seq, bundle_start, this->strand, start_offset, end_offset);
                    adjust_stop(this->mods.back(), this->strand);
                    adjust_start(this->mods.back(), this->strand);
                    extend_to_stop(this->mods.back(), this->exons, bundle_seq, bundle_start, strand, start_offset,
                                   bundle_len);

                    // should we do only if the original start codon is not found?
                    this->mods.back().start_codon = get_codon(this->mods.back(),0,this->strand);
                    if (std::get<2>(this->mods.back().start_codon)!=std::get<2>(orig_cds_mod.start_codon)) {
                        extend_to_start(this->mods.back(),orig_cds_mod, this->exons, bundle_seq, bundle_start, strand, start_offset,
                                        bundle_len);
                    }
                    compare_aa(orig_cds_mod, this->mods.back());
                }
            }
        }

        // should we try to extend it until the query start? and if not - extend as far as you can

        if(this->mods.back().new_chain.empty()){ // no matching chain - likely due to stop-codon
            this->mods.back().new_chain = CDS_CHAIN_TYPE{std::array<int,3>{0,0,0}};
        }

        // now we can take the cut piece and directly compare it to the desired CDS counting any changes
        int score = compare(orig_cds_mod.new_chain,this->mods.back().new_chain,this->mods.back());

        if(cut_len>0 &&
           !(std::get<0>(this->mods.back().new_chain.front())==0 &&
             std::get<1>(this->mods.back().new_chain.front())==0 &&
             this->mods.back().new_chain.size()==1)){ // no overlap found - create a dummy chain which will never overlap anything. This way the Mods will still get populated with data and will make it possible to filter
            assign_phase(this->mods.back().new_chain,this->strand,0);
        }

        return this->mods.back();
    }

    bool overlap(int s, int e) {
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
        int chain_i=0;
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
        int chain_i=0;
        std::string cur_tid = "";
        for(auto& m : this->mods){
            cur_tid = this->tid+":-:"+std::to_string(chain_i);
            gtf_str+=get_gtf(cur_tid,m)+"\n";
            chain_i++;
        }

        gtf_str.pop_back(); // remove trailing new-line

        return gtf_str;
    } // TODO: output fasta of nt and aa for each perfect and imperfect transcript to check against gffread and validate the results

private:
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
    std::string mods_gtf_str;

    std::vector<Mods> mods;

    std::map<std::string,std::string> attrs;
};

class Bundle{
public:
    Bundle()=default;
    ~Bundle(){}

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
        return true;
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

    std::string get_nt(CDS_CHAIN_TYPE& chain, const char* subseq,int start_offset,int end_offset){
        std::string cds_nt = "";
        if(this->strand=='+'){
            for(auto& c : chain){
                int cs = std::get<0>(c)-this->start; // start coordinate of the chain with respect to the bundle start
                int ce = std::get<1>(c)-this->start; // end coordinate of the chain with respect to the bundle start
                for(int i=cs;i<=ce;i++){
                    cds_nt+=subseq[i];
                }
            }
        }
        else{ // strand=='-'
            for(int ci=chain.size()-1;ci>=0;ci--){
                int cs = start_offset+(std::get<0>(chain[ci])-this->start); // -3 to adjust for the stop codon on -; start coordinate of the chain with respect to the bundle start
                int ce = start_offset+(std::get<1>(chain[ci])-this->start); // -3 to adjust for the stop codon on -; end coordinate of the chain with respect to the bundle start
                for (int i=ce;i>=cs;i--) {
                    cds_nt+=ntComplement(subseq[i]);
                }
            }
        }
        return cds_nt;
    }

    int process(Globals& globals){
        const char* bundle_seq = NULL;
        int start_offset = 0;
        int end_offset = 0;
        int bundle_len = 0;
        if(this->check_ref){
            this->seqid_seq=fastaSeqGet(gfasta, this->seqid.c_str());
            if(this->seqid_seq==NULL){
                std::cerr<<"Error: no genomic sequence available (check -r option!)."<<std::endl;
            }
            int bundle_start = this->start;
            start_offset = std::min(3,this->start-1); // if start < 3 - will compensate; -1 is to compensate for the 1-based used by gfaseqget
            end_offset = std::min(3,this->seqid_seq->getseqlen()-this->end);
            if(this->strand=='-'){
                bundle_start-=start_offset;
            }
            if(start_offset>3 || end_offset>3){
                std::cerr<<"shouldn't happen"<<std::endl;
                exit(-1);
            }
            bundle_len = this->end-this->start;
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
        std::tuple<std::string,char,std::vector<int>> nt_na_nc;
        CDS_CHAIN_TYPE cur_cds_chain;
        bool valid = true;
        for(auto& tx : this->txs){
//            if(std::strcmp(tx.get_tid().c_str(),"rna-XM_005244749.3")==0){
//                std::cout<<"found"<<std::endl;
//            }
            if(tx.has_cds()){
                cur_cds_chain.clear();
                tx.build_cds_chain(cur_cds_chain);
                if(chain_len(cur_cds_chain)%3!=0){
                    std::cout<<"discarding transcript: "<<tx.get_tid()<<" - len(CDS)%3!=0"<<std::endl;
                    continue;
                }
                cit = cur_cds_chains.insert(cur_cds_chain);
                if(cit.second){ // successfully inserted
                    Mods m;
                    m.orig_chain = cur_cds_chain;
                    m.new_chain = cur_cds_chain;
                    m.orig_cds_tid = tx.get_tid();
                    if(bundle_seq!=NULL){
                        evaluate_fasta(m,bundle_seq,this->start,this->strand,start_offset,end_offset);
                        int adjusted_end = trim_stop(m,this->strand); // trim stop codon if exists
                        tx.set_cds_end(adjusted_end);
                        CDS_CHAIN_TYPE exons = tx.get_exons();
                        nt_na_nc = next_codon(m,exons,bundle_seq,this->start,this->strand,start_offset,bundle_len);
                        m.next_codon = std::get<1>(nt_na_nc);

                        m.start_codon = get_codon(m,0,tx.get_strand());

                        valid = check_valid_aa(m);
                    }
                    if(valid){ // only add valid transcripts
                        cds_chains.push_back(m);
                    }
                }
            }
        }

        // iterate over each transcript-ORF pair to gauge compatibility - assign compatibility scores
        for(auto& tx : this->txs){
//            if(std::strcmp(tx.get_tid().c_str(),"rna-XM_005244749.3")==0){
//                std::cout<<"found"<<std::endl;
//            }
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
                Mods longest_perfect_mod,longest_imperfect_mod;
                int longest_perfect_mod_len = 0, longest_imperfect_mod_len = 0;
                for(auto& chain : cds_chains){
                    if(chain.new_chain.empty()){
                        std::cerr<<"reference chain is empty"<<std::endl;
                        exit(-1);
                    }
                    Mods new_mod = tx.fit_new(chain,bundle_seq,this->start,start_offset,end_offset,bundle_len);
                    int new_mod_len = chain_len(new_mod.new_chain);
                    if(new_mod.get_score()==0 && new_mod_len>longest_perfect_mod_len){
                        longest_perfect_mod = new_mod;
                        longest_perfect_mod_len = new_mod_len;
                    }
                    else if(new_mod_len>longest_imperfect_mod_len){
                        longest_imperfect_mod = new_mod;
                        longest_imperfect_mod_len = new_mod_len;
                    }
                    else{
                        continue;
                    }
                }
                globals.out_gtf_fp<<tx.get_all_gtf()<<std::endl;
                globals.out_stats_fp<<tx.get_stats_str()<<std::endl;
                if(longest_perfect_mod_len>0){ // found perfect fit
                    globals.out_gtf_perfect_fp<<tx.get_gtf(cur_tid,longest_perfect_mod)<<std::endl;
                }
                else{ // no perfect fits found
                    globals.out_gtf_imperfect_fp<<tx.get_gtf(cur_tid, longest_imperfect_mod)<<std::endl;
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


        for(int i=0;i<gffReader.gflst.Count();++i) {
            GffObj *pGffObj = gffReader.gflst.Get(i);
            TX tmp(pGffObj,tx_vec.size());
            tx_vec.push_back(tmp);
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

        for(int i=0;i<gffReader.gflst.Count();++i) {
            GffObj *pGffObj = gffReader.gflst.Get(i);
            TX tmp(pGffObj,tx_vec.size());
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

int run(std::vector<std::string> known_gtf_fnames, std::string novel_gtf_fname,std::string out_fname,std::string ref_fa_fname){
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

    if(!ref_fa_fname.empty()){
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
    FILTER    = 'f'};

int main(int argc, char** argv) {

    ArgParse args("orfanage");
    args.add_multi_string(Opt::CDS,"cds","","Comma-separated list of GTF filenames with known transcripts and CDS annotations",true);
    args.add_string(Opt::INPUT, "input", "", "Input GTF with transcripts to which CDSs are to be ported", true);
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

    // run
    std::string ref_fname = args.is_set(Opt::REFERENCE) ? args.get_string(Opt::REFERENCE) : "";
    run(knowns,args.get_string(Opt::INPUT),args.get_string(Opt::OUTPUT),ref_fname);

    return 0;
}