#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <numeric>

#include "arg_parse.h"
#include "TrackingStats.h"
#include "TrackingTree.h"

enum INTRON_TYPE {FALSE = -1,
                  AMB   = 0,
                  GOOD  = 1,
                  KNOWN = 2
};

struct SuppData{
public:
    SuppData(const std::string& introns_fname){
        std::ifstream int_ss;
        int_ss.open(introns_fname.c_str(),std::ios::in);
        if (!int_ss.good()){
            std::cerr<<"@ERROR::Couldn't open the intron file: "<<introns_fname<<std::endl;
            exit(1);
        }
        std::ios::sync_with_stdio(false);

        std::string aline;
        while (std::getline(int_ss,aline)) {
            if(aline.front()=='#'){ // comment
                continue;
            }
            std::stringstream *col_stream = new std::stringstream(aline);
            std::string col;

            std::getline(*col_stream,col,'\t'); // get chr
            std::string chrid = col;

            std::getline(*col_stream,col,'\t');
            int strand = col[0];

            std::getline(*col_stream,col,'\t');
            int start = std::stoi(col);

            std::getline(*col_stream,col,'\t');
            int end = std::stoi(col);

            std::getline(*col_stream,col,'\t');
            float cov = std::stof(col);

            std::getline(*col_stream,col,'\t');
            int type = std::stoi(col); // -1(fake);1(good);0(ambiguous);2(known)

            this->iit = this->introns.insert(std::make_pair(std::make_tuple(chrid,strand,start,end),Data(cov)));
            if((type==INTRON_TYPE::FALSE ||
                    type==INTRON_TYPE::AMB ||
                    type==INTRON_TYPE::GOOD ||
                    type==INTRON_TYPE::KNOWN)){
                this->iit.first->second.set_type(type);
            }
            else{
                std::cerr<<"wrong type in the introns file"<<std::endl;
                exit(-1);
            }
            if(!this->iit.second){
                std::cerr<<"duplicate introns?"<<std::endl;
                exit(-1);
            }
        }
        int_ss.close();
    }
    ~SuppData() = default;

    int get_cov(std::string& chr,int strand,int donor,int acceptor){
        this->iit.first = this->introns.find(std::make_tuple(chr,strand,donor,acceptor));
        if(this->iit.first == this->introns.end()){
            return -1;
        }
        return this->iit.first->second.cov;
    }

    int get_type(std::string& chr,int strand,int donor,int acceptor) {
        this->iit.first = this->introns.find(std::make_tuple(chr, strand, donor, acceptor));
        if (this->iit.first == this->introns.end()) {
            return false;
        }
        return this->iit.first->second.type;
    }

    bool is_type(std::string& chr,int strand,int donor,int acceptor,int type){
        this->iit.first = this->introns.find(std::make_tuple(chr,strand,donor,acceptor));
        if(this->iit.first == this->introns.end()){
            return false;
        }
        return this->iit.first->second.is_type(type);
    }

private:
    struct Data{
        explicit Data(float cov):cov(cov){};
        void add_cov(float cov){this->cov+=cov;}
        void set_type(int type){this->type=type;}
        bool is_type(int type){return this->type==type;}
        float cov;
        int type = 0;
    };
    std::map<std::tuple<std::string,int,int,int>,Data> introns; // key: chrid,strand,start,end
    std::pair<std::map<std::tuple<std::string,int,int,int>,Data>::iterator,bool> iit;
};

struct Transcript{
public:
    Transcript(std::string tid,int strand,float tpm):tid(tid),strand(strand),tpm(tpm){};
    ~Transcript()=default;

    void add_intron(int strand,int start,int end){
        this->introns.insert(std::make_pair(std::make_tuple(strand,start,end),false));
    }

    bool set_intron_pass(std::tuple<int,int,int> intron){
        this->iit.first = this->introns.find(intron);
        this->iit.first->second=true;
    }

    void set_pass(){
        this->pass=true;
    }

    void set_pass_order(){
        this->pass_order=true;
    }

    void set_pass_order_unq(){
        this->pass_order_unq=true;
    }

    bool all_introns_set(){
        for(auto& intron : this->introns){
            if(!intron.second){
                return false;
            }
        }
        return true;
    }

    void set_false(){
        this->fake = true;
    }
    bool is_false(){
        return this->fake;
    }

    float get_tpm(){
        return this->tpm;
    }

    std::string get_tid(){
        return this->tid;
    }

    int get_order(){
        return this->order;
    }

    int get_strand(){
        return this->strand;
    }

    bool passes(){return this->pass;}
    bool passes_order(){return this->pass_order;}
    bool passes_order_unq(){return this->pass_order_unq;}

    void set_order(int pos){
        this->order=pos;
    }

//    std::map<std::pair<int,int>,bool>* get_introns(){return &this->introns;}

    typedef std::map<std::tuple<int,int,int>,bool>::iterator iterator;
    typedef std::map<std::tuple<int,int,int>,bool>::const_iterator const_iterator;
    typedef std::map<std::tuple<int,int,int>,bool>::reference reference;
    iterator begin() { return introns.begin(); }
    const_iterator begin() const { return introns.cbegin(); }
    iterator end() { return introns.end(); }
    const_iterator end() const { return introns.cend(); }

private:
    std::string tid;
    int strand;
    float tpm;
    std::map<std::tuple<int,int,int>,bool> introns; // bool tells whether the intron passes or not
    std::pair<std::map<std::tuple<int,int,int>,bool>::iterator,bool> iit;
    bool pass=false; // has the transcript already been written
    bool pass_order=false; // has the transcript already been written based on the ordered output
    bool pass_order_unq = false; // has the transcript been written to the unique stream yet
    int order = -1; // expression order
    bool fake = false; // contains fake introns
};

struct Intron{
public:
    Intron(std::string seqid,int strand, int start,int end):seqid(seqid),strand(strand),start(start),end(end){};
    ~Intron()=default;

    std::string get_seqid(){return this->seqid;}
    int get_strand(){return this->strand;}
    int get_start(){return this->start;}
    int get_end(){return this->end;}

    void set_used(){
        this->used = true;
    }
    void set_false(){
        this->fake = true;
    }
    void set_good(){
        this->good = true;
    }

    bool is_used(){return this->used;}
    bool is_false(){return this->fake;}
    bool is_good(){return this->good;}

    void set_cov(float cov){
        this->cov = cov;
    }

    int get_cov(){return this->cov;}

    float tpm_sum(){
        return std::accumulate(this->tpms.begin(),this->tpms.end(),0.0);
    }

    void add_tpm(float tpm){
        this->tpms.push_back(tpm);
    }

    void add_tx(int tidx){
        this->tx_idxs.push_back(tidx);
    }

    void set_pass(){
        this->pass=true;
    }

    int get_txidx(int i){
        return this->tx_idxs[i];
    }

    int txCount(){
        return this->tx_idxs.size();
    }

    bool passes(){return this->pass;}

private:
    std::string seqid;
    int strand;
    int start;
    int end;
    float cov = -1;
    std::vector<float> tpms;
    std::vector<int> tx_idxs; // indices of transcripts in the vector which holds transcripts
    bool pass = false; // set to true when this intron is outputted from the sorted list. Used to check for compatible transcripts
    bool used = false; // intron has been used in a transcript
    bool fake = false; // intron is fake and should not be included in any transcript
    bool good = false; // intron is good and must be used in a transcript
};

bool intron_sort (std::map<std::tuple<int,int,int>,Intron>::iterator i,std::map<std::tuple<int,int,int>,Intron>::iterator j) { return (i->second.get_cov()>j->second.get_cov()); }

struct Gene{
public:
    Gene(std::string out_fname,SuppData* suppData){
        this->suppData = suppData;
        std::string out_intron_fname = out_fname+".introns";
        std::string out_tx_fname = out_fname+".txs";
        std::string out_tx_ordered_fname = out_fname+".txs_ord";
        std::string out_tx_ordered_uniq_fname = out_fname+".txs_ord_unq";

        this->out_intron_ss.open(out_intron_fname,std::ios::out);
        this->out_intron_ss<<"gid,seqid,strand,start,end,tx_count,frac_count,frac,cumulative"<<std::endl;

        this->out_tx_ss.open(out_tx_fname,std::ios::out);
        this->out_tx_ss<<"gid,tid,frac_ints,frac_int_quant,cum_frac_int_quant,frac_tx_tpm,cum_frac_tx_tpm"<<std::endl;

        this->out_tx_ss_ordered.open(out_tx_ordered_fname,std::ios::out);
        this->out_tx_ss_ordered<<"gid,tid,int_cov,frac_ints,frac_int_quant,cum_frac_int_quant,frac_tx_tpm,cum_frac_tx_tpm"<<std::endl;

        this->out_tx_ss_ordered_uniq.open(out_tx_ordered_uniq_fname,std::ios::out);
        this->out_tx_ss_ordered_uniq<<"gid,tid,frac_ints,frac_int_quant,cum_frac_int_quant,frac_tx_tpm,cum_frac_tx_tpm"<<std::endl;
    }
    ~Gene(){
        if(out_intron_ss.is_open()){
            this->out_intron_ss.close();
        }
        if(out_tx_ss.is_open()){
            this->out_tx_ss.close();
        }
        if(out_tx_ss_ordered.is_open()){
            this->out_tx_ss_ordered.close();
        }
        if(out_tx_ss_ordered_uniq.is_open()){
            this->out_tx_ss_ordered_uniq.close();
        }
    }
    std::string gid = "";
    float gene_tpm = 0.0;
    std::map<std::string,std::vector<float>> single_exons; // key tid - value is a list of TPMs
    std::pair<std::map<std::string,std::vector<float>>::iterator,bool> eit;
    std::map<std::tuple<int,int,int>,Intron> introns; // key is donor-acceptor pair - value is a list of TPMs
    std::pair<std::map<std::tuple<int,int,int>,Intron>::iterator,bool> iit;
    std::vector<Transcript> txs;
    std::vector<int> tx_order; // order of transcripts from most abundant to least abundant

    void add_tx(GffObj *p_gffObj){
        if(std::strcmp(p_gffObj->getGeneID(),this->gid.c_str())!=0){ // new gene found
            sort_transcripts();
            write_introns_sorted();
            this->clear();
            this->gid = p_gffObj->getGeneID();
        }
        _add_tx(p_gffObj);
    }
private:
    SuppData* suppData;
    std::fstream out_intron_ss,out_tx_ss,out_tx_ss_ordered,out_tx_ss_ordered_uniq;
    float cur_total_tx_tpm = 0.0; // keeps track of how much of the gene expression is explained by the transcripts written thus far
    float cur_total_tx_tpm_ordered = 0.0; // keeps track of how much of the gene expression is explained by the transcripts written thus far

    // sorts transcripts in descending order by expression by assigning each transcript a place in the order as opposed to re-organizing the array
    void sort_transcripts(){
        // gather all transcripts into an array of pairs (index,expression)(
        std::vector<std::pair<int,float>> tx_exp;
        for(int i=0;i<this->txs.size();i++){
            tx_exp.push_back(std::make_pair(i,this->txs[i].get_tpm()));
        }

        // sort by the second value in array (expression)
        std::sort(tx_exp.begin(),tx_exp.end(),[](std::pair<int,float> a, std::pair<int,float> b) {return a.second > b.second; });

        // set order of each transcript accordingly
        this->tx_order.clear();
        for(int i=0;i<tx_exp.size();i++){
            this->txs[tx_exp[i].first].set_order(i);
            this->tx_order.push_back(tx_exp[i].first);
        }
    }

    float get_tpm_sum(GffObj *p_gffObj){
        int tpm_attid = p_gffObj->names->attrs.getId("tpms");
        if(tpm_attid==-1){ // TPM atribute not found
            std::cerr<<"no attribute tpms found in file"<<std::endl;
            exit(-1);
        }

        std::string tpm_str = p_gffObj->attrs->getAttr(tpm_attid);
        std::string cur_tpm_str;
        std::istringstream tpm_ss(tpm_str);

        float total_tpm = 0.0;
        while (getline(tpm_ss,cur_tpm_str,',')) {
            total_tpm += std::stof(cur_tpm_str);
        }
        return total_tpm;
    }

    void _add_tx(GffObj *p_gffObj){
        float sum_tpm = get_tpm_sum(p_gffObj);

        if(p_gffObj->exons.Count()==1){ // single exon
            this->eit = this->single_exons.insert(std::make_pair(p_gffObj->getID(),std::vector<float>{}));
            this->eit.first->second.push_back(sum_tpm);
            return;
        }

        this->gene_tpm+=sum_tpm;
        this->txs.push_back(Transcript(p_gffObj->getID(),p_gffObj->strand,sum_tpm));
        int tidx = this->txs.size()-1;

        int donor=p_gffObj->exons.Get(0)->end;
        int acceptor=0;
        for(int i=1;i<p_gffObj->exons.Count();i++){
            acceptor = p_gffObj->exons.Get(i)->start;
            this->txs.back().add_intron(p_gffObj->strand,donor,acceptor);

            // find intron in the suppdata
            std::string chrid = p_gffObj->getGSeqName();

            this->iit = this->introns.insert(std::make_pair(std::make_tuple(p_gffObj->strand,donor,acceptor),Intron(p_gffObj->getGSeqName(),p_gffObj->strand,donor,acceptor)));
            float cov = suppData->get_cov(chrid,p_gffObj->strand,donor,acceptor); // return of -1 means the intron was not found
            if(cov==-1){
                cov=0;
            }
            this->iit.first->second.set_cov(cov);

            int cur_int_type = suppData->get_type(chrid,p_gffObj->strand,donor,acceptor);
            if(cur_int_type == INTRON_TYPE::KNOWN || cur_int_type == INTRON_TYPE::AMB){ // both ambiguous and known can be used in forming transcripts but only if a good novel intron is present
                this->iit.first->second.set_used();
            }
            else if(cur_int_type == INTRON_TYPE::FALSE){ // the intron is added to the transcript but not to the introns - so it will never be written out
                this->iit.first->second.set_false();
                this->txs.back().set_false();
            }
            else if(cur_int_type ==  INTRON_TYPE::GOOD){
                this->iit.first->second.set_good();
            }
            this->iit.first->second.add_tpm(sum_tpm);
            this->iit.first->second.add_tx(tidx);
            donor = p_gffObj->exons.Get(i)->end;
        }
    }

    void clear(){
        this->gid = "";
        this->single_exons.clear();
        this->introns.clear();
        this->txs.clear();
        this->cur_total_tx_tpm=0.0;
        this->gene_tpm=0.0;
        this->tx_order.clear();
        this->cur_total_tx_tpm_ordered=0.0;
    }

    void write_introns_sorted(){
        if(this->introns.empty()){
            return;
        }
        // compute total quantity
        std::vector<std::map<std::tuple<int,int,int>,Intron>::iterator> sorted_introns;
        float total_cov = 0;
//        for(auto& exon : this->single_exons){
//            total_tpm += std::accumulate(exon.second.begin(),exon.second.end(),0.0);
//        }
        for(auto intron_it=this->introns.begin();intron_it!=this->introns.end();intron_it++){
            if(intron_it->second.is_false()){
                continue; // don't use false introns
            }
            // get total expression of the intron across all transcripts where it occurs
            float cur_cov = intron_it->second.get_cov();
            if(cur_cov==-1){ // intron not found
                continue;
            }
            total_cov += cur_cov;
            sorted_introns.push_back(intron_it);
        }
        // sort the introns by quantity in descending order
        std::sort(sorted_introns.begin(),sorted_introns.end(),intron_sort);

        // output introns
        float cumulative_frac = 0;
        int int_count = 1;
        for(std::map<std::tuple<int,int,int>,Intron>::iterator& intron : sorted_introns){
            intron->second.set_pass();
            // set pass for this intron in all related transcripts
            for(int i=0;i<intron->second.txCount();i++){
                int tidx = intron->second.get_txidx(i);
                this->txs[tidx].set_intron_pass(std::make_tuple(intron->second.get_strand(),intron->second.get_start(),intron->second.get_end()));
            }
            cumulative_frac += (intron->second.get_cov()/total_cov)*100.0;
            this->out_intron_ss<<this->gid<<","
                               <<intron->second.get_seqid()<<","
                               <<(char)intron->second.get_strand()<<","
                               <<intron->second.get_start()<<","
                               <<intron->second.get_end()<<","
                               <<(float)intron->second.txCount()/(float)this->txs.size()<<","
                               <<((float)int_count/(float)sorted_introns.size())*100<<","
                               <<(intron->second.get_cov()/total_cov)*100.0<<","
                               <<cumulative_frac<<std::endl;

            // now check all transcripts and write those which can be written
            write_txs(((float)int_count/(float)sorted_introns.size())*100,(intron->second.get_cov()/total_cov)*100.0,cumulative_frac);

            // now check all highest order transcripts
            write_exp_ordered_txs(((float)int_count/(float)sorted_introns.size())*100,(intron->second.get_cov()/total_cov)*100.0,cumulative_frac,intron);

            int_count++;
        }
    }

    bool all_introns_used(Transcript& tx){
        for(auto& tx_int : tx){
            this->iit.first = this->introns.find(tx_int.first);
            if(this->iit.first==this->introns.end()){
                std::cerr<<"intron not found"<<std::endl;
                exit(-1);
            }
            if(!this->iit.first->second.is_used()){
                return false;
            }
        }
        return true;
    }

    void write_exp_ordered_txs(float frac_ints,float frac_int_quant,float cum_frac_int_quant,std::map<std::tuple<int,int,int>,Intron>::iterator intron){ // search for new transcripts which fit the current set of introns but consider only transcripts in descending order by expression (even if a new transcripts is satisfied by the set of introns, it is not outputted if it next in order)
        for(int i=0;i<this->tx_order.size();i++) {
            if(this->txs[this->tx_order[i]].is_false()){ // has fake introns - needs to be skipped
                continue;
            }
            if(this->txs[this->tx_order[i]].all_introns_set()) {
                if(!this->txs[this->tx_order[i]].passes_order()){
                    cur_total_tx_tpm_ordered+=this->txs[this->tx_order[i]].get_tpm();
                    this->out_tx_ss_ordered<<this->gid<<","
                                           <<this->txs[this->tx_order[i]].get_tid()<<","
                                           <<intron->second.get_cov()<<","
                                           <<frac_ints<<","
                                           <<frac_int_quant<<","
                                           <<cum_frac_int_quant<<","
                                           <<(this->txs[this->tx_order[i]].get_tpm()/this->gene_tpm)*100<<","
                                           <<(cur_total_tx_tpm_ordered/this->gene_tpm)*100
                                           <<std::endl;
                    this->txs[this->tx_order[i]].set_pass_order();
                }
                if(!this->txs[this->tx_order[i]].passes_order_unq() &&
                   !this->all_introns_used(this->txs[this->tx_order[i]])){ // has not been previously written out and has unused introns
                    this->out_tx_ss_ordered_uniq<<this->gid<<","
                                                <<this->txs[this->tx_order[i]].get_tid()<<","
                                                <<frac_ints<<","
                                                <<frac_int_quant<<","
                                                <<cum_frac_int_quant<<","
                                                <<(this->txs[this->tx_order[i]].get_tpm()/this->gene_tpm)*100<<","
                                                <<(cur_total_tx_tpm_ordered/this->gene_tpm)*100
                                                <<std::endl;
                    for(auto& tx_intron : this->txs[this->tx_order[i]]){
                        this->iit.first = this->introns.find(tx_intron.first);
                        this->iit.first->second.set_used();
                    }
                    this->txs[this->tx_order[i]].set_pass_order_unq();
                }
            }
            else{ // does not pass - no need to check the other transcripts until this one is explained
                break;
            }
        }
    }

//    void write_exp_ordered_txs(float frac_ints,float frac_int_quant,float cum_frac_int_quant,std::map<std::tuple<int,int,int>,Intron>::iterator intron){ // search for new transcripts which fit the current set of introns but consider only transcripts in descending order by expression (even if a new transcripts is satisfied by the set of introns, it is not outputted if it next in order)
//        for(int i=0;i<this->tx_order.size();i++) {
//            if (this->txs[this->tx_order[i]].passes_order()) {
//                continue;
//            }
//            if (this->txs[this->tx_order[i]].all_introns_set()) {
//                cur_total_tx_tpm_ordered+=this->txs[this->tx_order[i]].get_tpm();
//                this->out_tx_ss_ordered<<this->gid<<","
//                                       <<this->txs[this->tx_order[i]].get_tid()<<","
//                                       <<intron->second.get_cov()<<","
//                                       <<frac_ints<<","
//                                       <<frac_int_quant<<","
//                                       <<cum_frac_int_quant<<","
//                                       <<(this->txs[this->tx_order[i]].get_tpm()/this->gene_tpm)*100<<","
//                                       <<(cur_total_tx_tpm_ordered/this->gene_tpm)*100
//                                       <<std::endl;
//                this->txs[this->tx_order[i]].set_pass_order();
//
//                // now also check if any of the introns have not been used before
//                if(!intron->second.is_used()){
//                    this->out_tx_ss_ordered_uniq<<this->gid<<","
//                                                <<this->txs[this->tx_order[i]].get_tid()<<","
//                                                <<intron->second.get_cov()<<","
//                                                <<frac_ints<<","
//                                                <<frac_int_quant<<","
//                                                <<cum_frac_int_quant<<","
//                                                <<(this->txs[this->tx_order[i]].get_tpm()/this->gene_tpm)*100<<","
//                                                <<(cur_total_tx_tpm_ordered/this->gene_tpm)*100
//                                                <<std::endl;
//                    for(auto& tx_intron : this->txs[this->tx_order[i]]){
//                        this->iit.first = this->introns.find(tx_intron.first);
//                        this->iit.first->second.set_used();
//                    }
//                    continue;
//                }
//                else{ // all introns have already been used
//                    continue;
//                }
//            }
//            else{ // does not pass - no need to check the other transcripts until this one is explained
//                break;
//            }
//        }
//    }

    void write_txs(float frac_ints,float frac_int_quant,float cum_frac_int_quant){ // search for new transcripts which fit the current set of introns
        for(auto& tx : this->txs){
            if(tx.passes()){
                continue;
            }
            if(tx.is_false()){
                continue;
            }
            // for all introns in the transcript check if all of them satisfy
            if(tx.all_introns_set()){
                // write
                cur_total_tx_tpm+=tx.get_tpm();
                this->out_tx_ss<<this->gid<<","
                               <<tx.get_tid()<<","
                               <<frac_ints<<","
                               <<frac_int_quant<<","
                               <<cum_frac_int_quant<<","
                               <<(tx.get_tpm()/this->gene_tpm)*100<<","
                               <<(cur_total_tx_tpm/this->gene_tpm)*100
                               <<std::endl;
                tx.set_pass();
            }
        }
    }
};

enum Opt {GTF       = 'g',
        INTRONS   = 'i',
        OUTPUT    = 'o'};

int main(int argc, char** argv) {

    ArgParse args("novel_introns");
    args.add_string(Opt::GTF, "gtf", "", "Annotation file produced by assembly stats", true);
    args.add_string(Opt::INTRONS, "introns", "", "Introns in BED format with 4th column indicating coverage by spliced reads and 5th column indicating: -1 is fake intron and no transcript with that intron will be considered; 1 is good intron and will be used to select which transcripts to include; 0 is ambiguous intron and can be included if good introns exist in the same transcript", true);
    args.add_string(Opt::OUTPUT,"out","","Basename for the output files",true);

    if (argc <= 1 || strcmp(argv[1], "--help") == 0) {
        std::cerr << args.get_help() << std::endl;
        exit(1);
    }

    args.parse_args(argc, argv);

    // first create the execution string
    std::string cl = "novel_introns ";
    for (int i = 0; i < argc; i++) {
        if (i == 0) {
            cl += argv[i];
        } else {
            cl += " ";
            cl += argv[i];
        }
    }

    SuppData suppData(args.get_string(Opt::INTRONS));

    Gene gene(args.get_string(Opt::OUTPUT),&suppData);

    FILE* gff_file = fopen(args.get_string(Opt::GTF).c_str(), "r");
    if (gff_file == nullptr){
        std::cerr << "@ERROR::Couldn't open annotation: " << args.get_string(Opt::GTF)<< std::endl;
        exit(1);
    }

    GffReader gffReader;
    gffReader.init(gff_file,true);
    gffReader.readAll(true);

    GffObj *p_gffObj;

    std::string curGeneID = "";

    for (int i = 0; i < gffReader.gflst.Count(); ++i){
        p_gffObj = gffReader.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0){
            continue;
        }
        p_gffObj->store_elen();

        gene.add_tx(p_gffObj);
    }

    return 0;
}