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
    SuppData(const std::string& introns_fname,float thresh=0.5){
        std::ifstream int_ss;
        int_ss.open(introns_fname.c_str(),std::ios::in);
        if (!int_ss.good()){
            std::cerr<<"@ERROR::Couldn't open the intron file: "<<introns_fname<<std::endl;
            exit(1);
        }
        std::ios::sync_with_stdio(false);

        std::string aline;
        // skip header
        std::getline(int_ss,aline);
        while (std::getline(int_ss,aline)) {
            if(aline.front()=='#'){ // comment
                continue;
            }
            std::stringstream *col_stream = new std::stringstream(aline);
            std::string col;

            std::getline(*col_stream,col,','); // probability of FALSE
            float prob0 = std::stof(col);

            std::getline(*col_stream,col,','); // probability of TRUE
            float prob1 = std::stof(col);

            std::getline(*col_stream,col,','); // prior type
            int prior_type = std::stoi(col);

            if(prob1<thresh){ // if doesn't hold the threshold - do not load
                continue;
            }

            std::getline(*col_stream,col,','); // get chr
            std::string chrid = col;

            std::getline(*col_stream,col,',');
            int start = std::stoi(col);

            std::getline(*col_stream,col,',');
            int end = std::stoi(col);

            std::getline(*col_stream,col,',');
            int strand = col[0];

            std::getline(*col_stream,col,','); // skip column with med_num_samples

            std::getline(*col_stream,col,',');
            float cov = std::stof(col);

            this->iit = this->introns.insert(std::make_pair(std::make_tuple(chrid,strand,start,end),Data()));
            if(!this->iit.second){
                std::cerr<<"duplicate introns?"<<std::endl;
                exit(-1);
            }
            this->iit.first->second.add_cov(cov);
            this->iit.first->second.add_prob0(prob0);
            this->iit.first->second.add_prob1(prob1);
            this->iit.first->second.set_prior_type(prior_type);
            this->iit.first->second.set_prob1_thresh(thresh);
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

    bool remove(std::string chr,int strand,int donor,int acceptor) { // removes an intron and return true if found
        this->iit.first = this->introns.find(std::make_tuple(chr, strand, donor, acceptor));
        if(this->iit.first!=this->introns.end()){
            this->introns.erase(this->iit.first);
            return true;
        }
        return false;
    }

    bool exists(std::string chr,int strand,int donor,int acceptor) { // removes an intron and return true if found
        this->iit.first = this->introns.find(std::make_tuple(chr, strand, donor, acceptor));
        if (this->iit.first != this->introns.end()) {
            return true;
        }
        return false;
    }

private:
    struct Data{
        Data() = default;
        void add_prob0(float prob0){this->prob0 = prob0;}
        void add_prob1(float prob1){this->prob1 = prob1;}
        void add_cov(float cov){this->cov+=cov;}
        void set_prior_type(int type){this->type=type;}
        void set_prob1_thresh(float thresh){this->thresh=thresh;}
        bool is_type(int type){return this->type==type;}
        float cov;
        int type = 0;
        float prob0,prob1;
        float thresh = 0.5;
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

    bool all_introns_set(){
        for(auto& intron : this->introns){
            if(!intron.second){
                return false;
            }
        }
        return true;
    }

    float get_tpm(){
        return this->tpm;
    }

    std::string get_tid(){
        return this->tid;
    }

    void set_order(int pos){
        this->order=pos;
    }

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
    int order = -1; // expression order
};

enum SELECTION_TYPES {MIN       = 0,
                      ALLG      = 1,
                      ALL       = 2};

struct Gene{
public:
    Gene(std::string out_fname,SuppData* suppData,int selection_type){
        this->suppData = suppData;
        this->out_ss.open(out_fname,std::ios::out);
        this->out_ss<<"tid,tpm_sum,introns"<<std::endl;
        this->selection_type = selection_type;
    }
    ~Gene(){
        if(out_ss.is_open()){
            this->out_ss.close();
        }
    }
    std::string gid = "";
    std::string chrid = "";
    std::vector<Transcript> txs;
    std::vector<int> tx_order; // order of transcripts from most abundant to least abundant

    void add_tx(GffObj *p_gffObj){
        if(std::strcmp(p_gffObj->getGeneID(),this->gid.c_str())!=0){ // new gene found
            sort_transcripts();
            if(selection_type==SELECTION_TYPES::MIN){
                find_transcripts_min();
            }
            else if(selection_type==SELECTION_TYPES::ALLG){
                find_transcripts_allg();
            }
            else if(selection_type==SELECTION_TYPES::ALL){
                find_transcripts_all();
            }
            else{
                std::cerr<<"incorrect selection type used"<<std::endl;
                exit(1);
            }
            this->clear();
            this->gid = p_gffObj->getGeneID();
        }
        _add_tx(p_gffObj);
    }
private:
    SuppData* suppData;
    std::fstream out_ss;
    int selection_type = SELECTION_TYPES::MIN;

    void find_transcripts_allg(){ // select all transcripts for which all introns are good
        bool pass = false;
        for(auto& i : this->tx_order) {
            auto tx = this->txs[i];
            for(auto& intron : this->txs[i]){
                if(!suppData->exists(this->chrid,std::get<0>(intron.first),std::get<1>(intron.first),std::get<2>(intron.first))){
                    pass=false;
                    break;
                }
                else{
                    pass=true;
                }
            }
            if(pass){
                out_ss<<tx.get_tid()<<std::endl;
            }
            pass = false;
        }
    }

    void find_transcripts_all(){ // select all transcripts with at least one good intron
        bool pass = false;
        for(auto& i : this->tx_order) {
            auto tx = this->txs[i];
            for(auto& intron : this->txs[i]){
                if(suppData->exists(this->chrid,std::get<0>(intron.first),std::get<1>(intron.first),std::get<2>(intron.first))){
                    pass=true;
                    break;
                }
            }
            if(pass){
                out_ss<<tx.get_tid()<<std::endl;
            }
            pass = false;
        }
    }

    void find_transcripts_min(){
        std::vector<std::tuple<int,int,int>> passing_introns; // introns which a transcript satisfies
        bool pass = false;
        for(auto& i : this->tx_order) {
            auto tx = this->txs[i];
            for(auto& intron : this->txs[i]){
                if(suppData->remove(this->chrid,std::get<0>(intron.first),std::get<1>(intron.first),std::get<2>(intron.first))){
                    passing_introns.push_back(intron.first);
                    pass = true;
                }
            }
            if(pass){
                out_ss<<tx.get_tid()<<std::endl;
            }
            pass = false;
            passing_introns.clear();
            // check if any introns in this transcript are in the loaded list
        }
    }

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
            return;
        }

        this->chrid = p_gffObj->getGSeqName();

        this->txs.push_back(Transcript(p_gffObj->getID(),p_gffObj->strand,sum_tpm));
        int tidx = this->txs.size()-1;

        int donor=p_gffObj->exons.Get(0)->end;
        int acceptor=0;
        for(int i=1;i<p_gffObj->exons.Count();i++){
            acceptor = p_gffObj->exons.Get(i)->start;
            this->txs.back().add_intron(p_gffObj->strand,donor,acceptor);
            donor = p_gffObj->exons.Get(i)->end;
        }
    }

    void clear(){
        this->gid = "";
        this->chrid = "";
        this->txs.clear();
        this->tx_order.clear();
    }
};

enum Opt {GTF       = 'g',
          INTRONS   = 'i',
          OUTPUT    = 'o',
          SELECTION_TYPE = 't'};

int main(int argc, char** argv) {

    ArgParse args("novel_introns");
    args.add_string(Opt::GTF, "gtf", "", "Annotation file produced by assembly stats", true);
    args.add_string(Opt::INTRONS, "introns", "", "Introns in BED format with 4th column indicating coverage by spliced reads and 5th column indicating: -1 is fake intron and no transcript with that intron will be considered; 1 is good intron and will be used to select which transcripts to include; 0 is ambiguous intron and can be included if good introns exist in the same transcript", true);
    args.add_string(Opt::OUTPUT,"out","","Basename for the output files",true);
    args.add_string(Opt::SELECTION_TYPE,"type","min","Type of selection to perform on the transcripts. Choices: min (default) - will select the minimum intronically non-redundant set of transcripts with maximum coverage; allg - will select all transcripts that select only good introns; all - will select all transcripts which include good introns.",false);

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

    std::string tmp = args.get_string(Opt::SELECTION_TYPE);
    int selection_type = 0;
    if(std::strcmp(tmp.c_str(),"min")==0){
        selection_type=SELECTION_TYPES::MIN;
    }
    else if(std::strcmp(tmp.c_str(),"allg")==0){
        selection_type=SELECTION_TYPES::ALLG;
    }
    else if(std::strcmp(tmp.c_str(),"all")==0){
        selection_type=SELECTION_TYPES::ALL;
    }
    else{
        std::cerr<<"incorrect selection method requested. Please select from: min,allg,all "<<std::endl;
        exit(1);
    }

    Gene gene(args.get_string(Opt::OUTPUT),&suppData,selection_type);

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