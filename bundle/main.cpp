#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "arg_parse.h"

enum Opt {GTF       = 'g',
          OUTPUT    = 'o'};

/* define scores
 1. number of bases of overlap - needs to be normalized somehow since variable length overlaps will
    ultimately dominate over any other score - could do percentage of total effective length normalized
    then to min max (0-2 where max defines the maximum weight the overlap feature may have)
 2. number of exons shared
 3. number of introns shared
 4. number of donors shared
 5. number of acceptors shared

 Each feature needs a weight assigned (much like blossum etc.)

 Then we can perform hierarchical clustering to get a distance matrix

 Alternative solution (more complicated than hierarchical clustering) is to use splice graphs
    - build a splice graph of all transcripts that overlap
    - iterate over all transcripts and add weights to corresponding features
    - assign edges between components of the graph which are increased in weight by the new transcript
        (add transcript IDs to these new edges)
    - when all is done - traverse the graph and identify clusters of weights */

struct Transcript{
public:
    Transcript()=default;
    ~Transcript()=default;

    int add_exon(){

    }
private:
};

struct Bundle{
public:
    Bundle() = default;
    ~Bundle() = default;

    int add_transcript(){

    }
private:
};

int main(int argc, char** argv) {

    ArgParse args("bundle");
    args.add_string(Opt::GTF, "gtf", "", "File containing the merged gtf of all tissues. The file is expected to have classification codes", true);
    args.add_string(Opt::OUTPUT, "output", "", "Basename for the output files", true);

    if (argc <= 1 || strcmp(argv[1], "--help") == 0) {
        std::cerr << args.get_help() << std::endl;
        exit(1);
    }

    args.parse_args(argc, argv);

    // first create the execution string
    std::string cl = "bundle ";
    for (int i = 0; i < argc; i++) {
        if (i == 0) {
            cl += argv[i];
        } else {
            cl += " ";
            cl += argv[i];
        }
    }



    return 0;
}