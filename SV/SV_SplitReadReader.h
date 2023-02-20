

using namespace std;



using namespace BamTools;

struct split_read_parameters {
	string bam_file,
		sample_name;
	unsigned int min_non_overlap,
				 back_distance,
				 min_mapping_threshold,
                 min_clip;
	int weight;
	vector<string> read_group;
};

class SV_SplitReadReader : public SV_EvidenceReader
{
    public:
        string bam_file;
        unsigned int min_non_overlap,
                     back_distance,
                     min_mapping_threshold,
                 min_clip;
        int weight;
        vector<string> read_group;
        bool is_open,
        have_next_alignment;

        BamAlignment bam;
        BamReader reader;
        map<string, BamAlignment> mapped_splits;
        string header;
        RefVector refs;
        bool inited;

        ~SV_SplitReadReader();
        SV_SplitReadReader();
        string check_params();
        void initialize();
        void set_statics();
        void unset_statics();

        void process_input( BamAlignment &_bam,
        RefVector &_ref,

        void process_input(BamAlignment &_bam,
                           RefVector &_ref,
                           BamWriter &inter_chrom_reads,

        void terminate();
        string get_curr_chr();
        CHR_POS get_curr_pos();
        bool has_next();
        string get_source_file_name();
};

