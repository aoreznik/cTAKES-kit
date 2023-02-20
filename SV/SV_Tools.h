


using namespace std;


int read_histo_file(string file_name,

int read_distro_file(string file_name,

bool sort_inter_chrom_bam(string in_file_name,
			  string out_file_name);

bool create_sorted_temp_file(vector<BamAlignment>& buffer,
                             string out_file_name,
                             int num_runs,
                             string header_text,
                             RefVector &ref);

bool merge_sorted_files(string out_file_name,
			int buff_count,
                        string header_text,
                        RefVector &ref);

bool write_temp_file(vector<BamAlignment>& buffer,
		     string temp_file_name,
                     string header_text,
                     RefVector &ref);



void parse_exclude_file(string exclude_bed_file,
                        UCSCBins<int> &exclude_regions);
uint32_t
count_clipped(vector< CigarOp > cigar_data);

