

using namespace BamTools;



using namespace std;

class SV_Pair: public SV_Evidence
{
    friend ostream& operator<<(ostream& out, const SV_Pair& p);

    private:
                                              unsigned int back_distance,
                                              unsigned int distro_size);

	public:
		static double insert_Z;

		unsigned int min_mapping_quality;
		struct interval read_l;
		struct interval read_r;
		bool read_l_is_split, read_r_is_split;
		std::string read_id;


		SV_Pair(const BamAlignment &bam_a,
				const BamAlignment &bam_b,
				const RefVector &refs,
				int weight,
				int ev_id,

		static void process_pair(const BamAlignment &curr,
								const RefVector refs,
								map<string, BamAlignment> &mapped_pairs,
								int weight,
								int ev_id,

		static void process_intra_chrom_pair(
								 const BamAlignment &curr,
								 const RefVector refs,
								 BamWriter &inter_chrom_reads,
								 map<string, BamAlignment> &mapped_pairs,
								 int weight,
								 int ev_id,

													  int distro_size,

		bool is_aberrant();
		bool is_sane();
		bool is_interchromosomal();

		static void set_distro_from_histo ();
		static int set_distro_from_histo(int back_distance,
										 int histo_start,
										 int histo_end,

		void print_evidence();

		void print_bedpe(int score);
		string evidence_type();
};

