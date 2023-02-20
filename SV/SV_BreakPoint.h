
class SV_Evidence;




using namespace std;

struct interval {
	CHR_POS start, end;
	CHR_POS start_clip, end_clip;
	char strand;
	string chr;
};

struct breakpoint_interval
{
	struct interval i;
};

class SV_BreakPoint
{
	friend ostream& operator<<(ostream& out, const SV_BreakPoint& b);

	private:
			                  unsigned int size,

	public:
		static const int DELETION = 1;
		static const int DUPLICATION = 2;
		static const int INVERSION = 3;
		static const int TRANSLOCATION = 4;
		static double p_trim_threshold;
		static double p_merge_threshold;
		static string ascii_interval_prob(
		int type;
		int weight;
		map<int, int> ev_ids;
		struct breakpoint_interval interval_l, interval_r;
		SV_BreakPoint();
		~SV_BreakPoint();
		void free_evidence();
                                           bool check_strand);
                static bool test_interval_merge(
                void print_evidence(string pre);
                int trim_intervals();
                void init_interval_probabilities();
                void free_interval_probabilities();
                void print_bedpe(int id, int print_prob);
                void print_interval_probabilities();
                vector<int> get_evidence_ids();
                void do_it();
                static pair<CHR_POS,CHR_POS> min_pair(
                        vector< vector< pair<CHR_POS,CHR_POS> > > &m);
		int get_max_sample_weight();
};
