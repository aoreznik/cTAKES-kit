class SV_BreakPoint;



using namespace std;


class SV_Evidence
{

public:
        static map<int, int> distros_size;
        static UCSCBins<int> exclude_regions;
        int weight;
        int ev_id;
        int type;


        virtual string evidence_type();

        virtual void print_evidence();
        virtual ~SV_Evidence();

};
