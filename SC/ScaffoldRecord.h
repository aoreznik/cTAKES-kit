

struct ResolveStats
{
    ResolveStats()
    {
        numGapsResolved = 0;
        numGapsAttempted = 0;
        numScaffolds = 0;
        
        graphWalkFound = 0;
        graphWalkTooMany = 0;
        graphWalkNoPath = 0;

        overlapFound = 0;
        overlapFailed = 0;
    }

    ~ResolveStats() { print(); }
    
    void print() const
    {
        printf("Num scaffolds: %d\n", numScaffolds);
        printf("Num gaps attempted: %d\n", numGapsAttempted);
        printf("Num gaps resolved: %d\n", numGapsResolved);

        printf("Num gaps resolved by graph walk: %d\n", graphWalkFound);
        printf("Num graph walks failed because of too many solutions: %d\n", graphWalkTooMany);
        printf("Num graph walks failed because of no path: %d\n", graphWalkNoPath);

        printf("Num gaps resolved by overlap: %d\n", overlapFound);
        printf("Num overlaps failed: %d\n", overlapFailed);
    }

    int numGapsResolved;
    int numGapsAttempted;
    int numScaffolds;

    int graphWalkFound;
    int graphWalkTooMany;
    int graphWalkNoPath;

    int overlapFound;
    int overlapFailed;
};

struct ResolveParams
{

    int minOverlap;
    int maxOverlap;
    double maxErrorRate;
    int resolveMask;
    int minGapLength;
    double distanceFactor;
};

const int RESOLVE_GRAPH_UNIQUE = 1;
const int RESOLVE_GRAPH_BEST = 2;
const int RESOLVE_OVERLAP = 4;

class ScaffoldRecord
{
    public:
        ScaffoldRecord();

        void setRoot(const std::string& root);
        void addLink(const ScaffoldLink& link);
        size_t getNumComponents() const;

        std::string generateString(const ResolveParams& params, StringVector& ids) const;

        bool graphResolve(const ResolveParams& params, const std::string& startID, 
                          const ScaffoldLink& link, std::string& extensionString) const;

        bool overlapResolve(const ResolveParams& params, const std::string& s1, const std::string& s2, 
                            const ScaffoldLink& link, std::string& outString) const;

        bool introduceGap(int minGapLength, const std::string& contigString, const ScaffoldLink& link, std::string& outString) const;

        void parse(const std::string& text);

    private:
        
        typedef std::vector<ScaffoldLink> LinkVector;
        
        std::string m_rootID;
        LinkVector m_links;
};

