


struct LinkVertexPair
{
    ScaffoldLink link;
};

typedef std::list<LinkVertexPair> LinkList;
typedef LinkList::iterator LinkListIterator;
typedef std::vector<ScaffoldLink> LinkVector;
typedef std::vector<LinkVertexPair> LinkPairVector;
typedef LinkPairVector::iterator LinkVectorPairIterator;
typedef LinkPairVector::const_iterator LinkVectorPairConstIterator;

class ScaffoldGroup
{
    public:


        bool isOrderAmbiguous();
        bool markPolymorphic(double p_cutoff, double cn_cutoff);

        bool hasConsistentLayout();

        bool areLinksAmbiguous(const ScaffoldLink& linkA,
                               const ScaffoldLink& linkB,
                               double p);

        void getSecondaryLinks();

        void computeBestOrdering();
        std::string getBestOrderingString() const;

        int calculateLongestOverlap();
        double calculateProbACloserThanB(const ScaffoldLink& linkA,
                                         const ScaffoldLink& linkB);

        
        void getLinearLinks(LinkVector& outLinks);

    private:

        double normCDF(double x, double mean, double sd);
        int scoreLinkPlacement(const ScaffoldLink& link, const LinkList& unplacedList);

        int m_maxOverlap;
        LinkPairVector m_links;

        bool m_isOrdered;

};

