

struct ScaffoldDistanceFunction
{
    {
        return pEdge->getDistance();
    }
};

typedef GraphSearchTree<ScaffoldVertex, ScaffoldEdge, ScaffoldDistanceFunction> ScaffoldSearchTree;

struct ScaffoldWalkBuilder
{
    public:
        ScaffoldWalkBuilder(ScaffoldWalkVector& outWalks);
        ~ScaffoldWalkBuilder();

        void finishCurrentWalk();

    private:
        ScaffoldWalkVector& m_outWalks;
};

namespace ScaffoldSearch
{
                          EdgeDir initialDir, 
                          int maxDistance,
                          size_t maxWalks, 
                          ScaffoldWalkVector& outWalks);

    
    int findCoveringWalk(const ScaffoldWalkVector& allWalks, 
                         const ScaffoldVertexPtrVector& coverVector);

                          EdgeDir intialDir,
                          int maxDistance,
                          size_t maxNodes, 
                          ScaffoldWalkVector& outWalks);

    void printWalks(const ScaffoldWalkVector& walkVector);

};

