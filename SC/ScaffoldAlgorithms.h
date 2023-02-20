

typedef std::vector<ScaffoldVertexPtrVector> ScaffoldConnectedComponents;

namespace ScaffoldAlgorithms
{




    void computeTerminalsForConnectedComponent(const ScaffoldVertexPtrVector& component, 
                                               ScaffoldVertexPtrVector& terminals);

    
    struct LayoutNode
    {
    };

    typedef std::queue<LayoutNode> LayoutQueue;


};

