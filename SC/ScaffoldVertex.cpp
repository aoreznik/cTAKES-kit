
ScaffoldVertex::ScaffoldVertex(VertexID id, size_t seqLen) : m_id(id), 
                                                             m_seqLen(seqLen), 
                                                             m_AStatistic(0.0f),
                                                             m_estCopyNumber(-1.0f),
                                                             m_classification(SVC_UNKNOWN), 
                                                             m_color(GC_WHITE),
                                                             m_hasConflictingLink(false)
{

}

ScaffoldVertex::~ScaffoldVertex()
{
    for(ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
          iter != m_edges.end(); ++iter)
    {
    }
}

{
    m_edges.push_back(pEdge);
}

{
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
    }

    return NULL;
}

{
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
    }

    return NULL;
}

ScaffoldEdgePtrVector ScaffoldVertex::getEdges()
{
    return m_edges;
}

ScaffoldEdgePtrVector ScaffoldVertex::getEdges(EdgeDir dir)
{
    ScaffoldEdgePtrVector out;
    for(ScaffoldEdgePtrVector::iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
    }
    return out;
}

void ScaffoldVertex::deleteEdges()
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        ++iter;
    }
    m_edges.clear();
}

void ScaffoldVertex::deleteEdgesAndTwins()
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        pEdge->getEnd()->deleteEdge(pEdge->getTwin());
        delete pEdge;
        ++iter;
    }
    m_edges.clear();
}

void ScaffoldVertex::deleteEdgesAndTwins(EdgeDir dir)
{
    ScaffoldEdgePtrVector out;

    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {

        if(pEdge->getDir() == dir || (pEdge->getEnd() == this && pEdge->getTwin()->getDir() == dir))
        {
            if(pEdge->getEnd() != this)
            {
                pEdge->getEnd()->deleteEdge(pEdge->getTwin());
            }
            delete pEdge;
        }
        else
        {
        }
        ++iter;
    }
    
    m_edges = out;
}

{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
            break;
        ++iter;
    }
    assert(iter != m_edges.end());
    m_edges.erase(iter);
    delete pEdge;
}

{
    assert(pEdge != pEdge->getTwin());
    pEdge->getEnd()->deleteEdge(pEdge->getTwin());
    deleteEdge(pEdge);
}

void ScaffoldVertex::deleteEdgesAndTwinsByColor(GraphColor c)
{
    ScaffoldEdgePtrVector delVec;
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        ++iter;
    }

    iter = delVec.begin();
    for(; iter != delVec.end(); ++iter)
    {
    }
}

void ScaffoldVertex::markEdgesInDir(EdgeDir dir, GraphColor c)
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        if(pEdge->getDir() == dir)
            pEdge->setColor(c);
        ++iter;
    }
}

void ScaffoldVertex::setEdgeColors(GraphColor c)
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        ++iter;
    }
}

void ScaffoldVertex::setAStatistic(double v)
{
    m_AStatistic = v;
}

void ScaffoldVertex::setEstCopyNumber(double v)
{
    m_estCopyNumber = v;
}

void ScaffoldVertex::setClassification(ScaffoldVertexClassification classification)
{
    m_classification = classification;
}

void ScaffoldVertex::setColor(GraphColor c)
{
    m_color = c;
}

void ScaffoldVertex::setConflictingFlag(bool b)
{
    m_hasConflictingLink = b;
}

VertexID ScaffoldVertex::getID() const
{
    return m_id;
}

bool ScaffoldVertex::isRepeat() const
{
    return m_classification == SVC_REPEAT;
}
size_t ScaffoldVertex::getNumEdges() const
{
    return m_edges.size();
}

size_t ScaffoldVertex::getSeqLen() const
{
    return m_seqLen;
}

double ScaffoldVertex::getAStatistic() const
{
    return m_AStatistic;
}

double ScaffoldVertex::getEstCopyNumber() const
{
    return m_estCopyNumber;
}

ScaffoldVertexClassification ScaffoldVertex::getClassification() const
{
    return m_classification;
}

GraphColor ScaffoldVertex::getColor() const
{
    return m_color;
}

bool ScaffoldVertex::hasConflictingLink() const
{
    return m_hasConflictingLink;
}

std::string ScaffoldVertex::getColorString() const
{
    switch(m_classification)
    {
        case SVC_UNKNOWN:
            return "gray";
        case SVC_UNIQUE:
            return "white";
        case SVC_REPEAT:
            return "red";
        default:
            return "white";
    }
}

{
   VertexID id = getID();
   writeEdgesDot(pWriter);
}

{
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
    }
}
