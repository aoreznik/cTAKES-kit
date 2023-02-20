
ScaffoldGraph::ScaffoldGraph()
{
    m_vertices.set_deleted_key("");
}

ScaffoldGraph::~ScaffoldGraph()
{
    for(ScaffoldVertexMap::iterator iter = m_vertices.begin();
         iter != m_vertices.end(); ++iter)
    {
        iter->second->deleteEdges();
        delete iter->second;
        iter->second = NULL;
    }
}

{
    m_vertices.insert(std::make_pair(pVertex->getID(), pVertex));
}

{
    assert(pVertex != NULL);
    pVertex->addEdge(pEdge);
}

ScaffoldVertexPtrVector ScaffoldGraph::getAllVertices() const
{
    ScaffoldVertexPtrVector outVertices;
    ScaffoldVertexMap::const_iterator iter = m_vertices.begin();
    for(; iter != m_vertices.end(); ++iter)
        outVertices.push_back(iter->second);
    return outVertices;
}

{
    ScaffoldVertexMap::const_iterator iter = m_vertices.find(id);
    if(iter == m_vertices.end())
        return NULL;
    return iter->second;
}

void ScaffoldGraph::deleteVertices(ScaffoldVertexClassification classification)
{
    ScaffoldVertexMap::iterator iter = m_vertices.begin(); 
    while(iter != m_vertices.end())
    {
        if(iter->second->getClassification() == classification)
        {
            iter->second->deleteEdgesAndTwins();
            delete iter->second;
            iter->second = NULL;
            m_vertices.erase(iter++);
        }
        else
        {
            ++iter;
        }
    }
}

void ScaffoldGraph::deleteEdgesByColor(GraphColor c)
{
    ScaffoldVertexMap::iterator iter = m_vertices.begin(); 
    while(iter != m_vertices.end())
    {
        iter->second->deleteEdgesAndTwinsByColor(c);
        ++iter;
    }
}



void ScaffoldGraph::setVertexColors(GraphColor c)
{
    ScaffoldVertexMap::iterator iter = m_vertices.begin(); 
    while(iter != m_vertices.end())
    {
        iter->second->setColor(c);
        ++iter;
    }
}

void ScaffoldGraph::setEdgeColors(GraphColor c)
{
    ScaffoldVertexMap::iterator iter = m_vertices.begin(); 
    while(iter != m_vertices.end())
    {
        iter->second->setEdgeColors(c);
        ++iter;
    }
}

void ScaffoldGraph::loadVertices(const std::string& filename, int minLength)
{
    SeqReader reader(filename, SRF_NO_VALIDATION | SRF_KEEP_CASE);
    SeqRecord sr;
    while(reader.get(sr))
    {
        int contigLength = sr.seq.length();
        if(contigLength >= minLength)
        {
            addVertex(pVertex);
        }
    }    
}

void ScaffoldGraph::loadDistanceEstimateEdges(const std::string& filename, bool isMatePair, int verbose)
{
    std::cout << "Reading distance estimates from " << filename << "\n";
    std::string line;

    {
        assert(line.substr(0,4) != "Mate");
        StringVector fields = split(line, ' ');
        assert(fields.size() >= 1);

        std::string rootID = fields[0];

        for(size_t i = 1; i < fields.size(); ++i)
        {
            std::string record = fields[i];
            if(record == ";")
            {
                currDir = !currDir;
                continue;
            }

            std::string id;
            EdgeComp comp;
            int distance;
            int numPairs;
            double stdDev;
            parseDERecord(record, id, comp, distance, numPairs, stdDev);


            if(pVertex1 != NULL && pVertex2 != NULL)
            {
                if(pVertex1 == pVertex2)
                {
                    std::cout << "Self-edges not allowed\n";
                    assert(false);
                    continue;
                }

                ScaffoldLink link1(id, currDir, comp, distance, stdDev, numPairs, pVertex2->getSeqLen(), SLT_DISTANCEEST);
                ScaffoldLink link2(rootID, !correctDir(currDir, comp), comp, distance, stdDev, numPairs, pVertex1->getSeqLen(), SLT_DISTANCEEST);

                if(pEdge != NULL)
                {
                    if(!isMatePair && pEdge->getLink().stdDev < stdDev)
                    {
                        pEdge->setLink(link1);
                        pEdge->getTwin()->setLink(link2);
                    }
                    else
                    {
                        if(abs(pEdge->getDistance() - link1.distance) > 100)
                        {
                            if(verbose >= 1)
                            {
                                printf("LL skipped from %s to %s. Distance1: %d Distance2: %d\n", pVertex1->getID().c_str(), 
                                                                                                  link1.endpointID.c_str(), 
                                                                                                  link1.distance, 
                                                                                                  pEdge->getDistance());
                            }
                            pVertex1->setConflictingFlag(true);
                            pVertex2->setConflictingFlag(true);
                        }
                    }
                }
                else
                {

                    pEdge1->setTwin(pEdge2);
                    pEdge2->setTwin(pEdge1);

                    addEdge(pVertex1, pEdge1);
                    addEdge(pVertex2, pEdge2);
                }
            }
        }
    }

    delete pReader;
}

void ScaffoldGraph::loadAStatistic(const std::string& filename)
{
    std::string line;

    {
        StringVector fields = split(line, '\t');
        assert(fields.size() == 6);

        VertexID id = fields[0];

        std::stringstream cn_parser(fields[4]);
        double cn;
        cn_parser >> cn;

        std::stringstream as_parser(fields[5]);
        double as;
        as_parser >> as;

        if(pVertex != NULL)
        {
            pVertex->setAStatistic(as);
            pVertex->setEstCopyNumber(cn);
        }
    }
    delete pReader;
}


void ScaffoldGraph::parseDERecord(const std::string& record, std::string& id, 
                                  EdgeComp& comp, int& distance, int& numPairs, double& stdDev)
{
    StringVector fields = split(record, ',');
    if(fields.size() != 4)
    {
        std::cerr << "Distance Estimate record is not formatted correctly: " << record << "\n";
        exit(1);
    }

    id = fields[0].substr(0, fields[0].size() - 1);
    comp = (fields[0][fields[0].size() - 1] == '+' ? EC_SAME : EC_REVERSE);

    std::stringstream d_parser(fields[1]);
    d_parser >> distance;
    
    std::stringstream np_parser(fields[2]);
    np_parser >> numPairs;

    std::stringstream sd_parser(fields[3]);
    sd_parser >> stdDev;
}

void ScaffoldGraph::writeDot(const std::string& outFile) const
{
    
    std::string graphType = "digraph";

    ScaffoldVertexMap::const_iterator iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        iter->second->writeDot(pWriter);
    }
    delete pWriter;
}
