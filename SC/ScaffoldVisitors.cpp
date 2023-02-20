
{
    const double AMBIGUOUS_TOLERANCE = 2.0f;
    assert(pXY->getDistance() <= pXZ->getDistance());
    if(translated < pXY->getDistance())
        return true;
    else
        return false;
}

                                             int& dist, double& sd, EdgeDir& dir_yz, 
                                             EdgeDir& dir_zy, EdgeComp& comp)
{
    comp = (pXY->getComp() == pXZ->getComp() ? EC_SAME : EC_REVERSE);
    dist = pXZ->getDistance() - (pXY->getDistance() + pXY->getEnd()->getSeqLen());
    sd = sqrt(pow(pXY->getStdDev(), 2.0) + pow(pXZ->getStdDev(), 2.0));
}

{
    m_numVertices = 0;
    m_numEdges = 0;
}

{
    ++m_numVertices;
    m_numEdges += pVertex->getNumEdges();

    return false;
}

{
    printf("Scaffold Stats -- Num vertices: %zu Num edges: %zu\n", m_numVertices, m_numEdges);
}

ScaffoldAStatisticVisitor::ScaffoldAStatisticVisitor(double uniqueThreshold, 
                                                     double minCopyNumber) : m_uniqueThreshold(uniqueThreshold),
                                                                             m_minCopyNumber(minCopyNumber)

{
                                                                               
}

{
    m_sumUnique = 0;
    m_sumRepeat = 0;
    m_sumLowCN = 0;
    m_numUnique = 0;
    m_numRepeat = 0;
    m_numLowCN = 0;
}

{
    if(pVertex->getClassification() != SVC_REPEAT)
    {
        if(pVertex->getAStatistic() > m_uniqueThreshold)
        {
            pVertex->setClassification(SVC_UNIQUE);
            ++m_numUnique;
            m_sumUnique += pVertex->getSeqLen();
        }
        else
        {
            pVertex->setClassification(SVC_REPEAT);
            ++m_numRepeat;
            m_sumRepeat += pVertex->getSeqLen();
        }

        if(pVertex->getEstCopyNumber() < m_minCopyNumber)
        {
            pVertex->setClassification(SVC_REPEAT);
            ++m_numLowCN;
            m_sumLowCN += pVertex->getSeqLen();
        }
    }
    return false;   
}

{
}

ScaffoldLinkValidator::ScaffoldLinkValidator(int maxOverlap, 
                                             double threshold,
                                             int verbose) : m_maxOverlap(maxOverlap),
                                                            m_threshold(threshold),
                                                            m_verbose(verbose)
{

}

{
    m_numUnique = 0;
    m_numRepeat = 0;
    m_numCut = 0;
    pGraph->setEdgeColors(GC_WHITE);
}

{
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            ScaffoldGroup group(pVertex, m_maxOverlap);
            for(size_t i = 0; i < edgeVec.size(); ++i)
                group.addLink(edgeVec[i]->getLink(), edgeVec[i]->getEnd());
        
            bool isAmbiguous = group.isOrderAmbiguous();
            int longestOverlap = group.calculateLongestOverlap();
            group.computeBestOrdering();
            bool overlapCheck = longestOverlap < 400;
            bool result = overlapCheck;

            if(m_verbose >= 1)
            {
                std::string ambiStr = isAmbiguous ? "ambiguous" : "unambiguous";
                std::string overStr = overlapCheck ? "good-ordering" : "no-ordering";
                std::string resultStr = result ? "PASS" : "FAIL";
                std::string orderStr = group.getBestOrderingString();
                printf("LV %s %d CL:%d EC:%2.2lf AS:%2.2lf %s %s %s LO:%d BO:%s\n", pVertex->getID().c_str(), 
                                                                                  (int)idx, 
                                                                                  (int)pVertex->getSeqLen(), 
                                                                                  pVertex->getEstCopyNumber(),
                                                                                  pVertex->getAStatistic(),
                                                                                  ambiStr.c_str(), 
                                                                                  overStr.c_str(), 
                                                                                  resultStr.c_str(), 
                                                                                  longestOverlap,
                                                                                  orderStr.c_str());
            }

            if(!result)
            {
                for(size_t i = 0; i < edgeVec.size(); ++i)
                {
                    edgeVec[i]->setColor(GC_BLACK);
                    edgeVec[i]->getEnd()->markEdgesInDir(edgeVec[i]->getTwin()->getDir(), GC_BLACK);
                }
                m_numCut += 1;
            }
        }
    }
    return false;   
}

{
    printf("ScaffoldLinkValidator removed edges for %zu vertices\n", m_numCut);
    pGraph->deleteEdgesByColor(GC_BLACK);
}

ScaffoldChainVisitor::ScaffoldChainVisitor(int maxOverlap) : m_maxOverlap(maxOverlap)
{
}

{
}

{
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    bool changed_graph = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            
            std::sort(edgeVec.begin(), edgeVec.end(), ScaffoldEdgePtrDistanceCompare);
            assert(pXY->getDistance() <= pXZ->getDistance());

            if(pY->isRepeat() || pZ->isRepeat())
            {
                std::cerr << "Warning, skipping repeat\n";
                continue;
            }



            bool isAmbiguous = ScaffoldAlgorithms::areEdgesAmbiguous(pXY, pXZ);
            if(isAmbiguous)
            {
                              " are ambiguously ordered\n";
            }

            int dist;
            double sd;
            EdgeDir dir_yz;
            EdgeDir dir_zy; 
            EdgeComp comp;
            ScaffoldAlgorithms::inferScaffoldEdgeYZ(pXY, pXZ, dist, sd, dir_yz, dir_zy, comp);

            if(!isConsistent)
            {
                std::cout << "\tEdge is not consistent with max overlap: " << dist << "\n";
            }

            if(pCheckEdge == NULL)
            {
                ScaffoldLink linkYZ(pZ->getID(), dir_yz, comp, dist, sd, 0, pZ->getSeqLen(), SLT_INFERRED);
                ScaffoldLink linkZY(pY->getID(), dir_zy, comp, dist, sd, 0, pY->getSeqLen(), SLT_INFERRED);
                pYZ->setTwin(pZY);
                pZY->setTwin(pYZ);

                pY->addEdge(pYZ);
                pZ->addEdge(pZY);
            }
            else
            {
            }


            
            pXZ = 0;
            pZX = 0;
            changed_graph = true;
        }
    }

    return changed_graph;
}

{

}

{
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            pVertex->deleteEdgesAndTwins(dir);
        }
    }
    return false;
}

ScaffoldTransitiveReductionVisitor::ScaffoldTransitiveReductionVisitor()
{

}

{
    pGraph->setEdgeColors(GC_WHITE);
}

{
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    bool changed_graph = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() <= 1)
            continue;

        ScaffoldWalkVector walkVector;
        ScaffoldSearch::findVariantWalks(pVertex, dir, 1000000, 100, walkVector);

        int walkIdx = -1;
        int lowestIdxInWalk = std::numeric_limits<int>::max();
        int lowestIdxInVec = -1;
        for(size_t i = 0; i < walkVector.size(); ++i)
        {
            bool allContained = true;
            lowestIdxInWalk = std::numeric_limits<int>::max();
            lowestIdxInVec = -1;

            for(size_t j = 0; j < edgeVec.size(); ++j)
            {
                EdgeComp expected_orientation = edgeVec[j]->getComp();
                int idxInWalk = walkVector[i].findVertex(edgeVec[j]->getEnd());

                if(idxInWalk == -1)
                {
                    allContained = false;
                    break;
                }

                EdgeComp walk_orientation = walkVector[i].findOrientation(edgeVec[j]->getEnd());
                if(walk_orientation != expected_orientation)
                {
                    allContained = false;
                    break;
                }
                
                if(idxInWalk < lowestIdxInWalk)
                {
                    lowestIdxInWalk = idxInWalk;
                    lowestIdxInVec = j;
                }
            }

            if(allContained)
            {
                walkIdx = i;
                break;
            }
        }

        if(walkIdx == -1)
        {
            continue;
        }

        std::cout << "Walk " << walkIdx << " contains all links ";
        walkVector[walkIdx].print();
        std::cout << "keeping link: " << lowestIdxInVec << " " << edgeVec[lowestIdxInVec]->getLink() << "\n";

        for(int i = 0; i < (int)edgeVec.size(); ++i)
        {
            if(i != lowestIdxInVec)
            {
                edgeVec[i]->setColor(GC_BLACK);
                edgeVec[i]->getTwin()->setColor(GC_BLACK);
                changed_graph = true;
            }
        }
    }
    return changed_graph;
}

{
    pGraph->deleteEdgesByColor(GC_BLACK);
}

ScaffoldPolymorphismVisitor::ScaffoldPolymorphismVisitor(int maxOverlap) : m_maxOverlap(maxOverlap)
{

}

{
    m_numMarked = 0;
}

{
    (void)pGraph;
    (void)pVertex;
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            ScaffoldGroup group(pVertex, m_maxOverlap);

            for(size_t i = 0; i < edgeVec.size(); ++i)
            {
                group.addLink(edgeVec[i]->getLink(), edgeVec[i]->getEnd());
            }

            double p_cutoff = 0.01;
            double cn_cutoff = 1.5f;
            bool nodeMarked = group.markPolymorphic(p_cutoff, cn_cutoff);
            if(nodeMarked)
                m_numMarked++;
            return nodeMarked;
        }
    }    
    return false;
}

{
    printf("Marked %d nodes as polymorphic\n", m_numMarked);
    pGraph->deleteVertices(SVC_POLYMORPHIC);
}


ScaffoldSVVisitor::ScaffoldSVVisitor(int maxSize) : m_maxSVSize(maxSize)
{

}

{
    pGraph->setEdgeColors(GC_WHITE);
    m_numMarked = 0;
}

{
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    bool changed_graph = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() <= 1)
            continue;

        ScaffoldWalkVector walkVector;
        ScaffoldSearch::findVariantWalks(pVertex, dir, 50000, 100, walkVector);

        ScaffoldVertexPtrVector linkedVertices;
        for(size_t i = 0; i < edgeVec.size(); ++i)
            linkedVertices.push_back(edgeVec[i]->getEnd());

        int coveringWalkIdx = ScaffoldSearch::findCoveringWalk(walkVector, linkedVertices);
        if(coveringWalkIdx == -1)
        {
            continue;
        }

        
        int maxDiff = 0;
        int closestIdx = 0;
        int closestDistance = std::numeric_limits<int>::max();
        for(int i = 0; i < (int)edgeVec.size(); ++i)
        {
            int distance = walkVector[coveringWalkIdx].getDistanceToVertex(edgeVec[i]->getEnd());
            int diff = abs(distance - edgeVec[i]->getDistance());
          
            printf("SVW: %s -> %s DE: %d WD: %d DF: %d TL: %.2lf\n", edgeVec[i]->getStartID().c_str(), 
                                                                     edgeVec[i]->getEndID().c_str(),
                                                                     edgeVec[i]->getDistance(), distance, 
                                                                     diff, 

            if(diff > maxDiff)
                maxDiff = diff;

            if(distance < closestDistance)
            {
                closestDistance = distance;
                closestIdx = i;
            }
        }

        if(maxDiff > m_maxSVSize)
            continue;

        for(int i = 0; i < (int)edgeVec.size(); ++i)
        {
            if(i != closestIdx)
            {
                edgeVec[i]->setColor(GC_BLACK);
                edgeVec[i]->getTwin()->setColor(GC_BLACK);
                m_numMarked += 1;
            }
        }
    }
    return changed_graph;
}

{
    pGraph->deleteEdgesByColor(GC_BLACK);
    std::cout << "Structural variation resolver marked " << m_numMarked << " edges\n";
}

{
    m_numMarked = 0;
}

{
    if(pVertex->hasConflictingLink())
    {
        pVertex->setClassification(SVC_REPEAT);
        m_numMarked += 1;
    }
    return false;
}

{
    std::cout << "[conflict] marked: " << m_numMarked << "\n";
}

ScaffoldLayoutVisitor::ScaffoldLayoutVisitor()
{

}
  
{
    pGraph->setEdgeColors(GC_WHITE);
}

{
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() <= 1)
            continue;

        ScaffoldGroup group(pVertex, 100);
        for(size_t i = 0; i < edgeVec.size(); ++i)
            group.addLink(edgeVec[i]->getLink(), edgeVec[i]->getEnd());

        group.computeBestOrdering();
        int longestOverlap = group.calculateLongestOverlap();

        if(longestOverlap > 150)
        {
            std::cout << "cannot scaffold from " << pVertex->getID() << " OL: " << longestOverlap << "\n";
            continue;
        }

        std::cout << "Computing layout for " << pVertex->getID() << "\n";
        LinkVector linearLinks;
        group.getLinearLinks(linearLinks);

        std::cout << "Linearized links:\n";

        bool allResolved = true;
        for(size_t i = 0; i < linearLinks.size(); ++i)
        {
            std::cout << "\t" << linearLinks[i] << "\n";

            if(pRealEdge != NULL)
            {
                std::cout << "\t\tresolved:" << pRealEdge->getLink() << "\n";
            }
            else
            {
                pCurrStartVertex->setClassification(SVC_REPEAT);
                allResolved = false;
                std::cout << "\t\tnotfound\n";
            }

            assert(pCurrEndVertex != NULL);
            pCurrStartVertex = pCurrEndVertex;
        }

        group.getSecondaryLinks();

        if(allResolved)
        {
            std::string firstID = linearLinks[0].endpointID;
            for(size_t i = 0; i < edgeVec.size(); ++i)
            {
                if(edgeVec[i]->getEndID() == firstID)
                {
                    edgeVec[i]->setColor(GC_BLACK);
                    edgeVec[i]->getTwin()->setColor(GC_BLACK);
                }
            }
        }
    }
    return false;
}

{
    pGraph->deleteEdgesByColor(GC_BLACK);
}


{

}

{
}

{
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);

        for(size_t i = 0; i < edgeVec.size(); ++i)
        {
            assert(pX != NULL && pY != NULL);

            SGWalkVector walks;
            SGSearch::findWalks(pX, pY, edgeVec[i]->getDir(), edgeVec[i]->getDistance() + pY->getSeqLen() + 1000, true, 100000, walks); 

            if(walks.size() > 0)
            {
                int closest = std::numeric_limits<int>::max();
                int est = edgeVec[i]->getDistance();

                size_t idx = -1;
                for(size_t j = 0; j < walks.size(); ++j)
                {
                    int diff = abs(walks[j].getEndToStartDistance() - est);
                    if(diff < closest)
                    {
                        closest = diff;
                        idx = j;
                    }
                }

                printf("%s -> %s\t%d\t%d\t%s\n", pX->getID().c_str(), 
                                                 pY->getID().c_str(), 
                                                 edgeVec[i]->getDistance(), 
                                                 walks[idx].getEndToStartDistance(),
                                                 walks[idx].pathSignature().c_str());
            }
            else
            {
            }

        }
    }    
    return false;
}

{
}

ScaffoldWriterVisitor::ScaffoldWriterVisitor(const std::string& filename)
{
    m_pWriter = createWriter(filename);
}

ScaffoldWriterVisitor::~ScaffoldWriterVisitor()
{
    delete m_pWriter;
}

{
    pGraph->setVertexColors(GC_WHITE);
    m_statsVector.clear();
}

{
    if(pVertex->getColor() == GC_RED)

    ScaffoldEdgePtrVector edges = pVertex->getEdges();
    
    if(edges.size() <= 1)
    {
        size_t num_contigs = 0;

        pVertex->setColor(GC_RED);
    
        ScaffoldRecord record;
        record.setRoot(pVertex->getID());
        
        bases += pVertex->getSeqLen();
        span += pVertex->getSeqLen();
        num_contigs += 1;
        
        if(edges.size() == 1)
        {
            while(1)
            {
                record.addLink(pXY->getLink());
                if(pY->getColor() == GC_RED)

                pY->setColor(GC_RED);
                bases += pY->getSeqLen();
                span += pY->getSeqLen() + pXY->getDistance();
                num_contigs += 1;

                EdgeDir nextDir = !pXY->getTwin()->getDir();
                ScaffoldEdgePtrVector nextEdges = pY->getEdges(nextDir);
                
                if(nextEdges.size() == 1)
                    pXY = nextEdges[0];
                else
                    break;
            }
        }
        record.writeScaf(m_pWriter);

        ScaffoldStats stats;
        stats.numContigs = num_contigs;
        stats.bases = bases;
        stats.span = span;
        m_statsVector.push_back(stats);
    }
    return false;
}

{
    size_t totalScaffolds = m_statsVector.size();
    size_t singleton = 0;
    size_t sumSpan = 0;
    size_t sumBases = 0;
    size_t totalContigs = 0;
    std::vector<ScaffoldStats>::iterator iter = m_statsVector.begin();
    for(; iter != m_statsVector.end(); ++iter)
    {
        sumSpan += iter->span;
        sumBases += iter->bases;       
        totalContigs += iter->numContigs;
        if(iter->numContigs == 1)
            ++singleton;
    }


    std::sort(m_statsVector.begin(), m_statsVector.end(), ScaffoldStats::sortSpanDesc);

    size_t spanN50 = -1;
    size_t maxSpan = 0;
    size_t runningSum = 0;
    iter = m_statsVector.begin();
    for(; iter != m_statsVector.end(); ++iter)
    {
        if((size_t)iter->span > maxSpan)
            maxSpan = iter->span;
        runningSum += iter->span;
        if(runningSum > targetSpan)
        {
            spanN50 = iter->span;
            break;
        }
    }

    printf("\n======\n");
    printf("Constructed %zu scaffolds from %zu contigs (%zu singleton scaffolds)\n", totalScaffolds, totalContigs, singleton);
    
    printf("Max scaffold span: %zu bp\n", maxSpan);
    printf("N50 scaffold span: %zu bp\n", spanN50);
}
