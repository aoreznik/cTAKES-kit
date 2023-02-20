
{
    pGraph->setEdgeColors(GC_BLACK);

    ScaffoldConnectedComponents connectedComponents;
    ScaffoldAlgorithms::connectedComponents(pGraph, connectedComponents);
 
    for(size_t i = 0; i < connectedComponents.size(); ++i)
    {
        ScaffoldVertexPtrVector& component = connectedComponents[i];

        if(component.size() == 1) 

        ScaffoldVertexPtrVector terminalVertices;
        computeTerminalsForConnectedComponent(component, terminalVertices);

        if(terminalVertices.empty())
        {
            std::cerr << "Warning: scaffold component of size " << component.size() << 
                         " does not have a terminal vertex. Skipping\n";
            continue;
        }
        size_t bestLayoutBases = 0;
        ScaffoldWalk bestWalk(NULL);

        for(size_t j = 0; j < terminalVertices.size(); j++)
        {
            ScaffoldWalkVector layouts;
            computeLayout(terminalVertices[j], layouts);

            for(size_t k = 0; k < layouts.size(); ++k)
            {
                size_t layoutSum = layouts[k].getContigLengthSum();

                if(layoutSum > bestLayoutBases)
                {
                    bestLayoutBases = layoutSum;
                    bestWalk = layouts[k];
                }
            }
        }

        assert(bestLayoutBases > 0 && bestWalk.getStartVertex() != NULL);

        ScaffoldEdgePtrVector keptEdges = bestWalk.getEdges();
        for(size_t j = 0; j < keptEdges.size(); ++j)
        {
            keptEdges[j]->setColor(GC_WHITE);
            keptEdges[j]->getTwin()->setColor(GC_WHITE);
        }
    }

    pGraph->deleteEdgesByColor(GC_BLACK);
}

{
    bool done = false;
    while(!done)
    {
        bool cycleFound = false;
        pGraph->setVertexColors(GC_WHITE);

        ScaffoldConnectedComponents connectedComponents;
        ScaffoldAlgorithms::connectedComponents(pGraph, connectedComponents);

        for(size_t i = 0; i < connectedComponents.size(); ++i)
        {
            ScaffoldVertexPtrVector& component = connectedComponents[i];

            if(component.size() == 1) 

            ScaffoldVertexPtrVector terminalVertices;
            computeTerminalsForConnectedComponent(component, terminalVertices);     

            for(size_t j = 0; j < terminalVertices.size(); j++)
            {
                if(pBackEdge != NULL)
                {
                    std::cout << "Internal cycle found between: " << pBackEdge->getStartID() << " and " << pBackEdge->getEndID() << "\n";

                    pBackEdge->getStart()->setClassification(SVC_REPEAT);
                    pBackEdge->getEnd()->setClassification(SVC_REPEAT);
                    cycleFound = true;
                    break;
                }
            }
        }

        if(cycleFound)
        {
            pGraph->deleteVertices(SVC_REPEAT);
        }
        else
        {
            done = true;
        }
    }
}

{
    std::ofstream cycle_writer(out_filename.c_str());

    bool done = false;
    while(!done)
    {
        done = true;
        pGraph->setVertexColors(GC_WHITE);
        pGraph->setEdgeColors(GC_WHITE);

        ScaffoldConnectedComponents connectedComponents;
        ScaffoldAlgorithms::connectedComponents(pGraph, connectedComponents);

        for(size_t i = 0; i < connectedComponents.size(); ++i)
        {
            ScaffoldVertexPtrVector& component = connectedComponents[i];

            if(component.size() == 1)
                continue;
            
            assert(!component.empty());

            ScaffoldVertexVector cycle_vertices = checkForStrictCycle(test_vertex);
            if(!cycle_vertices.empty())
            {

                for(size_t j = 0; j < cycle_vertices.size(); ++j)
                {
                    cycle_vertices[j]->deleteEdgesAndTwins();
                    cycle_writer << cycle_vertices[j]->getID() << "\n";
                }
                
                done = false;
            }
        }
    }
}

{
    ScaffoldVertexPtrVector allVertices = pGraph->getAllVertices();
    ScaffoldSearchTree::connectedComponents(allVertices, outComponents);
}

void ScaffoldAlgorithms::computeTerminalsForConnectedComponent(const ScaffoldVertexPtrVector& component, 
                                                               ScaffoldVertexPtrVector& terminals)
{
    for(size_t i = 0; i < component.size(); ++i)
    {
        size_t asCount = pVertex->getEdges(ED_ANTISENSE).size();
        size_t sCount = pVertex->getEdges(ED_SENSE).size();

        if(asCount == 0 || sCount == 0)
            terminals.push_back(pVertex);
    }
}


{
    size_t asCount = pStartVertex->getEdges(ED_ANTISENSE).size();
    size_t sCount = pStartVertex->getEdges(ED_SENSE).size();
    assert(asCount == 0 || sCount == 0);
    if(asCount == 0 && sCount == 0)
    EdgeDir dir = (asCount > 0) ? ED_ANTISENSE : ED_SENSE;

    LayoutDistanceMap distanceMap;
    LayoutEdgeMap edgeMap;
    LayoutTerminalSet terminalSet;

    distanceMap[pStartVertex] = 0;
    edgeMap[pStartVertex] = NULL;

    LayoutQueue visitQueue;
    ScaffoldEdgePtrVector startEdges = pStartVertex->getEdges(dir);
    for(size_t i = 0; i < startEdges.size(); ++i)
    {
        LayoutNode initial;
        initial.pEdge = startEdges[i];
        initial.distance = startEdges[i]->getDistance();

        distanceMap[pZ] = initial.distance;
        edgeMap[pZ] = startEdges[i];

        visitQueue.push(initial);
    }

    while(!visitQueue.empty())
    {
        LayoutNode node = visitQueue.front();
        visitQueue.pop();

        EdgeDir yDir = !pXY->getTwin()->getDir();
        ScaffoldEdgePtrVector yEdges = pY->getEdges(yDir);
        if(yEdges.empty())
        {
            terminalSet.insert(pY);
        }

        for(size_t i = 0; i < yEdges.size(); ++i)
        {

            int zDistance = node.distance + pYZ->getDistance();

            LayoutDistanceMap::iterator distIter = distanceMap.find(pZ);
            if(distIter == distanceMap.end() || distIter->second > zDistance)
            {
                distanceMap[pZ] = zDistance;
                edgeMap[pZ] = pYZ;
                LayoutNode updateNode;
                updateNode.pEdge = pYZ;
                updateNode.distance = zDistance;
                visitQueue.push(updateNode);
            }
        }
    }

    LayoutTerminalSet::iterator terminalIter = terminalSet.begin();
    for(; terminalIter != terminalSet.end(); ++terminalIter)
    {
        ScaffoldEdgePtrVector reverseEdges;
        while(pCurrent != pStartVertex)
        {
            reverseEdges.push_back(pEdge);
            assert(pEdge != NULL);
            pCurrent = pEdge->getStart();
        }

        ScaffoldWalk walk(pStartVertex);
        ScaffoldEdgePtrVector::reverse_iterator rIter = reverseEdges.rbegin();
        for(; rIter != reverseEdges.rend(); ++rIter)
        outWalks.push_back(walk);
    }
}

{
    size_t asCount = pVertex->getEdges(ED_ANTISENSE).size();
    size_t sCount = pVertex->getEdges(ED_SENSE).size();
    assert(asCount == 0 || sCount == 0);
    if(asCount == 0 && sCount == 0)
    EdgeDir dir = (asCount > 0) ? ED_ANTISENSE : ED_SENSE;

    assert(pVertex->getColor() == GC_WHITE);
    LayoutEdgeMap predMap;
    predMap[pVertex] = NULL;

    for(LayoutEdgeMap::iterator iter = predMap.begin(); iter != predMap.end(); ++iter)
        iter->first->setColor(GC_WHITE);
    return pBackEdge;
}

{
    pVertex->setColor(GC_GRAY);
    ScaffoldEdgePtrVector edges = pVertex->getEdges(dir);
    for(size_t i = 0; i < edges.size(); ++i)
    {
        if(pEnd->getColor() == GC_GRAY)
        {
            return pEdge;
        }

        if(pEnd->getColor() == GC_WHITE)
        {
            predMap[pEnd] = pEdge;
            EdgeDir continueDir = !pEdge->getTwin()->getDir();
            if(pBackEdge != NULL)
                return pBackEdge;
        }
    }
    pVertex->setColor(GC_BLACK);
    return NULL;
}

{
    assert(pVertex->getColor() == GC_WHITE);
    
    LayoutEdgeMap predMap;
    predMap[pVertex] = NULL;

    ScaffoldVertexVector out;
    if(pCycleEdge != NULL)
    {
        while(current != NULL)
        {
            out.push_back(current->getStart());
            current = predMap[current->getStart()];
        }
    }

    for(LayoutEdgeMap::iterator iter = predMap.begin(); iter != predMap.end(); ++iter)
    {
        iter->first->setColor(GC_WHITE);

        if(iter->second != NULL)
        {
            iter->second->setColor(GC_WHITE);
            iter->second->getTwin()->setColor(GC_WHITE);
        }
    }

    return out;
}

{
    pVertex->setColor(GC_GRAY);
    ScaffoldEdgePtrVector edges = pVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
    {

        if(pEdge->getColor() == GC_RED)
            continue;
        
        if(pEnd->getColor() == GC_GRAY)
            return pEdge;

        if(pEnd->getColor() == GC_WHITE)
        {
            pEdge->setColor(GC_RED);
            pEdge->getTwin()->setColor(GC_RED);
            predMap[pEnd] = pEdge;
            
            if(pBackEdge != NULL)
                return pBackEdge;
        }
    }
    pVertex->setColor(GC_BLACK);
    return NULL;
}
