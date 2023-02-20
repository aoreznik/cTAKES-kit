

ScaffoldRecord::ScaffoldRecord() 
{

}

void ScaffoldRecord::setRoot(const std::string& root)
{
    m_rootID = root;
}

void ScaffoldRecord::addLink(const ScaffoldLink& link)
{
    m_links.push_back(link);
}

size_t ScaffoldRecord::getNumComponents() const
{
    if(m_rootID.empty())
        return 0;
    
    return 1 + m_links.size();
}

std::string ScaffoldRecord::generateString(const ResolveParams& params, StringVector& ids) const
{
    assert(params.pSequenceCollection != NULL);

    params.pStats->numScaffolds += 1;

    std::string sequence = params.pSequenceCollection->getSequence(m_rootID);
    params.pSequenceCollection->setPlaced(m_rootID);
    ids.push_back(m_rootID + "+");
 
    if(m_links.empty())
        return sequence;

    EdgeDir rootDir = m_links[0].getDir();
    EdgeComp relativeComp = EC_SAME;
    EdgeComp prevComp = EC_SAME;
    std::string currID = m_rootID;

    bool reverseAll = (rootDir == ED_ANTISENSE);
    if(reverseAll)
        sequence = reverse(sequence);

    for(size_t i = 0; i < m_links.size(); ++i)
    {   
        params.pStats->numGapsAttempted += 1;

        const ScaffoldLink& link = m_links[i];
        params.pSequenceCollection->setPlaced(link.endpointID);

        if(link.getComp() == EC_REVERSE)
            relativeComp = !relativeComp;

        std::string resolvedSequence;
        
        bool resolved = false;
        if(params.resolveMask & RESOLVE_GRAPH_BEST || params.resolveMask & RESOLVE_GRAPH_UNIQUE)
        {
            resolved = graphResolve(params, currID, link, resolvedSequence);
            if(resolved)
            {
                if(prevComp == EC_REVERSE)
                    resolvedSequence = reverseComplementIUPAC(resolvedSequence);

                if(reverseAll)
                    resolvedSequence = reverse(resolvedSequence);
                params.pStats->numGapsResolved += 1;
            }
        }

        if(!resolved)
        {
            std::string toAppend = params.pSequenceCollection->getSequence(link.endpointID);
            if(relativeComp == EC_REVERSE)
                toAppend = reverseComplementIUPAC(toAppend);
            if(reverseAll)
                toAppend = reverse(toAppend);
            
            if(link.distance < 0 && params.resolveMask & RESOLVE_OVERLAP)
            {
                resolved = overlapResolve(params, sequence, toAppend, link, resolvedSequence);
                if(resolved)
                {
                    params.pStats->numGapsResolved += 1;
                    params.pStats->overlapFound += 1;
                }
                else
                {
                    params.pStats->overlapFailed += 1;
                }
            }

            if(!resolved)
                introduceGap(params.minGapLength, toAppend, link, resolvedSequence);

        }

        sequence.append(resolvedSequence);
        currID = link.endpointID;

        std::string outID = currID;
        outID.append(relativeComp == EC_SAME ? "+" : "-");
        ids.push_back(outID);
        prevComp = relativeComp;
    }

    if(reverseAll) 
    {
        sequence = reverse(sequence);
        std::reverse(ids.begin(), ids.end());
    }
    return sequence;
}

bool ScaffoldRecord::graphResolve(const ResolveParams& params, const std::string& startID, 
                                  const ScaffoldLink& link, std::string& outExtensionString) const
{
    assert(params.pGraph != NULL);

    assert(pStartVertex != NULL && pEndVertex != NULL);

    int maxDistance = link.distance + threshold;
    int maxExtensionDistance = maxDistance + pEndVertex->getSeqLen();
    SGWalkVector walks;
    SGSearch::findWalks(pStartVertex, pEndVertex, link.getDir(), maxExtensionDistance, 10000, true, walks);

    int numWalksValid = 0;
    int numWalksClosest = 0;
    int selectedIdx = -1;
    int closestDist = std::numeric_limits<int>::max();

            std::cout << "Attempting graph resolve of link " << startID << " -- " << link.endpointID << " expected distance: " << link.distance << " orientation: " << link.edgeData.getComp() << "\n";
    
    for(size_t i = 0; i < walks.size(); ++i)
    {
        std::vector<EdgeComp> vertexOrientations = walks[i].getOrientationsToStart();
        assert(walks[i].getLastEdge()->getEndID() == link.endpointID);

        if(vertexOrientations.back() != link.edgeData.getComp())
        {
            std::cout << "SKIPPING WALK OF THE WRONG ORIENTATION\n";
            continue;
        }

        int walkDistance = walks[i].getEndToStartDistance();
        int diff = abs(abs(link.distance - walkDistance));
        if(diff <= threshold)
        {

            std::cout << "  Walk distance: " << walkDistance << " diff: " << diff << " threshold: " << threshold << " close: " << closestDist << "\n";
            ++numWalksValid;
            if(diff < closestDist)
            {
                selectedIdx = i;
                closestDist = diff;
                numWalksClosest = 1;
            }
            else if(diff == closestDist)
            {
                numWalksClosest += 1;
            }

        }
    }

    bool useWalk = false;

    if(numWalksValid > 0)
    {
        if(params.resolveMask & RESOLVE_GRAPH_BEST)
        {
            if(!(params.resolveMask & RESOLVE_GRAPH_UNIQUE) || numWalksClosest == 1)
                useWalk = true;
            else if((params.resolveMask & RESOLVE_GRAPH_UNIQUE) && numWalksClosest > 1)
                params.pStats->graphWalkTooMany += 1;
        }
        else
        {
            if(numWalksValid == 1)
                useWalk = true;
            else if(numWalksValid > 1)
                params.pStats->graphWalkTooMany += 1;
        }
    }

    std::cout << "  Num walks: " << walks.size() << " Num valid: " << numWalksValid << " Num closest: " << numWalksClosest << " using: " << useWalk << "\n";

    if(useWalk)
    {
        assert(selectedIdx != -1);
        outExtensionString = walks[selectedIdx].getString(SGWT_EXTENSION);
        params.pStats->graphWalkFound += 1;

        VertexPtrVec vertexPtrVector = walks[selectedIdx].getVertices();
        for(size_t i = 0; i < vertexPtrVector.size(); ++i)
            params.pSequenceCollection->setPlaced(vertexPtrVector[i]->getID());
        return true;
    }
    else
    {
        if(numWalksValid == 0)
            params.pStats->graphWalkNoPath += 1;
        assert(outExtensionString.empty());
        return false;
    }
}

bool ScaffoldRecord::overlapResolve(const ResolveParams& params, const std::string& s1, const std::string& s2, 
                                    const ScaffoldLink& link, std::string& outString) const
{

    std::cout << "Attempting overlap resolve of link to " << link.endpointID << " expected distance: " << link.distance << " orientation: " << link.edgeData.getComp() << "\n";


    int upperBound = 0;
    if(params.maxOverlap == -1)
    else
        upperBound = params.maxOverlap;
    
    Match match;
    bool overlapFound = OverlapTools::boundedOverlapDP(s1, s2, params.minOverlap, upperBound, params.maxErrorRate, match);
    if(overlapFound)
    {
        std::cout << "Overlap found, length: " << match.coord[1].length() << "\n";
        SeqCoord overlapCoord = match.coord[1];
        SeqCoord overhangCoord = overlapCoord.complement();
        outString = overhangCoord.getSubstring(s2);
        return true;
    }
    else
    {
        return false;
    }
}

bool ScaffoldRecord::introduceGap(int minGapLength, const std::string& contigString, const ScaffoldLink& link, std::string& out) const
{
    assert(out.empty());
    if(link.distance < 0)
    {
        out.append(minGapLength, 'N');
        assert(expectedOverlap < (int)contigString.length());
        out.append(contigString.substr(expectedOverlap));
    }
    else
    {
        int gap = std::max(link.distance, minGapLength);
        out.append(gap, 'N');
        out.append(contigString);
    }
    return true;
}

void ScaffoldRecord::parse(const std::string& text)
{
    StringVector fields = split(text, '\t');
    assert(fields.size() >= 1);

    m_rootID = fields[0];
    for(size_t i = 1; i < fields.size(); ++i)
    {
        ScaffoldLink link;
        link.parse(fields[i]);
        m_links.push_back(link);
    }
}

{
    for(size_t i = 0; i < m_links.size(); ++i)
}
