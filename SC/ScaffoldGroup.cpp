
                             int maxOverlap) : m_pRootVertex(pRootVertex),
                                               m_maxOverlap(maxOverlap),
                                               m_isOrdered(false)
{

}

{
    LinkVertexPair pair = {link, pVertex};
    m_links.push_back(pair);
}

bool ScaffoldGroup::isOrderAmbiguous()
{
    double ambiguity_p = 0.01f;
    LinkVectorPairIterator i = m_links.begin();
    for(; i != m_links.end(); ++i)
    {
        LinkVectorPairIterator j = i + 1;
        for(; j != m_links.end(); ++j)
        {
            bool isAmbiguous = areLinksAmbiguous(i->link, j->link, ambiguity_p);
            if(isAmbiguous)
            {
                return true;
            }
        }
    }
    return false;
}

bool ScaffoldGroup::markPolymorphic(double p_cutoff, double cn_cutoff)
{
    LinkVectorPairIterator i = m_links.begin();
    for(; i != m_links.end(); ++i)
    {
        LinkVectorPairIterator j = i + 1;
        for(; j != m_links.end(); ++j)
        {
            bool isAmbiguous = areLinksAmbiguous(i->link, j->link, p_cutoff);
            if(isAmbiguous)
            {
                double sum = i->pEndpoint->getEstCopyNumber() + j->pEndpoint->getEstCopyNumber();
                std::cout << "\tLinks " << i->link << " and " << j->link << " are ambiguous\n";
                std::cout << "\tECN I: "<< i->pEndpoint->getEstCopyNumber() << " ECN J: " << j->pEndpoint->getEstCopyNumber() << "\n";
                std::cout << "\tsum: " << sum << "\n";
                if(sum < cn_cutoff)
                {
                    if(i->pEndpoint->getEstCopyNumber() < j->pEndpoint->getEstCopyNumber())
                        i->pEndpoint->setClassification(SVC_POLYMORPHIC);
                    else
                        j->pEndpoint->setClassification(SVC_POLYMORPHIC);
                    return true;
                }
            }
        }
    }
    
    return false;
}

bool ScaffoldGroup::areLinksAmbiguous(const ScaffoldLink& linkA,
                                      const ScaffoldLink& linkB,
                                      double p)
{
    double p_AB = calculateProbACloserThanB(linkA, linkB);
    double p_BA = 1.0f - p_AB;
    double best = std::max(p_AB, p_BA);
    double p_wrong = 1.0f - best;
    return p_wrong > p;
}

bool ScaffoldGroup::hasConsistentLayout()
{
    return true;
}

int ScaffoldGroup::calculateLongestOverlap()
{
    int longestOverlap = 0;
    LinkVectorPairIterator i = m_links.begin();
    for(; i != m_links.end(); ++i)
    {
        Interval interval_i(i->link.distance, i->link.getEndpoint() - 1);

        LinkVectorPairIterator j = i + 1;
        for(; j != m_links.end(); ++j)
        {
            Interval interval_j(j->link.distance, j->link.getEndpoint() - 1);
            if(Interval::isIntersecting(interval_i.start, interval_i.end,
                                        interval_j.start, interval_j.end)) {
                Interval intersection;
                Interval::intersect(interval_i.start, interval_i.end,
                                    interval_j.start, interval_j.end,
                                    intersection.start, intersection.end);

                int overlap = intersection.end - intersection.start + 1;
                if(overlap > longestOverlap)
                    longestOverlap = overlap;
            }
        }
    }

    return longestOverlap;
}

void ScaffoldGroup::computeBestOrdering()
{
    LinkList unplacedLinks(m_links.begin(), m_links.end());
    LinkPairVector placedLinks;

    int totalScore = 0;
    while(!unplacedLinks.empty())
    {
        int bestScore = std::numeric_limits<int>::max();
        LinkListIterator bestLink = unplacedLinks.end();

        for(LinkListIterator iter = unplacedLinks.begin();
                             iter != unplacedLinks.end();
                             ++iter)
        {
            int score = scoreLinkPlacement(iter->link, unplacedLinks);
            if(score < bestScore)
            {
                bestScore = score;
                bestLink = iter;
            }
        }

        assert(bestLink != unplacedLinks.end());

        unplacedLinks.erase(bestLink);
        totalScore += bestScore;
    }

    m_links = placedLinks;
    m_isOrdered = true;
}

std::string ScaffoldGroup::getBestOrderingString() const
{
    assert(m_isOrdered);
    std::stringstream ss;
    for(LinkVectorPairConstIterator iter = m_links.begin(); iter != m_links.end(); ++iter)
    {
        ss << iter->link.endpointID << ":" << iter->link.distance << "-" << iter->link.getEndpoint() << " ";
    }
    return ss.str();
}

void ScaffoldGroup::getLinearLinks(LinkVector& outLinks)
{
    assert(m_isOrdered);

    if(m_links.empty())
        return;


    LinkPairVector::const_iterator iter = m_links.begin();
    ScaffoldLink previousLink = iter->link;
    outLinks.push_back(previousLink);
    ++iter;
    for(; iter != m_links.end(); ++iter)
    {
        const ScaffoldLink& currentLink = iter->link;
        EdgeComp orientation = previousLink.getComp() == currentLink.getComp() ? EC_SAME : EC_REVERSE;

        EdgeDir dir;
        if(previousLink.getComp() == EC_SAME)
            dir = previousLink.getDir();
        else
            dir = !previousLink.getDir();

        int distance = currentLink.distance - previousLink.getEndpoint();
        double sd = sqrt(pow(currentLink.stdDev, 2.0) + pow(previousLink.stdDev, 2.0));
        ScaffoldLink inferred(currentLink.endpointID, dir, orientation, distance, sd, 0, currentLink.seqLen, SLT_INFERRED);
        outLinks.push_back(inferred);
        previousLink = currentLink;
    }
}

void ScaffoldGroup::getSecondaryLinks()
{
    assert(m_isOrdered);

    if(m_links.empty())
        return;

    LinkPairVector::const_iterator iter = m_links.begin();
    for(; iter != m_links.end(); ++iter)
    {
        ScaffoldLink xyLink = iter->link;
        ScaffoldEdgePtrVector yEdges = pY->getEdges();


        int xyDistance = xyLink.distance;
        EdgeDir yxDir = xyLink.getTwinDir();
        EdgeComp yxComp = xyLink.getComp();
        for(size_t i = 0; i < yEdges.size(); ++i)
        {
            EdgeDir yzDir = pYZ->getDir();
            EdgeComp yzComp = pYZ->getComp();
            int yzDistance = pYZ->getDistance();

            int xzDistance;
            if(yzDir == yxDir)
            {
                xzDistance = xyDistance - (pZ->getSeqLen() + yzDistance);
                printf("SAME xy %d - (%d + %d) = %d\n", xyDistance, (int)pZ->getSeqLen(), yzDistance, xzDistance);
            }
            else
            {
                xzDistance = xyDistance + pY->getSeqLen() + yzDistance;
                printf("DIFF xy %d + %d + %d = %d\n", xyDistance, (int)pY->getSeqLen(), yzDistance, xzDistance);
            }

            std::cout << "yzLink: " << pYZ->getLink() << "\n";

            WARN_ONCE("Recalculate xzDir");
            EdgeDir xzDir = xyLink.getDir();
            double xzSD = sqrt(pow(xyLink.stdDev, 2.0) + pow(pYZ->getStdDev(), 2.0));
            EdgeComp xzComp = (yxComp == yzComp) ? EC_SAME : EC_REVERSE;
            ScaffoldLink xzLink(pZ->getID(), xzDir, xzComp, xzDistance, xzSD, 0, pZ->getSeqLen(), SLT_INFERRED);
            std::cout << "Inferred secondary edge: " << xzLink << "\n";
        }
    }
}

int ScaffoldGroup::scoreLinkPlacement(const ScaffoldLink& link,
                                      const LinkList& unplacedList)
{
    int end = link.getEndpoint();

    int score = 0;
    for(LinkList::const_iterator iter = unplacedList.begin();
                                 iter != unplacedList.end();
                                 ++iter)
    {
        if(link.endpointID != iter->link.endpointID)
        {
            int next_valid_position = end - m_maxOverlap - 1;
            if(iter->link.distance < next_valid_position)
            {
                score += next_valid_position - iter->link.distance;
            }
        }
    }
    return score;
}

double ScaffoldGroup::calculateProbACloserThanB(const ScaffoldLink& linkA,
                                                const ScaffoldLink& linkB)
{
    double mean = linkA.distance - linkB.distance;
    double variance = pow(linkA.stdDev, 2.0f) + pow(linkB.stdDev, 2.0f);

    return normCDF(0.0f, mean, sqrt(variance));
}

double ScaffoldGroup::normCDF(double x, double mean, double sd)
{
}
