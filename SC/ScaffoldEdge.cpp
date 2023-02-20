
{
}

{
    m_pTwin = pTwin;
}

void ScaffoldEdge::setLink(ScaffoldLink link)
{
    m_link = link;
}

void ScaffoldEdge::setColor(GraphColor c)
{
    m_color = c;
}

VertexID ScaffoldEdge::getStartID() const
{
    assert(m_pTwin != NULL);
    return m_pTwin->getEndID();
}

VertexID ScaffoldEdge::getEndID() const
{
    return m_pEnd->getID();
}

{
    assert(m_pTwin != NULL);
    return m_pTwin->getEnd();
}

{
    return m_pEnd;
}

{
    return m_pTwin;
}

EdgeDir ScaffoldEdge::getDir() const
{
    return m_link.edgeData.getDir();
}

EdgeComp ScaffoldEdge::getComp() const
{
    return m_link.edgeData.getComp();
}

int ScaffoldEdge::getDistance() const
{
    return m_link.distance;
}

double ScaffoldEdge::getStdDev() const
{
    return m_link.stdDev;
}

GraphColor ScaffoldEdge::getColor() const
{
    return m_color;
}   

ScaffoldLinkType ScaffoldEdge::getType() const
{
    return m_link.type;
}

const ScaffoldLink& ScaffoldEdge::getLink() const
{
    return m_link;
}

std::string ScaffoldEdge::makeLinkString() const
{
    std::stringstream ss;
    ss << m_link;
    return ss.str();
}

std::ostream& operator<<(std::ostream& out, const ScaffoldEdge& edge)
{
    out << edge.getStartID() << " -- " << edge.getEndID() << "," << edge.getDistance() << "," << edge.getStdDev() 
        << "," << edge.getDir() << "," << edge.getComp() << "," << edge.m_link.getTypeCode();
    return out;
}

{
    return pXY->getDistance() < pXZ->getDistance();
}

