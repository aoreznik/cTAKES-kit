

class ScaffoldVertex;

class ScaffoldEdge
{
    public:


        void setLink(ScaffoldLink link);
        void setColor(GraphColor c);

        VertexID getStartID() const;
        VertexID getEndID() const;
        const ScaffoldLink& getLink() const;

        EdgeDir getDir() const;
        EdgeComp getComp() const;        
        ScaffoldLinkType getType() const;
        int getDistance() const;
        double getStdDev() const;
        char getTypeCode() const;
        GraphColor getColor() const;

        std::string makeLinkString() const;

        friend std::ostream& operator<<(std::ostream& out, const ScaffoldEdge& edge);

    private:
        GraphColor m_color;
};



