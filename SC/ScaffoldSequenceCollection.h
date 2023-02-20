

class ScaffoldSequenceCollection
{
    public:
        ScaffoldSequenceCollection() {}
        
        virtual ~ScaffoldSequenceCollection() {}

        virtual std::string getSequence(const std::string& id) const = 0;

        virtual void setPlaced(const std::string& id) = 0;
        
};

class GraphSequenceCollection : public ScaffoldSequenceCollection
{
    public:

        ~GraphSequenceCollection() {}
        
        std::string getSequence(const std::string& id) const;

        void setPlaced(const std::string& id);


    private:
        
};

class MapSequenceCollection : public ScaffoldSequenceCollection
{
    public:

        MapSequenceCollection(std::string filename);
        ~MapSequenceCollection() {}
        
        std::string getSequence(const std::string& id) const;

        void setPlaced(const std::string& id);


    private:
        
        struct SequenceMapData
        {
            std::string sequence;
            bool isPlaced;
        };
        typedef std::map<std::string, SequenceMapData> SMPMap;
        SMPMap m_map;
};

