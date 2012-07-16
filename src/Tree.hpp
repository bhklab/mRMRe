#ifndef ensemble_Tree_hpp
#define ensemble_Tree_hpp

class Tree
{
private:
    unsigned int* const mpChildrenCountPerLevel;
    unsigned int const mLevelCount;
    unsigned int* mpData;

public:
    Tree(unsigned int* const pChildrenCountPerLevel, unsigned int const levelCount);

    void
    build();

};

#endif /* ensemble_Tree_hpp */
