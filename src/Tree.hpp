#ifndef ensemble_Tree_hpp
#define ensemble_Tree_hpp

#include <iostream>

class Tree
{
private:
    unsigned int* const mpChildrenCountPerLevel;
    unsigned int const mLevelCount;
    unsigned int* mpStartingIndexPerLevel;
    unsigned int* mpIndexTree;
    float* mpInformativeContributionTree;
    float* mpRedundantContributionTree;

public:
    Tree(unsigned int* const pChildrenCountPerLevel, unsigned int const levelCount);

    ~Tree();

    void
    build();

};

#endif /* ensemble_Tree_hpp */
