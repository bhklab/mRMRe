#ifndef ensemble_Tree_hpp
#define ensemble_Tree_hpp

class Tree
{
private:
    unsigned int* const mpChildrenCountPerLevel;
    unsigned int const mLevelCount;
    float* const mpFeatureInformationMatrix;
    unsigned int const mFeatureCount;
    unsigned int* mpStartingIndexPerLevel;
    unsigned int* mpIndexTree;
    float* mpInformativeContributionTree;
    float* mpRedundantContributionTree;

public:
    Tree(unsigned int* const pChildrenCountPerLevel, unsigned int const levelCount,
            float* const pFeatureInformationMatrix, unsigned int const featureCount);

    ~Tree();

    void
    build();

    inline unsigned int const
    getParentAbsoluteIndex(unsigned int const absoluteIndex, unsigned int const level) const;
};

#endif /* ensemble_Tree_hpp */
