#ifndef ensemble_Tree_hpp
#define ensemble_Tree_hpp

#include <limits>
#include <vector>

#include <omp.h>

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
    unsigned int mTreeElementCount;

public:
    Tree(unsigned int* const pChildrenCountPerLevel, unsigned int const levelCount,
            float* const pFeatureInformationMatrix, unsigned int const featureCount,
            unsigned int const targetFeatureIndex);

    ~Tree();

    void
    build();

    unsigned int const
    selectBestFeature(unsigned int absoluteIndex, unsigned int level);

    bool
    isRedundantSolution(unsigned int absoluteIndex, unsigned int level);

    bool
    hasSameAncestry(unsigned int targetAbsoluteIndex, unsigned int consideredAbsoluteIndex,
            unsigned int level);

    float
    computeQualityScore(unsigned int absoluteIndex, unsigned int level);

    bool
    hasAncestorByIndex(unsigned int absoluteIndex, unsigned int consideredFeatureIndex, unsigned int level);

    bool
    hasSiblingByIndex(unsigned int absoluteIndex, unsigned int consideredFeatureIndex, unsigned int level);

    void
    getPaths(std::vector<unsigned int>* pPaths) const;

    inline unsigned int const
    getParentAbsoluteIndex(unsigned int const absoluteIndex, unsigned int const level) const;
};

#endif /* ensemble_Tree_hpp */
