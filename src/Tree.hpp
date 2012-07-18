#ifndef ensemble_Tree_hpp
#define ensemble_Tree_hpp

#include <cmath>
#include <limits>
#include <vector>

#include <omp.h>

#include "MatrixInterface.hpp"

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

    void const
    build();

    inline float const
    computeQualityScore(unsigned int const absoluteIndex, unsigned int const level) const;

    inline unsigned int const
    getParentAbsoluteIndex(unsigned int const absoluteIndex, unsigned int const level) const;

    void const
    getPaths(std::vector<unsigned int>* pPaths) const;

    bool const
    hasAncestorByFeatureIndex(unsigned int const absoluteIndex, unsigned int const featureIndex,
            unsigned int level) const;

    bool const
    hasSamePath(unsigned int const absoluteIndex1, unsigned int const absoluteIndex2,
            unsigned int const level) const;

    bool const
    hasSiblingByFeatureIndex(unsigned int const absoluteIndex, unsigned int const featureIndex,
            unsigned int const level) const;

    bool const
    isRedundantPath(unsigned int const absoluteIndex, unsigned int const featureIndex,
            unsigned int const level) const;

    void const
    placeElement(unsigned int const absoluteIndex, unsigned int const level);
};

#endif /* ensemble_Tree_hpp */
