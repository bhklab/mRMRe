#ifndef ensemble_Tree_h
#define ensemble_Tree_h

#include <cmath>
#include <limits>
#include <omp.h>
#include <vector>

#include "Matrix.h"

class Tree
{
private:
    unsigned int const* const mpChildrenCountPerLevel;
    unsigned int const mLevelCount;
    Matrix const* const mpFeatureInformationMatrix;
    unsigned int* const mpStartingIndexPerLevel;
    unsigned int* mpIndexTree;
    float* mpInformativeContributionTree;
    float* mpRedundantContributionTree;
    unsigned int mTreeElementCount;
    std::vector<unsigned int> mPaths;
    std::vector<float> mScores;

public:
    Tree(unsigned int const* const pChildrenCountPerLevel, unsigned int const levelCount,
            Matrix const* const pFeatureInformationMatrix, unsigned int const targetFeatureIndex);

    ~Tree();

    void const
    build();

    inline float const
    computeQualityScore(unsigned int const absoluteIndex, unsigned int const level) const;

    inline unsigned int const
    getParentAbsoluteIndex(unsigned int const absoluteIndex, unsigned int const level) const;

    std::vector<unsigned int> const
    getPaths() const;

    std::vector<float> const
    getScores() const;

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

#endif /* ensemble_Tree_h */
