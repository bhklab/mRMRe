#ifndef ensemble_Tree_h
#define ensemble_Tree_h

#include <cmath>
#include <limits>
#include <omp.h>
#include <vector>

#include "Math.h"
#include "Matrix.h"

class Tree
{
private:
    unsigned int const* const mpChildrenCountPerLevel;
    unsigned int const mLevelCount;
    Matrix* const mpFeatureInformationMatrix;
    unsigned int* const mpStartingIndexPerLevel;
    unsigned int* mpIndexTree;
    unsigned int mTreeElementCount;
    std::vector<unsigned int> mPaths;

    Tree(const Tree&);

    Tree&
    operator=(const Tree&);

public:
    Tree(unsigned int const* const pChildrenCountPerLevel, unsigned int const levelCount,
            Matrix* const pFeatureInformationMatrix, unsigned int const targetFeatureIndex);

    ~Tree();

    void const
    build();

    inline unsigned int const
    getParentAbsoluteIndex(unsigned int const absoluteIndex, unsigned int const level) const;

    bool const
    hasAncestorByFeatureIndex(unsigned int const absoluteIndex, unsigned int const featureIndex,
            unsigned int level) const;

    bool const
    hasSamePath(unsigned int const absoluteIndex1, unsigned int const absoluteIndex2,
            unsigned int const level) const;

    bool const
    isRedundantPath(unsigned int const absoluteIndex, unsigned int const featureIndex,
            unsigned int const level) const;

    operator std::vector<unsigned int>() const;

    void const
    placeElements(unsigned int const startingIndex, unsigned int childrenCount,
            unsigned int const level);
};

#endif /* ensemble_Tree_h */
