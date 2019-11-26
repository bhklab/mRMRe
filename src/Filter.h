#ifndef mRMRe_Filter_h
#define mRMRe_Filter_h

#include <cmath>
#include <limits>


#ifdef _OPENMP
#include <omp.h>
#endif 

#include <vector>

#include "Math.h"
#include "Matrix.h"

class Filter
{
private:
    Filter(const Filter&);

    Filter&
    operator=(const Filter&);

    int const* const mpChildrenCountPerLevel;
    unsigned int const mLevelCount;
    Matrix* const mpFeatureInformationMatrix;
    unsigned int* const mpStartingIndexPerLevel;
    unsigned int const mFixedFeatureCount;
    unsigned int* mpIndexTree;
    double* mpScoreTree;
    unsigned int mTreeElementCount;

public:
    Filter(int const* const pChildrenCountPerLevel, unsigned int const levelCount,
            Matrix* const pFeatureInformationMatrix, unsigned int const targetFeatureIndex,
            unsigned int const fixedFeatureCount);
    
    ~Filter();

    void const
    build();

    inline unsigned int const
    getParentAbsoluteIndex(unsigned int const absoluteIndex, unsigned int const level) const;

    void const
    getSolutions(int* const solutions) const;

    void const
    getScores(double* const scores) const;

    bool const
    hasAncestorByFeatureIndex(unsigned int const absolute);

    bool const
    hasAncestorByFeatureIndex(unsigned int const absoluteIndex, unsigned int const featureIndex,
            unsigned int level) const;

    bool const
    isRedundantPath(unsigned int const absoluteIndex, unsigned int const featureIndex,
            unsigned int const level) const;

    void const
    placeElements(unsigned int const startingIndex, unsigned int childrenCount,
            unsigned int const level);
};

#endif /* mRMRe_Filter_h */