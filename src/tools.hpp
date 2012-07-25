#ifndef ensemble_tools_hpp
#define ensemble_tools_hpp

#include <algorithm>
#include <cmath>

#include "Matrix.hpp"

class DataMatrixComparator
{
private:
    unsigned int const mFeatureIndex;
    Matrix* const mpDataMatrix;

public:
    DataMatrixComparator(unsigned int const featureIndex, Matrix* const pDataMatrix);

    bool const
    operator()(unsigned int const i, unsigned int const j) const;
};

void const
placeRanksByFeatureIndex(unsigned int const index, Matrix* const pRankedDataMatrix,
        Matrix* const pDataMatrix);

float const
computeSpearmanCorrelation(unsigned int const i, unsigned int const j,
        Matrix* const pRankedDataMatrix);

float const
computeConcordanceIndex(unsigned int const discreteFeatureIndex,
        unsigned int const continuousFeatureIndex, int const timeFeatureIndex,
        Matrix* const pDataMatrix, float const* const pSampleWeights,
        unsigned int const* const pSampleStrata, bool const outX);

float const
computeCramersV(unsigned int const featureIndex1, unsigned int const featureIndex2,
        Matrix* pDataMatrix, float const* const pSampleWeights);

#endif /* ensemble_tools_hpp */
