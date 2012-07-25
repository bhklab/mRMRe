#ifndef ensemble_tools_hpp
#define ensemble_tools_hpp

#include <algorithm>
#include <cmath>

#include "Matrix.hpp"

class DataMatrixComparator
{
private:
    unsigned int const mFeatureIndex;
    Matrix const* const mpDataMatrix;

public:
    DataMatrixComparator(unsigned int const featureIndex, Matrix const* const pDataMatrix);

    bool const
    operator()(unsigned int const i, unsigned int const j) const;
};

float const
computeConcordanceIndex(unsigned int const discreteFeatureIndex,
        unsigned int const continuousFeatureIndex, int const timeFeatureIndex,
        Matrix const* const pDataMatrix, float const* const pSampleWeights,
        unsigned int const* const pSampleStrata, bool const outX);

bool const
isComparablePair(unsigned int const i, unsigned int const j, int const timeFeatureIndex,
        unsigned int const discreteFeatureIndex, Matrix const* const pDataMatrix);

float const
computeCramersV(unsigned int const featureIndex1, unsigned int const featureIndex2,
        Matrix const* const pDataMatrix, float const* const pSampleWeights);

float const
computePearsonCorrelation(unsigned int const i, unsigned int const j,
        Matrix const* const pDataMatrix, float const* const pSampleWeights);

/*float const
computeSpearmanCorrelation(unsigned int const i, unsigned int const j,
        Matrix const* const pRankedDataMatrix, float const* const pSampleWeights);*/

float const
convertCorrelationToMi(float const correlation);

void const
placeRanksByFeatureIndex(unsigned int const index, Matrix* const pRankedDataMatrix,
        Matrix const* const pDataMatrix);

#endif /* ensemble_tools_hpp */
