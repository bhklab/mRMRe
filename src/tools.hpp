#ifndef ensemble_tools_hpp
#define ensemble_tools_hpp

#include <algorithm>
#include <cmath>

#include "Matrix.hpp"

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

#endif /* ensemble_tools_hpp */
