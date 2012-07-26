#ifndef ensemble_Math_hpp
#define ensemble_Math_hpp

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

// ***********

class Math
{
private:
    virtual ~Math() = 0; // Prevent instantiation of ::Math

public:
    class IndirectComparator
    {
    private:
        float const* const mpSamples;
        unsigned int const* const mpSampleIndices;

    public:
        IndirectComparator(float const* const pSamples, unsigned int const* const pSampleIndices);

        bool const
        operator()(unsigned int const i, unsigned int const j) const;
    };

    static float const
    computeMi(float const r);

    static float const
    computePearsonCorrelation(float const* const pSamplesX, float const* const pSamplesY,
            float const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);

    static float const
    computePearsonCorrelation(float const* const pSamplesX, float const* const pSamplesY,
            float const* const pSampleWeights, unsigned int const* const pSampleIndices,
            unsigned int const sampleCount);

    static void const
    placeRanksByFeatureIndex(float const* const pSamples, float* const pRanks,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);
};

#endif /* ensemble_Math_hpp */
