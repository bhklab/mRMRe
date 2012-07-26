#ifndef ensemble_Math_hpp
#define ensemble_Math_hpp

#include <algorithm>
#include <cmath>

//#include "Matrix.hpp"
//
//float const
//computeCramersV(unsigned int const featureIndex1, unsigned int const featureIndex2,
//        Matrix const* const pDataMatrix, float const* const pSampleWeights);

class Math
{
private:
    virtual
    ~Math() = 0; // Prevent instantiation of ::Math

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
    computeConcordanceIndex(float const* const pDiscreteSamples,
            float const* const pContinuousSamples, float const* const pTimeSamples,
            bool const timeSwitch, float const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            bool const outX);

    static float const
    computeFisherTransformation(float const r);

    static float const
    computeFisherTransformationReverse(float const z);

    static float const
    computeMi(float const r);

    static float const
    computePearsonCorrelation(float const* const pSamplesX, float const* const pSamplesY,
            float const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            float const* const pTotalWeightPerStratum,
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
