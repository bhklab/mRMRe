#ifndef ensemble_Math_h
#define ensemble_Math_h

#include <algorithm>
#include <cmath>
#include <ctime>
#include <limits>

#include "Matrix.h"

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
            float const* const pContinuousSamples, float const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            bool const outX, float* const pConcordantWeight = 0, float* const pDiscordantWeight = 0,
            float* const pUninformativeWeight = 0, float* const pRelevantWeight = 0);

    static float const
    computeConcordanceIndexWithTime(float const* const pDiscreteSamples,
            float const* const pContinuousSamples, float const* const pTimeSamples,
            float const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            bool const outX, float* const pConcordantWeight = 0, float* const pDiscordantWeight = 0,
            float* const pUninformativeWeight = 0, float* const pRelevantWeight = 0);

    static float const
    computeCramersV(float const* const pSamplesX, float const* const pSamplesY,
            float const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            float const* const pTotalWeightPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const bootstrapCount);

    static float const
    computeCramersV(float const* const pSamplesX, float const* const pSamplesY,
            float const* const pSampleWeights, unsigned int const* const pSampleIndices,
            unsigned int const sampleCount);

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
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const bootstrapCount);

    static float const
    computePearsonCorrelation(float const* const pSamplesX, float const* const pSamplesY,
            float const* const pSampleWeights, unsigned int const* const pSampleIndices,
            unsigned int const sampleCount);

    static int const
    computeRandomNumber(unsigned int* const seed);

    static float const
    computeSomersD(float const c);

    static float const
    computeVariance(float const* const pSamples, unsigned int const sampleCount);

    static void const
    placeOrders(float const* const pSamples, float* const pOrders,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);

    static void const
    placeRanksFromOrders(float const* const pSamplesX, float const* const pSamplesY,
            float const* const pOrdersX, float const* const pOrdersY, float* const pRanksX,
            float* const pRanksY, unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);

    static void const
    placeRanksFromSamples(float const* const pSamples, float* const pRanks,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);

    static void const
    placeStratificationData(unsigned int const* const pSampleStrata,
            float const* const pSampleWeights, unsigned int** const pSampleIndicesPerStratum,
            float* const pTotalWeightPerStratum, unsigned int* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount, unsigned int const sampleCount);
};

#endif /* ensemble_Math_h */
