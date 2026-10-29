#ifndef mRMRe_Math_h
#define mRMRe_Math_h

#include <algorithm>
#include <cmath>
#include <limits>

#include "Matrix.h"

namespace Math
{
    class IndirectComparator
    {
    private:
        double const* const mpSamples;
        unsigned int const* const mpSampleIndices;

    public:
        IndirectComparator(double const* const pSamples, unsigned int const* const pSampleIndices);

        bool
        operator()(unsigned int const i, unsigned int const j) const;
    };

    void
    computeCausality(double* const pCausalityArray, Matrix const* const pMiMatrix,
            int const* const pSolutions, unsigned int const solutionCount,
            unsigned int const featureCountPerSolution, unsigned int const featureCount,
            unsigned int const targetFeatureIndex);

    double
    computeCoInformationLattice(double const cor_ij, double const cor_ik, double const cor_jk);

    double
    computeConcordanceIndex(double const* const pDiscreteSamples,
            double const* const pContinuousSamples, double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            bool const outX, double* const pConcordantWeights = 0, double* const pDiscordantWeights =
                    0, double* const pUninformativeWeights = 0, double* const pRelevantWeights = 0);

    double
    computeConcordanceIndex(double const* const pDiscreteSamples,
            double const* const pContinuousSamples, double const* const pTimeSamples,
            double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            bool const outX, double* const pConcordantWeights = 0, double* const pDiscordantWeights =
                    0, double* const pUninformativeWeights = 0, double* const pRelevantWeights = 0);

    double
    computeConcordanceIndex(double const* const pDiscreteSamplesX,
            double const* const pDiscreteSamplesY, double const* const pTimeSamplesX,
            double const* const pTimeSamplesY, double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            bool const outX, double* const pConcordantWeights = 0, double* const pDiscordantWeights =
                    0, double* const pUninformativeWeights = 0, double* const pRelevantWeights = 0);

    double
    computeCramersV(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const bootstrapCount, unsigned int const seed);

    double
    computeCramersV(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights, unsigned int const* const pSampleIndices,
            unsigned int const sampleCount, double* const pTotalWeight = 0);

    double
    computeFrequency(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const bootstrapCount, unsigned int const seed);

    double
    computeFrequency(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights, unsigned int const* const pSamplesIndices,
            unsigned int const sampleCount, double* const pTotalWeight = 0);

    double
    computeMi(double const r);

    double
    computePearsonCorrelation(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const bootstrapCount, unsigned int const seed);

    double
    computePearsonCorrelation(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights, unsigned int const* const pSampleIndices,
            unsigned int const sampleCount, double* const pTotalWeight = 0);

    double
    computeSomersD(double const c);

    double
    computeVariance(double const* const pSamples, unsigned int const sampleCount);

    void
    placeOrders(double const* const pSamples, double* const pOrders,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);

    void
    placeRanksFromOrders(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pOrdersX, double const* const pOrdersY, double* const pRanksX,
            double* const pRanksY, unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);

    void
    placeStratificationData(int const* const pSampleStrata, double const* const pSampleWeights,
            unsigned int** const pSampleIndicesPerStratum,
            unsigned int* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const sampleCount);
}

#endif /* mRMRe_Math_h */
