#ifndef mRMRe_Math_h
#define mRMRe_Math_h

#include <algorithm>
#include <cmath>
#include <ctime>
#include <limits>

#include "Matrix.h"

class Math
{
private:
    Math();

    Math(const Math&);

    Math&
    operator=(const Math&);

    virtual
    ~Math() = 0;

public:
    class IndirectComparator
    {
    private:
        double const* const mpSamples;
        unsigned int const* const mpSampleIndices;

    public:
        IndirectComparator(double const* const pSamples, unsigned int const* const pSampleIndices);

        bool const
        operator()(unsigned int const i, unsigned int const j) const;
    };

    static void const
    computeCausality(double* const pCausalityArray, Matrix const* const pMiMatrix,
            int const* const pSolutions, unsigned int const solutionCount,
            unsigned int const featureCountPerSolution, unsigned int const featureCount,
            unsigned int const targetFeatureIndex);

    static double const
    computeCoInformationLattice(double const cor_ij, double const cor_ik, double const cor_jk);

    static double const
    computeConcordanceIndex(double const* const pDiscreteSamples,
            double const* const pContinuousSamples, double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            bool const outX, double* const pConcordantWeights = 0, double* const pDiscordantWeights =
                    0, double* const pUninformativeWeights = 0, double* const pRelevantWeights = 0);

    static double const
    computeConcordanceIndex(double const* const pDiscreteSamples,
            double const* const pContinuousSamples, double const* const pTimeSamples,
            double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            bool const outX, double* const pConcordantWeights = 0, double* const pDiscordantWeights =
                    0, double* const pUninformativeWeights = 0, double* const pRelevantWeights = 0);

    static double const
    computeConcordanceIndex(double const* const pDiscreteSamplesX,
            double const* const pDiscreteSamplesY, double const* const pTimeSamplesX,
            double const* const pTimeSamplesY, double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            bool const outX, double* const pConcordantWeights = 0, double* const pDiscordantWeights =
                    0, double* const pUninformativeWeights = 0, double* const pRelevantWeights = 0);

    static double const
    computeCramersV(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const bootstrapCount);

    static double const
    computeCramersV(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights, unsigned int const* const pSampleIndices,
            unsigned int const sampleCount, double* const pTotalWeight = 0);

    static double const
    computeFisherTransformation(double const r);

    static double const
    computeFisherTransformationReverse(double const z);

    static double const
    computeFrequency(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const bootstrapCount);

    static double const
    computeFrequency(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights, unsigned int const* const pSamplesIndices,
            unsigned int const sampleCount, double* const pTotalWeight = 0);

    static double const
    computeMi(double const r);

    static double const
    computePearsonCorrelation(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const bootstrapCount);

    static double const
    computePearsonCorrelation(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights, unsigned int const* const pSampleIndices,
            unsigned int const sampleCount, double* const pTotalWeight = 0);

    static int const
    computeRandomNumber(unsigned int* const seed);

    static double const
    computeSomersD(double const c);

    static double const
    computeSpearmanCorrelation(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pSampleWeights,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const bootstrapCount, unsigned int const sampleCount);

    static double const
    computeVariance(double const* const pSamples, unsigned int const sampleCount);

    static void const
    placeOrders(double const* const pSamples, double* const pOrders,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);

    static void const
    placeRanksFromOrders(double const* const pSamplesX, double const* const pSamplesY,
            double const* const pOrdersX, double const* const pOrdersY, double* const pRanksX,
            double* const pRanksY, unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);

    static void const
    placeRanksFromSamples(double const* const pSamples, double* const pRanks,
            unsigned int const* const * const pSampleIndicesPerStratum,
            unsigned int const* const pSampleCountPerStratum,
            unsigned int const sampleStratumCount);

    static void const
    placeStratificationData(int const* const pSampleStrata, double const* const pSampleWeights,
            unsigned int** const pSampleIndicesPerStratum,
            unsigned int* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
            unsigned int const sampleCount);
};

#endif /* mRMRe_Math_h */
