#ifndef ensemble_Data_h
#define ensemble_Data_h

#include <limits>

#include "Math.h"
#include "Matrix.h"

class Data
{
private:
    Matrix const* const mpDataMatrix;
    Matrix* const mpOrderMatrix;
    bool* const mpHasOrderCached;
    unsigned int const* const mpSampleStrata;
    float const* const mpSampleWeights;
    unsigned int const* const mpFeatureTypes;
    unsigned int const mSampleStratumCount;
    unsigned int** const mpSampleIndicesPerStratum;
    float* const mpTotalWeightPerStratum;
    unsigned int* const mpSampleCountPerStratum;
    bool const mUsesRanks;
    bool const mOutX;
    unsigned int const mBootstrapCount;

public:
    static unsigned int const FEATURE_CONTINUOUS = 0;
    static unsigned int const FEATURE_DISCRETE = 1;
    static unsigned int const FEATURE_SURVIVAL_EVENT = 2;
    static unsigned int const FEATURE_SURVIVAL_TIME = 3;

    Data(float* const pData, unsigned int const sampleCount, unsigned int const featureCount,
            unsigned int const* const pSampleStrata, float const* const pSampleWeights,
            unsigned int const* const pFeatureTypes, unsigned int const sampleStratumCount,
            bool const usesRanks, bool const outX, unsigned int const mBootstrapCount);

    ~Data();

    float const
    computeMiBetweenFeatures(unsigned int const i, unsigned int const j) const;

    float const
    computeCorrelationBetweenContinuousFeatures(unsigned int const i, unsigned int const j) const;

    unsigned int const
    getSampleCount() const;

    unsigned int const
    getFeatureCount() const;
};

#endif /* ensemble_Data_h */
