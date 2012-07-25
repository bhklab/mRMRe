#ifndef ensemble_Data_hpp
#define ensemble_Data_hpp

#include <limits>

#include "Matrix.hpp"
#include "tools.hpp"

class Data
{
private:
    Matrix const* const mpDataMatrix;
    Matrix* const mpRankedDataMatrix;
    bool* const mpHasFeatureRanksCached;
    unsigned int const* const mpSampleStrata;
    float const* const mpSampleWeights;
    unsigned int const* const mpFeatureTypes;
    unsigned int const mSampleStratumCount;
    unsigned int** const mpSampleIndicesPerStratum;
    unsigned int* const mpSampleCountPerStratum;

public:
    static unsigned int const FEATURE_CONTINUOUS = 0;
    static unsigned int const FEATURE_DISCRETE = 1;
    static unsigned int const FEATURE_SURVIVAL_EVENT = 2;
    static unsigned int const FEATURE_SURVIVAL_TIME = 3;

    Data(float* const pData, unsigned int const sampleCount, unsigned int const featureCount,
            unsigned int const* const pSampleStrata, float const* const pSampleWeights,
            unsigned int const* const pFeatureTypes, unsigned int const sampleStratumCount);

    ~Data();

    float const
    computeMiBetweenFeatures(unsigned int const i, unsigned int const j) const;

    unsigned int const
    getSampleCount() const;

    unsigned int const
    getFeatureCount() const;
};

#endif /* ensemble_Data_hpp */
