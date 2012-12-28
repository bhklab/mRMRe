#ifndef mRMRe_Data_h
#define mRMRe_Data_h

#include <limits>

#include "Math.h"
#include "Matrix.h"

class Data
{
private:
    Data(const Data&);

    Data&
    operator=(const Data&);

    Matrix const* const mpDataMatrix;
    Matrix* const mpOrderMatrix;
    Matrix const* const mpPriorsMatrix;
    bool* const mpHasOrderCached;
    int const* const mpSampleStrata;
    double const* const mpSampleWeights;
    int const* const mpFeatureTypes;
    unsigned int const mSampleStratumCount;
    unsigned int** const mpSampleIndicesPerStratum;
    unsigned int** const mpMasterSampleIndicesPerStratum;
    unsigned int* const mpSampleCountPerStratum;
    unsigned int const mContinuousEstimator;
    bool const mOutX;
    unsigned int const mBootstrapCount;
    double const mPriorsWeight;

public:
    static int const FEATURE_CONTINUOUS = 0;
    static int const FEATURE_DISCRETE = 1;
    static int const FEATURE_SURVIVAL_EVENT = 2;
    static int const FEATURE_SURVIVAL_TIME = 3;

    static int const PEARSON_ESTIMATOR = 0;
    static int const SPEARMAN_ESTIMATOR = 1;
    static int const KENDALL_ESTIMATOR = 2;
    static int const FREQUENCY_ESTIMATOR = 3;

    Data(double* const pData, Matrix const* const pPriorsMatrix, double const priorsWeight,
            unsigned int const sampleCount, unsigned int const featureCount,
            int const* const pSampleStrata, double const* const pSampleWeights,
            int const* const pFeatureTypes, unsigned int const sampleStratumCount,
            unsigned int const continuousEstimator, bool const outX,
            unsigned int const bootstrapCount);

    ~Data();

    void const
    bootstrap();

    void const
    computeMiBetweenFeatures(unsigned int const i, unsigned int const j, double* const mi_ij,
            double* const mi_ji) const;

    unsigned int const
    getSampleCount() const;

    unsigned int const
    getFeatureCount() const;
};

#endif /* mRMRe_Data_h */
