#ifndef ensemble_MutualInformationMatrix_hpp
#define ensemble_MutualInformationMatrix_hpp

#include <limits>

#include "SymmetricMatrix.hpp"
#include "tools.hpp"

class MutualInformationMatrix : public SymmetricMatrix
{
protected:
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

    explicit
    MutualInformationMatrix(Matrix const* const pMatrix, unsigned int const* const pSampleStrata,
            float const* const pSampleWeights, unsigned int const* const pFeatureTypes,
            unsigned int const sampleStratumCount);

    virtual
    ~MutualInformationMatrix();

    virtual float&
    operator()(unsigned int const i, unsigned int const j);
};

#endif /* ensemble_MutualInformationMatrix_hpp */
