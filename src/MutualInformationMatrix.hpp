#ifndef ensemble_MutualInformationMatrix_hpp
#define ensemble_MutualInformationMatrix_hpp

#include <limits>

#include "SymmetricMatrix.hpp"
#include "tools.hpp"

class MutualInformationMatrix : public SymmetricMatrix
{
protected:
    Matrix* const mpDataMatrix;
    Matrix* const mpRankedDataMatrix;
    bool* const mpHasFeatureRanksCached;
    unsigned int const* const mpSampleStrata;
    float const* const mpSampleWeights;
    unsigned int const* const mpFeatureTypes;

public:
    static unsigned int const FEATURE_CONTINUOUS = 0;
    static unsigned int const FEATURE_DISCRETE = 1;
    static unsigned int const FEATURE_SURVIVAL_EVENT = 2;
    static unsigned int const FEATURE_SURVIVAL_TIME = 3;

    explicit
    MutualInformationMatrix(Matrix* const pMatrix, unsigned int const* const pSampleStrata,
            float const* const pSampleWeights, unsigned int const* const pFeatureTypes);

    virtual
    ~MutualInformationMatrix();

    virtual float&
    operator()(unsigned int const i, unsigned int const j);
};

#endif /* ensemble_MutualInformationMatrix_hpp */
