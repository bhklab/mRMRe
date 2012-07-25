#include "MutualInformationMatrix.hpp"

/* explicit */
MutualInformationMatrix::MutualInformationMatrix(Matrix* const pDataMatrix,
        unsigned int const* const pSampleStrata, float const* const pSampleWeights,
        unsigned int const* const pFeatureTypes) :
        SymmetricMatrix(pDataMatrix->getColumnCount()), mpDataMatrix(pDataMatrix), mpRankedDataMatrix(
                new Matrix(pDataMatrix->getRowCount(), pDataMatrix->getColumnCount())), mpHasFeatureRanksCached(
                new bool[pDataMatrix->getColumnCount()]), mpSampleStrata(pSampleStrata), mpSampleWeights(
                pSampleWeights), mpFeatureTypes(pFeatureTypes)
{
    for (unsigned int i = 0; i < mpDataMatrix->getColumnCount(); ++i)
        mpHasFeatureRanksCached[i] = false;

    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = i; j < mColumnCount; ++j)
            SymmetricMatrix::operator()(i, j) = std::numeric_limits<float>::quiet_NaN();
}

/* virtual */
MutualInformationMatrix::~MutualInformationMatrix()
{
    delete mpRankedDataMatrix;
    delete[] mpHasFeatureRanksCached;
}

/* virtual */float&
MutualInformationMatrix::operator()(unsigned int const i, unsigned int const j)
{
    if (SymmetricMatrix::operator()(i, j) != SymmetricMatrix::operator()(i, j))
    {
        float r = 0;
        // Check correlation type
        bool const A_is_continuous = mpFeatureTypes[i] == FEATURE_CONTINUOUS;
        bool const A_is_discrete = mpFeatureTypes[i] == FEATURE_DISCRETE;
        bool const A_is_survival_event = mpFeatureTypes[i] == FEATURE_SURVIVAL_EVENT;

        bool const B_is_continuous = mpFeatureTypes[j] == FEATURE_CONTINUOUS;
        bool const B_is_discrete = mpFeatureTypes[j] == FEATURE_DISCRETE;

        if (A_is_continuous && B_is_continuous)
        {
            if (!mpHasFeatureRanksCached[i])
            {
                placeRanksByFeatureIndex(i, mpRankedDataMatrix, mpDataMatrix);
                mpHasFeatureRanksCached[i] = true;
            }

            if (!mpHasFeatureRanksCached[j])
            {
                placeRanksByFeatureIndex(j, mpRankedDataMatrix, mpDataMatrix);
                mpHasFeatureRanksCached[j] = true;
            }

            r = computeSpearmanCorrelation(i, j, mpRankedDataMatrix);
        }

        else if (A_is_survival_event && B_is_continuous)
            r = computeConcordanceIndex(i, j, i + 1, mpDataMatrix, mpSampleWeights, mpSampleStrata,
                    true);
        else if (A_is_discrete && B_is_continuous)
            r = computeConcordanceIndex(i, j, -1, mpDataMatrix, mpSampleWeights, mpSampleStrata,
                    true);
        else if (A_is_continuous && B_is_discrete)
            r = computeConcordanceIndex(j, i, -1, mpDataMatrix, mpSampleWeights, mpSampleStrata,
                    true);
        else if (A_is_discrete && B_is_discrete)
            r = computeCramersV(i, j, mpDataMatrix, mpSampleWeights);
        else
            r = std::numeric_limits<double>::quiet_NaN();

        SymmetricMatrix::operator()(i, j) = r; //-0.5 * log(1 - (r * r));
    }

    return SymmetricMatrix::operator()(i, j);
}
