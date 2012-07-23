#include "MutualInformationMatrix.hpp"

/* explicit */
MutualInformationMatrix::MutualInformationMatrix(Matrix* const pDataMatrix) :
        SymmetricMatrix(pDataMatrix->getColumnCount()), mpDataMatrix(pDataMatrix), mpRankedDataMatrix(
                new Matrix(pDataMatrix->getRowCount(), pDataMatrix->getColumnCount())), mpHasFeatureRanksCached(
                new bool[pDataMatrix->getColumnCount()])
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
        if (mpHasFeatureRanksCached[i])
        {
            placeRanksByFeatureIndex(i, mpRankedDataMatrix, mpDataMatrix);
            mpHasFeatureRanksCached[i] = true;
        }

        if (mpHasFeatureRanksCached[j])
        {
            placeRanksByFeatureIndex(j, mpRankedDataMatrix, mpDataMatrix);
            mpHasFeatureRanksCached[j] = true;
        }

        SymmetricMatrix::operator()(i, j) = computeSpearmanCorrelation(i, j, mpRankedDataMatrix);
    }

    return SymmetricMatrix::operator()(i, j);
}
