#include "MutualInformationMatrix.hpp"

/* explicit */
MutualInformationMatrix::MutualInformationMatrix(Matrix* const pDataMatrix) :
        SymmetricMatrix(pDataMatrix->getColumnCount()), mpDataMatrix(pDataMatrix), mpRankedDataMatrix(
                new Matrix(pDataMatrix->getRowCount(), pDataMatrix->getColumnCount()))
{

    mpTestArray = new unsigned int[pDataMatrix->getColumnCount()];
    for (unsigned int i = 0; i < mpDataMatrix->getColumnCount(); ++i)
        (*mpRankedDataMatrix)(0, i) = std::numeric_limits<float>::quiet_NaN();

#pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = i; j < mColumnCount; ++j)
            SymmetricMatrix::operator()(i, j) = std::numeric_limits<float>::quiet_NaN();
}

/* virtual */
MutualInformationMatrix::~MutualInformationMatrix()
{
    delete mpRankedDataMatrix;
}

/* virtual */float&
MutualInformationMatrix::operator()(unsigned int const i, unsigned int const j)
{
    if (SymmetricMatrix::operator()(i, j) != SymmetricMatrix::operator()(i, j))
    {

        if ((*mpRankedDataMatrix)(0, i) != (*mpRankedDataMatrix)(0, i))
            placeRanksByFeatureIndex(i, mpRankedDataMatrix, mpDataMatrix);

        if ((*mpRankedDataMatrix)(0, j) != (*mpRankedDataMatrix)(0, j))
            placeRanksByFeatureIndex(j, mpRankedDataMatrix, mpDataMatrix);

        SymmetricMatrix::operator()(i, j) = computeSpearmanCorrelation(i, j, mpRankedDataMatrix);
    }

    return SymmetricMatrix::operator()(i, j);
}

void const
MutualInformationMatrix::build()
{
#pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = i; j < mColumnCount; ++j)
            operator()(i, j);
}
