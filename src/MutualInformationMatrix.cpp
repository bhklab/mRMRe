#include "MutualInformationMatrix.hpp"

/* explicit */
MutualInformationMatrix::MutualInformationMatrix(Matrix* const pDataMatrix) :
        SymmetricMatrix(pDataMatrix->getColumnCount()), mpDataMatrix(pDataMatrix), mpRankedDataMatrix(
                new Matrix(pDataMatrix->getRowCount(), pDataMatrix->getColumnCount()))
{
    for (unsigned int i = 0; i < mpDataMatrix->getColumnCount(); ++i)
        (*mpRankedDataMatrix)(0, i) = std::numeric_limits<float>::quiet_NaN();

    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = i; j < mColumnCount; ++j)
            SymmetricMatrix::operator()(i, j) = std::numeric_limits<float>::quiet_NaN();
}

/* virtual */
MutualInformationMatrix::~MutualInformationMatrix()
{
    delete mpRankedDataMatrix;
}

#include <Rcpp.h>

/* virtual */float&
MutualInformationMatrix::operator()(unsigned int const i, unsigned int const j)
{
    if (SymmetricMatrix::operator()(i, j) != SymmetricMatrix::operator()(i, j))
    {
        if ((*mpRankedDataMatrix)(0, i) != (*mpRankedDataMatrix)(0, i))
            placeRanksByFeatureIndex(i);

        if ((*mpRankedDataMatrix)(0, j) != (*mpRankedDataMatrix)(0, j))
            placeRanksByFeatureIndex(j);

        placeSpearmanCorrelation(i, j);
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

void const
MutualInformationMatrix::placeRanksByFeatureIndex(unsigned int const dataMatrixFeatureIndex)
{
    unsigned int const sample_count = mpRankedDataMatrix->getRowCount();
    unsigned int p_order[sample_count];

    for (unsigned int i = 0; i < sample_count; ++i)
        p_order[i] = i;

    std::sort(p_order, p_order + sample_count,
            DataMatrixComparator(dataMatrixFeatureIndex, mpDataMatrix));

    for (unsigned int i = 0; i < sample_count; ++i)
        (*mpRankedDataMatrix)(p_order[i], dataMatrixFeatureIndex) = i;
}

void const
MutualInformationMatrix::placeSpearmanCorrelation(unsigned int const i, unsigned int const j)
{
    float *a = &(*mpRankedDataMatrix)(0, i);
    float *b = &(*mpRankedDataMatrix)(0, j);
    unsigned int const sample_count = mpRankedDataMatrix->getRowCount();
    float sum = 0.;

    for (unsigned int n = 0; n < sample_count; ++n)
    {
        float const difference = a[n] - b[n];
        sum += difference * difference;
    }

    float const correlation = 1
            - ((6 * sum) / (sample_count * ((sample_count * sample_count) - 1)));

    SymmetricMatrix::operator()(i, j) = -0.5 * log(1 - (correlation * correlation));
}

MutualInformationMatrix::DataMatrixComparator::DataMatrixComparator(
        unsigned int const dataMatrixFeatureIndex, Matrix* const pDataMatrix) :
        mDataMatrixFeatureIndex(dataMatrixFeatureIndex), mpDataMatrix(pDataMatrix)
{

}

bool const
MutualInformationMatrix::DataMatrixComparator::operator()(unsigned int const i,
        unsigned int const j) const
{
    return (*mpDataMatrix)(i, mDataMatrixFeatureIndex) < (*mpDataMatrix)(j, mDataMatrixFeatureIndex);
}
