#include "MutualInformationMatrix.h"

MutualInformationMatrix::MutualInformationMatrix(Data const* const pData) :
        Matrix(pData->getFeatureCount() * pData->getFeatureCount(), pData->getFeatureCount(),
                pData->getFeatureCount()), mpData(pData)
{
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = 0; j < mColumnCount; ++j)
            Matrix::at(i, j) = std::numeric_limits<double>::quiet_NaN();
}

MutualInformationMatrix::MutualInformationMatrix(Data const* const pData,
        double* const pInternalData) :
        Matrix(pInternalData, pData->getFeatureCount(), pData->getFeatureCount()), mpData(pData)
{

}

/* virtual */
MutualInformationMatrix::~MutualInformationMatrix()
{

}

/* virtual */double&
MutualInformationMatrix::at(unsigned int const i, unsigned int const j)
{
    if (Matrix::at(i, j) != Matrix::at(i, j))
        mpData->computeMiBetweenFeatures(i, j, &Matrix::at(i, j), &Matrix::at(j, i));

    return Matrix::at(i, j);
}

/* virtual */double const&
MutualInformationMatrix::at(unsigned int const i, unsigned int const j) const
{
    return Matrix::at(i, j);
}

void const
MutualInformationMatrix::build()
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = 0; j < mColumnCount; ++j)
            at(i, j);
}
