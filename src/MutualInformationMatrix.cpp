#include "MutualInformationMatrix.h"

MutualInformationMatrix::MutualInformationMatrix(Data const* const pData) :
        Matrix(pData->getFeatureCount(), pData->getFeatureCount()), mpData(pData), mIsComputed(static_cast<std::size_t>(mColumnCount) * mColumnCount, 0)
{
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = 0; j < mColumnCount; ++j)
            Matrix::at(i, j) = std::numeric_limits<double>::quiet_NaN();
}

MutualInformationMatrix::MutualInformationMatrix(Data const* const pData,
        double* const pInternalData) :
        Matrix(pInternalData, pData->getFeatureCount(), pData->getFeatureCount()), mpData(pData),
        mIsComputed(static_cast<std::size_t>(mColumnCount) * mColumnCount, 0)
{
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = 0; j < mColumnCount; ++j)
            mIsComputed[index(i, j)] = Matrix::at(i, j) == Matrix::at(i, j);

}

/* virtual */
MutualInformationMatrix::~MutualInformationMatrix()
{

}

/* virtual */double&
MutualInformationMatrix::at(unsigned int const i, unsigned int const j)
{
    // The estimators can be asymmetric (the reverse cell may stay NaN), so (i, j) and (j, i) are
    // each attempted once. A cell that evaluates to NaN is flagged too, so it is not recomputed.
    // Callers in parallel regions must own the pair {i, j}: this writes both cells.
    if (!mIsComputed[index(i, j)])
    {
        mpData->computeMiBetweenFeatures(i, j, &Matrix::at(i, j), &Matrix::at(j, i));
        mIsComputed[index(i, j)] = 1;
        if (Matrix::at(j, i) == Matrix::at(j, i))
            mIsComputed[index(j, i)] = 1;
    }

    return Matrix::at(i, j);
}

/* virtual */double const&
MutualInformationMatrix::at(unsigned int const i, unsigned int const j) const
{
    return Matrix::at(i, j);
}

void
MutualInformationMatrix::build()
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = i; j < mColumnCount; ++j)
        {
            at(i, j);
            at(j, i);
        }
}
