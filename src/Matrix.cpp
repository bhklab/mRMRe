#include "Matrix.h"

Matrix::Matrix(unsigned int const rowCount, unsigned int const columnCount) :
        mpData(new float[rowCount * columnCount]), mRowCount(rowCount), mColumnCount(columnCount), mHasAllocation(
                true)
{

}

/* explicit */
Matrix::Matrix(unsigned int const size, unsigned int const rowCount, unsigned int const columnCount) :
        mpData(new float[size]), mRowCount(rowCount), mColumnCount(columnCount), mHasAllocation(
                true)
{

}

/* explicit */
Matrix::Matrix(float* const pData, unsigned int const rowCount, unsigned int const columnCount) :
        mpData(pData), mRowCount(rowCount), mColumnCount(columnCount), mHasAllocation(false)
{

}

/* virtual */
Matrix::~Matrix()
{
    if (mHasAllocation)
        delete[] mpData;
}

/* virtual */float&
Matrix::at(unsigned int const i, unsigned int const j)
{
    return mpData[(j * mRowCount) + i];
}

/* virtual */float const&
Matrix::at(unsigned int const i, unsigned int const j) const
{
    return mpData[(j * mRowCount) + i];
}

unsigned int const
Matrix::getColumnCount() const
{
    return mColumnCount;
}

unsigned int const
Matrix::getRowCount() const
{
    return mRowCount;
}

std::vector<float> const
Matrix::getVectorizedData() const
{
    std::vector<float> vector;
    vector.resize(mRowCount * mColumnCount);

#pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = 0; j < mRowCount; ++j)
            vector[(i * mRowCount) + j] = at(i, j);

    return vector;
}

Matrix const
Matrix::getTransposedMatrix() const
{
    Matrix transposed(mColumnCount, mRowCount);
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = 0; j < mRowCount; ++j)
            transposed.at(j, i) = at(i, j);
    return transposed;
}

Matrix const
Matrix::operator+(Matrix const& matrix) const
{
    Matrix result(mRowCount, mColumnCount);
    for (unsigned int i = 0; i < getRowCount(); ++i)
        for (unsigned int j = 0; j < matrix.getColumnCount(); ++j)
            result.at(i, j) = at(i, k) + matrix.at(k, j);
    return result;
}

Matrix const
Matrix::operator-(Matrix const& matrix) const
{
    Matrix result(mRowCount, mColumnCount);
    for (unsigned int i = 0; i < getRowCount(); ++i)
        for (unsigned int j = 0; j < matrix.getColumnCount(); ++j)
            result.at(i, j) = at(i, k) - matrix.at(k, j);
    return result;
}

Matrix const
Matrix::operator*(Matrix const& matrix) const
{
    Matrix result(mRowCount, matrix.mColumnCount);
    for (unsigned int i = 0; i < getRowCount(); ++i)
        for (unsigned int j = 0; j < matrix.getColumnCount(); ++j)
        {
            result.at(i, j) = 0;
            for (unsigned int k = 0; k < matrix.getRowCount(); ++k)
                result.at(i, j) += at(i, k) * matrix.at(k, j);
        }
    return result;
}

Matrix const
Matrix::operator/(Matrix const& matrix) const
{
    Matrix result(mRowCount, mColumnCount);
    for (unsigned int i = 0; i < getRowCount(); ++i)
        for (unsigned int j = 0; j < matrix.getColumnCount(); ++j)
            result.at(i, j) = at(i, k) / matrix.at(k, j);
    return result;
}
