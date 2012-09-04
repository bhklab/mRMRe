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

Matrix::operator std::vector<float>() const
{
    std::vector<float> elements;
    elements.reserve(mRowCount * mColumnCount);

    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = 0; j < mRowCount; ++j)
            elements.push_back(at(i, j));

    return elements;
}
