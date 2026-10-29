#include "Matrix.h"

Matrix::Matrix(unsigned int const rowCount, unsigned int const columnCount) :
        mpData(new double[static_cast<std::size_t>(rowCount) * columnCount]), mRowCount(rowCount), mColumnCount(columnCount), mHasAllocation(
                true)
{

}

/* explicit */
Matrix::Matrix(double* const pData, unsigned int const rowCount, unsigned int const columnCount) :
        mpData(pData), mRowCount(rowCount), mColumnCount(columnCount), mHasAllocation(false)
{

}

/* virtual */
Matrix::~Matrix()
{
    if (mHasAllocation)
        delete[] mpData;
}

/* virtual */double&
Matrix::at(unsigned int const i, unsigned int const j)
{
    return mpData[(static_cast<std::size_t>(j) * mRowCount) + i];
}

/* virtual */double const&
Matrix::at(unsigned int const i, unsigned int const j) const
{
    return mpData[(static_cast<std::size_t>(j) * mRowCount) + i];
}

unsigned int
Matrix::getColumnCount() const
{
    return mColumnCount;
}

unsigned int
Matrix::getRowCount() const
{
    return mRowCount;
}
