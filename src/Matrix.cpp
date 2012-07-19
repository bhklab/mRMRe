#include "Matrix.hpp"

Matrix::Matrix(unsigned int const rowCount, unsigned int const columnCount) :
        MatrixInterface(rowCount * columnCount, rowCount, mColumnCount)
{

}

Matrix::Matrix(float* const data, unsigned int const rowCount, unsigned int const columnCount) :
        MatrixInterface(data, rowCount, mColumnCount)
{

}

Matrix::~Matrix()
{
}

/* inline */float&
Matrix::operator()(unsigned int const i, unsigned int const j)
{
    return mpData[(j * mRowCount) + i];
}

/* inline */unsigned int const
Matrix::getRowCount() const
{
    return mRowCount;
}

/* inline */unsigned int const
Matrix::getColumnCount() const
{
    return mColumnCount;
}
