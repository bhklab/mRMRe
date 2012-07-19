#include "MatrixInterface.hpp"

MatrixInterface::MatrixInterface(unsigned int const size, unsigned int const rowCount,
        unsigned int const columnCount) :
        mpData(new float[size]), mRowCount(rowCount), mColumnCount(rowCount), mHasAllocation(true)
{

}

MatrixInterface::MatrixInterface(float* const data, unsigned int const rowCount,
        unsigned int const columnCount) :
        mpData(data), mRowCount(rowCount), mColumnCount(rowCount), mHasAllocation(false)
{

}

MatrixInterface::~MatrixInterface()
{
    if (mHasAllocation)
        delete[] mpData;
}
