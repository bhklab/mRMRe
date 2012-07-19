#include "SymmetricMatrix.hpp"

SymmetricMatrix::SymmetricMatrix(unsigned int const rowCount) :
        mpData(new float[rowCount * (rowCount + 1) / 2]), mRowCount(rowCount), mHasAllocation(true)
{
}

SymmetricMatrix::SymmetricMatrix(float* const data, unsigned int const rowCount) :
        mpData(data), mRowCount(rowCount), mHasAllocation(false)
{
}

SymmetricMatrix::~SymmetricMatrix()
{
    if (mHasAllocation)
        delete[] mpData;
}

/* inline */float&
SymmetricMatrix::operator()(unsigned int const i, unsigned int const j)
{
    unsigned int const x = (i >= j) ? i : j;
    unsigned int const y = (i >= j) ? j : i;

    return mpData[((x * (x + 1) / 2)) + y];
}

/* inline */unsigned int const
SymmetricMatrix::getRowCount() const
{
    return mRowCount;
}

/* inline */unsigned int const
SymmetricMatrix::getColumnCount() const
{
    return mRowCount;
}
