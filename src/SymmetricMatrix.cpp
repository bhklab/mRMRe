#include "SymmetricMatrix.hpp"

SymmetricMatrix::SymmetricMatrix(unsigned int const rowCount) :
        MatrixInterface(rowCount * (rowCount + 1) / 2, rowCount, rowCount)
{

}

SymmetricMatrix::SymmetricMatrix(float* const data, unsigned int const rowCount) :
        MatrixInterface(data, rowCount, rowCount)
{

}

SymmetricMatrix::~SymmetricMatrix()
{
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
    return mColumnCount;
}
