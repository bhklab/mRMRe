#include "SymmetricMatrix.h"

SymmetricMatrix::SymmetricMatrix(unsigned int const rowCount) :
        Matrix(rowCount * (rowCount + 1) / 2, rowCount, rowCount)
{

}

SymmetricMatrix::SymmetricMatrix(double* const pData, unsigned int const rowCount) :
        Matrix(pData, rowCount, rowCount)
{

}

/* virtual */
SymmetricMatrix::~SymmetricMatrix()
{

}

/* virtual */double&
SymmetricMatrix::at(unsigned int const i, unsigned int const j)
{
    unsigned int const x = (i >= j) ? i : j;
    unsigned int const y = (i >= j) ? j : i;

    return mpData[((x * (x + 1) / 2)) + y];
}

/* virtual */double const&
SymmetricMatrix::at(unsigned int const i, unsigned int const j) const
{
    unsigned int const x = (i >= j) ? i : j;
    unsigned int const y = (i >= j) ? j : i;

    return mpData[((x * (x + 1) / 2)) + y];
}
