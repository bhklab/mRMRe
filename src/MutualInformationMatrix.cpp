#include "MutualInformationMatrix.hpp"

MutualInformationMatrix::MutualInformationMatrix(Matrix* const pMatrix) :
        SymmetricMatrix(pMatrix->getColumnCount())
{

}

/* virtual */
MutualInformationMatrix::~MutualInformationMatrix()
{

}

/* virtual */float&
MutualInformationMatrix::operator()(unsigned int const i, unsigned int const j)
{
    return SymmetricMatrix::operator()(i, j);
}
