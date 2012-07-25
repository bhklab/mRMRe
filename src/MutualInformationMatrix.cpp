#include "MutualInformationMatrix.hpp"

/* explicit */
MutualInformationMatrix::MutualInformationMatrix(Data const* const pData) :
        SymmetricMatrix(pData->getFeatureCount()), mpData(pData)
{
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = i; j < mColumnCount; ++j)
            SymmetricMatrix::operator()(i, j) = std::numeric_limits<float>::quiet_NaN();
}

/* virtual */
MutualInformationMatrix::~MutualInformationMatrix()
{

}

/* virtual */float&
MutualInformationMatrix::operator()(unsigned int const i, unsigned int const j)
{
    if (SymmetricMatrix::operator()(i, j) != SymmetricMatrix::operator()(i, j))
        SymmetricMatrix::operator()(i, j) = mpData->computeMiBetweenFeatures(i, j);

    return SymmetricMatrix::operator()(i, j);
}
