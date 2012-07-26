#include "MutualInformationMatrix.hpp"

MutualInformationMatrix::MutualInformationMatrix(Data const* const pData) :
        SymmetricMatrix(pData->getFeatureCount()), mpData(pData)
{
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = i; j < mColumnCount; ++j)
            SymmetricMatrix::at(i, j) = std::numeric_limits<float>::quiet_NaN();
}

/* virtual */
MutualInformationMatrix::~MutualInformationMatrix()
{

}

/* virtual */float&
MutualInformationMatrix::at(unsigned int const i, unsigned int const j)
{
    if (SymmetricMatrix::at(i, j) != SymmetricMatrix::at(i, j))
        SymmetricMatrix::at(i, j) = mpData->computeMiBetweenFeatures(i, j);

    return SymmetricMatrix::at(i, j);
}

void const
MutualInformationMatrix::build()
{
#pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < mColumnCount; ++i)
        for (unsigned int j = i; j < mColumnCount; ++j)
            at(i, j);
}
