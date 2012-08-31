#ifndef ensemble_MutualInformationMatrix_h
#define ensemble_MutualInformationMatrix_h

#include <omp.h>

#include "Data.h"
#include "SymmetricMatrix.h"

class MutualInformationMatrix : public SymmetricMatrix
{
private:
    MutualInformationMatrix(const MutualInformationMatrix&);

    MutualInformationMatrix&
    operator=(const MutualInformationMatrix&);

protected:
    Data const* const mpData;

public:
    MutualInformationMatrix(Data const* const pData);

    MutualInformationMatrix(Data const* const pData, float* const pInternalData);

    virtual
    ~MutualInformationMatrix();

    virtual float&
    at(unsigned int const i, unsigned int const j);

    virtual float const&
    at(unsigned int const i, unsigned int const j) const;

    void const
    build();
};

#endif /* ensemble_MutualInformationMatrix_h */
