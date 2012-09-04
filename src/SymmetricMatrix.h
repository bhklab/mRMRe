#ifndef mRMRe_SymmetricMatrix_h
#define mRMRe_SymmetricMatrix_h

#include "Matrix.h"

class SymmetricMatrix : public Matrix
{
private:
    SymmetricMatrix(const SymmetricMatrix&);

    SymmetricMatrix&
    operator=(const SymmetricMatrix&);

public:
    SymmetricMatrix(unsigned int const rowCount);

    SymmetricMatrix(float* const pData, unsigned int const rowCount);

    virtual
    ~SymmetricMatrix();

    virtual float&
    at(unsigned int const i, unsigned int const j);

    virtual float const&
    at(unsigned int const i, unsigned int const j) const;
};

#endif /* mRMRe_SymmetricMatrix_h */
