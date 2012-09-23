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

    SymmetricMatrix(double* const pData, unsigned int const rowCount);

    virtual
    ~SymmetricMatrix();

    virtual double&
    at(unsigned int const i, unsigned int const j);

    virtual double const&
    at(unsigned int const i, unsigned int const j) const;
};

#endif /* mRMRe_SymmetricMatrix_h */
