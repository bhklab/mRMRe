#ifndef ensemble_SymmetricMatrix_h
#define ensemble_SymmetricMatrix_h

#include "Matrix.h"

class SymmetricMatrix : public Matrix
{
public:
    SymmetricMatrix(unsigned int const rowCount);

    SymmetricMatrix(float* const pData, unsigned int const rowCount);

    virtual
    ~SymmetricMatrix();

    virtual float&
    at(unsigned int const i, unsigned int const j);
};

#endif /* ensemble_SymmetricMatrix_h */
