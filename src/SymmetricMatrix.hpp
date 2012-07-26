#ifndef ensemble_SymmetricMatrix_hpp
#define ensemble_SymmetricMatrix_hpp

#include "Matrix.hpp"

class SymmetricMatrix : public Matrix
{
public:
    SymmetricMatrix(unsigned int const rowCount);

    virtual
    ~SymmetricMatrix();

    virtual float&
    at(unsigned int const i, unsigned int const j);
};

#endif /* ensemble_SymmetricMatrix_hpp */
