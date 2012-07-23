#ifndef ensemble_MutualInformationMatrix_hpp
#define ensemble_MutualInformationMatrix_hpp

#include <limits>
#include <omp.h>

#include "SymmetricMatrix.hpp"
#include "tools.hpp"

class MutualInformationMatrix : public SymmetricMatrix
{
protected:
    Matrix* const mpDataMatrix;
    Matrix* const mpRankedDataMatrix;
    unsigned int* mpTestArray;

public:
    explicit
    MutualInformationMatrix(Matrix* const pMatrix);

    virtual
    ~MutualInformationMatrix();

    virtual float&
    operator()(unsigned int const i, unsigned int const j);

    void const
    build();
};

#endif /* ensemble_MutualInformationMatrix_hpp */
