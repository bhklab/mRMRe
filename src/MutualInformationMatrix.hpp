#ifndef ensemble_MutualInformationMatrix_hpp
#define ensemble_MutualInformationMatrix_hpp

#include "SymmetricMatrix.hpp"

class MutualInformationMatrix : public SymmetricMatrix
{
public:
    MutualInformationMatrix(Matrix* const pMatrix);

    virtual
    ~MutualInformationMatrix();

    virtual float&
    operator()(unsigned int const i, unsigned int const j);
};

#endif /* ensemble_MutualInformationMatrix_hpp */
