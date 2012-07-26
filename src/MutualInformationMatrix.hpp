#ifndef ensemble_MutualInformationMatrix_hpp
#define ensemble_MutualInformationMatrix_hpp

#include "Data.hpp"
#include "SymmetricMatrix.hpp"

class MutualInformationMatrix : public SymmetricMatrix
{
protected:
    Data const* const mpData;

public:
    MutualInformationMatrix(Data const* const pData);

    virtual
    ~MutualInformationMatrix();

    virtual float&
    operator()(unsigned int const i, unsigned int const j);
};

#endif /* ensemble_MutualInformationMatrix_hpp */
