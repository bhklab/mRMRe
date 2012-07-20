#ifndef ensemble_MutualInformationMatrix_hpp
#define ensemble_MutualInformationMatrix_hpp

#include <algorithm>
#include <cmath>
#include <limits>

#include "SymmetricMatrix.hpp"

class MutualInformationMatrix : public SymmetricMatrix
{
protected:
    Matrix* const mpDataMatrix;
    Matrix* mpRankedDataMatrix;

    class DataMatrixComparator
    {
    private:
        unsigned int const mDataMatrixFeatureIndex;
        Matrix* const mpDataMatrix;

    public:
        DataMatrixComparator(unsigned int const dataMatrixFeatureIndex, Matrix* const pDataMatrix);

        bool const
        operator()(unsigned int const i, unsigned int const j) const;
    };

public:
    MutualInformationMatrix(Matrix* const pMatrix);

    virtual
    ~MutualInformationMatrix();

    virtual float&
    operator()(unsigned int const i, unsigned int const j);

    void const
    placeRanksByFeatureIndex(unsigned int const dataMatrixFeatureIndex);

    void const
    placeSpearmanCorrelation(unsigned int const i, unsigned int const j);
};

#endif /* ensemble_MutualInformationMatrix_hpp */
