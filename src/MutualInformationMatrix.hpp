#ifndef ensemble_MutualInformationMatrix_hpp
#define ensemble_MutualInformationMatrix_hpp

#include <algorithm>
#include <cmath>
#include <math.h>
#include <limits>
#include <omp.h>

#include "SymmetricMatrix.hpp"

class MutualInformationMatrix : public SymmetricMatrix
{
protected:
    Matrix* const mpDataMatrix;
    Matrix* const mpRankedDataMatrix;

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
    explicit
    MutualInformationMatrix(Matrix* const pMatrix);

    virtual
    ~MutualInformationMatrix();

    virtual float&
    operator()(unsigned int const i, unsigned int const j);

    void const
    build();

    void const
    placeRanksByFeatureIndex(unsigned int const dataMatrixFeatureIndex);

    void const
    placeSpearmanCorrelation(unsigned int const i, unsigned int const j);
};

#endif /* ensemble_MutualInformationMatrix_hpp */
