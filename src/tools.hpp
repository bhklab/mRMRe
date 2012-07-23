#ifndef ensemble_tools_hpp
#define ensemble_tools_hpp

#include <algorithm>
#include <cmath>

#include "Matrix.hpp"

class DataMatrixComparator
{
private:
    unsigned int const mFeatureIndex;
    Matrix* const mpDataMatrix;

public:
    DataMatrixComparator(unsigned int const featureIndex, Matrix* const pDataMatrix);

    bool const
    operator()(unsigned int const i, unsigned int const j) const;
};

void const
placeRanksByFeatureIndex(unsigned int const index, Matrix* const pRankedDataMatrix,
        Matrix* const pDataMatrix);

float const
computeSpearmanCorrelation(unsigned int const i, unsigned int const j,
        Matrix* const pRankedDataMatrix);

#endif /* ensemble_tools_hpp */