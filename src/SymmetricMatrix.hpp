#ifndef ensemble_SymmetricMatrix_hpp
#define ensemble_SymmetricMatrix_hpp

#include "MatrixInterface.hpp"

class SymmetricMatrix : public MatrixInterface
{
public:
    SymmetricMatrix(unsigned int const rowCount);

    SymmetricMatrix(float* const data, unsigned int const rowCount);

    ~SymmetricMatrix();

    inline float&
    operator()(unsigned int const i, unsigned int const j);

    inline unsigned int const
    getRowCount() const;

    inline unsigned int const
    getColumnCount() const;
};

#endif /* ensemble_SymmetricMatrix_hpp */
