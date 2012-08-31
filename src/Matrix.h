#ifndef ensemble_Matrix_h
#define ensemble_Matrix_h

#include <omp.h>
#include <vector>

class Matrix
{
private:
    Matrix(const Matrix&);

    Matrix&
    operator=(const Matrix&);

protected:
    float* const mpData;
    unsigned int const mRowCount;
    unsigned int const mColumnCount;
    bool const mHasAllocation;

public:
    Matrix(unsigned int const rowCount, unsigned int const columnCount);

    explicit
    Matrix(unsigned int const size, unsigned int const rowCount, unsigned int const columnCount);

    explicit
    Matrix(float* const pData, unsigned int const rowCount, unsigned int const columnCount);

    virtual
    ~Matrix();

    virtual float&
    at(unsigned int const i, unsigned int const j);

    virtual float const&
    at(unsigned int const i, unsigned int const j) const;

    unsigned int const
    getColumnCount() const;

    unsigned int const
    getRowCount() const;

    operator std::vector<float>() const;
};

#endif /* ensemble_Matrix_h */
