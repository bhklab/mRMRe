#ifndef mRMRe_Matrix_h
#define mRMRe_Matrix_h

#include <vector>

class Matrix
{
private:
    Matrix(const Matrix&);

    Matrix&
    operator=(const Matrix&);

protected:
    double* const mpData;
    unsigned int const mRowCount;
    unsigned int const mColumnCount;
    bool const mHasAllocation;

public:
    Matrix(unsigned int const rowCount, unsigned int const columnCount);

    explicit
    Matrix(unsigned int const size, unsigned int const rowCount, unsigned int const columnCount);

    explicit
    Matrix(double* const pData, unsigned int const rowCount, unsigned int const columnCount);

    virtual
    ~Matrix();

    virtual double&
    at(unsigned int const i, unsigned int const j);

    virtual double const&
    at(unsigned int const i, unsigned int const j) const;

    unsigned int const
    getColumnCount() const;

    unsigned int const
    getRowCount() const;
};

#endif /* mRMRe_Matrix_h */
