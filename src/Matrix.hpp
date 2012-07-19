#ifndef ensemble_Matrix_hpp
#define ensemble_Matrix_hpp

class Matrix
{
private:
    float* const mpData;
    unsigned int const mRowCount;
    unsigned int const mColumnCount;
    bool const mHasAllocation;

protected:
    Matrix();

public:
    Matrix(unsigned int const rowCount, unsigned int const columnCount);

    Matrix(float* const data, unsigned int const rowCount, unsigned int const columnCount);

    ~Matrix();

    virtual inline float&
    operator()(unsigned int const i, unsigned int const j);

    inline unsigned int const
    getRowCount() const;

    inline unsigned int const
    getColumnCount() const;
};

#endif /* ensemble_Matrix_hpp */
