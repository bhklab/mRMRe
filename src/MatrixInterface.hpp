#ifndef ensemble_MatrixInterface_hpp
#define ensemble_MatrixInterface_hpp

class MatrixInterface
{
protected:
    float* const mpData;
    unsigned int const mRowCount;
    unsigned int const mColumnCount;
    bool const mHasAllocation;

public:
    MatrixInterface(unsigned int const size, unsigned int const rowCount,
            unsigned const columnCount);

    MatrixInterface(float* const data, unsigned int const rowCount, unsigned const columnCount);

    ~MatrixInterface();

    virtual float&
    operator()(unsigned int const i, unsigned int const j) = 0;

    virtual unsigned int const
    getRowCount() const = 0;

    virtual unsigned int const
    getColumnCount() const = 0;
};

#endif /* ensemble_MatrixInterface_hpp */
