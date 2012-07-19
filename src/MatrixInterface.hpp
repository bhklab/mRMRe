#ifndef ensemble_MatrixInterface_hpp
#define ensemble_MatrixInterface_hpp

class MatrixInterface
{
private:
    float* const mpData;
    unsigned int const mRowCount;
    unsigned int const mColumnCount;
    bool const mHasAllocation;

public:
    virtual
    ~MatrixInterface();

    virtual float&
    operator()(unsigned int const i, unsigned int const j) = 0;

    virtual unsigned int const
    getRowCount() const = 0;

    virtual unsigned int const
    getColumnCount() const = 0;
};

#endif /* ensemble_MatrixInterface_hpp */
