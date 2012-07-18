#ifndef ensemble_MatrixInterface_hpp
#define ensemble_MatrixInterface_hpp

class MatrixInterface
{
public:
    virtual
    ~MatrixInterface() {};

    virtual float&
    operator()(unsigned int const i, unsigned int const j) = 0;

    virtual unsigned int const
    getRowCount() const = 0;

    virtual unsigned int const
    getColumnCount() const = 0;
};

#endif /* ensemble_MatrixInterface_hpp */
