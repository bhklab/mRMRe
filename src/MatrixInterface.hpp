#ifndef ensemble_MatrixInterface_hpp
#define ensemble_MatrixInterface_hpp

template<typename T>
class MatrixInterface
{
public:
    virtual
    ~MatrixInterface() {};

    virtual inline T&
    operator()(unsigned int const i, unsigned int const j);
};

#endif /* ensemble_MatrixInterface_hpp */
