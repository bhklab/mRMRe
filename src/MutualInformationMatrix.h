#ifndef mRMRe_MutualInformationMatrix_h
#define mRMRe_MutualInformationMatrix_h

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>

#include "Data.h"
#include "Matrix.h"

class MutualInformationMatrix : public Matrix
{
private:
    MutualInformationMatrix(const MutualInformationMatrix&);

    MutualInformationMatrix&
    operator=(const MutualInformationMatrix&);

protected:
    Data const* const mpData;
    std::vector<char> mIsComputed; // char, not bool: adjacent flags must be writable from different threads

public:
    MutualInformationMatrix(Data const* const pData);

    MutualInformationMatrix(Data const* const pData, double* const pInternalData);

    virtual
    ~MutualInformationMatrix();

    virtual double&
    at(unsigned int const i, unsigned int const j);

    virtual double const&
    at(unsigned int const i, unsigned int const j) const;

    void const
    build();
};

#endif /* mRMRe_MutualInformationMatrix_h */
