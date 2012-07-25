#ifndef ensemble_exports_hpp
#define ensemble_exports_hpp

#include <Rcpp.h>
#include <vector>

#include "MutualInformationMatrix.hpp"
#include "Tree.hpp"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_SampleStrata, SEXP R_SampleWeights, SEXP R_FeatureType,
        SEXP R_RowCount, SEXP R_ColumnCount, SEXP R_SampleStratumCount);

#endif /* ensemble_exports_hpp */
