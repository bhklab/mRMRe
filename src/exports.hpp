#ifndef ensemble_exports_hpp
#define ensemble_exports_hpp

#include <cmath>
#include <omp.h>
#include <Rcpp.h>
#include <vector>

#include "MutualInformationMatrix.hpp"
#include "Tree.hpp"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_Strata, SEXP R_Weights, SEXP R_FeatureType, SEXP R_RowCount,
        SEXP R_ColumnCount);

/*
 extern "C" SEXP
 mRMR_filter_with_mim(SEXP R_FeatureInformationMatrix, SEXP R_ChildrenCountPerLevel,
 SEXP R_TargetFeatureIndex);

 extern "C" SEXP
 mRMR_filter_with_data(SEXP R_DataMatrix, SEXP R_Strata, SEXP R_Weights, SEXP R_FeatureType,
 SEXP R_RowCount, SEXP R_ColumnCount, SEXP R_ChildrenCountPerLevel,
 SEXP R_TargetFeatureIndex);
 */

#endif /* ensemble_exports_hpp */
