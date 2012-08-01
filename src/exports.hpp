#ifndef ensemble_exports_hpp
#define ensemble_exports_hpp

#include <Rcpp.h>
#include <vector>

#include "MutualInformationMatrix.hpp"
#include "Tree.hpp"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_SampleStrata, SEXP R_SampleWeights, SEXP R_FeatureTypes,
        SEXP R_SampleCount, SEXP R_FeatureCount, SEXP R_SampleStratumCount);

extern "C" SEXP
build_mRMR_tree_from_data(SEXP R_ChildrenCountPerLevel, SEXP R_DataMatrix, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_FeatureTypes, SEXP R_SampleCount, SEXP R_FeatureCount,
        SEXP R_SampleStratumCount, SEXP R_TargetFeatureIndex);

extern "C" SEXP
build_mRMR_tree_from_mim(SEXP R_ChildrenCountPerLevel, SEXP R_MiMatrix, SEXP R_FeatureCount,
        SEXP R_TargetFeatureIndex);

#endif /* ensemble_exports_hpp */
