#ifndef ensemble_exports_hpp
#define ensemble_exports_hpp

#include <Rcpp.h>
#include <vector>

#include "Math.hpp"
#include "MutualInformationMatrix.hpp"
#include "Tree.hpp"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_SampleStrata, SEXP R_SampleWeights, SEXP R_FeatureTypes,
        SEXP R_SampleCount, SEXP R_FeatureCount, SEXP R_SampleStratumCount, SEXP R_UsesRanks,
        SEXP R_OutX);

extern "C" SEXP
build_mRMR_tree_from_data(SEXP R_ChildrenCountPerLevel, SEXP R_DataMatrix, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_FeatureTypes, SEXP R_SampleCount, SEXP R_FeatureCount,
        SEXP R_SampleStratumCount, SEXP R_TargetFeatureIndex, SEXP R_UsesRanks, SEXP R_OutX);

extern "C" SEXP
build_mRMR_tree_from_mim(SEXP R_ChildrenCountPerLevel, SEXP R_MiMatrix, SEXP R_FeatureCount,
        SEXP R_TargetFeatureIndex);

extern "C" SEXP
compute_concordance_index(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_SampleIndicesPerStratum, SEXP R_outX);

extern "C" SEXP
compute_concordance_index_with_time(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_Time,
        SEXP R_SampleWeights, SEXP R_SampleStrata, SEXP R_SampleStratumCount, SEXP R_outX);

extern "C" SEXP
compute_cramers_v(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights, SEXP R_SampleStrata,
        SEXP R_SampleStratumCount);

extern "C" SEXP
compute_pearson_correlation(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights,
        SEXP R_SampleStrata, SEXP R_SampleStratumCount);

extern "C" SEXP
compute_spearman_correlation(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights,
        SEXP R_SampleStrata, SEXP R_SampleStratumCount);

#endif /* ensemble_exports_hpp */
