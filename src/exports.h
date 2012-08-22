#ifndef ensemble_exports_h
#define ensemble_exports_h

#include <Rcpp.h>
#include <vector>

#include "Math.h"
#include "MutualInformationMatrix.h"
#include "Tree.h"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_SampleStrata, SEXP R_SampleWeights, SEXP R_FeatureTypes,
        SEXP R_SampleCount, SEXP R_FeatureCount, SEXP R_SampleStratumCount, SEXP R_UsesRanks,
        SEXP R_OutX, SEXP R_BootstrapCount);

extern "C" SEXP
build_mRMR_tree(SEXP R_ChildrenCountPerLevel, SEXP R_DataMatrix, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_FeatureTypes, SEXP R_SampleCount, SEXP R_FeatureCount,
        SEXP R_SampleStratumCount, SEXP R_TargetFeatureIndex, SEXP R_UsesRanks, SEXP R_OutX,
        SEXP R_BootstrapCount);

extern "C" SEXP
compute_concordance_index(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_SampleIndicesPerStratum, SEXP R_OutX);

extern "C" SEXP
compute_concordance_index_with_time(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_Time,
        SEXP R_SampleWeights, SEXP R_SampleStrata, SEXP R_SampleStratumCount, SEXP R_OutX);

extern "C" SEXP
compute_cramers_v(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights, SEXP R_SampleStrata,
        SEXP R_SampleStratumCount, SEXP R_BootstrapCount);

extern "C" SEXP
compute_pearson_correlation(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights,
        SEXP R_SampleStrata, SEXP R_SampleStratumCount, SEXP R_BootstrapCount);

extern "C" SEXP
compute_spearman_correlation(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights,
        SEXP R_SampleStrata, SEXP R_SampleStratumCount, SEXP R_BootstrapCount);

extern "C" SEXP
set_thread_count(SEXP R_ThreadCount);

#endif /* ensemble_exports_h */
