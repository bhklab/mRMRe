#ifndef mRMRe_exports_h
#define mRMRe_exports_h

#define USE_RINTERNALS

#include <cstdlib>
#include <vector>

#include "Filter.h"
#include "Math.h"
#include "MutualInformationMatrix.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

extern "C" void
R_init_mRMRe(DllInfo* info);

extern "C" SEXP
export_concordance_index(SEXP samplesA, SEXP samplesB, SEXP samplesC, SEXP samplesD,
        SEXP sampleStrata, SEXP sampleWeights, SEXP sampleStratumCount, SEXP outX, SEXP ratio,
        SEXP concordantWeights, SEXP discordantWeights, SEXP uninformativeWeights,
        SEXP relevantWeights);

extern "C" SEXP
export_filters(SEXP childrenCountPerLevel, SEXP dataMatrix, SEXP priorsMatrix, SEXP priorsWeight,
        SEXP sampleStrata, SEXP sampleWeights, SEXP featureTypes, SEXP sampleCount,
        SEXP featureCount, SEXP sampleStratumCount, SEXP targetFeatureIndices, SEXP fixedFeatureCount,
        SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount, SEXP miMatrix);

extern "C" SEXP
export_filters_bootstrap(SEXP solutionCount, SEXP solutionLength, SEXP dataMatrix,
        SEXP priorsMatrix, SEXP priorsWeight, SEXP sampleStrata, SEXP sampleWeights,
        SEXP featureTypes, SEXP sampleCount, SEXP featureCount, SEXP sampleStratumCount,
        SEXP targetFeatureIndices, SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount, SEXP fixedFeatureCount,
        SEXP miMatrix);

extern "C" SEXP
export_mim(SEXP dataMatrix, SEXP priorsMatrix, SEXP priorsWeight, SEXP sampleStrata,
        SEXP sampleWeights, SEXP featureTypes, SEXP sampleCount, SEXP featureCount,
        SEXP sampleStratumCount, SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount,
        SEXP miMatrix);

extern "C" SEXP
get_thread_count(SEXP threadCount);

extern "C" SEXP
set_thread_count(SEXP threadCount);

#endif /* mRMRe_exports_h */
