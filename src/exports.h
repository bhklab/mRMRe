#ifndef mRMRe_exports_h
#define mRMRe_exports_h

#include <Rcpp.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <string>
#include <vector>

#include "Filter.h"
#include "Math.h"
#include "MutualInformationMatrix.h"

/*
 extern "C" SEXP
 export_association(SEXP R_SamplesA, SEXP R_SamplesB, SEXP R_SamplesC, SEXP R_SampleStrata,
 SEXP R_SampleWeights, SEXP R_SampleStratumCount, SEXP R_OutX, SEXP R_BootstrapCount,
 SEXP R_Method);

 extern "C" SEXP
 export_filter(SEXP R_ChildrenCountPerLevel, SEXP R_DataMatrix, SEXP R_PriorsMatrix,
 SEXP R_PriorsWeight, SEXP R_SampleStrata, SEXP R_SampleWeights, SEXP R_FeatureTypes,
 SEXP R_SampleCount, SEXP R_FeatureCount, SEXP R_SampleStratumCount,
 SEXP R_TargetFeatureIndex, SEXP R_UsesRanks, SEXP R_OutX, SEXP R_BootstrapCount);
 */

extern "C" SEXP
get_thread_count(SEXP threadCount);

extern "C" SEXP
set_thread_count(SEXP threadCount);

#endif /* mRMRe_exports_h */
