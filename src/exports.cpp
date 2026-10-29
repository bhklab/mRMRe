#include "exports.h"
#include <algorithm>
#include <climits>
#include <cstdio>
#include <exception>
#include <vector>
#include <R.h>
// borrowed from Matrix/rcpp
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef callEntries[] = {
    CALLDEF(export_concordance_index, 8),
    CALLDEF(export_filters, 15),
    CALLDEF(export_filters_bootstrap, 16),
    CALLDEF(export_mim, 12),
    CALLDEF(set_thread_count, 1),
    CALLDEF(get_thread_count, 0),
    {NULL, NULL, 0}
};

extern "C" void
R_init_mRMRe(DllInfo* info) {
  R_registerRoutines(info, NULL, callEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

// Sizes handed to Rf_allocVector are checked here, before any allocation and outside parallel
// regions; Rf_error is not safe to call from an OpenMP thread.
static R_xlen_t
checkedProduct(R_xlen_t const a, R_xlen_t const b)
{
    if (a < 0 || b < 0 || (b != 0 && a > R_XLEN_T_MAX / b))
        Rf_error("requested result is too large");
    return a * b;
}

// Number of solutions of a tree with the given children per level. Filter indexes its tree with
// unsigned int, so the total node count must fit in one.
static unsigned int
checkedSolutionCount(int const* const pChildrenCountPerLevel, R_xlen_t const levelCount)
{
    unsigned long long per_level = 1, total = 1;
    for (R_xlen_t i = 0; i < levelCount; ++i)
    {
        unsigned long long const c = pChildrenCountPerLevel[i] < 0 ? UINT_MAX : pChildrenCountPerLevel[i];
        if (c != 0 && per_level > UINT_MAX / c)
            Rf_error("the requested number of solutions is too large");
        per_level *= c;
        total += per_level;
        if (total > UINT_MAX)
            Rf_error("the requested number of solutions is too large");
    }
    return static_cast<unsigned int>(per_level);
}

// A C++ exception must not unwind into R: it is turned into an R error after leaving the handler
// (so the exception object and every local are destroyed first). Rf_error never returns.
#define MRMRE_TRY try {
#define MRMRE_CATCH \
    } catch (std::exception const& e) { std::snprintf(err, sizeof err, "%s", e.what()); } \
      catch (...) { std::snprintf(err, sizeof err, "unknown C++ exception"); } \
    Rf_error("%s", err);

extern "C" SEXP
export_concordance_index(SEXP samplesA, SEXP samplesB, SEXP samplesC, SEXP samplesD,
        SEXP sampleStrata, SEXP sampleWeights, SEXP sampleStratumCount, SEXP outX)
{
    char err[512] = "";
    MRMRE_TRY
    R_xlen_t const sample_count = XLENGTH(samplesA);
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    char const* const result_names[5] = { "ratio", "concordant", "discordant", "uninformative",
            "relevant" };
    for (int i = 0; i < 5; ++i)
    {
        SET_STRING_ELT(names, i, Rf_mkChar(result_names[i]));
        SET_VECTOR_ELT(result, i, Rf_allocVector(REALSXP, i == 0 ? 1 : sample_count));
        std::fill(REAL(VECTOR_ELT(result, i)), REAL(VECTOR_ELT(result, i)) + XLENGTH(VECTOR_ELT(result, i)), 0.0);
    }
    Rf_setAttrib(result, R_NamesSymbol, names);
    double* const ratio = REAL(VECTOR_ELT(result, 0));
    double* const concordantWeights = REAL(VECTOR_ELT(result, 1));
    double* const discordantWeights = REAL(VECTOR_ELT(result, 2));
    double* const uninformativeWeights = REAL(VECTOR_ELT(result, 3));
    double* const relevantWeights = REAL(VECTOR_ELT(result, 4));
    unsigned int const stratum_count = Rf_asInteger(sampleStratumCount);
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[stratum_count];
    unsigned int* const p_sample_count_per_stratum = new unsigned int[stratum_count];
    Math::placeStratificationData(INTEGER(sampleStrata), REAL(sampleWeights),
            p_sample_indices_per_stratum, p_sample_count_per_stratum,
            stratum_count, sample_count);

    if (LENGTH(samplesD) != 0 && LENGTH(samplesC) != 0)
        ratio[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                REAL(samplesC), REAL(samplesD), REAL(sampleWeights), p_sample_indices_per_stratum,
                p_sample_count_per_stratum, stratum_count, Rf_asLogical(outX) != 0,
                concordantWeights, discordantWeights, uninformativeWeights,
                relevantWeights);
    else if (LENGTH(samplesC) != 0)
        ratio[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                REAL(samplesC), REAL(sampleWeights), p_sample_indices_per_stratum,
                p_sample_count_per_stratum, stratum_count, Rf_asLogical(outX) != 0,
                concordantWeights, discordantWeights, uninformativeWeights,
                relevantWeights);
    else
        ratio[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                REAL(sampleWeights), p_sample_indices_per_stratum, p_sample_count_per_stratum,
                stratum_count, Rf_asLogical(outX) != 0, concordantWeights,
                discordantWeights, uninformativeWeights, relevantWeights);

    delete[] p_sample_count_per_stratum;
    for (unsigned int i = 0; i < stratum_count; ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;

    UNPROTECT(2);
    return result;
    MRMRE_CATCH
}

extern "C" SEXP
export_filters(SEXP childrenCountPerLevel, SEXP dataMatrix, SEXP priorsMatrix, SEXP priorsWeight,
        SEXP sampleStrata, SEXP sampleWeights, SEXP featureTypes, SEXP sampleCount,
        SEXP featureCount, SEXP sampleStratumCount, SEXP targetFeatureIndices, SEXP fixedFeatureCount,
        SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount)
{
    char err[512] = "";
    MRMRE_TRY
    unsigned int const feature_count = Rf_asInteger(featureCount);
    R_xlen_t const target_count = XLENGTH(targetFeatureIndices);
    R_xlen_t const mi_size = checkedProduct(feature_count, feature_count);
    R_xlen_t const feature_count_per_solution = XLENGTH(childrenCountPerLevel);
    unsigned int const solution_count =
            checkedSolutionCount(INTEGER(childrenCountPerLevel), feature_count_per_solution);
    R_xlen_t const chunk_size = checkedProduct(solution_count, feature_count_per_solution);

    // Every R object is allocated before the (large) C++ objects below.
    SEXP miMatrix = PROTECT(Rf_allocVector(REALSXP, mi_size));
    std::fill(REAL(miMatrix), REAL(miMatrix) + mi_size, NA_REAL);

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result, 0, Rf_allocVector(VECSXP, target_count));
    SET_VECTOR_ELT(result, 1, Rf_allocVector(VECSXP, target_count));
    SET_VECTOR_ELT(result, 2, Rf_allocVector(VECSXP, target_count));
    SET_VECTOR_ELT(result, 3, miMatrix);

    for (R_xlen_t i = 0; i < target_count; ++i)
    {
        SET_VECTOR_ELT(VECTOR_ELT(result, 0), i, Rf_allocVector(INTSXP, chunk_size));
        SET_VECTOR_ELT(VECTOR_ELT(result, 1), i, Rf_allocVector(REALSXP, feature_count));
        SET_VECTOR_ELT(VECTOR_ELT(result, 2), i, Rf_allocVector(REALSXP, chunk_size));
    }

    Matrix const priors_matrix(REAL(priorsMatrix), feature_count, feature_count);
    Matrix const* const p_priors_matrix =
            XLENGTH(priorsMatrix) == mi_size ? &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, Rf_asReal(priorsWeight), Rf_asInteger(sampleCount),
            feature_count, INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), Rf_asInteger(sampleStratumCount), Rf_asInteger(continuousEstimator),
            Rf_asLogical(outX) != 0, Rf_asInteger(bootstrapCount));
    MutualInformationMatrix mi_matrix(&data, REAL(miMatrix));

    for (R_xlen_t i = 0; i < target_count; ++i)
    {
        Filter filter(INTEGER(childrenCountPerLevel), feature_count_per_solution, &mi_matrix,
                INTEGER(targetFeatureIndices)[i], Rf_asInteger(fixedFeatureCount));
        filter.build();

        filter.getSolutions(INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), i)));
        filter.getScores(REAL(VECTOR_ELT(VECTOR_ELT(result, 2), i)));

        std::fill(REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i)),
                REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i)) + feature_count,
                std::numeric_limits<double>::quiet_NaN());

        Math::computeCausality(REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i)), &mi_matrix,
                INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), i)), solution_count,
                feature_count_per_solution, feature_count,
                INTEGER(targetFeatureIndices)[i]);
    }
    UNPROTECT(2);

    return result;
    MRMRE_CATCH
}

extern "C" SEXP
export_filters_bootstrap(SEXP solutionCount, SEXP solutionLength, SEXP dataMatrix,
        SEXP priorsMatrix, SEXP priorsWeight, SEXP sampleStrata, SEXP sampleWeights,
        SEXP featureTypes, SEXP sampleCount, SEXP featureCount, SEXP sampleStratumCount,
        SEXP targetFeatureIndices, SEXP fixedFeatureCount, SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount)
{
    char err[512] = "";
    MRMRE_TRY
    unsigned int const feature_count = Rf_asInteger(featureCount);
    R_xlen_t const target_count = XLENGTH(targetFeatureIndices);
    R_xlen_t const solution_count = Rf_asInteger(solutionCount);
    R_xlen_t const feature_count_per_solution = Rf_asInteger(solutionLength);
    R_xlen_t const chunk_size = checkedProduct(solution_count, feature_count_per_solution);
    R_xlen_t const mi_size = checkedProduct(feature_count, feature_count);

    // Every R object is allocated before the (large) C++ objects below.
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(result, 0, Rf_allocVector(VECSXP, target_count));
    SET_VECTOR_ELT(result, 1, Rf_allocVector(VECSXP, target_count));
    SET_VECTOR_ELT(result, 2, Rf_allocVector(VECSXP, target_count));

    for (R_xlen_t i = 0; i < target_count; ++i)
    {
        SET_VECTOR_ELT(VECTOR_ELT(result, 0), i, Rf_allocVector(INTSXP, chunk_size));
        SET_VECTOR_ELT(VECTOR_ELT(result, 1), i, Rf_allocVector(REALSXP, feature_count));
        SET_VECTOR_ELT(VECTOR_ELT(result, 2), i, Rf_allocVector(REALSXP, chunk_size));

        std::fill(REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i)),
                REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i)) + feature_count,
                std::numeric_limits<double>::quiet_NaN());
    }

    Matrix const priors_matrix(REAL(priorsMatrix), feature_count, feature_count);
    Matrix const* const p_priors_matrix =
            XLENGTH(priorsMatrix) == mi_size ? &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, Rf_asReal(priorsWeight), Rf_asInteger(sampleCount),
            feature_count, INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), Rf_asInteger(sampleStratumCount), Rf_asInteger(continuousEstimator),
            Rf_asLogical(outX) != 0, Rf_asInteger(bootstrapCount));

    std::vector<int> const children_count_per_level(feature_count_per_solution, 1);

    for (R_xlen_t i = 0; i < solution_count; ++i)
    {
        MutualInformationMatrix mi_matrix(&data);

        for (R_xlen_t j = 0; j < target_count; ++j)
        {
            Filter filter(children_count_per_level.data(), feature_count_per_solution, &mi_matrix,
                    INTEGER(targetFeatureIndices)[j], Rf_asInteger(fixedFeatureCount));
            filter.build();
            filter.getSolutions(
                    INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), j))
                            + (i * feature_count_per_solution));
            filter.getScores(
                    REAL(VECTOR_ELT(VECTOR_ELT(result, 2), j)) + (i * feature_count_per_solution));
        }

        data.bootstrap();
    }

    UNPROTECT(1);
    return result;
    MRMRE_CATCH
}

extern "C" SEXP
export_mim(SEXP dataMatrix, SEXP priorsMatrix, SEXP priorsWeight, SEXP sampleStrata,
        SEXP sampleWeights, SEXP featureTypes, SEXP sampleCount, SEXP featureCount,
        SEXP sampleStratumCount, SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount)
{
    char err[512] = "";
    MRMRE_TRY
    unsigned int const feature_count = Rf_asInteger(featureCount);
    R_xlen_t const mi_size = checkedProduct(feature_count, feature_count);

    SEXP miMatrix = PROTECT(Rf_allocVector(REALSXP, mi_size));
    std::fill(REAL(miMatrix), REAL(miMatrix) + mi_size, NA_REAL);

    Matrix const priors_matrix(REAL(priorsMatrix), feature_count, feature_count);
    Matrix const* const p_priors_matrix =
            XLENGTH(priorsMatrix) == mi_size ? &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, Rf_asReal(priorsWeight), Rf_asInteger(sampleCount),
            feature_count, INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), Rf_asInteger(sampleStratumCount), Rf_asInteger(continuousEstimator),
            Rf_asLogical(outX) != 0, Rf_asInteger(bootstrapCount));
    MutualInformationMatrix mi_matrix(&data, REAL(miMatrix));
    mi_matrix.build();

    UNPROTECT(1);
    return miMatrix;
    MRMRE_CATCH
}

extern "C" SEXP
get_thread_count()
{
    int count = 0;
#ifdef _OPENMP
    count = omp_get_max_threads();
#endif

    return Rf_ScalarInteger(count);
}

extern "C" SEXP
set_thread_count(SEXP threadCount)
{
    int count = 0;
#ifdef _OPENMP
    count = std::max(1, std::min(Rf_asInteger(threadCount), omp_get_num_procs()));
    omp_set_num_threads(count);
#endif

    return Rf_ScalarInteger(count);
}
