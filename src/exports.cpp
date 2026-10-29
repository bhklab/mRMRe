#include "exports.h"
#include <algorithm>
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

extern "C" SEXP
export_concordance_index(SEXP samplesA, SEXP samplesB, SEXP samplesC, SEXP samplesD,
        SEXP sampleStrata, SEXP sampleWeights, SEXP sampleStratumCount, SEXP outX)
{
    unsigned int const sample_count = LENGTH(samplesA);
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    char const* const result_names[5] = { "ratio", "concordant", "discordant", "uninformative",
            "relevant" };
    for (int i = 0; i < 5; ++i)
    {
        SET_STRING_ELT(names, i, Rf_mkChar(result_names[i]));
        SET_VECTOR_ELT(result, i, Rf_allocVector(REALSXP, i == 0 ? 1 : sample_count));
        std::fill(REAL(VECTOR_ELT(result, i)), REAL(VECTOR_ELT(result, i)) + LENGTH(VECTOR_ELT(result, i)), 0.0);
    }
    Rf_setAttrib(result, R_NamesSymbol, names);
    double* const ratio = REAL(VECTOR_ELT(result, 0));
    double* const concordantWeights = REAL(VECTOR_ELT(result, 1));
    double* const discordantWeights = REAL(VECTOR_ELT(result, 2));
    double* const uninformativeWeights = REAL(VECTOR_ELT(result, 3));
    double* const relevantWeights = REAL(VECTOR_ELT(result, 4));
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[Rf_asInteger(sampleStratumCount)];
    unsigned int* const p_sample_count_per_stratum =
            new unsigned int[Rf_asInteger(sampleStratumCount)];
    Math::placeStratificationData(INTEGER(sampleStrata), REAL(sampleWeights),
            p_sample_indices_per_stratum, p_sample_count_per_stratum,
            Rf_asInteger(sampleStratumCount), sample_count);

    if (LENGTH(samplesD) != 0 && LENGTH(samplesC) != 0)
        ratio[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                REAL(samplesC), REAL(samplesD), REAL(sampleWeights), p_sample_indices_per_stratum,
                p_sample_count_per_stratum, Rf_asInteger(sampleStratumCount), Rf_asLogical(outX) != 0,
                concordantWeights, discordantWeights, uninformativeWeights,
                relevantWeights);
    else if (LENGTH(samplesC) != 0)
        ratio[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                REAL(samplesC), REAL(sampleWeights), p_sample_indices_per_stratum,
                p_sample_count_per_stratum, Rf_asInteger(sampleStratumCount), Rf_asLogical(outX) != 0,
                concordantWeights, discordantWeights, uninformativeWeights,
                relevantWeights);
    else
        ratio[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                REAL(sampleWeights), p_sample_indices_per_stratum, p_sample_count_per_stratum,
                Rf_asInteger(sampleStratumCount), Rf_asLogical(outX) != 0, concordantWeights,
                discordantWeights, uninformativeWeights, relevantWeights);

    delete[] p_sample_count_per_stratum;
    for (unsigned int i = 0; i < Rf_asInteger(sampleStratumCount); ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;

    UNPROTECT(2);
    return result;
}

extern "C" SEXP
export_filters(SEXP childrenCountPerLevel, SEXP dataMatrix, SEXP priorsMatrix, SEXP priorsWeight,
        SEXP sampleStrata, SEXP sampleWeights, SEXP featureTypes, SEXP sampleCount,
        SEXP featureCount, SEXP sampleStratumCount, SEXP targetFeatureIndices, SEXP fixedFeatureCount,
        SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount)
{
    Matrix const priors_matrix(REAL(priorsMatrix), Rf_asInteger(featureCount),
            Rf_asInteger(featureCount));
    Matrix const* const p_priors_matrix =
            LENGTH(priorsMatrix) == Rf_asInteger(featureCount) * Rf_asInteger(featureCount) ?
                    &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, Rf_asReal(priorsWeight), Rf_asInteger(sampleCount),
            Rf_asInteger(featureCount), INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), Rf_asInteger(sampleStratumCount), Rf_asInteger(continuousEstimator),
            Rf_asLogical(outX) != 0, Rf_asInteger(bootstrapCount));
    SEXP miMatrix = PROTECT(Rf_allocVector(REALSXP,
            (R_xlen_t) Rf_asInteger(featureCount) * Rf_asInteger(featureCount)));
    std::fill(REAL(miMatrix), REAL(miMatrix) + XLENGTH(miMatrix), NA_REAL);
    MutualInformationMatrix mi_matrix(&data, REAL(miMatrix));

    unsigned int solution_count = 1;
    for (unsigned int i = 0; i < LENGTH(childrenCountPerLevel); ++i)
        solution_count *= INTEGER(childrenCountPerLevel)[i];
    unsigned int const feature_count_per_solution = LENGTH(childrenCountPerLevel);
    unsigned int const chunk_size = solution_count * feature_count_per_solution;

    SEXP result;
    PROTECT(result = Rf_allocVector(VECSXP, 4));

    SET_VECTOR_ELT(result, 0, Rf_allocVector(VECSXP, LENGTH(targetFeatureIndices)));
    SET_VECTOR_ELT(result, 1, Rf_allocVector(VECSXP, LENGTH(targetFeatureIndices)));
    SET_VECTOR_ELT(result, 2, Rf_allocVector(VECSXP, LENGTH(targetFeatureIndices)));
    SET_VECTOR_ELT(result, 3, miMatrix);

    for (unsigned int i = 0; i < LENGTH(targetFeatureIndices); ++i)
    {
        Filter filter(INTEGER(childrenCountPerLevel), LENGTH(childrenCountPerLevel), &mi_matrix,
                INTEGER(targetFeatureIndices)[i], Rf_asInteger(fixedFeatureCount));
        filter.build();

        SET_VECTOR_ELT(VECTOR_ELT(result, 0), i, Rf_allocVector(INTSXP, chunk_size));
        SET_VECTOR_ELT(VECTOR_ELT(result, 1), i, Rf_allocVector(REALSXP, Rf_asInteger(featureCount)));
        SET_VECTOR_ELT(VECTOR_ELT(result, 2), i, Rf_allocVector(REALSXP, chunk_size));

        filter.getSolutions(INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), i)));
        filter.getScores(REAL(VECTOR_ELT(VECTOR_ELT(result, 2), i)));

        for (unsigned int k = 0; k < Rf_asInteger(featureCount); ++k)
            REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i))[k] =
                    std::numeric_limits<double>::quiet_NaN();

        Math::computeCausality(REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i)), &mi_matrix,
                INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), i)), solution_count,
                feature_count_per_solution, Rf_asInteger(featureCount),
                INTEGER(targetFeatureIndices)[i]);
    }
    UNPROTECT(2);

    return result;
}

extern "C" SEXP
export_filters_bootstrap(SEXP solutionCount, SEXP solutionLength, SEXP dataMatrix,
        SEXP priorsMatrix, SEXP priorsWeight, SEXP sampleStrata, SEXP sampleWeights,
        SEXP featureTypes, SEXP sampleCount, SEXP featureCount, SEXP sampleStratumCount,
        SEXP targetFeatureIndices, SEXP fixedFeatureCount, SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount)
{
    Matrix const priors_matrix(REAL(priorsMatrix), Rf_asInteger(featureCount),
            Rf_asInteger(featureCount));
    Matrix const* const p_priors_matrix =
            LENGTH(priorsMatrix) == Rf_asInteger(featureCount) * Rf_asInteger(featureCount) ?
                    &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, Rf_asReal(priorsWeight), Rf_asInteger(sampleCount),
            Rf_asInteger(featureCount), INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), Rf_asInteger(sampleStratumCount), Rf_asInteger(continuousEstimator),
            Rf_asLogical(outX) != 0, Rf_asInteger(bootstrapCount));

    unsigned int solution_count = Rf_asInteger(solutionCount);
    unsigned int const feature_count_per_solution = Rf_asInteger(solutionLength);
    unsigned int const chunk_size = solution_count * feature_count_per_solution;

    int* const p_children_count_per_level = new int[feature_count_per_solution];
    for (unsigned int i = 0; i < feature_count_per_solution; ++i)
        p_children_count_per_level[i] = 1;

    SEXP result;
    PROTECT(result = Rf_allocVector(VECSXP, 3));

    SET_VECTOR_ELT(result, 0, Rf_allocVector(VECSXP, LENGTH(targetFeatureIndices)));
    SET_VECTOR_ELT(result, 1, Rf_allocVector(VECSXP, LENGTH(targetFeatureIndices)));
    SET_VECTOR_ELT(result, 2, Rf_allocVector(VECSXP, LENGTH(targetFeatureIndices)));

    for (unsigned int i = 0; i < LENGTH(targetFeatureIndices); ++i)
    {
        SET_VECTOR_ELT(VECTOR_ELT(result, 0), i, Rf_allocVector(INTSXP, chunk_size));
        SET_VECTOR_ELT(VECTOR_ELT(result, 1), i, Rf_allocVector(REALSXP, Rf_asInteger(featureCount)));
        SET_VECTOR_ELT(VECTOR_ELT(result, 2), i, Rf_allocVector(REALSXP, chunk_size));

        for (unsigned int k = 0; k < Rf_asInteger(featureCount); ++k)
            REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i))[k] =
                    std::numeric_limits<double>::quiet_NaN();
    }

    for (unsigned int i = 0; i < solution_count; ++i)
    {
        MutualInformationMatrix mi_matrix(&data);

        for (unsigned int j = 0; j < LENGTH(targetFeatureIndices); ++j)
        {
            Filter filter(p_children_count_per_level, feature_count_per_solution, &mi_matrix,
                    INTEGER(targetFeatureIndices)[j], Rf_asInteger(fixedFeatureCount));
            filter.build();
            filter.getSolutions(
                    INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), j))
                            + (i * feature_count_per_solution));
            filter.getScores(
                    REAL(VECTOR_ELT(VECTOR_ELT(result, 2), j)) + (i * feature_count_per_solution));

            /*            Math::computeCausality(REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i)), &mi_matrix,
             INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), i)) + (i * chunk_size), 1,
             feature_count_per_solution, Rf_asInteger(featureCount),
             INTEGER(targetFeatureIndices)[i]);*/
        }

        data.bootstrap();
    }

    UNPROTECT(1);
    delete[] p_children_count_per_level;
    return result;
}

extern "C" SEXP
export_mim(SEXP dataMatrix, SEXP priorsMatrix, SEXP priorsWeight, SEXP sampleStrata,
        SEXP sampleWeights, SEXP featureTypes, SEXP sampleCount, SEXP featureCount,
        SEXP sampleStratumCount, SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount)
{
    Matrix const priors_matrix(REAL(priorsMatrix), Rf_asInteger(featureCount),
            Rf_asInteger(featureCount));
    Matrix const* const p_priors_matrix =
            LENGTH(priorsMatrix) == Rf_asInteger(featureCount) * Rf_asInteger(featureCount) ?
                    &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, Rf_asReal(priorsWeight), Rf_asInteger(sampleCount),
            Rf_asInteger(featureCount), INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), Rf_asInteger(sampleStratumCount), Rf_asInteger(continuousEstimator),
            Rf_asLogical(outX) != 0, Rf_asInteger(bootstrapCount));
    SEXP miMatrix = PROTECT(Rf_allocVector(REALSXP,
            (R_xlen_t) Rf_asInteger(featureCount) * Rf_asInteger(featureCount)));
    std::fill(REAL(miMatrix), REAL(miMatrix) + XLENGTH(miMatrix), NA_REAL);
    MutualInformationMatrix mi_matrix(&data, REAL(miMatrix));
    mi_matrix.build();

    UNPROTECT(1);
    return miMatrix;
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
