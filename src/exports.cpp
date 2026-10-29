#include "exports.h"
#include <R.h>
// borrowed from Matrix/rcpp
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef callEntries[] = {
    CALLDEF(export_concordance_index, 13),
    CALLDEF(export_filters, 16),
    CALLDEF(export_filters_bootstrap, 17),
    CALLDEF(export_mim, 13),
    CALLDEF(set_thread_count, 1),
    CALLDEF(get_thread_count, 1),
    {NULL, NULL, 0}
};

extern "C" void
R_init_mRMRe(DllInfo* info) {
  R_registerRoutines(info, NULL, callEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

extern "C" SEXP
export_concordance_index(SEXP samplesA, SEXP samplesB, SEXP samplesC, SEXP samplesD,
        SEXP sampleStrata, SEXP sampleWeights, SEXP sampleStratumCount, SEXP outX, SEXP ratio,
        SEXP concordantWeights, SEXP discordantWeights, SEXP uninformativeWeights,
        SEXP relevantWeights)
{
    unsigned int const sample_count = LENGTH(samplesA);
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[INTEGER(sampleStratumCount)[0]];
    unsigned int* const p_sample_count_per_stratum =
            new unsigned int[INTEGER(sampleStratumCount)[0]];
    Math::placeStratificationData(INTEGER(sampleStrata), REAL(sampleWeights),
            p_sample_indices_per_stratum, p_sample_count_per_stratum,
            INTEGER(sampleStratumCount)[0], sample_count);

    if (LENGTH(samplesD) != 0 && LENGTH(samplesC) != 0)
        REAL(ratio)[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                REAL(samplesC), REAL(samplesD), REAL(sampleWeights), p_sample_indices_per_stratum,
                p_sample_count_per_stratum, INTEGER(sampleStratumCount)[0], INTEGER(outX)[0] != 0,
                REAL(concordantWeights), REAL(discordantWeights), REAL(uninformativeWeights),
                REAL(relevantWeights));
    else if (LENGTH(samplesC) != 0)
        REAL(ratio)[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                REAL(samplesC), REAL(sampleWeights), p_sample_indices_per_stratum,
                p_sample_count_per_stratum, INTEGER(sampleStratumCount)[0], INTEGER(outX)[0] != 0,
                REAL(concordantWeights), REAL(discordantWeights), REAL(uninformativeWeights),
                REAL(relevantWeights));
    else
        REAL(ratio)[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                REAL(sampleWeights), p_sample_indices_per_stratum, p_sample_count_per_stratum,
                INTEGER(sampleStratumCount)[0], INTEGER(outX)[0] != 0, REAL(concordantWeights),
                REAL(discordantWeights), REAL(uninformativeWeights), REAL(relevantWeights));

    delete[] p_sample_count_per_stratum;
    for (unsigned int i = 0; i < INTEGER(sampleStratumCount)[0]; ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;

    return R_NilValue;
}

extern "C" SEXP
export_filters(SEXP childrenCountPerLevel, SEXP dataMatrix, SEXP priorsMatrix, SEXP priorsWeight,
        SEXP sampleStrata, SEXP sampleWeights, SEXP featureTypes, SEXP sampleCount,
        SEXP featureCount, SEXP sampleStratumCount, SEXP targetFeatureIndices, SEXP fixedFeatureCount,
        SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount, SEXP miMatrix)
{
    Matrix const priors_matrix(REAL(priorsMatrix), INTEGER(featureCount)[0],
            INTEGER(featureCount)[0]);
    Matrix const* const p_priors_matrix =
            LENGTH(priorsMatrix) == INTEGER(featureCount)[0] * INTEGER(featureCount)[0] ?
                    &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, REAL(priorsWeight)[0], INTEGER(sampleCount)[0],
            INTEGER(featureCount)[0], INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), INTEGER(sampleStratumCount)[0], INTEGER(continuousEstimator)[0],
            INTEGER(outX)[0] != 0, INTEGER(bootstrapCount)[0]);
    MutualInformationMatrix mi_matrix(&data, REAL(miMatrix));

    unsigned int solution_count = 1;
    for (unsigned int i = 0; i < LENGTH(childrenCountPerLevel); ++i)
        solution_count *= INTEGER(childrenCountPerLevel)[i];
    unsigned int const feature_count_per_solution = LENGTH(childrenCountPerLevel);
    unsigned int const chunk_size = solution_count * feature_count_per_solution;

    SEXP result;
    PROTECT(result = Rf_allocVector(VECSXP, 3));

    SET_VECTOR_ELT(result, 0, Rf_allocVector(VECSXP, LENGTH(targetFeatureIndices)));
    SET_VECTOR_ELT(result, 1, Rf_allocVector(VECSXP, LENGTH(targetFeatureIndices)));
    SET_VECTOR_ELT(result, 2, Rf_allocVector(VECSXP, LENGTH(targetFeatureIndices)));

    for (unsigned int i = 0; i < LENGTH(targetFeatureIndices); ++i)
    {
        Filter filter(INTEGER(childrenCountPerLevel), LENGTH(childrenCountPerLevel), &mi_matrix,
                INTEGER(targetFeatureIndices)[i], INTEGER(fixedFeatureCount)[0]);
        filter.build();

        SET_VECTOR_ELT(VECTOR_ELT(result, 0), i, Rf_allocVector(INTSXP, chunk_size));
        SET_VECTOR_ELT(VECTOR_ELT(result, 1), i, Rf_allocVector(REALSXP, INTEGER(featureCount)[0]));
        SET_VECTOR_ELT(VECTOR_ELT(result, 2), i, Rf_allocVector(REALSXP, chunk_size));

        filter.getSolutions(INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), i)));
        filter.getScores(REAL(VECTOR_ELT(VECTOR_ELT(result, 2), i)));

        for (unsigned int k = 0; k < INTEGER(featureCount)[0]; ++k)
            REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i))[k] =
                    std::numeric_limits<double>::quiet_NaN();

        Math::computeCausality(REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i)), &mi_matrix,
                INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), i)), solution_count,
                feature_count_per_solution, INTEGER(featureCount)[0],
                INTEGER(targetFeatureIndices)[i]);
    }
    //PrintValue(result);
    UNPROTECT(1);

    return result;
}

extern "C" SEXP
export_filters_bootstrap(SEXP solutionCount, SEXP solutionLength, SEXP dataMatrix,
        SEXP priorsMatrix, SEXP priorsWeight, SEXP sampleStrata, SEXP sampleWeights,
        SEXP featureTypes, SEXP sampleCount, SEXP featureCount, SEXP sampleStratumCount,
        SEXP targetFeatureIndices, SEXP fixedFeatureCount, SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount,
        SEXP miMatrix)
{
    Matrix const priors_matrix(REAL(priorsMatrix), INTEGER(featureCount)[0],
            INTEGER(featureCount)[0]);
    Matrix const* const p_priors_matrix =
            LENGTH(priorsMatrix) == INTEGER(featureCount)[0] * INTEGER(featureCount)[0] ?
                    &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, REAL(priorsWeight)[0], INTEGER(sampleCount)[0],
            INTEGER(featureCount)[0], INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), INTEGER(sampleStratumCount)[0], INTEGER(continuousEstimator)[0],
            INTEGER(outX)[0] != 0, INTEGER(bootstrapCount)[0]);
    //MutualInformationMatrix mi_matrix(&data, REAL(miMatrix));

    unsigned int solution_count = INTEGER(solutionCount)[0];
    unsigned int const feature_count_per_solution = INTEGER(solutionLength)[0];
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
        SET_VECTOR_ELT(VECTOR_ELT(result, 1), i, Rf_allocVector(REALSXP, INTEGER(featureCount)[0]));
        SET_VECTOR_ELT(VECTOR_ELT(result, 2), i, Rf_allocVector(REALSXP, chunk_size));

        for (unsigned int k = 0; k < INTEGER(featureCount)[0]; ++k)
            REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i))[k] =
                    std::numeric_limits<double>::quiet_NaN();
    }

    for (unsigned int i = 0; i < solution_count; ++i)
    {
        MutualInformationMatrix mi_matrix(&data);

        for (unsigned int j = 0; j < LENGTH(targetFeatureIndices); ++j)
        {
            Filter filter(p_children_count_per_level, feature_count_per_solution, &mi_matrix,
                    INTEGER(targetFeatureIndices)[j], INTEGER(fixedFeatureCount)[0]);
            filter.build();
            filter.getSolutions(
                    INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), j))
                            + (i * feature_count_per_solution));
            //filter.getScores(
             //       REAL(VECTOR_ELT(VECTOR_ELT(result, 2), i)) + (i * feature_count_per_solution));

            /*            Math::computeCausality(REAL(VECTOR_ELT(VECTOR_ELT(result, 1), i)), &mi_matrix,
             INTEGER(VECTOR_ELT(VECTOR_ELT(result, 0), i)) + (i * chunk_size), 1,
             feature_count_per_solution, INTEGER(featureCount)[0],
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
        SEXP sampleStratumCount, SEXP continuousEstimator, SEXP outX, SEXP bootstrapCount,
        SEXP miMatrix)
{
    Matrix const priors_matrix(REAL(priorsMatrix), INTEGER(featureCount)[0],
            INTEGER(featureCount)[0]);
    Matrix const* const p_priors_matrix =
            LENGTH(priorsMatrix) == INTEGER(featureCount)[0] * INTEGER(featureCount)[0] ?
                    &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, REAL(priorsWeight)[0], INTEGER(sampleCount)[0],
            INTEGER(featureCount)[0], INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), INTEGER(sampleStratumCount)[0], INTEGER(continuousEstimator)[0],
            INTEGER(outX)[0] != 0, INTEGER(bootstrapCount)[0]);
    MutualInformationMatrix mi_matrix(&data, REAL(miMatrix));
    mi_matrix.build();

    return R_NilValue;
}

extern "C" SEXP
get_thread_count(SEXP threadCount)
{
#ifdef _OPENMP
    INTEGER(threadCount)[0] = omp_get_max_threads();
#endif

    return R_NilValue;
}

extern "C" SEXP
set_thread_count(SEXP threadCount)
{
#ifdef _OPENMP
    omp_set_num_threads(INTEGER(threadCount)[0]);
#endif

    return R_NilValue;
}
