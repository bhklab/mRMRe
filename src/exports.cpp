#include "exports.h"

extern "C" SEXP
export_association(SEXP samplesA, SEXP samplesB, SEXP samplesC, SEXP sampleStrata,
        SEXP sampleWeights, SEXP sampleStratumCount, SEXP outX, SEXP bootstrapCount, SEXP method,
        SEXP out)
{
    unsigned int const sample_count = LENGTH(samplesA);
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[INTEGER(sampleStratumCount)[0]];
    unsigned int* const p_sample_count_per_stratum =
            new unsigned int[INTEGER(sampleStratumCount)[0]];
    Math::placeStratificationData(INTEGER(sampleStrata), REAL(sampleWeights),
            p_sample_indices_per_stratum, p_sample_count_per_stratum,
            INTEGER(sampleStratumCount)[0], sample_count);

    bool const is_pearson = INTEGER(method)[0] == 0;
    bool const is_spearman = INTEGER(method)[0] == 1;
    bool const is_cramers_v = INTEGER(method)[0] == 2;
    bool const is_concordance_index = INTEGER(method)[0] == 3;
    bool const is_concordance_index_with_time = INTEGER(method)[0] == 4;

    if (is_pearson || is_spearman || is_cramers_v)
    {
        if (is_pearson)
            REAL(out)[0] = Math::computePearsonCorrelation(REAL(samplesA), REAL(samplesB),
                    REAL(sampleWeights), p_sample_indices_per_stratum, p_sample_count_per_stratum,
                    INTEGER(sampleStratumCount)[0], INTEGER(bootstrapCount)[0]);
        else if (is_spearman)
            REAL(out)[0] = Math::computeSpearmanCorrelation(REAL(samplesA), REAL(samplesB),
                    REAL(sampleWeights), p_sample_indices_per_stratum, p_sample_count_per_stratum,
                    INTEGER(sampleStratumCount)[0], INTEGER(bootstrapCount)[0], sample_count);
        else
            REAL(out)[0] = Math::computeCramersV(REAL(samplesA), REAL(samplesB),
                    REAL(sampleWeights), p_sample_indices_per_stratum, p_sample_count_per_stratum,
                    INTEGER(sampleStratumCount)[0], INTEGER(bootstrapCount)[0]);
    }
    else
    {
        if (LENGTH(samplesC) == 0)
            REAL(out)[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                    REAL(sampleWeights), p_sample_indices_per_stratum, p_sample_count_per_stratum,
                    INTEGER(sampleStratumCount)[0], INTEGER(outX)[0] != 0, &(REAL(out)[1]),
                    &(REAL(out)[2]), &(REAL(out)[3]), &(REAL(out)[4]));
        else
            REAL(out)[0] = Math::computeConcordanceIndex(REAL(samplesA), REAL(samplesB),
                    REAL(samplesC), REAL(sampleWeights), p_sample_indices_per_stratum,
                    p_sample_count_per_stratum, INTEGER(sampleStratumCount)[0],
                    INTEGER(outX)[0] != 0, &(REAL(out)[1]), &(REAL(out)[2]), &(REAL(out)[3]),
                    &(REAL(out)[4]));
    }

    delete[] p_sample_count_per_stratum;
    for (unsigned int i = 0; i < INTEGER(sampleStratumCount)[0]; ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;

    return R_NilValue;
}

extern "C" SEXP
export_filter(SEXP childrenCountPerLevel, SEXP dataMatrix, SEXP priorsMatrix, SEXP priorsWeight,
        SEXP sampleStrata, SEXP sampleWeights, SEXP featureTypes, SEXP sampleCount,
        SEXP featureCount, SEXP sampleStratumCount, SEXP targetFeatureIndex, SEXP usesRanks,
        SEXP outX, SEXP bootstrapCount, SEXP miMatrix)
{
    Matrix const priors_matrix(REAL(priorsMatrix), INTEGER(featureCount)[0],
            INTEGER(featureCount)[0]);
    Matrix const* const p_priors_matrix =
            LENGTH(priorsMatrix) == INTEGER(featureCount)[0] * INTEGER(featureCount)[0] ?
                    &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, REAL(priorsWeight)[0], INTEGER(sampleCount)[0],
            INTEGER(featureCount)[0], INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), INTEGER(sampleStratumCount)[0], INTEGER(usesRanks)[0] != 0,
            INTEGER(outX)[0] != 0, INTEGER(bootstrapCount)[0]);
    MutualInformationMatrix mi_matrix(&data, REAL(miMatrix));

    Filter filter(INTEGER(childrenCountPerLevel), LENGTH(childrenCountPerLevel), &mi_matrix,
            INTEGER(targetFeatureIndex)[0]);
    filter.build();
    std::vector<int> filter_solutions;
    filter.getSolutions(&filter_solutions);
    SEXP solutions;
    PROTECT(solutions = allocVector(INTSXP, filter_solutions.size()));
    memcpy(INTEGER(solutions), &filter_solutions[0], sizeof(int) * filter_solutions.size());
    UNPROTECT(1);

    return solutions;
}

extern "C" SEXP
export_mim(SEXP dataMatrix, SEXP priorsMatrix, SEXP priorsWeight, SEXP sampleStrata,
        SEXP sampleWeights, SEXP featureTypes, SEXP sampleCount, SEXP featureCount,
        SEXP sampleStratumCount, SEXP usesRanks, SEXP outX, SEXP bootstrapCount, SEXP miMatrix)
{
    Matrix const priors_matrix(REAL(priorsMatrix), INTEGER(featureCount)[0],
            INTEGER(featureCount)[0]);
    Matrix const* const p_priors_matrix =
            LENGTH(priorsMatrix) == INTEGER(featureCount)[0] * INTEGER(featureCount)[0] ?
                    &priors_matrix : 0;
    Data data(REAL(dataMatrix), p_priors_matrix, REAL(priorsWeight)[0], INTEGER(sampleCount)[0],
            INTEGER(featureCount)[0], INTEGER(sampleStrata), REAL(sampleWeights),
            INTEGER(featureTypes), INTEGER(sampleStratumCount)[0], INTEGER(usesRanks)[0] != 0,
            INTEGER(outX)[0] != 0, INTEGER(bootstrapCount)[0]);
    MutualInformationMatrix mi_matrix(&data, REAL(miMatrix));
    mi_matrix.build();

    return R_NilValue;
}

extern "C" SEXP
get_thread_count(SEXP threadCount)
{
    INTEGER(threadCount)[0] = omp_get_max_threads();

    return R_NilValue;
}

extern "C" SEXP
set_thread_count(SEXP threadCount)
{
    omp_set_num_threads(INTEGER(threadCount)[0]);

    return R_NilValue;
}
