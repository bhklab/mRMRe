#include "exports.h"

extern "C" SEXP
export_association(SEXP R_SamplesA, SEXP R_SamplesB, SEXP R_SamplesC, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_SampleStratumCount, SEXP R_OutX, SEXP R_BootstrapCount,
        SEXP R_Method)
{
//    std::vector<double> S_SamplesA = Rcpp::as < std::vector<double> > (R_SamplesA);
//    std::vector<double> S_SamplesB = Rcpp::as < std::vector<double> > (R_SamplesB);
//    std::vector<double> S_SamplesC = Rcpp::as < std::vector<double> > (R_SamplesC);
//    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
//            > (R_SampleStrata);
//    std::vector<double> S_SampleWeights = Rcpp::as < std::vector<double> > (R_SampleWeights);
//    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
//    bool const outX = Rcpp::as<bool>(R_OutX);
//    unsigned int const bootstrap_count = Rcpp::as<unsigned int>(R_BootstrapCount);
//    std::string S_Method = Rcpp::as < std::string > (R_Method);
//
//    unsigned int const sample_count = S_SamplesA.size();
//    unsigned int** p_sample_indices_per_stratum = new unsigned int*[sample_stratum_count];
//    unsigned int* const p_sample_count_per_stratum = new unsigned int[sample_stratum_count];
//    Math::placeStratificationData(&S_SampleStrata[0], &S_SampleWeights[0],
//            p_sample_indices_per_stratum, p_sample_count_per_stratum, sample_stratum_count,
//            sample_count);
//
//    bool const is_pearson = S_Method.compare("pearson") == 0;
//    bool const is_spearman = S_Method.compare("spearman") == 0;
//    bool const is_cramers_v = S_Method.compare("cramer") == 0;
//    bool const is_concordance_index = S_Method.compare("cindex") == 0;
//    bool const is_concordance_index_with_time = S_Method.compare("cindex_with_time") == 0;
//
//    Rcpp::List result;
//
//    if (is_pearson || is_spearman || is_cramers_v)
//    {
//        double statistic;
//
//        if (is_pearson)
//            statistic = Math::computePearsonCorrelation(&S_SamplesA[0], &S_SamplesB[0],
//                    &S_SampleWeights[0], p_sample_indices_per_stratum, p_sample_count_per_stratum,
//                    sample_stratum_count, bootstrap_count);
//        else if (is_spearman)
//            statistic = Math::computeSpearmanCorrelation(&S_SamplesA[0], &S_SamplesB[0],
//                    &S_SampleWeights[0], p_sample_indices_per_stratum, p_sample_count_per_stratum,
//                    sample_stratum_count, bootstrap_count, sample_count);
//        else
//            statistic = Math::computeCramersV(&S_SamplesA[0], &S_SamplesB[0], &S_SampleWeights[0],
//                    p_sample_indices_per_stratum, p_sample_count_per_stratum, sample_stratum_count,
//                    bootstrap_count);
//
//        result = Rcpp::List::create(Rcpp::Named("statistic") = Rcpp::wrap<double>(statistic));
//    }
//    else if (is_concordance_index)
//    {
//        double statistic;
//        double concordant_weight;
//        double discordant_weight;
//        double uninformative_weight;
//        double relevant_weight;
//
//        if (S_SamplesC.size() == 0)
//            statistic = Math::computeConcordanceIndex(&S_SamplesA[0], &S_SamplesB[0],
//                    &S_SampleWeights[0], p_sample_indices_per_stratum, p_sample_count_per_stratum,
//                    sample_stratum_count, outX, &concordant_weight, &discordant_weight,
//                    &uninformative_weight, &relevant_weight);
//        else
//            statistic = Math::computeConcordanceIndex(&S_SamplesA[0], &S_SamplesB[0],
//                    &S_SamplesC[0], &S_SampleWeights[0], p_sample_indices_per_stratum,
//                    p_sample_count_per_stratum, sample_stratum_count, outX, &concordant_weight,
//                    &discordant_weight, &uninformative_weight, &relevant_weight);
//
//        result = Rcpp::List::create(Rcpp::Named("statistic") = Rcpp::wrap<double>(statistic),
//                Rcpp::Named("concordant_weight") = Rcpp::wrap<double>(concordant_weight),
//                Rcpp::Named("discordant_weight") = Rcpp::wrap<double>(discordant_weight),
//                Rcpp::Named("uninformative_weight") = Rcpp::wrap<double>(uninformative_weight),
//                Rcpp::Named("relevant_weight") = Rcpp::wrap<double>(relevant_weight));
//    }
//
//    delete[] p_sample_count_per_stratum;
//    for (unsigned int i = 0; i < sample_stratum_count; ++i)
//        delete[] p_sample_indices_per_stratum[i];
//    delete[] p_sample_indices_per_stratum;
//
//    return result;

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
    memcpy(&filter_solutions[0], INTEGER(solutions), sizeof(int) * filter_solutions.size());
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
    INTEGER(threadCount)[0] = omp_get_num_threads();

    return R_NilValue;
}

extern "C" SEXP
set_thread_count(SEXP threadCount)
{
    omp_set_num_threads(INTEGER(threadCount)[0]);

    return R_NilValue;
}
