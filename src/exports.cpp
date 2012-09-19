#include "exports.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

extern "C" SEXP
rcpp_man(SEXP bigman)
{
    std::vector<float> bigman2 = Rcpp::as < std::vector<float> > (bigman);

    for (unsigned int i = 0; i < bigman2.size(); i++)
        bigman2[i]++;

    return R_NilValue;
}

extern "C" SEXP
r_man(SEXP bigman)
{
    for (unsigned int i = 0; i < LENGTH(bigman); i++)
        REAL(bigman)[i] = REAL(bigman)[i] + 1;

    return R_NilValue;
}

extern "C" SEXP
export_association(SEXP R_SamplesA, SEXP R_SamplesB, SEXP R_SamplesC, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_SampleStratumCount, SEXP R_OutX, SEXP R_BootstrapCount,
        SEXP R_Method)
{
    std::vector<float> S_SamplesA = Rcpp::as < std::vector<float> > (R_SamplesA);
    std::vector<float> S_SamplesB = Rcpp::as < std::vector<float> > (R_SamplesB);
    std::vector<float> S_SamplesC = Rcpp::as < std::vector<float> > (R_SamplesC);
    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
            > (R_SampleStrata);
    std::vector<float> S_SampleWeights = Rcpp::as < std::vector<float> > (R_SampleWeights);
    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
    bool const outX = Rcpp::as<bool>(R_OutX);
    unsigned int const bootstrap_count = Rcpp::as<unsigned int>(R_BootstrapCount);
    std::string S_Method = Rcpp::as < std::string > (R_Method);

    unsigned int const sample_count = S_SamplesA.size();
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[sample_stratum_count];
    unsigned int* const p_sample_count_per_stratum = new unsigned int[sample_stratum_count];
    Math::placeStratificationData(&S_SampleStrata[0], &S_SampleWeights[0],
            p_sample_indices_per_stratum, p_sample_count_per_stratum, sample_stratum_count,
            sample_count);

    bool const is_pearson = S_Method.compare("pearson") == 0;
    bool const is_spearman = S_Method.compare("spearman") == 0;
    bool const is_cramers_v = S_Method.compare("cramer") == 0;
    bool const is_concordance_index = S_Method.compare("cindex") == 0;
    bool const is_concordance_index_with_time = S_Method.compare("cindex_with_time") == 0;

    Rcpp::List result;

    if (is_pearson || is_spearman || is_cramers_v)
    {
        float statistic;

        if (is_pearson)
            statistic = Math::computePearsonCorrelation(&S_SamplesA[0], &S_SamplesB[0],
                    &S_SampleWeights[0], p_sample_indices_per_stratum, p_sample_count_per_stratum,
                    sample_stratum_count, bootstrap_count);
        else if (is_spearman)
            statistic = Math::computeSpearmanCorrelation(&S_SamplesA[0], &S_SamplesB[0],
                    &S_SampleWeights[0], p_sample_indices_per_stratum, p_sample_count_per_stratum,
                    sample_stratum_count, bootstrap_count, sample_count);
        else
            statistic = Math::computeCramersV(&S_SamplesA[0], &S_SamplesB[0], &S_SampleWeights[0],
                    p_sample_indices_per_stratum, p_sample_count_per_stratum, sample_stratum_count,
                    bootstrap_count);

        result = Rcpp::List::create(Rcpp::Named("statistic") = Rcpp::wrap<float>(statistic));
    }
    else if (is_concordance_index)
    {
        float statistic;
        float concordant_weight;
        float discordant_weight;
        float uninformative_weight;
        float relevant_weight;

        if (S_SamplesC.size() == 0)
            statistic = Math::computeConcordanceIndex(&S_SamplesA[0], &S_SamplesB[0],
                    &S_SampleWeights[0], p_sample_indices_per_stratum, p_sample_count_per_stratum,
                    sample_stratum_count, outX, &concordant_weight, &discordant_weight,
                    &uninformative_weight, &relevant_weight);
        else
            statistic = Math::computeConcordanceIndex(&S_SamplesA[0], &S_SamplesB[0],
                    &S_SamplesC[0], &S_SampleWeights[0], p_sample_indices_per_stratum,
                    p_sample_count_per_stratum, sample_stratum_count, outX, &concordant_weight,
                    &discordant_weight, &uninformative_weight, &relevant_weight);

        result = Rcpp::List::create(Rcpp::Named("statistic") = Rcpp::wrap<float>(statistic),
                Rcpp::Named("concordant_weight") = Rcpp::wrap<float>(concordant_weight),
                Rcpp::Named("discordant_weight") = Rcpp::wrap<float>(discordant_weight),
                Rcpp::Named("uninformative_weight") = Rcpp::wrap<float>(uninformative_weight),
                Rcpp::Named("relevant_weight") = Rcpp::wrap<float>(relevant_weight));
    }

    delete[] p_sample_count_per_stratum;
    for (unsigned int i = 0; i < sample_stratum_count; ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;

    return result;
}

extern "C" SEXP
export_filter(SEXP R_ChildrenCountPerLevel, SEXP R_DataMatrix, SEXP R_PriorsMatrix,
        SEXP R_PriorsWeight, SEXP R_SampleStrata, SEXP R_SampleWeights, SEXP R_FeatureTypes,
        SEXP R_SampleCount, SEXP R_FeatureCount, SEXP R_SampleStratumCount,
        SEXP R_TargetFeatureIndex, SEXP R_UsesRanks, SEXP R_OutX, SEXP R_BootstrapCount)
{
    std::vector<unsigned int> S_ChildrenCountPerLevel = Rcpp::as < std::vector<unsigned int>
            > (R_ChildrenCountPerLevel);
    std::vector<float> S_DataMatrix = Rcpp::as < std::vector<float> > (R_DataMatrix);
    std::vector<float> S_PriorsMatrix = Rcpp::as < std::vector<float> > (R_PriorsMatrix);
    float const priors_weight = Rcpp::as<float>(R_PriorsWeight);
    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
            > (R_SampleStrata);
    std::vector<float> S_SampleWeights = Rcpp::as < std::vector<float> > (R_SampleWeights);
    std::vector<unsigned int> S_FeatureTypes = Rcpp::as < std::vector<unsigned int>
            > (R_FeatureTypes);
    unsigned int const sample_count = Rcpp::as<unsigned int>(R_SampleCount);
    unsigned int const feature_count = Rcpp::as<unsigned int>(R_FeatureCount);
    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
    unsigned int const bootstrap_count = Rcpp::as<unsigned int>(R_BootstrapCount);
    bool const uses_ranks = Rcpp::as<bool>(R_UsesRanks);
    bool const outX = Rcpp::as<bool>(R_OutX);
    Matrix const priors_matrix(&S_PriorsMatrix[0], feature_count, feature_count);
    Matrix const* const p_priors_matrix =
            (S_PriorsMatrix.size() == feature_count * feature_count) ? &priors_matrix : 0;
    Data data(&S_DataMatrix[0], p_priors_matrix, priors_weight, sample_count, feature_count,
            &S_SampleStrata[0], &S_SampleWeights[0], &S_FeatureTypes[0], sample_stratum_count,
            uses_ranks, outX, bootstrap_count);
    MutualInformationMatrix const mi_matrix(&data);
    unsigned int const target_feature_index = Rcpp::as<unsigned int>(R_TargetFeatureIndex);
    Filter filter(&S_ChildrenCountPerLevel[0], S_ChildrenCountPerLevel.size(),
            const_cast<MutualInformationMatrix* const >(&mi_matrix), target_feature_index);
    filter.build();
    return Rcpp::List::create(
            Rcpp::Named("solutions") = Rcpp::wrap < std::vector<unsigned int>
                    > (static_cast<std::vector<unsigned int> >(filter)),
            Rcpp::Named("mi_matrix") = Rcpp::wrap < std::vector<float>
                    > (static_cast<std::vector<float> >(mi_matrix)));
}

extern "C" SEXP
export_mim(SEXP R_DataMatrix, SEXP R_PriorsMatrix, SEXP R_PriorsWeight, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_FeatureTypes, SEXP R_SampleCount, SEXP R_FeatureCount,
        SEXP R_SampleStratumCount, SEXP R_UsesRanks, SEXP R_OutX, SEXP R_BootstrapCount)
{
    std::vector<float> S_DataMatrix = Rcpp::as < std::vector<float> > (R_DataMatrix);
    std::vector<float> S_PriorsMatrix = Rcpp::as < std::vector<float> > (R_PriorsMatrix);
    float const priors_weight = Rcpp::as<float>(R_PriorsWeight);
    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
            > (R_SampleStrata);
    std::vector<float> S_SampleWeights = Rcpp::as < std::vector<float> > (R_SampleWeights);
    std::vector<unsigned int> S_FeatureTypes = Rcpp::as < std::vector<unsigned int>
            > (R_FeatureTypes);
    unsigned int const sample_count = Rcpp::as<unsigned int>(R_SampleCount);
    unsigned int const feature_count = Rcpp::as<unsigned int>(R_FeatureCount);
    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
    bool const uses_ranks = Rcpp::as<bool>(R_UsesRanks);
    bool const outX = Rcpp::as<bool>(R_OutX);
    unsigned int const bootstrap_count = Rcpp::as<unsigned int>(R_BootstrapCount);
    Matrix const priors_matrix(&S_PriorsMatrix[0], feature_count, feature_count);
    Matrix const* const p_priors_matrix =
            (S_PriorsMatrix.size() == feature_count * feature_count) ? &priors_matrix : 0;
    Data data(&S_DataMatrix[0], p_priors_matrix, priors_weight, sample_count, feature_count,
            &S_SampleStrata[0], &S_SampleWeights[0], &S_FeatureTypes[0], sample_stratum_count,
            uses_ranks, outX, bootstrap_count);
    MutualInformationMatrix mi_matrix(&data);
    mi_matrix.build();
    return Rcpp::wrap < std::vector<float> > (static_cast<std::vector<float> >(mi_matrix));
}

extern "C" SEXP
set_thread_count(SEXP R_ThreadCount)
{
    unsigned int const thread_count = Rcpp::as<unsigned int>(R_ThreadCount);
    omp_set_num_threads(thread_count);
    return R_NilValue;
}
