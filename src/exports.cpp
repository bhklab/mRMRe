#include "exports.h"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_SampleStrata, SEXP R_SampleWeights, SEXP R_FeatureTypes,
        SEXP R_SampleCount, SEXP R_FeatureCount, SEXP R_SampleStratumCount, SEXP R_UsesRanks,
        SEXP R_OutX, SEXP R_BootstrapCount)
{
    std::vector<float> S_DataMatrix = Rcpp::as < std::vector<float> > (R_DataMatrix);
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
    Data data(&S_DataMatrix[0], sample_count, feature_count, &S_SampleStrata[0],
            &S_SampleWeights[0], &S_FeatureTypes[0], sample_stratum_count, uses_ranks, outX,
            bootstrap_count);
    MutualInformationMatrix mi_matrix(&data);
    mi_matrix.build();
    std::vector<float> S_MiMatrix = mi_matrix.getVectorizedData();
    return Rcpp::wrap < std::vector<float> > (S_MiMatrix);
}

extern "C" SEXP
build_mRMR_tree(SEXP R_ChildrenCountPerLevel, SEXP R_DataMatrix, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_FeatureTypes, SEXP R_SampleCount, SEXP R_FeatureCount,
        SEXP R_SampleStratumCount, SEXP R_TargetFeatureIndex, SEXP R_UsesRanks, SEXP R_OutX,
        SEXP R_BootstrapCount)
{
    std::vector<unsigned int> S_ChildrenCountPerLevel = Rcpp::as < std::vector<unsigned int>
            > (R_ChildrenCountPerLevel);
    std::vector<float> S_DataMatrix = Rcpp::as < std::vector<float> > (R_DataMatrix);
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
    Data data(&S_DataMatrix[0], sample_count, feature_count, &S_SampleStrata[0],
            &S_SampleWeights[0], &S_FeatureTypes[0], sample_stratum_count, uses_ranks, outX,
            bootstrap_count);
    MutualInformationMatrix mi_matrix(&data);
    unsigned int const target_feature_index = Rcpp::as<unsigned int>(R_TargetFeatureIndex);
    Tree mRMR_tree(&S_ChildrenCountPerLevel[0], S_ChildrenCountPerLevel.size(), &mi_matrix,
            target_feature_index);
    mRMR_tree.build();
    std::vector<unsigned int> S_Paths = mRMR_tree.getPaths();
    std::vector<float> S_Scores = mRMR_tree.getScores();
    std::vector<float> S_MiMatrix;
    S_MiMatrix.resize(feature_count * feature_count);
#pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < feature_count; ++i)
        for (unsigned int j = 0; j < feature_count; ++j)
            S_MiMatrix[(i * feature_count) + j] = mi_matrix.SymmetricMatrix::at(i, j);
    return Rcpp::List::create(
            Rcpp::Named("paths") = Rcpp::wrap < std::vector<unsigned int> > (S_Paths),
            Rcpp::Named("scores") = Rcpp::wrap < std::vector<float> > (S_Scores),
            Rcpp::Named("mim") = Rcpp::wrap < std::vector<float> > (S_MiMatrix));
}

extern "C" SEXP
compute_concordance_index(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights,
        SEXP R_SampleStrata, SEXP R_SampleStratumCount, SEXP R_OutX)
{
    std::vector<float> S_SamplesX = Rcpp::as < std::vector<float> > (R_SamplesX);
    std::vector<float> S_SamplesY = Rcpp::as < std::vector<float> > (R_SamplesY);
    std::vector<float> S_SampleWeights = Rcpp::as < std::vector<float> > (R_SampleWeights);
    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
            > (R_SampleStrata);
    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
    bool const outX = Rcpp::as<bool>(R_OutX);
    unsigned int const sample_count = S_SamplesX.size();
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[sample_stratum_count];
    float* const p_total_weight_per_stratum = new float[sample_stratum_count];
    unsigned int* const p_sample_count_per_stratum = new unsigned int[sample_stratum_count];
    Math::placeStratificationData(&S_SampleStrata[0], &S_SampleWeights[0],
            p_sample_indices_per_stratum, p_total_weight_per_stratum, p_sample_count_per_stratum,
            sample_stratum_count, sample_count);
    float const r = Math::computeConcordanceIndex(&S_SamplesX[0], &S_SamplesY[0],
            &S_SampleWeights[0], p_sample_indices_per_stratum, p_sample_count_per_stratum,
            sample_stratum_count, outX);
    delete[] p_sample_count_per_stratum;
    delete[] p_total_weight_per_stratum;
    for (unsigned int i = 0; i < sample_stratum_count; ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;
    return Rcpp::wrap<float>(r);
}

extern "C" SEXP
compute_concordance_index_with_time(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_Time,
        SEXP R_SampleWeights, SEXP R_SampleStrata, SEXP R_SampleStratumCount, SEXP R_OutX)
{
    std::vector<float> S_SamplesX = Rcpp::as < std::vector<float> > (R_SamplesX);
    std::vector<float> S_SamplesY = Rcpp::as < std::vector<float> > (R_SamplesY);
    std::vector<float> S_Time = Rcpp::as < std::vector<float> > (R_Time);
    std::vector<float> S_SampleWeights = Rcpp::as < std::vector<float> > (R_SampleWeights);
    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
            > (R_SampleStrata);
    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
    bool const outX = Rcpp::as<bool>(R_OutX);
    unsigned int const sample_count = S_SamplesX.size();
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[sample_stratum_count];
    float* const p_total_weight_per_stratum = new float[sample_stratum_count];
    unsigned int* const p_sample_count_per_stratum = new unsigned int[sample_stratum_count];
    Math::placeStratificationData(&S_SampleStrata[0], &S_SampleWeights[0],
            p_sample_indices_per_stratum, p_total_weight_per_stratum, p_sample_count_per_stratum,
            sample_stratum_count, sample_count);
    float const r = Math::computeConcordanceIndexWithTime(&S_SamplesX[0], &S_SamplesY[0],
            &S_Time[0], &S_SampleWeights[0], p_sample_indices_per_stratum,
            p_sample_count_per_stratum, sample_stratum_count, outX);
    delete[] p_sample_count_per_stratum;
    delete[] p_total_weight_per_stratum;
    for (unsigned int i = 0; i < sample_stratum_count; ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;
    return Rcpp::wrap<float>(r);
}

extern "C" SEXP
compute_cramers_v(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights, SEXP R_SampleStrata,
        SEXP R_SampleStratumCount, SEXP R_BootstrapCount)
{
    std::vector<float> S_SamplesX = Rcpp::as < std::vector<float> > (R_SamplesX);
    std::vector<float> S_SamplesY = Rcpp::as < std::vector<float> > (R_SamplesY);
    std::vector<float> S_SampleWeights = Rcpp::as < std::vector<float> > (R_SampleWeights);
    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
            > (R_SampleStrata);
    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
    unsigned int const sample_count = S_SamplesX.size();
    unsigned int const bootstrap_count = Rcpp::as<unsigned int>(R_BootstrapCount);
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[sample_stratum_count];
    float* const p_total_weight_per_stratum = new float[sample_stratum_count];
    unsigned int* const p_sample_count_per_stratum = new unsigned int[sample_stratum_count];
    Math::placeStratificationData(&S_SampleStrata[0], &S_SampleWeights[0],
            p_sample_indices_per_stratum, p_total_weight_per_stratum, p_sample_count_per_stratum,
            sample_stratum_count, sample_count);
    float const r = Math::computeCramersV(&S_SamplesX[0], &S_SamplesY[0], &S_SampleWeights[0],
            p_sample_indices_per_stratum, p_total_weight_per_stratum, p_sample_count_per_stratum,
            sample_stratum_count, bootstrap_count);
    delete[] p_sample_count_per_stratum;
    delete[] p_total_weight_per_stratum;
    for (unsigned int i = 0; i < sample_stratum_count; ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;
    return Rcpp::wrap<float>(r);
}

extern "C" SEXP
compute_pearson_correlation(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights,
        SEXP R_SampleStrata, SEXP R_SampleStratumCount, SEXP R_BootstrapCount)
{
    std::vector<float> S_SamplesX = Rcpp::as < std::vector<float> > (R_SamplesX);
    std::vector<float> S_SamplesY = Rcpp::as < std::vector<float> > (R_SamplesY);
    std::vector<float> S_SampleWeights = Rcpp::as < std::vector<float> > (R_SampleWeights);
    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
            > (R_SampleStrata);
    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
    unsigned int const sample_count = S_SamplesX.size();
    unsigned int const bootstrap_count = Rcpp::as<unsigned int>(R_BootstrapCount);
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[sample_stratum_count];
    float* const p_total_weight_per_stratum = new float[sample_stratum_count];
    unsigned int* const p_sample_count_per_stratum = new unsigned int[sample_stratum_count];
    Math::placeStratificationData(&S_SampleStrata[0], &S_SampleWeights[0],
            p_sample_indices_per_stratum, p_total_weight_per_stratum, p_sample_count_per_stratum,
            sample_stratum_count, sample_count);
    float const r = Math::computePearsonCorrelation(&S_SamplesX[0], &S_SamplesY[0],
            &S_SampleWeights[0], p_sample_indices_per_stratum, p_total_weight_per_stratum,
            p_sample_count_per_stratum, sample_stratum_count, bootstrap_count);
    delete[] p_sample_count_per_stratum;
    delete[] p_total_weight_per_stratum;
    for (unsigned int i = 0; i < sample_stratum_count; ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;
    return Rcpp::wrap<float>(r);
}

extern "C" SEXP
compute_spearman_correlation(SEXP R_SamplesX, SEXP R_SamplesY, SEXP R_SampleWeights,
        SEXP R_SampleStrata, SEXP R_SampleStratumCount, SEXP R_BootstrapCount)
{
    std::vector<float> S_SamplesX = Rcpp::as < std::vector<float> > (R_SamplesX);
    std::vector<float> S_SamplesY = Rcpp::as < std::vector<float> > (R_SamplesY);
    std::vector<float> S_SampleWeights = Rcpp::as < std::vector<float> > (R_SampleWeights);
    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
            > (R_SampleStrata);
    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
    unsigned int const sample_count = S_SamplesX.size();
    unsigned int const bootstrap_count = Rcpp::as<unsigned int>(R_BootstrapCount);
    unsigned int** p_sample_indices_per_stratum = new unsigned int*[sample_stratum_count];
    float* const p_total_weight_per_stratum = new float[sample_stratum_count];
    unsigned int* const p_sample_count_per_stratum = new unsigned int[sample_stratum_count];
    Math::placeStratificationData(&S_SampleStrata[0], &S_SampleWeights[0],
            p_sample_indices_per_stratum, p_total_weight_per_stratum, p_sample_count_per_stratum,
            sample_stratum_count, sample_count);
    float* const p_ordered_samples_x = new float[sample_count];
    float* const p_ordered_samples_y = new float[sample_count];
    Math::placeOrders(&S_SamplesX[0], p_ordered_samples_x, p_sample_indices_per_stratum,
            p_sample_count_per_stratum, sample_stratum_count);
    Math::placeOrders(&S_SamplesY[0], p_ordered_samples_y, p_sample_indices_per_stratum,
            p_sample_count_per_stratum, sample_stratum_count);
    float* const p_ranked_samples_x = new float[sample_count];
    float* const p_ranked_samples_y = new float[sample_count];
    Math::placeRanksFromOrders(&S_SamplesX[0], &S_SamplesY[0], p_ordered_samples_x,
            p_ordered_samples_y, p_ranked_samples_x, p_ranked_samples_y,
            p_sample_indices_per_stratum, p_sample_count_per_stratum, sample_stratum_count);
    float const r = Math::computePearsonCorrelation(p_ranked_samples_x, p_ranked_samples_y,
            &S_SampleWeights[0], p_sample_indices_per_stratum, p_total_weight_per_stratum,
            p_sample_count_per_stratum, sample_stratum_count, bootstrap_count);
    delete[] p_ordered_samples_x;
    delete[] p_ordered_samples_y;
    delete[] p_ranked_samples_x;
    delete[] p_ranked_samples_y;
    delete[] p_sample_count_per_stratum;
    delete[] p_total_weight_per_stratum;
    for (unsigned int i = 0; i < sample_stratum_count; ++i)
        delete[] p_sample_indices_per_stratum[i];
    delete[] p_sample_indices_per_stratum;
    return Rcpp::wrap<float>(r);
}

extern "C" SEXP
set_thread_count(SEXP R_ThreadCount)
{
    unsigned int const thread_count = Rcpp::as<unsigned int>(R_ThreadCount);
    omp_set_num_threads(thread_count);
    return R_NilValue;
}
