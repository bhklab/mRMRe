#include "exports.hpp"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_SampleStrata, SEXP R_SampleWeights, SEXP R_FeatureTypes,
        SEXP R_SampleCount, SEXP R_FeatureCount, SEXP R_SampleStratumCount)
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
    Data data(&S_DataMatrix[0], sample_count, feature_count, &S_SampleStrata[0],
            &S_SampleWeights[0], &S_FeatureTypes[0], sample_stratum_count);
    MutualInformationMatrix mi_matrix(&data);
    mi_matrix.build();
    std::vector<float> S_MiMatrix = mi_matrix.getVectorizedData();
    return Rcpp::wrap < std::vector<float> > (S_MiMatrix);
}

extern "C" SEXP
filter_mRMR_with_data(SEXP R_ChildrenCountPerLevel, SEXP R_DataMatrix, SEXP R_SampleStrata,
        SEXP R_SampleWeights, SEXP R_FeatureTypes, SEXP R_SampleCount, SEXP R_FeatureCount,
        SEXP R_SampleStratumCount, SEXP R_TargetFeatureIndex)
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
    Data data(&S_DataMatrix[0], sample_count, feature_count, &S_SampleStrata[0],
            &S_SampleWeights[0], &S_FeatureTypes[0], sample_stratum_count);
    MutualInformationMatrix mi_matrix(&data);
    unsigned int const target_feature_index = Rcpp::as<unsigned int>(R_TargetFeatureIndex);
    Tree mRMR_tree(&S_ChildrenCountPerLevel[0], S_ChildrenCountPerLevel.size(), &mi_matrix,
            target_feature_index);
    mRMR_tree.build();
    std::vector<unsigned int> S_Paths = mRMR_tree.getPaths();
    std::vector<float> S_Scores = mRMR_tree.getScores();
    return Rcpp::wrap < std::vector<unsigned int> > (S_Paths); // TODO: Add Scores to return
}

extern "C" SEXP
filter_mRMR_with_mim(SEXP R_ChildrenCountPerLevel, SEXP R_MiMatrix, SEXP R_FeatureCount,
        SEXP R_TargetFeatureIndex)
{
    std::vector<unsigned int> S_ChildrenCountPerLevel = Rcpp::as < std::vector<unsigned int>
            > (R_ChildrenCountPerLevel);
    std::vector<float> S_MiMatrix = Rcpp::as < std::vector<float> > (R_MiMatrix);
    unsigned int const feature_count = Rcpp::as<unsigned int>(R_FeatureCount);
    Matrix mi_matrix(&S_MiMatrix[0], feature_count, feature_count);
    unsigned int const target_feature_index = Rcpp::as<unsigned int>(R_TargetFeatureIndex);
    Tree mRMR_tree(&S_ChildrenCountPerLevel[0], S_ChildrenCountPerLevel.size(), &mi_matrix,
            target_feature_index);
    mRMR_tree.build();
    std::vector<unsigned int> S_Paths = mRMR_tree.getPaths();
    std::vector<float> S_Scores = mRMR_tree.getScores();
    return Rcpp::wrap < std::vector<unsigned int> > (S_Paths); // TODO: Add Scores to return
}
