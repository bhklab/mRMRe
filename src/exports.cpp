#include "exports.hpp"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_RowCount, SEXP R_ColumnCount)
{
    std::vector<float> S_DataMatrix = Rcpp::as < std::vector<float> > (R_DataMatrix);
    unsigned int const row_count = Rcpp::as<unsigned int>(R_RowCount);
    unsigned int const column_count = Rcpp::as<unsigned int>(R_ColumnCount);
    Matrix data_matrix(&S_DataMatrix[0], row_count, column_count);

    MutualInformationMatrix mi_matrix(&data_matrix);
    mi_matrix.build();
    std::vector<float> S_MiMatrix = mi_matrix.getVectorizedData();

    return Rcpp::wrap < std::vector<float> > (S_MiMatrix);
}

extern "C" SEXP
mRMR_filter_with_mim(SEXP R_FeatureInformationMatrix, SEXP R_ChildrenCountPerLevel,
        SEXP R_TargetFeatureIndex)
{
    std::vector<unsigned int> S_ChildrenCountPerLevel = Rcpp::as < std::vector<unsigned int>
            > (R_ChildrenCountPerLevel);
    std::vector<float> S_FeatureInformationVector = Rcpp::as < std::vector<float>
            > (R_FeatureInformationMatrix);
    unsigned int const target_feature_index = Rcpp::as<unsigned int>(R_TargetFeatureIndex);
    unsigned int const feature_count = sqrt(S_FeatureInformationVector.size());
    Matrix feature_information_matrix(&S_FeatureInformationVector[0], feature_count, feature_count);

    Tree tree(&S_ChildrenCountPerLevel, &feature_information_matrix, target_feature_index);
    tree.build();
    std::vector<unsigned int> S_Paths = tree.getPaths();

    return Rcpp::wrap < std::vector<unsigned int> > (S_Paths);
}

extern "C" SEXP
mRMR_filter_with_data(SEXP R_DataMatrix, SEXP R_RowCount, SEXP R_ColumnCount,
        SEXP R_ChildrenCountPerLevel, SEXP R_TargetFeatureIndex)
{
    std::vector<float> S_DataMatrix = Rcpp::as < std::vector<float> > (R_DataMatrix);
    unsigned int const row_count = Rcpp::as<unsigned int>(R_RowCount);
    unsigned int const column_count = Rcpp::as<unsigned int>(R_ColumnCount);
    Matrix data_matrix(&S_DataMatrix[0], row_count, column_count);

    MutualInformationMatrix mi_matrix(&data_matrix);

    std::vector<unsigned int> S_ChildrenCountPerLevel = Rcpp::as < std::vector<unsigned int>
            > (R_ChildrenCountPerLevel);
    unsigned int const target_feature_index = Rcpp::as<unsigned int>(R_TargetFeatureIndex);

    Tree tree(&S_ChildrenCountPerLevel, &mi_matrix, target_feature_index);
    tree.build();
    std::vector<unsigned int> S_Paths = tree.getPaths();

    return Rcpp::wrap < std::vector<unsigned int> > (S_Paths);
}
