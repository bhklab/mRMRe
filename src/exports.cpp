#include "exports.hpp"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_RowCount, SEXP R_ColumnCount)
{
    std::vector<float> S_DataMatrix = Rcpp::as < std::vector<float> > (R_DataMatrix);
    unsigned int const row_count = Rcpp::as<unsigned int>(R_RowCount);
    unsigned int const column_count = Rcpp::as<unsigned int>(R_ColumnCount);
    Matrix data_matrix(&S_DataMatrix[0], row_count, column_count);
    MutualInformationMatrix mi_matrix(&data_matrix);
    std::vector<float> S_MiMatrix;
    S_MiMatrix.reserve(column_count * column_count);

//#pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < column_count; ++i)
        for (unsigned int j = 0; j < column_count; ++j)
        {
            float const mi = mi_matrix(i, j);
            Rprintf("%d:%d -> %f \n", i, j, mi);
            S_MiMatrix[(i * column_count) + j] = mi;
        }

    return Rcpp::wrap < std::vector<float> > (S_MiMatrix);
}

extern "C" SEXP
mRMR_filter(SEXP R_ChildrenCountPerLevel, SEXP R_FeatureInformationMatrix,
        SEXP R_TargetFeatureIndex)
{
    std::vector<unsigned int> C_ChildrenCountPerLevel = Rcpp::as < std::vector<unsigned int>
            > (R_ChildrenCountPerLevel);

    std::vector<float> C_FeatureInformationVector = Rcpp::as < std::vector<float>
            > (R_FeatureInformationMatrix);
    unsigned int const feature_count = sqrt(C_FeatureInformationVector.size());
    Matrix feature_information_matrix(&C_FeatureInformationVector[0], feature_count, feature_count);

    unsigned int target_feature_index = Rcpp::as<unsigned int>(R_TargetFeatureIndex);

    Tree tree(&C_ChildrenCountPerLevel, &feature_information_matrix, target_feature_index);

    tree.build();

    std::vector<unsigned int> paths;
    tree.getPaths(&paths);

    return Rcpp::wrap < std::vector<unsigned int> > (paths);
}
