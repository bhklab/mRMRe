#include "exports.hpp"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_Strata, SEXP R_Weights, SEXP R_FeatureType, SEXP R_RowCount,
        SEXP R_ColumnCount)
{
    std::vector<float> S_DataMatrix = Rcpp::as < std::vector<float> > (R_DataMatrix);
    std::vector<unsigned int> S_Strata = Rcpp::as < std::vector<unsigned int> > (R_Strata);
    std::vector<unsigned int> S_FeatureType = Rcpp::as < std::vector<unsigned int>
            > (R_FeatureType);
    std::vector<float> S_Weights = Rcpp::as < std::vector<float> > (R_Weights);

    unsigned int const row_count = Rcpp::as<unsigned int>(R_RowCount);
    unsigned int const column_count = Rcpp::as<unsigned int>(R_ColumnCount);
    Matrix data_matrix(&S_DataMatrix[0], row_count, column_count);

    MutualInformationMatrix mi_matrix(&data_matrix, &S_Strata[0], &S_Weights[0], &S_FeatureType[0]);
    std::vector<float> S_MiMatrix = mi_matrix.getVectorizedData();

    return Rcpp::wrap < std::vector<float> > (S_MiMatrix);
}
