#include "exports.hpp"

extern "C" SEXP
build_mim(SEXP R_DataMatrix, SEXP R_SampleStrata, SEXP R_SampleWeights, SEXP R_FeatureType,
        SEXP R_RowCount, SEXP R_ColumnCount, SEXP R_SampleStratumCount)
{
    std::vector<float> S_DataMatrix = Rcpp::as < std::vector<float> > (R_DataMatrix);
    std::vector<unsigned int> S_SampleStrata = Rcpp::as < std::vector<unsigned int>
            > (R_SampleStrata);
    std::vector<float> S_SampleWeights = Rcpp::as < std::vector<float> > (R_SampleWeights);
    std::vector<unsigned int> S_FeatureType = Rcpp::as < std::vector<unsigned int>
            > (R_FeatureType);
    unsigned int const row_count = Rcpp::as<unsigned int>(R_RowCount);
    unsigned int const column_count = Rcpp::as<unsigned int>(R_ColumnCount);
    unsigned int const sample_stratum_count = Rcpp::as<unsigned int>(R_SampleStratumCount);
    Matrix data_matrix(&S_DataMatrix[0], row_count, column_count);

    MutualInformationMatrix mi_matrix(&data_matrix, &S_SampleStrata[0], &S_SampleWeights[0],
            &S_FeatureType[0], sample_stratum_count);
    std::vector<float> S_MiMatrix = mi_matrix.getVectorizedData();

    return Rcpp::wrap < std::vector<float> > (S_MiMatrix);
}
