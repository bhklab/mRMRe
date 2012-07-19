#include "exports.hpp"

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
