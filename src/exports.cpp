#include "exports.hpp"

extern "C" SEXP
mRMR_filter(SEXP R_ChildrenCountPerLevel, SEXP R_FeatureInformationMatrix,
        SEXP R_TargetFeatureIndex)
{
    std::vector<unsigned int> children_count_per_level = Rcpp::as < std::vector<unsigned int>
            > (R_ChildrenCountPerLevel);
    std::vector<float> feature_information_vector = Rcpp::as < std::vector<float>
            > (R_FeatureInformationMatrix);
    std::vector<unsigned int> paths;

    unsigned int const feature_count = sqrt(feature_information_vector.size());
    Matrix feature_information_matrix(&feature_information_vector[0], feature_count, feature_count);

    unsigned int target_feature_index = Rcpp::as<unsigned int>(R_TargetFeatureIndex);

    Tree tree(&(children_count_per_level[0]), children_count_per_level.size(),
            &feature_information_matrix,
            target_feature_index);
    tree.build();
    tree.getPaths(&paths);

    return Rcpp::wrap < std::vector<unsigned int> > (paths);
}
