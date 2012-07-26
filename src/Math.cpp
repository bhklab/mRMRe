#include "Math.hpp"

//float const
//computeConcordanceIndex(unsigned int const discreteFeatureIndex,
//        unsigned int const continuousFeatureIndex, int const timeFeatureIndex,
//        Matrix const* const pDataMatrix, float const* const pSampleWeights,
//        unsigned int const* const pSampleStrata, bool const outX)
//{
//    float concordant_weight = 0.;
//    float discordant_weight = 0.;
//    float uninformative_weight = 0.;
//    float relevant_weight = 0.;
//
//    unsigned int const sample_count = pDataMatrix->getRowCount();
//    for (unsigned int i = 0; i < sample_count; ++i)
//    {
//        if ((*pDataMatrix)(i, discreteFeatureIndex) != (*pDataMatrix)(i, discreteFeatureIndex))
//            continue;
//        if (timeFeatureIndex >= 0
//                && (*pDataMatrix)(i, timeFeatureIndex) != (*pDataMatrix)(i, timeFeatureIndex))
//            continue;
//        for (unsigned int j = 0; j < sample_count; ++j)
//        {
//            if ((*pDataMatrix)(j, discreteFeatureIndex) != (*pDataMatrix)(j, discreteFeatureIndex))
//                continue;
//            if (timeFeatureIndex >= 0
//                    && (*pDataMatrix)(j, timeFeatureIndex) != (*pDataMatrix)(j, timeFeatureIndex))
//                continue;
//            if (pSampleStrata[i] != pSampleStrata[j])
//                continue;
//
//            float pair_weight = pSampleWeights[i] * pSampleWeights[j];
//            if (isComparablePair(i, j, timeFeatureIndex, discreteFeatureIndex, pDataMatrix))
//            {
//                relevant_weight += pair_weight;
//                if ((*pDataMatrix)(i, continuousFeatureIndex)
//                        > (*pDataMatrix)(j, continuousFeatureIndex))
//                    concordant_weight += pair_weight;
//                else if ((*pDataMatrix)(i, continuousFeatureIndex)
//                        < (*pDataMatrix)(j, continuousFeatureIndex))
//                    discordant_weight += pair_weight;
//                else if (outX)
//                    uninformative_weight += pair_weight;
//                else
//                    discordant_weight += pair_weight;
//            }
//
//            else if (isComparablePair(j, i, timeFeatureIndex, discreteFeatureIndex, pDataMatrix))
//            {
//                relevant_weight += pair_weight;
//                if ((*pDataMatrix)(i, continuousFeatureIndex)
//                        < (*pDataMatrix)(j, continuousFeatureIndex))
//                    concordant_weight += pair_weight;
//                else if ((*pDataMatrix)(i, continuousFeatureIndex)
//                        > (*pDataMatrix)(j, continuousFeatureIndex))
//                    discordant_weight += pair_weight;
//                else if (outX)
//                    uninformative_weight += pair_weight;
//                else
//                    discordant_weight += pair_weight;
//            }
//
//        }
//    }
//    return concordant_weight / relevant_weight;
//}
//
//bool const
//isComparablePair(unsigned int const i, unsigned int const j, int const timeFeatureIndex,
//        unsigned int const discreteFeatureIndex, Matrix const* const pDataMatrix)
//{
//    return (timeFeatureIndex >= 0
//            && ((*pDataMatrix)(i, timeFeatureIndex) < (*pDataMatrix)(j, timeFeatureIndex)
//                    && (*pDataMatrix)(i, discreteFeatureIndex) == 1))
//            || (timeFeatureIndex < 0
//                    && (*pDataMatrix)(i, discreteFeatureIndex)
//                            > (*pDataMatrix)(j, discreteFeatureIndex));
//}
//
//float const
//computeCramersV(unsigned int const featureIndex1, unsigned int const featureIndex2,
//        Matrix const* const pDataMatrix, float const* const pSampleWeights)
//{
//
//    unsigned int const sample_count = pDataMatrix->getRowCount();
//    unsigned int x_class_count = 0;
//    unsigned int y_class_count = 0;
//
//    for (unsigned int i = 0; i < sample_count; ++i)
//    {
//        if (x_class_count <= (*pDataMatrix)(i, featureIndex1))
//            x_class_count = (*pDataMatrix)(i, featureIndex1) + 1;
//        if (y_class_count <= (*pDataMatrix)(i, featureIndex2))
//            y_class_count = (*pDataMatrix)(i, featureIndex2) + 1;
//    }
//
//    Matrix* const contingency_table = new Matrix(x_class_count + 1, y_class_count + 1);
//    for (unsigned int i = 0; i < x_class_count; ++i)
//        for (unsigned int j = 0; j < y_class_count; ++j)
//            (*contingency_table)(i, j) = 0;
//
//    for (unsigned int i = 0; i < sample_count; ++i)
//    {
//        float const sample_weight = pSampleWeights[i];
//        (*contingency_table)((*pDataMatrix)(i, featureIndex1), (*pDataMatrix)(i, featureIndex2)) +=
//                sample_weight;
//        (*contingency_table)(x_class_count, (*pDataMatrix)(i, featureIndex2)) += sample_weight;
//        (*contingency_table)((*pDataMatrix)(i, featureIndex1), y_class_count) += sample_weight;
//        (*contingency_table)(x_class_count, y_class_count) += sample_weight;
//    }
//    float chi_square = 0.;
//
//    for (unsigned int i = 0; i < x_class_count; ++i)
//        for (unsigned int j = 0; j < y_class_count; ++j)
//        {
//            float expected_value = (*contingency_table)(i, y_class_count)
//                    * (*contingency_table)(x_class_count, j)
//                    / (*contingency_table)(x_class_count, y_class_count);
//
//            chi_square += std::pow(((*contingency_table)(i, j) - expected_value), 2)
//                    / expected_value;
//        }
//
//    unsigned int min_classes = (x_class_count < y_class_count) ? x_class_count : y_class_count;
//    return std::sqrt(
//            chi_square / ((*contingency_table)(x_class_count, y_class_count) * (min_classes - 1)));
//}

Math::IndirectComparator::IndirectComparator(float const* const pSamples,
        unsigned int const* const pSampleIndices) :
        mpSamples(pSamples), mpSampleIndices(pSampleIndices)
{

}

bool const
Math::IndirectComparator::operator()(unsigned int const i, unsigned int const j) const
{
    return mpSamples[mpSampleIndices[i]] < mpSamples[mpSampleIndices[j]];
}

/* static */float const
Math::computeMi(float const r)
{
    return -0.5 * std::log(1 - (r * r));
}

/* static */float const
Math::computePearsonCorrelation(float const* const pSamplesX, float const* const pSamplesY,
        float const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        float const* const pTotalWeightPerStratum, unsigned int const* const pSampleCountPerStratum,
        unsigned int const sampleStratumCount)
{
    float r = 0.;
    float total_weight = 0.;

    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        r += pTotalWeightPerStratum[i]
                * computePearsonCorrelation(pSamplesX, pSamplesY, pSampleWeights,
                        pSampleIndicesPerStratum[i], pSampleCountPerStratum[i]);
        total_weight += pTotalWeightPerStratum[i];
    }

    r /= total_weight;

    return r;
}

/* static */float const
Math::computePearsonCorrelation(float const* const pSamplesX, float const* const pSamplesY,
        float const* const pSampleWeights, unsigned int const* const pSampleIndices,
        unsigned int const sampleCount)
{
    float sum_of_x = 0.;
    float sum_of_x_x = 0.;
    float sum_of_y = 0.;
    float sum_of_y_y = 0.;
    float sum_of_x_y = 0.;
    float sum_of_weights = 0.;

    for (unsigned int i = 0; i < sampleCount; ++i)
    {
        float const my_weight = pSampleWeights[pSampleIndices[i]];
        float const my_x = pSamplesX[pSampleIndices[i]];
        sum_of_x += my_x * my_weight;
        sum_of_x_x += my_x * my_x * my_weight;
        float const my_y = pSamplesY[pSampleIndices[i]];
        sum_of_y += my_y * my_weight;
        sum_of_y_y += my_y * my_y * my_weight;
        sum_of_x_y += my_x * my_y * my_weight;
        sum_of_weights += my_weight;
    }

    float const r = (sum_of_x_y - ((sum_of_x * sum_of_y) / sum_of_weights))
            / std::sqrt(
                    (sum_of_x_x - ((sum_of_x * sum_of_x) / sum_of_weights))
                            * (sum_of_y_y - ((sum_of_y * sum_of_y) / sum_of_weights)));

    return r;
}

/* static */void const
Math::placeRanksByFeatureIndex(float const* const pSamples, float* const pRanks,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount)
{
    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        unsigned int const* const p_sample_indices = pSampleIndicesPerStratum[i];
        unsigned int const sample_count = pSampleCountPerStratum[i];

        unsigned int p_order[sample_count];

        for (unsigned int j = 0; j < sample_count; ++j)
            p_order[j] = j;

        std::sort(p_order, p_order + sample_count,
                Math::IndirectComparator(pSamples, p_sample_indices));

        for (unsigned int j = 0; j < sample_count; ++j)
            pRanks[p_sample_indices[p_order[j]]] = j;
    }
}
