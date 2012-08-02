#include "Math.hpp"

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
Math::computeConcordanceIndex(float const* const pDiscreteSamples,
        float const* const pContinuousSamples, float const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        bool const outX)
{
    float concordant_weight = 0.;
    float discordant_weight = 0.;
    float uninformative_weight = 0.;
    float relevant_weight = 0.;

    for (unsigned int stratum = 0; stratum < sampleStratumCount; ++stratum)
    {
        for (unsigned int a = 0; a < pSampleCountPerStratum[stratum]; ++a)
        {
            unsigned int const i = pSampleIndicesPerStratum[stratum][a];

            if (pDiscreteSamples[i] != pDiscreteSamples[i])
                continue;

            for (unsigned int b = 0; b < pSampleCountPerStratum[stratum]; ++b)
            {
                unsigned int const j = pSampleIndicesPerStratum[stratum][b];

                if (pDiscreteSamples[j] != pDiscreteSamples[j])
                    continue;

                float pair_weight = pSampleWeights[i] * pSampleWeights[j];

                if (pDiscreteSamples[i] > pDiscreteSamples[j])
                {
                    relevant_weight += pair_weight;

                    if (pContinuousSamples[i] > pContinuousSamples[j])
                        concordant_weight += pair_weight;
                    else if (pContinuousSamples[i] < pContinuousSamples[j])
                        discordant_weight += pair_weight;
                    else if (outX)
                        uninformative_weight += pair_weight;
                    else
                        discordant_weight += pair_weight;
                }
                else if (pDiscreteSamples[i] < pDiscreteSamples[j])
                {
                    relevant_weight += pair_weight;

                    if (pContinuousSamples[i] < pContinuousSamples[j])
                        concordant_weight += pair_weight;
                    else if (pContinuousSamples[i] > pContinuousSamples[j])
                        discordant_weight += pair_weight;
                    else if (outX)
                        uninformative_weight += pair_weight;
                    else
                        discordant_weight += pair_weight;
                }

            }
        }
    }

    return concordant_weight / relevant_weight;
}

/*static*/float const
Math::computeConcordanceIndexWithTime(float const* const pDiscreteSamples,
        float const* const pContinuousSamples, float const* const pTimeSamples,
        float const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        bool const outX)
{
    float concordant_weight = 0.;
    float discordant_weight = 0.;
    float uninformative_weight = 0.;
    float relevant_weight = 0.;

    for (unsigned int stratum = 0; stratum < sampleStratumCount; ++stratum)
    {
        for (unsigned int a = 0; a < pSampleCountPerStratum[stratum]; ++a)
        {
            unsigned int const i = pSampleIndicesPerStratum[stratum][a];

            if (pDiscreteSamples[i] != pDiscreteSamples[i])
                continue;
            if (pTimeSamples[i] != pTimeSamples[i])
                continue;

            for (unsigned int b = 0; b < pSampleCountPerStratum[stratum]; ++b)
            {
                unsigned int const j = pSampleIndicesPerStratum[stratum][b];

                if (pDiscreteSamples[j] != pDiscreteSamples[j])
                    continue;
                if (pTimeSamples[j] != pTimeSamples[j])
                    continue;

                float pair_weight = pSampleWeights[i] * pSampleWeights[j];

                if (pTimeSamples[i] < pTimeSamples[j] && pDiscreteSamples[i] == 1)
                {
                    relevant_weight += pair_weight;

                    if (pContinuousSamples[i] > pContinuousSamples[j])
                        concordant_weight += pair_weight;
                    else if (pContinuousSamples[i] < pContinuousSamples[j])
                        discordant_weight += pair_weight;
                    else if (outX)
                        uninformative_weight += pair_weight;
                    else
                        discordant_weight += pair_weight;
                }
                else if (pTimeSamples[i] > pTimeSamples[j] && pDiscreteSamples[j] == 1)
                {
                    relevant_weight += pair_weight;

                    if (pContinuousSamples[i] < pContinuousSamples[j])
                        concordant_weight += pair_weight;
                    else if (pContinuousSamples[i] > pContinuousSamples[j])
                        discordant_weight += pair_weight;
                    else if (outX)
                        uninformative_weight += pair_weight;
                    else
                        discordant_weight += pair_weight;
                }
            }
        }
    }

    return concordant_weight / relevant_weight;
}

/* static */float const
Math::computeCramersV(float const* const pSamplesX, float const* const pSamplesY,
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
                * computeCramersV(pSamplesX, pSamplesY, pSampleWeights, pSampleIndicesPerStratum[i],
                        pSampleCountPerStratum[i]);
        total_weight += pTotalWeightPerStratum[i];
    }

    r /= total_weight;

    return r;
}

/* static */float const
Math::computeCramersV(float const* const pSamplesX, float const* const pSamplesY,
        float const* const pSampleWeights, unsigned int const* const pSampleIndices,
        unsigned int const sampleCount)
{
    unsigned int x_class_count = 0;
    unsigned int y_class_count = 0;

    for (unsigned int i = 0; i < sampleCount; ++i)
    {
        unsigned int const index = pSampleIndices[i];
        if (x_class_count <= pSamplesX[index])
            x_class_count = pSamplesX[index] + 1;
        if (y_class_count <= pSamplesY[index])
            y_class_count = pSamplesY[index] + 1;
    }

    Matrix contingency_table(x_class_count + 1, y_class_count + 1);
    for (unsigned int i = 0; i <= x_class_count; ++i)
        for (unsigned int j = 0; j <= y_class_count; ++j)
            contingency_table.at(i, j) = 0;

    for (unsigned int i = 0; i < sampleCount; ++i)
    {
        unsigned int const index = pSampleIndices[i];
        float const sample_weight = pSampleWeights[index];
        contingency_table.at(pSamplesX[index], pSamplesY[index]) += sample_weight;
        contingency_table.at(x_class_count, pSamplesY[index]) += sample_weight;
        contingency_table.at(pSamplesX[index], y_class_count) += sample_weight;
        contingency_table.at(x_class_count, y_class_count) += sample_weight;
    }

    float chi_square = 0.;

    for (unsigned int i = 0; i < x_class_count; ++i)
        for (unsigned int j = 0; j < y_class_count; ++j)
        {
            float expected_value = contingency_table.at(i, y_class_count)
                    * contingency_table.at(x_class_count, j)
                    / contingency_table.at(x_class_count, y_class_count);

            chi_square += std::pow((contingency_table.at(i, j) - expected_value), 2)
                    / expected_value;
        }

    unsigned int const min_classes =
            (x_class_count < y_class_count) ? x_class_count : y_class_count;
    float const v = std::sqrt(
            chi_square / (contingency_table.at(x_class_count, y_class_count) * (min_classes - 1)));

    return v;
}

/* static */float const
Math::computeFisherTransformation(float const r)
{
    return 0.5 * std::log((1 + r) / (1 - r));
}

/* static */float const
Math::computeFisherTransformationReverse(float const z)
{
    float const exp = std::exp(2 * z);
    return (exp - 1) / (exp + 1);
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

/* static */void const
Math::placeStratificationData(unsigned int const* const pSampleStrata,
        float const* const pSampleWeights, unsigned int** const pSampleIndicesPerStratum,
        float* const pTotalWeightPerStratum, unsigned int* const pSampleCountPerStratum,
        unsigned int const sampleStratumCount, unsigned int const sampleCount)
{
    unsigned int p_iterator_per_stratum[sampleStratumCount];
    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        pTotalWeightPerStratum[i] = 0.;
        p_iterator_per_stratum[i] = 0;
        pSampleCountPerStratum[i] = 0;
    }

    for (unsigned int i = 0; i < sampleCount; ++i)
        ++pSampleCountPerStratum[pSampleStrata[i]];

    for (unsigned int i = 0; i < sampleStratumCount; ++i)
        pSampleIndicesPerStratum[i] = new unsigned int[pSampleCountPerStratum[i]];

    for (unsigned int i = 0; i < sampleCount; ++i)
    {
        unsigned int const p_sample_stratum = pSampleStrata[i];
        pSampleIndicesPerStratum[p_sample_stratum][p_iterator_per_stratum[p_sample_stratum]++] = i;
        pTotalWeightPerStratum[p_sample_stratum] += pSampleWeights[i];
    }
}