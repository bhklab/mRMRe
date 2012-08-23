#include "Math.h"

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
        bool const outX, float* const pConcordantWeight, float* const pDiscordantWeight,
        float* const pUninformativeWeight, float* const pRelevantWeight)
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

    if (pConcordantWeight)
        *pConcordantWeight = concordant_weight;
    if (pDiscordantWeight)
        *pDiscordantWeight = discordant_weight;
    if (pUninformativeWeight)
        *pUninformativeWeight = uninformative_weight;
    if (pRelevantWeight)
        *pRelevantWeight = relevant_weight;

    return concordant_weight / relevant_weight;
}

/*static*/float const
Math::computeConcordanceIndexWithTime(float const* const pDiscreteSamples,
        float const* const pContinuousSamples, float const* const pTimeSamples,
        float const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        bool const outX, float* const pConcordantWeight, float* const pDiscordantWeight,
        float* const pUninformativeWeight, float* const pRelevantWeight)
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

    if (pConcordantWeight)
        *pConcordantWeight = concordant_weight;
    if (pDiscordantWeight)
        *pDiscordantWeight = discordant_weight;
    if (pUninformativeWeight)
        *pUninformativeWeight = uninformative_weight;
    if (pRelevantWeight)
        *pRelevantWeight = relevant_weight;

    return concordant_weight / relevant_weight;
}

/* static */float const
Math::computeCramersV(float const* const pSamplesX, float const* const pSamplesY,
        float const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        float const* const pTotalWeightPerStratum, unsigned int const* const pSampleCountPerStratum,
        unsigned int const sampleStratumCount, unsigned int const bootstrapCount)
{
    float* const p_error_per_stratum = new float[sampleStratumCount];
    float const* p_weight_per_stratum = pTotalWeightPerStratum;

    if (bootstrapCount > 3 && sampleStratumCount > 0)
    {
        unsigned int seed = std::time(NULL);
        Matrix bootstraps(bootstrapCount, sampleStratumCount);

#pragma omp parallel for schedule(dynamic) firstprivate(seed)
        for (unsigned int i = 0; i < bootstrapCount; ++i)
        {
            for (unsigned int j = 0; j < sampleStratumCount; ++j)
            {
                unsigned int const sample_count = pSampleCountPerStratum[j];
                unsigned int* const p_samples = new unsigned int[sample_count];

                for (unsigned int k = 0; k < sample_count; ++k)
                    p_samples[k] = pSampleIndicesPerStratum[j][Math::computeRandomNumber(&seed) % sample_count];

                float const correlation = computeCramersV(pSamplesX, pSamplesY, pSampleWeights,
                        p_samples, sample_count);
                bootstraps.at(i, j) = correlation;

                delete[] p_samples;
            }
        }

        for (unsigned int i = 0; i < sampleStratumCount; ++i)
            p_error_per_stratum[i] = 1.
                    / Math::computeVariance(&(bootstraps.at(0, i)), bootstrapCount);

        if (p_error_per_stratum[0] != 0.)
            p_weight_per_stratum = p_error_per_stratum;
    }

    float r = 0.;
    float total_weight = 0.;

    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        float const correlation = computeCramersV(pSamplesX, pSamplesY, pSampleWeights,
                pSampleIndicesPerStratum[i], pSampleCountPerStratum[i]);
        r += p_weight_per_stratum[i] * correlation;
        total_weight += p_weight_per_stratum[i];
    }

    r /= total_weight;

    delete[] p_error_per_stratum;

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
        unsigned int const sampleStratumCount, unsigned int const bootstrapCount)
{
    float* const p_error_per_stratum = new float[sampleStratumCount];
    float const* p_weight_per_stratum = pTotalWeightPerStratum;

    if (bootstrapCount > 3 && sampleStratumCount > 0)
    {
        unsigned int seed = std::time(NULL);
        Matrix bootstraps(bootstrapCount, sampleStratumCount);

#pragma omp parallel for schedule(dynamic) firstprivate(seed)
        for (unsigned int i = 0; i < bootstrapCount; ++i)
        {
            for (unsigned int j = 0; j < sampleStratumCount; ++j)
            {
                unsigned int const sample_count = pSampleCountPerStratum[j];
                unsigned int* const p_samples = new unsigned int[sample_count];

                for (unsigned int k = 0; k < sample_count; ++k)
                    p_samples[k] = pSampleIndicesPerStratum[j][Math::computeRandomNumber(&seed) % sample_count];

                float const correlation = computePearsonCorrelation(pSamplesX, pSamplesY,
                        pSampleWeights, p_samples, sample_count);
                bootstraps.at(i, j) = correlation;

                delete[] p_samples;
            }
        }

        for (unsigned int i = 0; i < sampleStratumCount; ++i)
            p_error_per_stratum[i] = 1.
                    / Math::computeVariance(&(bootstraps.at(0, i)), bootstrapCount);

        if (p_error_per_stratum[0] != 0.)
            p_weight_per_stratum = p_error_per_stratum;
    }

    float r = 0.;
    float total_weight = 0.;

    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        float const correlation = computePearsonCorrelation(pSamplesX, pSamplesY, pSampleWeights,
                pSampleIndicesPerStratum[i], pSampleCountPerStratum[i]);
        r += p_weight_per_stratum[i] * correlation;
        total_weight += p_weight_per_stratum[i];
    }

    r /= total_weight;

    delete[] p_error_per_stratum;

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
        float const my_x = pSamplesX[pSampleIndices[i]];
        float const my_y = pSamplesY[pSampleIndices[i]];

        if (my_x == my_x && my_y == my_y)
        {
            float const my_weight = pSampleWeights[pSampleIndices[i]];
            sum_of_x += my_x * my_weight;
            sum_of_x_x += my_x * my_x * my_weight;
            sum_of_y += my_y * my_weight;
            sum_of_y_y += my_y * my_y * my_weight;
            sum_of_x_y += my_x * my_y * my_weight;
            sum_of_weights += my_weight;
        }
    }

    float const r = (sum_of_x_y - ((sum_of_x * sum_of_y) / sum_of_weights))
            / std::sqrt(
                    (sum_of_x_x - ((sum_of_x * sum_of_x) / sum_of_weights))
                            * (sum_of_y_y - ((sum_of_y * sum_of_y) / sum_of_weights)));

    return r;
}

/* static */int const
Math::computeRandomNumber(unsigned int* const seed)
{
    unsigned int next = *seed;
    int result;

    next *= 1103515245;
    next += 12345;
    result = (unsigned int) (next / 65536) % 2048;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;

    *seed = next;
    return result;
}

/* static */float const
Math::computeSomersD(float const c)
{
    return (c - 0.5) * 2;
}

/* static */float const
Math::computeVariance(float const* const pSamples, unsigned int const sampleCount)
{
    if (sampleCount == 0)
        return 0.;

    float sum_for_mean = pSamples[0];
    float sum_for_error = 0.;

    for (unsigned int i = 1; i < sampleCount; ++i)
    {
        float const my_sum = pSamples[i] - sum_for_mean;
        float const my_mean = ((i - 1) * my_sum) / i;
        sum_for_mean += my_mean;
        sum_for_error += my_mean * my_sum;
    }

    return sum_for_error / (sampleCount - 1);
}

/* static */void const
Math::placeOrders(float const* const pSamples, float* const pOrders,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount)
{
    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        unsigned int const* const p_sample_indices = pSampleIndicesPerStratum[i];
        unsigned int const sample_count = pSampleCountPerStratum[i];
        unsigned int* const p_order = new unsigned int[sample_count];

        unsigned int offset = 0;
        for (unsigned int j = 0; j < sample_count; ++j)
        {
            if (pSamples[p_sample_indices[j]] == pSamples[p_sample_indices[j]])
                p_order[j - offset] = j;
            else
                p_order[sample_count - 1 - offset++] = j;
        }

        std::sort(p_order, p_order + sample_count - offset,
                Math::IndirectComparator(pSamples, p_sample_indices));

        for (unsigned int j = 0; j < sample_count; ++j)
            pOrders[p_sample_indices[j]] = p_order[j];

        delete[] p_order;
    }
}

/* static */void const
Math::placeRanksFromOrders(float const* const pSamplesX, float const* const pSamplesY,
        float const* const pOrdersX, float const* const pOrdersY, float* const pRanksX,
        float* const pRanksY, unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount)
{
    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        unsigned int const* const p_sample_indices = pSampleIndicesPerStratum[i];
        unsigned int const stratum_sample_count = pSampleCountPerStratum[i];

        unsigned int offset_x = 0;
        unsigned int offset_y = 0;

        for (unsigned int j = 0; j < stratum_sample_count; ++j)
        {
            unsigned int const order_x =
                    p_sample_indices[static_cast<unsigned int>(pOrdersX[p_sample_indices[j]])];
            unsigned int const order_y =
                    p_sample_indices[static_cast<unsigned int>(pOrdersY[p_sample_indices[j]])];

            bool const NA_x = pSamplesX[order_x] != pSamplesX[order_x]
                    || pSamplesY[order_x] != pSamplesY[order_x];
            bool const NA_y = pSamplesY[order_y] != pSamplesY[order_y]
                    || pSamplesX[order_y] != pSamplesX[order_y];

            if (NA_x)
                pRanksX[order_x] = std::numeric_limits<float>::quiet_NaN();
            else
                pRanksX[order_x] = offset_x++;

            if (NA_y)
                pRanksY[order_y] = std::numeric_limits<float>::quiet_NaN();
            else
                pRanksY[order_y] = offset_y++;
        }
    }
}

/* static */void const
Math::placeRanksFromSamples(float const* const pSamples, float* const pRanks,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount)
{
    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        unsigned int const* const p_sample_indices = pSampleIndicesPerStratum[i];
        unsigned int const sample_count = pSampleCountPerStratum[i];
        unsigned int* const p_order = new unsigned int[sample_count];

        unsigned int offset = 0;
        for (unsigned int j = 0; j < sample_count; ++j)
        {
            unsigned int const my_index = p_sample_indices[j];

            if (pSamples[my_index] == pSamples[my_index])
                p_order[j - offset] = j;
            else
                ++offset;
        }

        std::sort(p_order, p_order + sample_count - offset,
                Math::IndirectComparator(pSamples, p_sample_indices));

        for (unsigned int j = 0; j < sample_count; ++j)
            pRanks[j] = std::numeric_limits<float>::quiet_NaN();

        for (unsigned int j = 0; j < sample_count - offset; ++j)
            pRanks[p_sample_indices[p_order[j]]] = j;

        delete[] p_order;
    }
}

/* static */void const
Math::placeStratificationData(unsigned int const* const pSampleStrata,
        float const* const pSampleWeights, unsigned int** const pSampleIndicesPerStratum,
        float* const pTotalWeightPerStratum, unsigned int* const pSampleCountPerStratum,
        unsigned int const sampleStratumCount, unsigned int const sampleCount)
{
    unsigned int* const p_iterator_per_stratum = new unsigned int[sampleStratumCount];

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

    delete[] p_iterator_per_stratum;
}
