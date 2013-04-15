#include "Math.h"

Math::IndirectComparator::IndirectComparator(double const* const pSamples,
        unsigned int const* const pSampleIndices) :
        mpSamples(pSamples), mpSampleIndices(pSampleIndices)
{

}

bool const
Math::IndirectComparator::operator()(unsigned int const i, unsigned int const j) const
{
    return mpSamples[mpSampleIndices[i]] < mpSamples[mpSampleIndices[j]];
}

/* static */void const
Math::computeCausality(double* const pCausalityArray, Matrix const* const pMiMatrix,
        int const* const pSolutions, unsigned int const solutionCount,
        unsigned int const featureCountPerSolution, unsigned int const featureCount,
        unsigned int const targetFeatureIndex)
{
    for (unsigned int s = 0; s < solutionCount; ++s)
    {
        for (unsigned int i = 0; i < featureCountPerSolution - 1; ++i)
        {
            for (unsigned int j = i + 1; j < featureCountPerSolution; ++j)
            {
                int const a = pSolutions[(featureCountPerSolution * s) + i];
                int const b = pSolutions[(featureCountPerSolution * s) + j];

                double const cor_ij =
                        (std::fabs(pMiMatrix->at(a, b)) > std::fabs(pMiMatrix->at(b, a))) ?
                                pMiMatrix->at(a, b) : pMiMatrix->at(b, a);

                double const cor_ik = pMiMatrix->at(a, targetFeatureIndex);
                double const cor_jk = pMiMatrix->at(b, targetFeatureIndex);

                double const coefficient = Math::computeCoInformationLattice(cor_ij, cor_ik,
                        cor_jk);

                if (pCausalityArray[a] != pCausalityArray[a] || pCausalityArray[a] > coefficient)
                    pCausalityArray[a] = coefficient;

                if (pCausalityArray[b] != pCausalityArray[b] || pCausalityArray[b] > coefficient)
                    pCausalityArray[b] = coefficient;
            }
        }
    }
}

/* static */double const
Math::computeCoInformationLattice(double const cor_ij, double const cor_ik, double const cor_jk)
{
    double const cor_ij_sq = cor_ij * cor_ij;
    double const cor_jk_sq = cor_jk * cor_jk;
    double const cor_ik_sq = cor_ik * cor_ik;

    return -.5
            * std::log(
                    ((1 - cor_ij_sq) * (1 - cor_ik_sq) * (1 - cor_jk_sq))
                            / (1 + 2 * cor_ij * cor_ik * cor_jk - cor_ij_sq - cor_ik_sq - cor_jk_sq));
}

/* static */double const
Math::computeConcordanceIndex(double const* const pDiscreteSamples,
        double const* const pContinuousSamples, double const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        bool const outX, double* const pConcordantWeights, double* const pDiscordantWeights,
        double* const pUninformativeWeights, double* const pRelevantWeights)
{
    double sum_concordant_weight = 0.;
    double sum_relevant_weight = 0.;

    for (unsigned int stratum = 0; stratum < sampleStratumCount; ++stratum)
    {
        for (unsigned int a = 0; a < pSampleCountPerStratum[stratum]; ++a)
        {
            unsigned int const i = pSampleIndicesPerStratum[stratum][a];

            double concordant_weight = 0.;
            double discordant_weight = 0.;
            double uninformative_weight = 0.;
            double relevant_weight = 0.; 

            if (pDiscreteSamples[i] != pDiscreteSamples[i]
                    || pContinuousSamples[i] != pContinuousSamples[i])
                continue;

            for (unsigned int b = 0; b < pSampleCountPerStratum[stratum]; ++b)
            {
                unsigned int const j = pSampleIndicesPerStratum[stratum][b];

                if (pDiscreteSamples[j] != pDiscreteSamples[j]
                        || pContinuousSamples[j] != pContinuousSamples[j])
                    continue;

                double pair_weight = pSampleWeights[i] * pSampleWeights[j];

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

            sum_concordant_weight += concordant_weight;
            sum_relevant_weight += relevant_weight;

            if (pConcordantWeights != 0) // Implicity, the other similar vectors
            {                            // should also match this condition.
                pConcordantWeights[i] = concordant_weight;
                pDiscordantWeights[i] = discordant_weight;
                pUninformativeWeights[i] = uninformative_weight;
                pRelevantWeights[i] = relevant_weight;
            }
        }
    }

    return sum_concordant_weight / sum_relevant_weight;
}

/*static*/double const
Math::computeConcordanceIndex(double const* const pDiscreteSamples,
        double const* const pContinuousSamples, double const* const pTimeSamples,
        double const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        bool const outX, double* const pConcordantWeights, double* const pDiscordantWeights,
        double* const pUninformativeWeights, double* const pRelevantWeights)
{
    double sum_concordant_weight = 0.;
    double sum_relevant_weight = 0.;

    for (unsigned int stratum = 0; stratum < sampleStratumCount; ++stratum)
    {
        for (unsigned int a = 0; a < pSampleCountPerStratum[stratum]; ++a)
        {
            unsigned int const i = pSampleIndicesPerStratum[stratum][a];

            double concordant_weight = 0.;
            double discordant_weight = 0.;
            double uninformative_weight = 0.;
            double relevant_weight = 0.; 

            if (pDiscreteSamples[i] != pDiscreteSamples[i] || pTimeSamples[i] != pTimeSamples[i]
                    || pContinuousSamples[i] != pContinuousSamples[i])
                continue;

            for (unsigned int b = 0; b < pSampleCountPerStratum[stratum]; ++b)
            {
                unsigned int const j = pSampleIndicesPerStratum[stratum][b];

                if (pDiscreteSamples[j] != pDiscreteSamples[j] || pTimeSamples[j] != pTimeSamples[j]
                        || pContinuousSamples[j] != pContinuousSamples[j])
                    continue;

                double pair_weight = pSampleWeights[i] * pSampleWeights[j];

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

            sum_concordant_weight += concordant_weight;
            sum_relevant_weight += relevant_weight;

            if (pConcordantWeights != 0) // Implicity, the other similar vectors
            {                            // should also match this condition.
                pConcordantWeights[i] = concordant_weight;
                pDiscordantWeights[i] = discordant_weight;
                pUninformativeWeights[i] = uninformative_weight;
                pRelevantWeights[i] = relevant_weight;
            }
        }
    }

    return sum_concordant_weight / sum_relevant_weight;
}

/*static*/double const
Math::computeConcordanceIndex(double const* const pDiscreteSamplesX,
        double const* const pDiscreteSamplesY, double const* const pTimeSamplesX,
        double const* const pTimeSamplesY, double const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        bool const outX, double* const pConcordantWeights, double* const pDiscordantWeights,
        double* const pUninformativeWeights, double* const pRelevantWeights)
{
    double sum_concordant_weight = 0.;
    double sum_relevant_weight = 0.;

    for (unsigned int stratum = 0; stratum < sampleStratumCount; ++stratum)
    {
        for (unsigned int a = 0; a < pSampleCountPerStratum[stratum]; ++a)
        {
            unsigned int const i = pSampleIndicesPerStratum[stratum][a];

            double concordant_weight = 0.;
            double discordant_weight = 0.;
            double uninformative_weight = 0.;
            double relevant_weight = 0.; 

            if (pDiscreteSamplesX[i] != pDiscreteSamplesX[i]
                    || pDiscreteSamplesY[i] != pDiscreteSamplesY[i]
                    || pTimeSamplesX[i] != pTimeSamplesX[i] || pTimeSamplesY[i] != pTimeSamplesY[i])
                continue;

            for (unsigned int b = 0; b < pSampleCountPerStratum[stratum]; ++b)
            {
                unsigned int const j = pSampleIndicesPerStratum[stratum][b];

                if (pDiscreteSamplesX[j] != pDiscreteSamplesX[j]
                        || pDiscreteSamplesY[j] != pDiscreteSamplesY[j]
                        || pTimeSamplesX[j] != pTimeSamplesX[j]
                        || pTimeSamplesY[j] != pTimeSamplesY[j])
                    continue;

                double pair_weight = pSampleWeights[i] * pSampleWeights[j];

                if (pTimeSamplesX[i] < pTimeSamplesX[j] && pDiscreteSamplesX[i] == 1)
                {
                    relevant_weight += pair_weight;

                    if (pTimeSamplesY[i] > pTimeSamplesY[j] && pDiscreteSamplesY[j] == 1)
                        concordant_weight += pair_weight;
                    else if (pTimeSamplesY[i] < pTimeSamplesY[j] && pDiscreteSamplesY[j] == 1)
                        discordant_weight += pair_weight;
                    else if (outX)
                        uninformative_weight += pair_weight;
                    else
                        discordant_weight += pair_weight;
                }
                else if (pTimeSamplesX[i] > pTimeSamplesX[j] && pDiscreteSamplesX[j] == 1)
                {
                    relevant_weight += pair_weight;

                    if (pTimeSamplesY[i] > pTimeSamplesY[j] && pDiscreteSamplesY[j] == 1)
                        concordant_weight += pair_weight;
                    else if (pTimeSamplesY[i] < pTimeSamplesY[j] && pDiscreteSamplesY[j] == 1)
                        discordant_weight += pair_weight;
                    else if (outX)
                        uninformative_weight += pair_weight;
                    else
                        discordant_weight += pair_weight;
                }
            }

           sum_concordant_weight += concordant_weight;
            sum_relevant_weight += relevant_weight;

            if (pConcordantWeights != 0) // Implicity, the other similar vectors
            {                            // should also match this condition.
                pConcordantWeights[i] = concordant_weight;
                pDiscordantWeights[i] = discordant_weight;
                pUninformativeWeights[i] = uninformative_weight;
                pRelevantWeights[i] = relevant_weight;
            }
        }
    }

    return sum_concordant_weight / sum_relevant_weight;
}

/* static */double const
Math::computeCramersV(double const* const pSamplesX, double const* const pSamplesY,
        double const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        unsigned int const bootstrapCount)
{
    bool const useBootstrap = bootstrapCount > 3 && sampleStratumCount > 0;
    double* p_error_per_stratum = 0;

    if (useBootstrap)
    {
        p_error_per_stratum = new double[sampleStratumCount];
        unsigned int seed = std::time(NULL);
        Matrix bootstraps(bootstrapCount, sampleStratumCount);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) firstprivate(seed)
#endif
        for (unsigned int i = 0; i < bootstrapCount; ++i)
        {
            for (unsigned int j = 0; j < sampleStratumCount; ++j)
            {
                unsigned int const sample_count = pSampleCountPerStratum[j];
                unsigned int* const p_samples = new unsigned int[sample_count];

                for (unsigned int k = 0; k < sample_count; ++k)
                    p_samples[k] = pSampleIndicesPerStratum[j][Math::computeRandomNumber(&seed)
                            % sample_count];

                double const correlation = computeCramersV(pSamplesX, pSamplesY, pSampleWeights,
                        p_samples, sample_count);
                bootstraps.at(i, j) = correlation;

                delete[] p_samples;
            }
        }

        for (unsigned int i = 0; i < sampleStratumCount; ++i)
            p_error_per_stratum[i] = 1.
                    / Math::computeVariance(&(bootstraps.at(0, i)), bootstrapCount);
    }

    double r = 0.;
    double total_weight = 0.;

    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        double weight = 0.;
        double const correlation = computeCramersV(pSamplesX, pSamplesY, pSampleWeights,
                pSampleIndicesPerStratum[i], pSampleCountPerStratum[i], &weight);

        if (useBootstrap)
            weight = p_error_per_stratum[i];

        r += weight * correlation;
        total_weight += weight;
    }

    r /= total_weight;

    delete[] p_error_per_stratum;

    return r;
}

/* static */double const
Math::computeCramersV(double const* const pSamplesX, double const* const pSamplesY,
        double const* const pSampleWeights, unsigned int const* const pSampleIndices,
        unsigned int const sampleCount, double* const pTotalWeight)
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

        if (pSamplesX[index] != pSamplesX[index] || pSamplesY[index] != pSamplesY[index])
            continue;

        double const sample_weight = pSampleWeights[index];
        contingency_table.at(pSamplesX[index], pSamplesY[index]) += sample_weight;
        contingency_table.at(x_class_count, pSamplesY[index]) += sample_weight;
        contingency_table.at(pSamplesX[index], y_class_count) += sample_weight;
        contingency_table.at(x_class_count, y_class_count) += sample_weight;
    }

    double chi_square = 0.;

    for (unsigned int i = 0; i < x_class_count; ++i)
        for (unsigned int j = 0; j < y_class_count; ++j)
        {
            double expected_value = contingency_table.at(i, y_class_count)
                    * contingency_table.at(x_class_count, j)
                    / contingency_table.at(x_class_count, y_class_count);

            chi_square += std::pow((contingency_table.at(i, j) - expected_value), 2)
                    / expected_value;
        }

    unsigned int const min_classes =
            (x_class_count < y_class_count) ? x_class_count : y_class_count;

    *pTotalWeight = contingency_table.at(x_class_count, y_class_count);

    double const v = std::sqrt(chi_square / ((*pTotalWeight) * (min_classes - 1)));

    return v;
}

/* static */double const
Math::computeFrequency(double const* const pSamplesX, double const* const pSamplesY,
        double const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        unsigned int const bootstrapCount)
{
    bool const useBootstrap = bootstrapCount > 3 && sampleStratumCount > 0;
    double* p_error_per_stratum = 0;

    if (useBootstrap)
    {
        p_error_per_stratum = new double[sampleStratumCount];
        unsigned int seed = std::time(NULL);
        Matrix bootstraps(bootstrapCount, sampleStratumCount);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) firstprivate(seed)
#endif
        for (unsigned int i = 0; i < bootstrapCount; ++i)
        {
            for (unsigned int j = 0; j < sampleStratumCount; ++j)
            {
                unsigned int const sample_count = pSampleCountPerStratum[j];
                unsigned int* const p_samples = new unsigned int[sample_count];

                for (unsigned int k = 0; k < sample_count; ++k)
                    p_samples[k] = pSampleIndicesPerStratum[j][Math::computeRandomNumber(&seed)
                            % sample_count];

                double const correlation = computeFrequency(pSamplesX, pSamplesY, pSampleWeights,
                        p_samples, sample_count);
                bootstraps.at(i, j) = correlation;

                delete[] p_samples;
            }
        }

        for (unsigned int i = 0; i < sampleStratumCount; ++i)
            p_error_per_stratum[i] = 1.
                    / Math::computeVariance(&(bootstraps.at(0, i)), bootstrapCount);
    }

    double r = 0.;
    double total_weight = 0.;

    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        double weight = 0.;
        double const correlation = computeFrequency(pSamplesX, pSamplesY, pSampleWeights,
                pSampleIndicesPerStratum[i], pSampleCountPerStratum[i], &weight);

        if (useBootstrap)
            weight = p_error_per_stratum[i];

        r += weight * correlation;
        total_weight += weight;
    }

    r /= total_weight;

    delete[] p_error_per_stratum;

    return r;
}

/* static */double const
Math::computeFrequency(double const* const pSamplesX, double const* const pSamplesY,
        double const* const pSampleWeights, unsigned int const* const pSampleIndices,
        unsigned int const sampleCount, double* const pTotalWeight)
{
    double sum = 0.;
    double total_weight = 0.;
    double r = 0.;

    for (unsigned int i = 0; i < sampleCount; ++i)
    {
        unsigned int const sample_index = pSampleIndices[i];
        double const sample_weight = pSampleWeights[sample_index];

        if (pSamplesX[sample_index] == pSamplesX[sample_index]
                && pSamplesY[sample_index] == pSamplesY[sample_index])
        {
            total_weight += sample_weight;

            if (pSamplesX[sample_index] > pSamplesY[sample_index])
                sum += sample_weight;
        }
    }

    if (pTotalWeight != 0)
        *pTotalWeight = total_weight;
    
    r = sum / total_weight;
    
    return r;
}

/* static */double const
Math::computeFisherTransformation(double const r)
{
    return 0.5 * std::log((1 + r) / (1 - r));
}

/* static */double const
Math::computeFisherTransformationReverse(double const z)
{
    double const exp = std::exp(2 * z);
    return (exp - 1) / (exp + 1);
}

/* static */double const
Math::computeMi(double const r)
{
    return -0.5 * std::log(1 - (r * r));
}

/* static */double const
Math::computePearsonCorrelation(double const* const pSamplesX, double const* const pSamplesY,
        double const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        unsigned int const bootstrapCount)
{
    bool const useBootstrap = bootstrapCount > 3 && sampleStratumCount > 0;
    double* p_error_per_stratum = 0;

    if (useBootstrap)
    {
        p_error_per_stratum = new double[sampleStratumCount];
        unsigned int seed = std::time(NULL);
        Matrix bootstraps(bootstrapCount, sampleStratumCount);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) firstprivate(seed)
#endif
        for (unsigned int i = 0; i < bootstrapCount; ++i)
        {
            for (unsigned int j = 0; j < sampleStratumCount; ++j)
            {
                unsigned int const sample_count = pSampleCountPerStratum[j];
                unsigned int* const p_samples = new unsigned int[sample_count];

                for (unsigned int k = 0; k < sample_count; ++k)
                    p_samples[k] = pSampleIndicesPerStratum[j][Math::computeRandomNumber(&seed)
                            % sample_count];

                double const correlation = computeCramersV(pSamplesX, pSamplesY, pSampleWeights,
                        p_samples, sample_count);
                bootstraps.at(i, j) = correlation;

                delete[] p_samples;
            }
        }

        for (unsigned int i = 0; i < sampleStratumCount; ++i)
            p_error_per_stratum[i] = 1.
                    / Math::computeVariance(&(bootstraps.at(0, i)), bootstrapCount);
    }

    double r = 0.;
    double total_weight = 0.;

    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
        double weight = 0.;
        double const correlation = computePearsonCorrelation(pSamplesX, pSamplesY, pSampleWeights,
                pSampleIndicesPerStratum[i], pSampleCountPerStratum[i], &weight);

        if (useBootstrap)
            weight = p_error_per_stratum[i];

        r += weight * correlation;
        total_weight += weight;
    }

    r /= total_weight;

    delete[] p_error_per_stratum;

    return r;
}

/* static */double const
Math::computePearsonCorrelation(double const* const pSamplesX, double const* const pSamplesY,
        double const* const pSampleWeights, unsigned int const* const pSampleIndices,
        unsigned int const sampleCount, double* const pTotalWeight)
{
    double sum_of_x = 0.;
    double sum_of_x_x = 0.;
    double sum_of_y = 0.;
    double sum_of_y_y = 0.;
    double sum_of_x_y = 0.;
    double sum_of_weights = 0.;

    for (unsigned int i = 0; i < sampleCount; ++i)
    {
        double const my_x = pSamplesX[pSampleIndices[i]];
        double const my_y = pSamplesY[pSampleIndices[i]];

        if (my_x == my_x && my_y == my_y)
        {
            double const my_weight = pSampleWeights[pSampleIndices[i]];
            sum_of_x += my_x * my_weight;
            sum_of_x_x += my_x * my_x * my_weight;
            sum_of_y += my_y * my_weight;
            sum_of_y_y += my_y * my_y * my_weight;
            sum_of_x_y += my_x * my_y * my_weight;
            sum_of_weights += my_weight;
        }
    }

    double const r = (sum_of_x_y - ((sum_of_x * sum_of_y) / sum_of_weights))
            / std::sqrt(
                    (sum_of_x_x - ((sum_of_x * sum_of_x) / sum_of_weights))
                            * (sum_of_y_y - ((sum_of_y * sum_of_y) / sum_of_weights)));

    *pTotalWeight = sum_of_weights;

    return r;
}

/* static */int const
Math::computeRandomNumber(unsigned int* const seed)
{
    unsigned int next = *seed;
    int result;

    next *= 1103515245;
    next += 12345;
    result = static_cast<unsigned int>(next / 65536) % 2048;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= static_cast<unsigned int>(next / 65536) % 1024;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= static_cast<unsigned int>(next / 65536) % 1024;

    *seed = next;
    return result;
}

/* static */double const
Math::computeSomersD(double const c)
{
    return (c - 0.5) * 2;
}

/* static */double const
Math::computeSpearmanCorrelation(double const* const pSamplesX, double const* const pSamplesY,
        double const* const pSampleWeights,
        unsigned int const* const * const pSampleIndicesPerStratum,
        unsigned int const* const pSampleCountPerStratum, unsigned int const sampleStratumCount,
        unsigned int const bootstrapCount, unsigned int const sampleCount)
{
    double* const p_ordered_samples_x = new double[sampleCount];
    double* const p_ordered_samples_y = new double[sampleCount];

    Math::placeOrders(&pSamplesX[0], p_ordered_samples_x, pSampleIndicesPerStratum,
            pSampleCountPerStratum, sampleStratumCount);
    Math::placeOrders(&pSamplesY[0], p_ordered_samples_y, pSampleIndicesPerStratum,
            pSampleCountPerStratum, sampleStratumCount);

    double* const p_ranked_samples_x = new double[sampleCount];
    double* const p_ranked_samples_y = new double[sampleCount];

    Math::placeRanksFromOrders(&pSamplesX[0], &pSamplesY[0], p_ordered_samples_x,
            p_ordered_samples_y, p_ranked_samples_x, p_ranked_samples_y, pSampleIndicesPerStratum,
            pSampleCountPerStratum, sampleStratumCount);

    delete[] p_ordered_samples_x;
    delete[] p_ordered_samples_y;

    double const r = Math::computePearsonCorrelation(p_ranked_samples_x, p_ranked_samples_y,
            &pSampleWeights[0], pSampleIndicesPerStratum, pSampleCountPerStratum,
            sampleStratumCount, bootstrapCount);

    delete[] p_ranked_samples_x;
    delete[] p_ranked_samples_y;

    return r;
}

/* static */double const
Math::computeVariance(double const* const pSamples, unsigned int const sampleCount)
{
    if (sampleCount == 0)
        return 0.;

    double sum_for_mean = pSamples[0];
    double sum_for_error = 0.;

    for (unsigned int i = 1; i < sampleCount; ++i)
    {
        double const my_sum = pSamples[i] - sum_for_mean;
        double const my_mean = ((i - 1) * my_sum) / i;
        sum_for_mean += my_mean;
        sum_for_error += my_mean * my_sum;
    }

    return sum_for_error / (sampleCount - 1);
}

/* static */void const
Math::placeOrders(double const* const pSamples, double* const pOrders,
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
Math::placeRanksFromOrders(double const* const pSamplesX, double const* const pSamplesY,
        double const* const pOrdersX, double const* const pOrdersY, double* const pRanksX,
        double* const pRanksY, unsigned int const* const * const pSampleIndicesPerStratum,
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
                pRanksX[order_x] = std::numeric_limits<double>::quiet_NaN();
            else
                pRanksX[order_x] = offset_x++;

            if (NA_y)
                pRanksY[order_y] = std::numeric_limits<double>::quiet_NaN();
            else
                pRanksY[order_y] = offset_y++;
        }
    }
}

/* static */void const
Math::placeRanksFromSamples(double const* const pSamples, double* const pRanks,
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
            pRanks[j] = std::numeric_limits<double>::quiet_NaN();

        for (unsigned int j = 0; j < sample_count - offset; ++j)
            pRanks[p_sample_indices[p_order[j]]] = j;

        delete[] p_order;
    }
}

/* static */void const
Math::placeStratificationData(int const* const pSampleStrata, double const* const pSampleWeights,
        unsigned int** const pSampleIndicesPerStratum, unsigned int* const pSampleCountPerStratum,
        unsigned int const sampleStratumCount, unsigned int const sampleCount)
{
    unsigned int* const p_iterator_per_stratum = new unsigned int[sampleStratumCount];

    for (unsigned int i = 0; i < sampleStratumCount; ++i)
    {
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
    }

    delete[] p_iterator_per_stratum;
}
