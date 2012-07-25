#include "tools.hpp"

DataMatrixComparator::DataMatrixComparator(unsigned int const featureIndex,
        Matrix const* const pDataMatrix) :
        mFeatureIndex(featureIndex), mpDataMatrix(pDataMatrix)
{

}

bool const
DataMatrixComparator::operator()(unsigned int const i, unsigned int const j) const
{
    return (*mpDataMatrix)(i, mFeatureIndex) < (*mpDataMatrix)(j, mFeatureIndex);
}

float const
computeConcordanceIndex(unsigned int const discreteFeatureIndex,
        unsigned int const continuousFeatureIndex, int const timeFeatureIndex,
        Matrix const* const pDataMatrix, float const* const pSampleWeights,
        unsigned int const* const pSampleStrata, bool const outX)
{
    float concordant_weight = 0.;
    float discordant_weight = 0.;
    float uninformative_weight = 0.;
    float relevant_weight = 0.;

    unsigned int const sample_count = pDataMatrix->getRowCount();
    for (unsigned int i = 0; i < sample_count; ++i)
    {
        if ((*pDataMatrix)(i, discreteFeatureIndex) != (*pDataMatrix)(i, discreteFeatureIndex))
            continue;
        if (timeFeatureIndex >= 0
                && (*pDataMatrix)(i, timeFeatureIndex) != (*pDataMatrix)(i, timeFeatureIndex))
            continue;
        for (unsigned int j = 0; j < sample_count; ++j)
        {
            if ((*pDataMatrix)(j, discreteFeatureIndex) != (*pDataMatrix)(j, discreteFeatureIndex))
                continue;
            if (timeFeatureIndex >= 0
                    && (*pDataMatrix)(j, timeFeatureIndex) != (*pDataMatrix)(j, timeFeatureIndex))
                continue;
            if (pSampleStrata[i] != pSampleStrata[j])
                continue;

            float pair_weight = pSampleWeights[i] * pSampleWeights[j];
            if (isComparablePair(i, j, timeFeatureIndex, discreteFeatureIndex, pDataMatrix))
            {
                relevant_weight += pair_weight;
                if ((*pDataMatrix)(i, continuousFeatureIndex)
                        > (*pDataMatrix)(j, continuousFeatureIndex))
                    concordant_weight += pair_weight;
                else if ((*pDataMatrix)(i, continuousFeatureIndex)
                        < (*pDataMatrix)(j, continuousFeatureIndex))
                    discordant_weight += pair_weight;
                else if (outX)
                    uninformative_weight += pair_weight;
                else
                    discordant_weight += pair_weight;
            }

            else if (isComparablePair(j, i, timeFeatureIndex, discreteFeatureIndex, pDataMatrix))
            {
                relevant_weight += pair_weight;
                if ((*pDataMatrix)(i, continuousFeatureIndex)
                        < (*pDataMatrix)(j, continuousFeatureIndex))
                    concordant_weight += pair_weight;
                else if ((*pDataMatrix)(i, continuousFeatureIndex)
                        > (*pDataMatrix)(j, continuousFeatureIndex))
                    discordant_weight += pair_weight;
                else if (outX)
                    uninformative_weight += pair_weight;
                else
                    discordant_weight += pair_weight;
            }

        }
    }
    return concordant_weight / relevant_weight;
}

bool const
isComparablePair(unsigned int const i, unsigned int const j, int const timeFeatureIndex,
        unsigned int const discreteFeatureIndex, Matrix const* const pDataMatrix)
{
    return (timeFeatureIndex >= 0
            && ((*pDataMatrix)(i, timeFeatureIndex) < (*pDataMatrix)(j, timeFeatureIndex)
                    && (*pDataMatrix)(i, discreteFeatureIndex) == 1))
            || (timeFeatureIndex < 0
                    && (*pDataMatrix)(i, discreteFeatureIndex)
                            > (*pDataMatrix)(j, discreteFeatureIndex));
}

float const
computeCramersV(unsigned int const featureIndex1, unsigned int const featureIndex2,
        Matrix const* const pDataMatrix, float const* const pSampleWeights)
{

    unsigned int const sample_count = pDataMatrix->getRowCount();
    unsigned int pX_class_count = 0;
    unsigned int pY_class_count = 0;

    for (unsigned int i = 0; i < sample_count; ++i)
    {
        if (pX_class_count <= (*pDataMatrix)(i, featureIndex1))
            pX_class_count = (*pDataMatrix)(i, featureIndex1) + 1;
        if (pY_class_count <= (*pDataMatrix)(i, featureIndex2))
            pY_class_count = (*pDataMatrix)(i, featureIndex2) + 1;
    }

    Matrix contingency_table(pX_class_count + 1, pY_class_count + 1);
    for (unsigned int i = 0; i < pX_class_count; ++i)
        for (unsigned int j = 0; j < pY_class_count; ++j)
            contingency_table(i, j) = 0;

    for (unsigned int i = 0; i < sample_count; ++i)
    {
        float const sample_weight = pSampleWeights[i];
        contingency_table((*pDataMatrix)(i, featureIndex1), (*pDataMatrix)(i, featureIndex2)) +=
                sample_weight;
        contingency_table(pX_class_count, (*pDataMatrix)(i, featureIndex2)) += sample_weight;
        contingency_table((*pDataMatrix)(i, featureIndex1), pY_class_count) += sample_weight;
        contingency_table(pX_class_count, pY_class_count) += sample_weight;
    }

    float chi_square = computeChiSquare(contigency_table);
    unsigned int min_classes = (pX_class_count < pY_class_count) ? pX_class_count : pY_class_count;

    return std::sqrt(
            chi_square / (contingency_table(pX_class_count, pY_class_count) * (min_classes - 1)));
}

float const
computeChiSquare(Matrix const* const pContingencyTable)
{
    unsigned int const row_count = pContingencyTable->getRowCount();
    unsigned int const col_count = pContingencyTable->getColumnCount();

    for (unsigned int i = 0; i < row_count; ++i)
        for (unsigned int j = 0; j < col_count; ++j)
        {
            float expected_value = contingency_table(i, col_count) * contingency_table(row_count, j)
                    / contingency_table(row_count, col_count);
            chi_square += std::pow((contingency_table(i, j) - expected_value), 2) / expected_value;
        }
}

float const
computeSpearmanCorrelation(unsigned int const i, unsigned int const j,
        Matrix const* const pRankedDataMatrix, float const* const pSampleWeights)
{
    float const* const a = &(*pRankedDataMatrix)(0, i);
    float const* const b = &(*pRankedDataMatrix)(0, j);
    unsigned int const sample_count = pRankedDataMatrix->getRowCount();
    float sum = 0.;
    float total_weight = 0.;

    for (unsigned int n = 0; n < sample_count; ++n)
    {
        float const difference = a[n] - b[n];
        sum += pSampleWeights[n] * difference * difference;
        total_weight += pSampleWeights[n];
    }

    float const correlation = 1
            - ((6 * sum) / (total_weight * ((total_weight * total_weight) - 1)));

    return correlation;
}

void const
placeRanksByFeatureIndex(unsigned int const index, Matrix* const pRankedDataMatrix,
        Matrix const* const pDataMatrix)
{
    unsigned int const sample_count = pRankedDataMatrix->getRowCount();
    unsigned int p_order[sample_count];

    for (unsigned int i = 0; i < sample_count; ++i)
        p_order[i] = i;

    std::sort(p_order, p_order + sample_count, DataMatrixComparator(index, pDataMatrix));

    for (unsigned int i = 0; i < sample_count; ++i)
        (*pRankedDataMatrix)(p_order[i], index) = i;
}
