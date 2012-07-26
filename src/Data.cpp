#include "Data.hpp"

Data::Data(float* const pData, unsigned int const sampleCount, unsigned int const featureCount,
        unsigned int const* const pSampleStrata, float const* const pSampleWeights,
        unsigned int const* const pFeatureTypes, unsigned int const sampleStratumCount) :
        mpDataMatrix(new Matrix(pData, sampleCount, featureCount)), mpRankedDataMatrix(
                new Matrix(sampleCount, featureCount)), mpHasFeatureRanksCached(
                new bool[mpDataMatrix->getColumnCount()]), mpSampleStrata(pSampleStrata), mpSampleWeights(
                pSampleWeights), mpFeatureTypes(pFeatureTypes), mSampleStratumCount(
                sampleStratumCount), mpSampleIndicesPerStratum(
                new unsigned int*[sampleStratumCount]), mpSampleCountPerStratum(
                new unsigned int[sampleStratumCount])
{
    for (unsigned int i = 0; i < mpDataMatrix->getColumnCount(); ++i)
        mpHasFeatureRanksCached[i] = false;

    unsigned int p_iterator_per_stratum[mSampleStratumCount];
    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
    {
        mpSampleCountPerStratum[i] = 0;
        p_iterator_per_stratum[i] = 0;
    }

    for (unsigned int i = 0; i < mpDataMatrix->getRowCount(); ++i)
        ++mpSampleCountPerStratum[mpSampleStrata[i]];

    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
        mpSampleIndicesPerStratum[i] = new unsigned int[mpSampleCountPerStratum[i]];

    for (unsigned int i = 0; i < mpDataMatrix->getRowCount(); ++i)
        mpSampleIndicesPerStratum[mpSampleStrata[i]][p_iterator_per_stratum[mpSampleStrata[i]]++] =
                i;
}

Data::~Data()
{
    delete mpDataMatrix;
    delete mpRankedDataMatrix;
    delete[] mpHasFeatureRanksCached;
    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
        delete[] mpSampleIndicesPerStratum[i];
    delete[] mpSampleIndicesPerStratum;
    delete[] mpSampleCountPerStratum;
}

float const
Data::computeMiBetweenFeatures(unsigned int const i, unsigned int const j) const
{
    float r = std::numeric_limits<double>::quiet_NaN();

    bool const A_is_continuous = mpFeatureTypes[i] == FEATURE_CONTINUOUS;
    bool const A_is_discrete = mpFeatureTypes[i] == FEATURE_DISCRETE;
    bool const A_is_survival_event = mpFeatureTypes[i] == FEATURE_SURVIVAL_EVENT;

    bool const B_is_continuous = mpFeatureTypes[j] == FEATURE_CONTINUOUS;
    bool const B_is_discrete = mpFeatureTypes[j] == FEATURE_DISCRETE;

    if (A_is_continuous && B_is_continuous)
        r = computeSpearmanCorrelationBetweenFeatures(i, j);
    else if (A_is_survival_event && B_is_continuous)
        r = computeConcordanceIndex(i, j, i + 1, mpDataMatrix, mpSampleWeights, mpSampleStrata,
                true);
    else if (A_is_discrete && B_is_continuous)
        r = computeConcordanceIndex(i, j, -1, mpDataMatrix, mpSampleWeights, mpSampleStrata, true);
    else if (A_is_continuous && B_is_discrete)
        r = computeConcordanceIndex(j, i, -1, mpDataMatrix, mpSampleWeights, mpSampleStrata, true);
    else if (A_is_discrete && B_is_discrete)
        r = computeCramersV(i, j, mpDataMatrix, mpSampleWeights);

    return -0.5 * log(1 - (r * r));
}

float const
Data::computeSpearmanCorrelationBetweenFeatures(unsigned int const i, unsigned int const j) const
{
    if (!mpHasFeatureRanksCached[i])
    {
        placeRanksByFeatureIndex(i);
        mpHasFeatureRanksCached[i] = true;
    }

    if (!mpHasFeatureRanksCached[j])
    {
        placeRanksByFeatureIndex(j);
        mpHasFeatureRanksCached[j] = true;
    }

    return computePearsonCorrelation(i, j, mpRankedDataMatrix, mpSampleWeights);
}

unsigned int const
Data::getSampleCount() const
{
    return mpDataMatrix->getRowCount();
}

unsigned int const
Data::getFeatureCount() const
{
    return mpDataMatrix->getColumnCount();
}

void const
Data::placeRanksByFeatureIndex(unsigned int const index) const
{
    unsigned int const sample_count = getSampleCount();
    unsigned int p_order[sample_count];

    for (unsigned int i = 0; i < sample_count; ++i)
        p_order[i] = i;

    std::sort(p_order, p_order + sample_count, DataMatrixComparator(index, mpDataMatrix));

    for (unsigned int i = 0; i < sample_count; ++i)
        (*mpRankedDataMatrix)(p_order[i], index) = i;
}

Data::DataMatrixComparator::DataMatrixComparator(unsigned int const featureIndex,
        Matrix const* const pDataMatrix) :
        mFeatureIndex(featureIndex), mpDataMatrix(pDataMatrix)
{

}

bool const
Data::DataMatrixComparator::operator()(unsigned int const i, unsigned int const j) const
{
    return (*mpDataMatrix)(i, mFeatureIndex) < (*mpDataMatrix)(j, mFeatureIndex);
}
