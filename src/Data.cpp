#include "Data.hpp"

Data::Data(float* const pData, unsigned int const sampleCount, unsigned int const featureCount,
        unsigned int const* const pSampleStrata, float const* const pSampleWeights,
        unsigned int const* const pFeatureTypes, unsigned int const sampleStratumCount) :
        mpDataMatrix(new Matrix(pData, sampleCount, featureCount)), mpRankedDataMatrix(
                new Matrix(sampleCount, featureCount)), mpHasFeatureRanksCached(
                new bool[mpDataMatrix->getColumnCount()]), mpSampleStrata(pSampleStrata), mpSampleWeights(
                pSampleWeights), mpFeatureTypes(pFeatureTypes), mSampleStratumCount(
                sampleStratumCount), mpSampleIndicesPerStratum(
                new unsigned int*[sampleStratumCount]), mpTotalWeightPerStratum(
                new float[sampleStratumCount]), mpSampleCountPerStratum(
                new unsigned int[sampleStratumCount])
{
    for (unsigned int i = 0; i < mpDataMatrix->getColumnCount(); ++i)
        mpHasFeatureRanksCached[i] = false;

    unsigned int p_iterator_per_stratum[mSampleStratumCount];
    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
    {
        mpTotalWeightPerStratum[i] = 0.;
        p_iterator_per_stratum[i] = 0;
        mpSampleCountPerStratum[i] = 0;
    }

    for (unsigned int i = 0; i < mpDataMatrix->getRowCount(); ++i)
        ++mpSampleCountPerStratum[mpSampleStrata[i]];

    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
        mpSampleIndicesPerStratum[i] = new unsigned int[mpSampleCountPerStratum[i]];

    for (unsigned int i = 0; i < mpDataMatrix->getRowCount(); ++i)
    {
        unsigned int const p_sample_stratum = mpSampleStrata[i];
        mpSampleIndicesPerStratum[p_sample_stratum][p_iterator_per_stratum[p_sample_stratum]++] = i;
        mpTotalWeightPerStratum[p_sample_stratum] += pSampleWeights[i];
    }
}

Data::~Data()
{
    delete mpDataMatrix;
    delete mpRankedDataMatrix;
    delete[] mpHasFeatureRanksCached;
    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
        delete[] mpSampleIndicesPerStratum[i];
    delete[] mpSampleIndicesPerStratum;
    delete[] mpTotalWeightPerStratum;
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
    bool const B_is_survival_event = mpFeatureTypes[j] == FEATURE_SURVIVAL_EVENT;

    if (A_is_continuous && B_is_continuous)
        r = computeCorrelationBetweenContinuousFeatures(i, j);
    else if (A_is_discrete && B_is_continuous)
        r = Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                mSampleStratumCount, true);
    else if (A_is_continuous && B_is_discrete)
        r = Math::computeConcordanceIndex(&(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i)),
                mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                mSampleStratumCount, true);
    else if (A_is_discrete && B_is_discrete)
        r = Math::computeCramersV(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                mpSampleWeights, getSampleCount());
    else if (A_is_survival_event && B_is_continuous)
        r = Math::computeConcordanceIndexWithTime(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                &(mpDataMatrix->at(0, i + 1)), mpSampleWeights, mpSampleIndicesPerStratum,
                mpSampleCountPerStratum, mSampleStratumCount, true);
    else if (A_is_continuous && B_is_survival_event)
        r = Math::computeConcordanceIndexWithTime(&(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i)),
                &(mpDataMatrix->at(0, j + 1)), mpSampleWeights, mpSampleIndicesPerStratum,
                mpSampleCountPerStratum, mSampleStratumCount, true);

    return Math::computeMi(r);
}

float const
Data::computeCorrelationBetweenContinuousFeatures(unsigned int const i, unsigned int const j) const
{
    if (!mpHasFeatureRanksCached[i])
    {
        Math::placeRanksByFeatureIndex(&(mpDataMatrix->at(0, i)), &(mpRankedDataMatrix->at(0, i)),
                mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount);
        mpHasFeatureRanksCached[i] = true;
    }

    if (!mpHasFeatureRanksCached[j])
    {
        Math::placeRanksByFeatureIndex(&(mpDataMatrix->at(0, j)), &(mpRankedDataMatrix->at(0, j)),
                mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount);
        mpHasFeatureRanksCached[j] = true;
    }

    return Math::computePearsonCorrelation(&(mpRankedDataMatrix->at(0, i)),
            &(mpRankedDataMatrix->at(0, j)), mpSampleWeights, mpSampleIndicesPerStratum,
            mpTotalWeightPerStratum, mpSampleCountPerStratum, mSampleStratumCount);
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
