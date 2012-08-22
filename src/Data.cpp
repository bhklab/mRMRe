#include "Data.h"

Data::Data(float* const pData, unsigned int const sampleCount, unsigned int const featureCount,
        unsigned int const* const pSampleStrata, float const* const pSampleWeights,
        unsigned int const* const pFeatureTypes, unsigned int const sampleStratumCount,
        bool const usesRanks, bool const outX, unsigned int const bootstrapCount) :
        mpDataMatrix(new Matrix(pData, sampleCount, featureCount)), mpOrderMatrix(
                usesRanks ? new Matrix(sampleCount, featureCount) : 0), mpHasOrderCached(
                new bool[mpDataMatrix->getColumnCount()]), mpSampleStrata(pSampleStrata), mpSampleWeights(
                pSampleWeights), mpFeatureTypes(pFeatureTypes), mSampleStratumCount(
                sampleStratumCount), mpSampleIndicesPerStratum(
                new unsigned int*[sampleStratumCount]), mpTotalWeightPerStratum(
                new float[sampleStratumCount]), mpSampleCountPerStratum(
                new unsigned int[sampleStratumCount]), mUsesRanks(usesRanks), mOutX(outX), mBootstrapCount(
                bootstrapCount)
{
    for (unsigned int i = 0; i < mpDataMatrix->getColumnCount(); ++i)
        mpHasOrderCached[i] = false;

    Math::placeStratificationData(mpSampleStrata, mpSampleWeights, mpSampleIndicesPerStratum,
            mpTotalWeightPerStratum, mpSampleCountPerStratum, mSampleStratumCount, sampleCount);
}

Data::~Data()
{
    delete mpDataMatrix;
    if (mpOrderMatrix)
        delete mpOrderMatrix;
    delete[] mpHasOrderCached;
    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
        delete[] mpSampleIndicesPerStratum[i];
    delete[] mpSampleIndicesPerStratum;
    delete[] mpTotalWeightPerStratum;
    delete[] mpSampleCountPerStratum;
}

float const
Data::computeMiBetweenFeatures(unsigned int const i, unsigned int const j) const
{
    if (i == j)
        return std::numeric_limits<float>::infinity();

    float r = std::numeric_limits<float>::quiet_NaN();

    bool const A_is_continuous = mpFeatureTypes[i] == FEATURE_CONTINUOUS;
    bool const A_is_discrete = mpFeatureTypes[i] == FEATURE_DISCRETE;
    bool const A_is_survival_event = mpFeatureTypes[i] == FEATURE_SURVIVAL_EVENT;

    bool const B_is_continuous = mpFeatureTypes[j] == FEATURE_CONTINUOUS;
    bool const B_is_discrete = mpFeatureTypes[j] == FEATURE_DISCRETE;
    bool const B_is_survival_event = mpFeatureTypes[j] == FEATURE_SURVIVAL_EVENT;

    if (A_is_continuous && B_is_continuous)
        r = computeCorrelationBetweenContinuousFeatures(i, j);
    else if (A_is_discrete && B_is_continuous)
        r = Math::computeSomersD(
                Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                        mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                        mSampleStratumCount, true));
    else if (A_is_continuous && B_is_discrete)
        r = Math::computeSomersD(
                Math::computeConcordanceIndex(&(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i)),
                        mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                        mSampleStratumCount, true));
    else if (A_is_discrete && B_is_discrete)
        r = Math::computeCramersV(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                mpSampleWeights, mpSampleIndicesPerStratum, mpTotalWeightPerStratum,
                mpSampleCountPerStratum, mSampleStratumCount, mBootstrapCount);
    else if (A_is_survival_event && B_is_continuous)
        r = Math::computeSomersD(
                Math::computeConcordanceIndexWithTime(&(mpDataMatrix->at(0, i)),
                        &(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i + 1)), mpSampleWeights,
                        mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount,
                        mOutX));
    else if (A_is_continuous && B_is_survival_event)
        r = Math::computeSomersD(
                Math::computeConcordanceIndexWithTime(&(mpDataMatrix->at(0, j)),
                        &(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j + 1)), mpSampleWeights,
                        mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount,
                        mOutX));
    else if (A_is_survival_event && B_is_discrete)
        r = Math::computeSomersD(
                Math::computeConcordanceIndexWithTime(&(mpDataMatrix->at(0, i)),
                        &(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i + 1)), mpSampleWeights,
                        mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount,
                        mOutX));
    else if (A_is_discrete && B_is_survival_event)
        r = Math::computeSomersD(
                Math::computeConcordanceIndexWithTime(&(mpDataMatrix->at(0, j)),
                        &(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j + 1)), mpSampleWeights,
                        mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount,
                        mOutX));

    return Math::computeMi(r);
}

float const
Data::computeCorrelationBetweenContinuousFeatures(unsigned int const i, unsigned int const j) const
{
    if (mUsesRanks)
    {
        if (!mpHasOrderCached[i])
        {
            Math::placeOrders(&(mpDataMatrix->at(0, i)), &(mpOrderMatrix->at(0, i)),
                    mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount);
            mpHasOrderCached[i] = true;
        }

        if (!mpHasOrderCached[j])
        {
            Math::placeOrders(&(mpDataMatrix->at(0, j)), &(mpOrderMatrix->at(0, j)),
                    mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount);
            mpHasOrderCached[j] = true;
        }

        float* const p_ranked_samples_x = new float[mpDataMatrix->getRowCount()];
        float* const p_ranked_samples_y = new float[mpDataMatrix->getRowCount()];
        Math::placeRanksFromOrders(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                &(mpOrderMatrix->at(0, i)), &(mpOrderMatrix->at(0, j)), p_ranked_samples_x,
                p_ranked_samples_y, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                mSampleStratumCount);
        float const r = Math::computePearsonCorrelation(p_ranked_samples_x, p_ranked_samples_y,
                mpSampleWeights, mpSampleIndicesPerStratum, mpTotalWeightPerStratum,
                mpSampleCountPerStratum, mSampleStratumCount, mBootstrapCount);
        delete[] p_ranked_samples_x;
        delete[] p_ranked_samples_y;

        return r;
    }
    else
        return Math::computePearsonCorrelation(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                mpSampleWeights, mpSampleIndicesPerStratum, mpTotalWeightPerStratum,
                mpSampleCountPerStratum, mSampleStratumCount, mBootstrapCount);
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
