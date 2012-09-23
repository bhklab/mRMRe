#include "Data.h"

Data::Data(double* const pData, Matrix const* const pPriorsMatrix, double const priorsWeight,
        unsigned int const sampleCount, unsigned int const featureCount,
        int const* const pSampleStrata, double const* const pSampleWeights,
        int const* const pFeatureTypes, unsigned int const sampleStratumCount, bool const usesRanks,
        bool const outX, unsigned int const bootstrapCount) :
        mpDataMatrix(new Matrix(pData, sampleCount, featureCount)), mpOrderMatrix(
                usesRanks ? new Matrix(sampleCount, featureCount) : 0), mpPriorsMatrix(
                pPriorsMatrix), mpHasOrderCached(new bool[mpDataMatrix->getColumnCount()]), mpSampleStrata(
                pSampleStrata), mpSampleWeights(pSampleWeights), mpFeatureTypes(pFeatureTypes), mSampleStratumCount(
                sampleStratumCount), mpSampleIndicesPerStratum(
                new unsigned int*[sampleStratumCount]), mpSampleCountPerStratum(
                new unsigned int[sampleStratumCount]), mUsesRanks(usesRanks), mOutX(outX), mBootstrapCount(
                bootstrapCount), mPriorsWeight(priorsWeight)
{
    for (unsigned int i = 0; i < mpDataMatrix->getColumnCount(); ++i)
        mpHasOrderCached[i] = false;

    Math::placeStratificationData(mpSampleStrata, mpSampleWeights, mpSampleIndicesPerStratum,
            mpSampleCountPerStratum, mSampleStratumCount, sampleCount);
}

Data::~Data()
{
    delete mpDataMatrix;
    delete mpOrderMatrix;
    delete[] mpHasOrderCached;
    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
        delete[] mpSampleIndicesPerStratum[i];
    delete[] mpSampleIndicesPerStratum;
    delete[] mpSampleCountPerStratum;
}

double const
Data::computeMiBetweenFeatures(unsigned int const i, unsigned int const j) const
{
    double r = std::numeric_limits<double>::quiet_NaN();

    if (i == j)
        return 1.;
    else
    {
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
                    Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)),
                            &(mpDataMatrix->at(0, j)), mpSampleWeights, mpSampleIndicesPerStratum,
                            mpSampleCountPerStratum, mSampleStratumCount, true));
        else if (A_is_continuous && B_is_discrete)
            r = Math::computeSomersD(
                    Math::computeConcordanceIndex(&(mpDataMatrix->at(0, j)),
                            &(mpDataMatrix->at(0, i)), mpSampleWeights, mpSampleIndicesPerStratum,
                            mpSampleCountPerStratum, mSampleStratumCount, true));
        else if (A_is_discrete && B_is_discrete)
            r = Math::computeCramersV(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                    mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                    mSampleStratumCount, mBootstrapCount);
        else if (A_is_survival_event && B_is_continuous)
            r = Math::computeSomersD(
                    Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)),
                            &(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i + 1)),
                            mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                            mSampleStratumCount, mOutX));
        else if (A_is_continuous && B_is_survival_event)
            r = Math::computeSomersD(
                    Math::computeConcordanceIndex(&(mpDataMatrix->at(0, j)),
                            &(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j + 1)),
                            mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                            mSampleStratumCount, mOutX));
        else if (A_is_survival_event && B_is_discrete)
            r = Math::computeSomersD(
                    Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)),
                            &(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i + 1)),
                            mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                            mSampleStratumCount, mOutX));
        else if (A_is_discrete && B_is_survival_event)
            r = Math::computeSomersD(
                    Math::computeConcordanceIndex(&(mpDataMatrix->at(0, j)),
                            &(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j + 1)),
                            mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                            mSampleStratumCount, mOutX));
        else if (A_is_survival_event && B_is_survival_event)
            r = Math::computeSomersD(
                    Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)),
                            &(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i + 1)),
                            &(mpDataMatrix->at(0, j + 1)), mpSampleWeights,
                            mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount,
                            mOutX));
    }

    if (mpPriorsMatrix == 0)
        return r;
    else
    {
        double sign = r < 0 ? -1. : 1.;
        return (sign * (std::fabs(1.0 - mPriorsWeight) * r))
                + (mPriorsWeight * mpPriorsMatrix->at(i, j));
    }
}

double const
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

        double* const p_ranked_samples_x = new double[getSampleCount()];
        double* const p_ranked_samples_y = new double[getSampleCount()];
        Math::placeRanksFromOrders(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                &(mpOrderMatrix->at(0, i)), &(mpOrderMatrix->at(0, j)), p_ranked_samples_x,
                p_ranked_samples_y, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                mSampleStratumCount);
        double const r = Math::computePearsonCorrelation(p_ranked_samples_x, p_ranked_samples_y,
                mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                mSampleStratumCount, mBootstrapCount);
        delete[] p_ranked_samples_x;
        delete[] p_ranked_samples_y;

        return r;
    }
    else
        return Math::computePearsonCorrelation(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                mSampleStratumCount, mBootstrapCount);
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
