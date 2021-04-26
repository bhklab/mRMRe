#include "Data.h"
#include <sys/time.h>

Data::Data(double* const pData, Matrix const* const pPriorsMatrix, double const priorsWeight,
        unsigned int const sampleCount, unsigned int const featureCount,
        int const* const pSampleStrata, double const* const pSampleWeights,
        int const* const pFeatureTypes, unsigned int const sampleStratumCount,
        unsigned int const continuousEstimator, bool const outX, unsigned int const bootstrapCount) :
        mpDataMatrix(new Matrix(pData, sampleCount, featureCount)), mpOrderMatrix(
                continuousEstimator ? new Matrix(sampleCount, featureCount) : 0), mpPriorsMatrix(
                pPriorsMatrix), mpHasOrderCached(new bool[mpDataMatrix->getColumnCount()]), mpSampleStrata(
                pSampleStrata), mpSampleWeights(pSampleWeights), mpFeatureTypes(pFeatureTypes), mSampleStratumCount(
                sampleStratumCount), mpSampleIndicesPerStratum(
                new unsigned int*[sampleStratumCount]), mpMasterSampleIndicesPerStratum(
                new unsigned int*[sampleStratumCount]), mpSampleCountPerStratum(
                new unsigned int[sampleStratumCount]), mContinuousEstimator(continuousEstimator), mOutX(
                outX), mBootstrapCount(bootstrapCount), mPriorsWeight(priorsWeight)
{
    for (unsigned int i = 0; i < mpDataMatrix->getColumnCount(); ++i)
        mpHasOrderCached[i] = false;

    Math::placeStratificationData(mpSampleStrata, mpSampleWeights, mpSampleIndicesPerStratum,
            mpSampleCountPerStratum, mSampleStratumCount, sampleCount);

    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
    {
        mpMasterSampleIndicesPerStratum[i] = new unsigned int[mpSampleCountPerStratum[i]];
        for (unsigned int j = 0; j < mpSampleCountPerStratum[i]; ++j)
            mpMasterSampleIndicesPerStratum[i][j] = mpSampleIndicesPerStratum[i][j];
    }
}

Data::~Data()
{
    delete mpDataMatrix;
    delete mpOrderMatrix;
    delete[] mpHasOrderCached;
    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
    {
        delete[] mpSampleIndicesPerStratum[i];
        delete[] mpMasterSampleIndicesPerStratum[i];
    }
    delete[] mpSampleIndicesPerStratum;
    delete[] mpMasterSampleIndicesPerStratum;
    delete[] mpSampleCountPerStratum;
}

void const
Data::bootstrap()
{
    // unsigned int seed = std::time(NULL); too long for small datasets
    struct timeval start;
    gettimeofday(&start, NULL);
    unsigned int seed = start.tv_usec; //microseconds
    
    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
        for (unsigned int j = 0; j < mpSampleCountPerStratum[i]; ++j)
        {
            unsigned int index = Math::computeRandomNumber(&seed) % mpSampleCountPerStratum[i];
            mpSampleIndicesPerStratum[i][j] = mpMasterSampleIndicesPerStratum[i][index];
        }
}

void const
Data::computeMiBetweenFeatures(unsigned int const i, unsigned int const j, double* const mi_ij,
        double* const mi_ji) const
{
    double val_ij = std::numeric_limits<double>::quiet_NaN();
    double val_ji = std::numeric_limits<double>::quiet_NaN();

    bool const A_is_continuous = mpFeatureTypes[i] == FEATURE_CONTINUOUS;
    bool const A_is_discrete = mpFeatureTypes[i] == FEATURE_DISCRETE;
    bool const A_is_survival_event = mpFeatureTypes[i] == FEATURE_SURVIVAL_EVENT;

    bool const B_is_continuous = mpFeatureTypes[j] == FEATURE_CONTINUOUS;
    bool const B_is_discrete = mpFeatureTypes[j] == FEATURE_DISCRETE;
    bool const B_is_survival_event = mpFeatureTypes[j] == FEATURE_SURVIVAL_EVENT;

    if (A_is_continuous && B_is_continuous)
    {
        switch (mContinuousEstimator)
        {
        case PEARSON_ESTIMATOR:
        {
            val_ij = val_ji = Math::computePearsonCorrelation(&(mpDataMatrix->at(0, i)),
                    &(mpDataMatrix->at(0, j)), mpSampleWeights, mpSampleIndicesPerStratum,
                    mpSampleCountPerStratum, mSampleStratumCount, mBootstrapCount);
        }
            break;

        case SPEARMAN_ESTIMATOR:
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
            val_ij = val_ji = Math::computePearsonCorrelation(p_ranked_samples_x,
                    p_ranked_samples_y, mpSampleWeights, mpSampleIndicesPerStratum,
                    mpSampleCountPerStratum, mSampleStratumCount, mBootstrapCount);
            delete[] p_ranked_samples_x;
            delete[] p_ranked_samples_y;
        }
            break;

        case KENDALL_ESTIMATOR:
        {
            val_ij =  Math::computeSomersD(
              Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)),
                                            &(mpDataMatrix->at(0, j)), mpSampleWeights, mpSampleIndicesPerStratum,
                                            mpSampleCountPerStratum, mSampleStratumCount, mOutX));
            val_ji = Math::computeSomersD(
              Math::computeConcordanceIndex(&(mpDataMatrix->at(0, j)),
                                            &(mpDataMatrix->at(0, i)), mpSampleWeights, mpSampleIndicesPerStratum,
                                            mpSampleCountPerStratum, mSampleStratumCount, mOutX));
        }
            break;

        case FREQUENCY_ESTIMATOR:
        {
            val_ij = Math::computeFrequency(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                    mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                    mSampleStratumCount, mBootstrapCount);
            val_ji = 1 - val_ij;
        }
            break;
        }
    }
    else if (A_is_discrete && B_is_continuous) // Not symmetrical
        val_ij = Math::computeSomersD(
                Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                        mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                        mSampleStratumCount, mOutX));
    else if (A_is_continuous && B_is_discrete) // Not symmetrical
        val_ij = Math::computeSomersD(
                Math::computeConcordanceIndex(&(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i)),
                        mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                        mSampleStratumCount, mOutX));
    else if (A_is_discrete && B_is_discrete)
        val_ij = val_ji = Math::computeCramersV(&(mpDataMatrix->at(0, i)),
                &(mpDataMatrix->at(0, j)), mpSampleWeights, mpSampleIndicesPerStratum,
                mpSampleCountPerStratum, mSampleStratumCount, mBootstrapCount);
    else if (A_is_survival_event && B_is_continuous)
        val_ij = val_ji = Math::computeSomersD(
                Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                        &(mpDataMatrix->at(0, i + 1)), mpSampleWeights, mpSampleIndicesPerStratum,
                        mpSampleCountPerStratum, mSampleStratumCount, mOutX));
    else if (A_is_continuous && B_is_survival_event)
        val_ij = val_ji = Math::computeSomersD(
                Math::computeConcordanceIndex(&(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i)),
                        &(mpDataMatrix->at(0, j + 1)), mpSampleWeights, mpSampleIndicesPerStratum,
                        mpSampleCountPerStratum, mSampleStratumCount, mOutX));
    else if (A_is_survival_event && B_is_discrete)
        val_ij = val_ji = Math::computeSomersD(
                Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                        &(mpDataMatrix->at(0, i + 1)), mpSampleWeights, mpSampleIndicesPerStratum,
                        mpSampleCountPerStratum, mSampleStratumCount, mOutX));
    else if (A_is_discrete && B_is_survival_event)
        val_ij = val_ji = Math::computeSomersD(
                Math::computeConcordanceIndex(&(mpDataMatrix->at(0, j)), &(mpDataMatrix->at(0, i)),
                        &(mpDataMatrix->at(0, j + 1)), mpSampleWeights, mpSampleIndicesPerStratum,
                        mpSampleCountPerStratum, mSampleStratumCount, mOutX));
    else if (A_is_survival_event && B_is_survival_event) // Not symmetrical for some reason
        val_ij = Math::computeSomersD(
                Math::computeConcordanceIndex(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                        &(mpDataMatrix->at(0, i + 1)), &(mpDataMatrix->at(0, j + 1)),
                        mpSampleWeights, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                        mSampleStratumCount, mOutX));

    if (mpPriorsMatrix != 0)
    {
        val_ij = (std::fabs(1.0 - mPriorsWeight) * val_ij)
                	+ (mPriorsWeight * mpPriorsMatrix->at(i, j));

        val_ji = (std::fabs(1.0 - mPriorsWeight) * val_ji)
                	+ (mPriorsWeight * mpPriorsMatrix->at(j, i));
    }

    if (val_ij == val_ij)
        *mi_ij = val_ij;

    if (val_ji == val_ji)
        *mi_ji = val_ji;
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
