#include "Data.h"

#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>

// R's RNG is not thread-safe: draw from it serially, never inside an OpenMP region.
static unsigned int
drawSeed()
{
    GetRNGstate();
    unsigned int const seed = static_cast<unsigned int>(R_unif_index(4294967296.));
    PutRNGstate();
    return seed;
}

Data::Data(double* const pData, Matrix const* const pPriorsMatrix, double const priorsWeight,
        unsigned int const sampleCount, unsigned int const featureCount,
        int const* const pSampleStrata, double const* const pSampleWeights,
        int const* const pFeatureTypes, unsigned int const sampleStratumCount,
        unsigned int const continuousEstimator, bool const outX, unsigned int const bootstrapCount) :
        mpDataMatrix(new Matrix(pData, sampleCount, featureCount)), mpOrderMatrix(
                continuousEstimator ? new Matrix(sampleCount, featureCount) : 0), mpPriorsMatrix(
                pPriorsMatrix), mpSampleStrata(
                pSampleStrata), mpSampleWeights(pSampleWeights), mpFeatureTypes(pFeatureTypes), mSampleStratumCount(
                sampleStratumCount), mpSampleIndicesPerStratum(
                new unsigned int*[sampleStratumCount]), mpMasterSampleIndicesPerStratum(
                new unsigned int*[sampleStratumCount]), mpSampleCountPerStratum(
                new unsigned int[sampleStratumCount]), mContinuousEstimator(continuousEstimator), mOutX(
                outX), mBootstrapCount(bootstrapCount), mPriorsWeight(priorsWeight)
{
    // computeMiBetweenFeatures() reads the time column at i + 1 of every survival event column.
    for (unsigned int i = 0; i < featureCount; ++i)
    {
        bool const is_event = pFeatureTypes[i] == FEATURE_SURVIVAL_EVENT;
        bool const is_time = pFeatureTypes[i] == FEATURE_SURVIVAL_TIME;
        if (is_event && (i + 1 >= featureCount || pFeatureTypes[i + 1] != FEATURE_SURVIVAL_TIME))
            Rf_error("survival event column %u must be followed by a survival time column", i + 1);
        if (is_time && (i == 0 || pFeatureTypes[i - 1] != FEATURE_SURVIVAL_EVENT))
            Rf_error("survival time column %u must be preceded by a survival event column", i + 1);
    }

    mSeed = bootstrapCount > 0 ? drawSeed() : 0;

    Math::placeStratificationData(mpSampleStrata, mpSampleWeights, mpSampleIndicesPerStratum,
            mpSampleCountPerStratum, mSampleStratumCount, sampleCount);

    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
    {
        mpMasterSampleIndicesPerStratum[i] = new unsigned int[mpSampleCountPerStratum[i]];
        for (unsigned int j = 0; j < mpSampleCountPerStratum[i]; ++j)
            mpMasterSampleIndicesPerStratum[i][j] = mpSampleIndicesPerStratum[i][j];
    }

    placeOrders();
}

Data::~Data()
{
    delete mpDataMatrix;
    delete mpOrderMatrix;
    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
    {
        delete[] mpSampleIndicesPerStratum[i];
        delete[] mpMasterSampleIndicesPerStratum[i];
    }
    delete[] mpSampleIndicesPerStratum;
    delete[] mpMasterSampleIndicesPerStratum;
    delete[] mpSampleCountPerStratum;
}

void
Data::bootstrap()
{
    GetRNGstate();

    for (unsigned int i = 0; i < mSampleStratumCount; ++i)
        for (unsigned int j = 0; j < mpSampleCountPerStratum[i]; ++j)
        {
            unsigned int const index = static_cast<unsigned int>(R_unif_index(mpSampleCountPerStratum[i]));
            mpSampleIndicesPerStratum[i][j] = mpMasterSampleIndicesPerStratum[i][index];
        }

    PutRNGstate();

    if (mBootstrapCount > 0)
        mSeed = drawSeed();

    placeOrders();
}

// Spearman orders depend on the (bootstrapped) sample indices. Computing them all here, serially,
// keeps computeMiBetweenFeatures() read-only so it is safe inside OpenMP regions.
void
Data::placeOrders()
{
    if (mContinuousEstimator != SPEARMAN_ESTIMATOR)
        return;

    for (unsigned int i = 0; i < mpDataMatrix->getColumnCount(); ++i)
        if (mpFeatureTypes[i] == FEATURE_CONTINUOUS)
            Math::placeOrders(&(mpDataMatrix->at(0, i)), &(mpOrderMatrix->at(0, i)),
                    mpSampleIndicesPerStratum, mpSampleCountPerStratum, mSampleStratumCount);
}

void
Data::computeMiBetweenFeatures(unsigned int const i, unsigned int const j, double* const mi_ij,
        double* const mi_ji) const
{
    unsigned int const seed = mSeed + (i * getFeatureCount()) + j;
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
                    mpSampleCountPerStratum, mSampleStratumCount, mBootstrapCount, seed);
        }
            break;

        case SPEARMAN_ESTIMATOR:
        {
            double* const p_ranked_samples_x = new double[getSampleCount()];
            double* const p_ranked_samples_y = new double[getSampleCount()];
            Math::placeRanksFromOrders(&(mpDataMatrix->at(0, i)), &(mpDataMatrix->at(0, j)),
                    &(mpOrderMatrix->at(0, i)), &(mpOrderMatrix->at(0, j)), p_ranked_samples_x,
                    p_ranked_samples_y, mpSampleIndicesPerStratum, mpSampleCountPerStratum,
                    mSampleStratumCount);
            val_ij = val_ji = Math::computePearsonCorrelation(p_ranked_samples_x,
                    p_ranked_samples_y, mpSampleWeights, mpSampleIndicesPerStratum,
                    mpSampleCountPerStratum, mSampleStratumCount, mBootstrapCount, seed);
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
                    mSampleStratumCount, mBootstrapCount, seed);
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
                mpSampleCountPerStratum, mSampleStratumCount, mBootstrapCount, seed);
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

unsigned int
Data::getSampleCount() const
{
    return mpDataMatrix->getRowCount();
}

unsigned int
Data::getFeatureCount() const
{
    return mpDataMatrix->getColumnCount();
}
