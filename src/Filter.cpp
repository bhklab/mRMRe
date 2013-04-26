#include "Filter.h"

Filter::Filter(int const* const pChildrenCountPerLevel, unsigned int const levelCount,
        Matrix* const pFeatureInformationMatrix, unsigned int const targetFeatureIndex) :
        mpChildrenCountPerLevel(pChildrenCountPerLevel), mLevelCount(levelCount), mpFeatureInformationMatrix(
                pFeatureInformationMatrix), mpStartingIndexPerLevel(
                new unsigned int[mLevelCount + 2])
{
    unsigned int cumulative_element_count = 1;
    unsigned int children_per_level = 1;

    mpStartingIndexPerLevel[0] = 0;

    for (unsigned int level = 0; level < mLevelCount; ++level)
    {
        mpStartingIndexPerLevel[level + 1] = cumulative_element_count;
        children_per_level *= mpChildrenCountPerLevel[level];
        cumulative_element_count += children_per_level;
    }

    mpStartingIndexPerLevel[mLevelCount + 1] = cumulative_element_count;
    mTreeElementCount = cumulative_element_count;
    mpIndexTree = new unsigned int[cumulative_element_count];
    mpScoreTree = new double [cumulative_element_count];

    for (unsigned int i = 0; i < mTreeElementCount; ++i)
    {
        mpIndexTree[i] = targetFeatureIndex;
        mpScoreTree[i] = 0;
    }
}

Filter::~Filter()
{
    delete[] mpStartingIndexPerLevel;
    delete[] mpIndexTree;
    delete[] mpScoreTree;
}

void const
Filter::build()
{
    for (unsigned int level = 0; level < mLevelCount; ++level)
    {
        unsigned int const parent_count = mpStartingIndexPerLevel[level + 1]
                - mpStartingIndexPerLevel[level];

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (unsigned int parent = 0; parent < parent_count; ++parent)
            placeElements(
                    mpStartingIndexPerLevel[level + 1] + (parent * mpChildrenCountPerLevel[level]),
                    mpChildrenCountPerLevel[level], level + 1);
    }
}

/* inline */unsigned int const
Filter::getParentAbsoluteIndex(unsigned int const absoluteIndex, unsigned int const level) const
{
    return (absoluteIndex - mpStartingIndexPerLevel[level]) / mpChildrenCountPerLevel[level - 1]
            + mpStartingIndexPerLevel[level - 1];
}

void const
Filter::getSolutions(int* const solutions) const
{
    unsigned int counter = 0;

    for (unsigned int end_element_absolute_index = mTreeElementCount - 1;
            end_element_absolute_index >= mpStartingIndexPerLevel[mLevelCount];
            --end_element_absolute_index)
    {
        unsigned int element_absolute_index = end_element_absolute_index;

        for (unsigned int level = mLevelCount; level > 0; --level)
        {
            solutions[counter++] = mpIndexTree[element_absolute_index];
            element_absolute_index = getParentAbsoluteIndex(element_absolute_index, level);
        }
    }
}

void const
Filter::getScores(double* const scores) const
{
    unsigned int counter = 0;

    for (unsigned int end_element_absolute_index = mTreeElementCount - 1;
            end_element_absolute_index >= mpStartingIndexPerLevel[mLevelCount];
            --end_element_absolute_index)
    {
        unsigned int element_absolute_index = end_element_absolute_index;

        for (unsigned int level = mLevelCount; level > 0; --level)
        {
            scores[counter++] = mpScoreTree[element_absolute_index];
            element_absolute_index = getParentAbsoluteIndex(element_absolute_index, level);
        }
    }
}

bool const
Filter::hasAncestorByFeatureIndex(unsigned int const absoluteIndex, unsigned int const featureIndex,
        unsigned int level) const
{
    // This function only considers the ancestry of the putative absolute/featureIndex

    unsigned int parent_absolute_index = absoluteIndex;

    for (unsigned int i = level; i > 0; --i)
    {
        parent_absolute_index = getParentAbsoluteIndex(parent_absolute_index, i);

        if (mpIndexTree[parent_absolute_index] == featureIndex)
            return true;
    }

    return false;
}

bool const
Filter::isRedundantPath(unsigned int const absoluteIndex, unsigned int const featureIndex,
        unsigned int const level) const
{
    for (unsigned int i = mpStartingIndexPerLevel[level]; i < mpStartingIndexPerLevel[level + 1];
            ++i)
    {
        if (mpIndexTree[i] == mpIndexTree[0])
            continue;

        unsigned int candidate_absolute_index = absoluteIndex;
        unsigned int candidate_feature_index = featureIndex;

        bool solution_is_redundant = true;

        for (unsigned int j = level; j > 0; --j)
        {
            unsigned int parent_absolute_index = i;

            bool feature_is_redundant = false;

            for (unsigned int k = level; k > 0; --k)
            {
                if (mpIndexTree[parent_absolute_index] == candidate_feature_index)
                {
                    feature_is_redundant = true;
                    break;
                }

                parent_absolute_index = getParentAbsoluteIndex(parent_absolute_index, k);
            }

            if (!feature_is_redundant)
            {
                solution_is_redundant = false;
                break;
            }

            candidate_absolute_index = getParentAbsoluteIndex(candidate_absolute_index, j);
            candidate_feature_index = mpIndexTree[candidate_absolute_index];
        }

        if (solution_is_redundant)
            return true;
    }

    return false;
}

void const
Filter::placeElements(unsigned int const startingIndex, unsigned int childrenCount,
        unsigned int const level)
{
    unsigned int counter = 0;
    unsigned int const feature_count = mpFeatureInformationMatrix->getRowCount();
    unsigned int* const p_candidate_feature_indices = new unsigned int[feature_count];
    unsigned int* const p_order = new unsigned int[feature_count];
    unsigned int* const p_adaptor = new unsigned int[feature_count];
    double* const p_candidate_scores = new double[feature_count];

    for (unsigned int i = 0; i < feature_count; ++i)
    {
        if (hasAncestorByFeatureIndex(startingIndex, i, level))
            continue;

        double ancestry_score = 0.;

        if (level > 1)
        {
            unsigned int ancestor_absolute_index = startingIndex;

            for (unsigned int j = level; j > 0; --j)
            {
                ancestor_absolute_index = getParentAbsoluteIndex(ancestor_absolute_index, j);

                double ancestry_score_ij = Math::computeMi(
                        mpFeatureInformationMatrix->at(i, mpIndexTree[ancestor_absolute_index]));
                double ancestry_score_ji = Math::computeMi(
                        mpFeatureInformationMatrix->at(mpIndexTree[ancestor_absolute_index], i));

                ancestry_score += std::max(ancestry_score_ij, ancestry_score_ji);
            }
        }

        double const score = Math::computeMi(mpFeatureInformationMatrix->at(i, mpIndexTree[0]))
                - (ancestry_score / level);

        if (score == score)
        {
            p_order[counter] = counter;
            p_adaptor[counter] = counter;
            p_candidate_feature_indices[counter] = i;
            p_candidate_scores[counter] = score;

            ++counter;
        }
    }

    std::sort(p_order, p_order + counter, Math::IndirectComparator(p_candidate_scores, p_adaptor));

#pragma omp critical(filter_element_placement)
    {
        unsigned int children_counter = 0;
        unsigned int i = counter - 1;

        // Sorry about this
        //       ||
        //       \/
        // i < counter depends on the fact that i-- causes an overflow on the unsigned int
        while (i < counter && children_counter < childrenCount)
        {
            unsigned int const index = p_candidate_feature_indices[p_order[i--]];

            if (!isRedundantPath(startingIndex + children_counter, index, level))
            {
                mpIndexTree[startingIndex + children_counter++] = index;
                mpScoreTree[startingIndex + children_counter-1] = p_candidate_scores[p_order[i+1]];
            }
        }
    }

    delete[] p_order;
    delete[] p_adaptor;
    delete[] p_candidate_feature_indices;
    delete[] p_candidate_scores;
}
