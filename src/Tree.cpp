#include "Tree.hpp"

Tree::Tree(unsigned int* const pChildrenCountPerLevel, unsigned int const levelCount,
        float* const pFeatureInformationMatrix, unsigned int const featureCount,
        unsigned int const targetFeatureIndex) :
        mpChildrenCountPerLevel(pChildrenCountPerLevel), mLevelCount(levelCount), mpFeatureInformationMatrix(
                pFeatureInformationMatrix), mFeatureCount(featureCount)
{
    unsigned int cumulative_element_count = 1;
    unsigned int children_per_level = 1;

    mpStartingIndexPerLevel = new unsigned int[mLevelCount + 1];
    mpStartingIndexPerLevel[0] = 0;

    for (unsigned int level = 0; level < mLevelCount; ++level)
    {
        mpStartingIndexPerLevel[level + 1] = cumulative_element_count;
        children_per_level *= mpChildrenCountPerLevel[level];
        cumulative_element_count += children_per_level;
    }

    mpIndexTree = new unsigned int[cumulative_element_count];
    mpInformativeContributionTree = new float[cumulative_element_count];
    mpRedundantContributionTree = new float[cumulative_element_count];
    mTreeElementCount = cumulative_element_count;

    mpIndexTree[0] = targetFeatureIndex;
    mpInformativeContributionTree[0] = 0.;
    mpRedundantContributionTree[0] = 0.;
}

Tree::~Tree()
{
    delete[] mpStartingIndexPerLevel;
    delete[] mpIndexTree;
    delete[] mpInformativeContributionTree;
    delete[] mpRedundantContributionTree;
}

void
Tree::build()
{
    for (unsigned int level = 0; level < mLevelCount; ++level)
    {
        unsigned int const parent_count = mpStartingIndexPerLevel[level + 1]
                - mpStartingIndexPerLevel[level];

#pragma omp parallel for schedule(dynamic)
        for (unsigned int parent = 0; parent < parent_count; ++parent)
        {
            for (unsigned int child = 0; child < mpChildrenCountPerLevel[level]; ++child)
            {
                unsigned int const child_absolute_index = mpStartingIndexPerLevel[level + 1]
                        + (parent * mpChildrenCountPerLevel[level]) + child;

                // Selection criterion -> no two paths contain the same index set

                mpIndexTree[child_absolute_index] = mFeatureCount;// selectBestFeature(child_absolute_index, level+1); // OK
                // mpInformativeContributionTree[child_index]
                // mpRedundantContributionTree[child_index]
            }
        }
    }
}

// TODO: Add redundancy checking
//       Check if ancestry_score_sum is properly computed
//
unsigned int const
Tree::selectBestFeature(unsigned int absoluteIndex, unsigned int level)
{
    unsigned int max_candidate_index = 0;
    float max_candidate_value = -std::numeric_limits<double>::max();
    float ancestry_score_sum = 0;
    float max_candidate_ancestry_score_sum = 0;
    float max_candidate_target_mi = 0;

    for (unsigned int i = 0; i < mFeatureCount; ++i)
    {

        if (hasAncestorByIndex(absoluteIndex, i, level)
                || hasBrotherByIndex(absoluteIndex, i, level))
            continue;

        float target_score = mpFeatureInformationMatrix[mpIndexTree[0] * mFeatureCount + i];
        //float target_score = (*mpFeatureInformationMatrix)(mpIndexTree[0], i);

        // Compute the average redundancy of i with ancestry
        unsigned int absoluteAncestorIndex = absoluteIndex;
        if (level > 1)
            for (unsigned int j = level; j > 0; --j)
            {
                absoluteAncestorIndex = getParentAbsoluteIndex(absoluteAncestorIndex, j);
                ancestry_score_sum += mpFeatureInformationMatrix[absoluteAncestorIndex
                        * mFeatureCount + i];
                //ancestry_score_sum += (*mpFeatureInformationMatrix)(absoluteAncestorIndex, i);
            }

        float ancestry_mean_score = ancestry_score_sum / level;
        float score = target_score - ancestry_mean_score;

        // Pick the candidate only if it is a local maximum
        if (score > max_candidate_value)
        {
            max_candidate_index = i;
            max_candidate_value = score;
            max_candidate_ancestry_score_sum = ancestry_score_sum;
            max_candidate_target_mi = target_score;
        }
    }
    return max_candidate_index;
}

// TODO: Confirm that the method works
bool
Tree::isRedundantSolution(unsigned int absoluteIndex, unsigned int level)
{
    float target_score = computeQualityScore(absoluteIndex, level);
    for (unsigned int i = mpStartingIndexPerLevel[level]; i < absoluteIndex; ++i)
    {
        if (computeQualityScore(i, level) == target_score)
            if (hasSameAncestry(absoluteIndex, i, level))
                return true;
    }
    return false;
}

// TODO: Confirm that the method works
bool
Tree::hasSameAncestry(unsigned int targetAbsoluteIndex, unsigned int consideredAbsoluteIndex,
        unsigned int level)
{

    for (unsigned int i = 0; i < level; ++i)
    {
        targetAbsoluteIndex = getParentAbsoluteIndex(targetAbsoluteIndex, level);
        if (!hasAncestorByIndex(consideredAbsoluteIndex, mpIndexTree[targetAbsoluteIndex], level))
            return false;
    }
    return true;
}

// TODO: Confirm that the method works (divide by level -1 ?)
float
Tree::computeQualityScore(unsigned int absoluteIndex, unsigned int level)
{
    return mpInformativeContributionTree[absoluteIndex]
            - 2 * mpRedundantContributionTree[absoluteIndex] / level;
}

// TODO: Confirm that the method works
bool
Tree::hasAncestorByIndex(unsigned int absoluteIndex, unsigned int targetMimIndex,
        unsigned int level)
{
    for (unsigned int i = level; i > 0; --i)
    {
        absoluteIndex = getParentAbsoluteIndex(absoluteIndex, i);
        if (mpIndexTree[absoluteIndex] == targetMimIndex)
            return true;
    }
    return false;
}

// TODO: Confirm that the method works
bool
Tree::hasBrotherByIndex(unsigned int absoluteIndex, unsigned int targetMimIndex, unsigned int level)
{
    //  [2][3][4][5]...[6][7][8][9]...[10]
    //  [0][1][2][3]...[0][1][2][3]...[0]
    //  8: (8-2)%4 = 2
    //  9: (9-2)%4 = 3
    // 10: (10-2)%4 = 0
    unsigned int child_number = (absoluteIndex - mpStartingIndexPerLevel[level])
            % mpChildrenCountPerLevel[level];
    for (unsigned int i = absoluteIndex - child_number; i < absoluteIndex; ++i)
        if (mpIndexTree[i] == targetMimIndex)
            return true;
    return false;
}

void
Tree::getPaths(std::vector<unsigned int>* pPaths) const
{
    pPaths->reserve(mLevelCount * (mTreeElementCount - mpStartingIndexPerLevel[mLevelCount]));

    for (unsigned int end_element_absolute_index = mTreeElementCount - 1;
            end_element_absolute_index >= mpStartingIndexPerLevel[mLevelCount];
            --end_element_absolute_index)
    {
        unsigned int element_absolute_index = end_element_absolute_index;

        for (unsigned int level = mLevelCount; level > 0; --level)
        {
            pPaths->push_back(mpIndexTree[element_absolute_index]);
            element_absolute_index = getParentAbsoluteIndex(element_absolute_index, level);
        }
    }
}

/* inline */unsigned int const
Tree::getParentAbsoluteIndex(unsigned int const absoluteIndex, unsigned int const level) const
{
    return (absoluteIndex - mpStartingIndexPerLevel[level]) / mpChildrenCountPerLevel[level - 1]
            + mpStartingIndexPerLevel[level - 1];
}
