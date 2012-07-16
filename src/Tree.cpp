#include "Tree.hpp"

Tree::Tree(unsigned int* const pChildrenCountPerLevel, unsigned int const levelCount) :
        mpChildrenCountPerLevel(pChildrenCountPerLevel), mLevelCount(levelCount)
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
    mpIndexTree[0] = 0;

    for (unsigned int level = 0; level < mLevelCount; ++level)
    {
        unsigned int const parent_count = mpStartingIndexPerLevel[level + 1]
                - mpStartingIndexPerLevel[level];
        for (unsigned int parent = 0; parent < parent_count; ++parent)
        {
            for (unsigned int child = 0; child < mpChildrenCountPerLevel[level]; ++child)
            {
                unsigned int const child_index = mpStartingIndexPerLevel[level + 1]
                        + (parent * mpChildrenCountPerLevel[level]) + child;
                // create/select/validate child
                mpIndexTree[child_index] = child;
                // mpInformativeContributionTree[child_index]
                // mpRedundantContributionTree[child_index]
            }
        }
    }

    for (unsigned int i = 0; i < 33; ++i)
    {
        std::cout << mpIndexTree[i] << "\t";
    }
    std::cout << std::endl;
}
