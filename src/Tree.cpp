#include "Tree.hpp"

Tree::Tree(unsigned int* const pChildrenCountPerLevel, unsigned int const levelCount) :
        mpChildrenCountPerLevel(pChildrenCountPerLevel), mLevelCount(levelCount)
{
    unsigned int cumulative_element_count = 1;
    unsigned int children_per_level = 1;

    for (unsigned int i = 0; i < mLevelCount; ++i)
    {
        children_per_level *= mpChildrenCountPerLevel[i];
        cumulative_element_count += children_per_level;
    }

    mpData = new unsigned int[cumulative_element_count];

    std::cout << cumulative_element_count << std::endl;
}

Tree::~Tree()
{
    delete[] mpData;
}

void
Tree::build()
{
    mpData[0] = 0;

    unsigned int last_level_element_count = 1;
    unsigned int this_level_starting_index = 1;

    for (unsigned int level = 1; level < mLevelCount; ++level)
    {
        for (unsigned int parent = 0; parent < last_level_element_count; ++parent)
        {
            for (unsigned int child = 0; child < mpChildrenCountPerLevel[level - 1]; ++child)
            {
                // create/select/validate child
                mpData[last_level_ending_index + (parent * mpChildrenCountPerLevel[level - 1])
                        + child] = child;
            }
        }

        last_level_element_count *= mpChildrenCountPerLevel[level - 1];
        this_level_starting_index += last_level_element_count;
    }
}
