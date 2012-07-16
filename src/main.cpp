#include "Tree.hpp"

int
main()
{
    unsigned int levels[] = { 2, 3, 2, 1 };
    Tree tree(levels, 4);
    tree.build();
    return 0;
}
