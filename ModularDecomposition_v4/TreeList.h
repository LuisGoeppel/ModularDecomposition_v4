#pragma once
#include "MD_Tree.h"

class TreeList
{
    MD_Tree* firstElement;
    MD_Tree* lastElement;

public:
    void replaceElement(MD_Tree* toReplace, MD_Tree* leftNewElement, MD_Tree* rightNewElement);
    MD_Tree* getCorrespondingTree(const TreeNode* node) const;
    void insert(MD_Tree* tree);
    MD_Tree* getStart();
    void print();
};


