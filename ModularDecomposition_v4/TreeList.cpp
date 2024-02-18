#include "TreeList.h"


/**
* Replaces an element in the tree list with two new elements.
*
* @param toReplace The element to replace.
* @param leftNewElement The first new element.
* @param rightNewElement The second new element.
*/
void TreeList::replaceElement(MD_Tree* toReplace, MD_Tree* leftNewElement, MD_Tree* rightNewElement)
{
    setNeighbor(toReplace->left, leftNewElement);
    setNeighbor(leftNewElement, rightNewElement);
    setNeighbor(rightNewElement, toReplace->right);
    if (toReplace == firstElement) {
        firstElement = leftNewElement;
    }
    //delete toReplace;
}

/**
* Returns the MD_Tree in which a specific TreeNode can be found (in the forest)
*
* @param root The root MD_Tree, representing the entry point of the tree list.
* @param node The node, of which's tree the index is returned.
* @return The searched index.
*/
MD_Tree* TreeList::getCorrespondingTree(const TreeNode* node) const {
    if (node == nullptr) {
        cerr << "getRootIndex Was Called From Nullpointer" << endl;
    }
    else {
        while (node->parent != nullptr) {
            node = node->parent;
        }

        MD_Tree* currentTree = firstElement;

        while (currentTree != nullptr) {
            if (currentTree->root == node) {
                return currentTree;
            }
            currentTree = currentTree->right;
        }
    }
    return nullptr;
}

/**
* Inserts an MD_Tree at the end of the list of MD_Trees.
*
* @param tree The element to insert.
*/
void TreeList::insert(MD_Tree* tree)
{
    if (firstElement == nullptr) {
        firstElement = tree;
        lastElement = tree;
    }
    else {
        setNeighbor(lastElement, tree);
        lastElement = tree;
    }
}

/**
* Returns the start of the TreeList.
*/
MD_Tree* TreeList::getStart()
{
    return firstElement;
}

/**
* Prints the current TreeList to the command line.
*/
void TreeList::print()
{
    int treeCount = 1;
    MD_Tree* current = firstElement;
    while (current != nullptr) {
        cout << "MD Tree " << treeCount << ": " << endl;
        TreeNode* treeRoot = current->root;
        printTree(treeRoot);
        cout << endl;

        current = current->right;
        treeCount++;
    }
}
