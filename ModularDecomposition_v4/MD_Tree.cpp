#include "MD_Tree.h"

/**
* Constructor implementations
*/
TreeNode::TreeNode(int val) : value(val), timestamp(-1), nChildNodes(0), nMarkedChildNodes(0), includeNode(false),
                              label(LEAF), markedLeft(false), markedRight(false), child(nullptr), sibling(nullptr), parent(nullptr) {}

TreeNode::TreeNode(Label l) : value(-1), timestamp(-1), nChildNodes(0), nMarkedChildNodes(0), label(l), includeNode(false),
                              markedLeft(false), markedRight(false), child(nullptr), sibling(nullptr), parent(nullptr) {}

MD_Tree::MD_Tree(TreeNode* r) : root(r), leftIndex(-1), rightIndex(-1), isLeftOfPivot(false),
                                left(nullptr), right(nullptr) {}

MD_Tree::MD_Tree() : root(nullptr), leftIndex(-1), rightIndex(-1), isLeftOfPivot(false),
                     left(nullptr), right(nullptr) {}


/**
* Operator implementations
*/

bool compareTreeNodePointers(const TreeNode* lhs, const TreeNode* rhs)
{
    if (lhs->value == rhs->value) {
        return getMaxLeafValue(lhs) < getMaxLeafValue(rhs);
    }
    return lhs->value < rhs->value;
}
bool operator==(const MD_Tree& lhs, const TreeNode* rhs) {
    return lhs.root == rhs;
}

bool operator!=(const MD_Tree& lhs, const TreeNode* rhs) {
    return !(lhs == rhs);
}

bool operator==(const TreeNode& lhs, const TreeNode& rhs) {
    return lhs.value == rhs.value && lhs.label == rhs.label && lhs.child == rhs.child
           && lhs.sibling == rhs.sibling && lhs.parent == rhs.parent;
}

/**
* Sets the child of a given TreeNode.
*
* @param lhs The parent
* @param rhs The child
*/
void setChild(TreeNode* lhs, TreeNode* rhs)
{
    if (lhs != nullptr) {
        lhs->child = rhs;
    }
    if (rhs != nullptr) {
        rhs->parent = lhs;
    }
}

/**
* Sets the sibling of a given TreeNode.
*
* @param lhs The "left" sibling.
* @param rhs The "right" sibling.
*/
void setSibling(TreeNode* lhs, TreeNode* rhs)
{
    if (lhs != nullptr) {
        lhs->sibling = rhs;
    }
    if (rhs != nullptr && lhs != nullptr) {
        rhs->parent = lhs->parent;
    }
}

/**
* Sets the neighbor of a given MD_Tree.
*
* @param lhs The "left" neighbor.
* @param rhs The "right" neighbor.
*/
void setNeighbor(MD_Tree* lhs, MD_Tree* rhs)
{
    if (lhs != nullptr) {
        lhs->right = rhs;
    }
    if (rhs != nullptr) {
        rhs->left = lhs;
    }
}

/**
* Recursively prints a MD_Tree on the command line
*
* @param node The starting node for the recursive process.
* @param depth The current depth of the tree (default: 0)
*/
void printTree(const TreeNode* node, int depth) {
    for (int i = 0; i < depth; ++i) {
        cout << "  ";
    }

    string mark = "";
    if (node->markedLeft && node->markedRight) {
        mark = ", Mark: Left & Right";
    }
    else if (node->markedLeft) {
        mark = ", Mark: Left";
    }
    else if (node->markedRight) {
        mark = ", Mark: Right";
    }

    if (node->label == LEAF) {
        cout << "Leaf: " << node->value << mark << endl;
    }
    else {
        cout << "Label: " << node->label << mark << endl;
    }

    if (node->child) {
        printTree(node->child, depth + 1);
    }
    if (depth > 0) {
        if (node->sibling) {
            printTree(node->sibling, depth);
        }
    }
}

string generateTreeString(const TreeNode* node, int depth) {
    string result;

    for (int i = 0; i < depth; ++i) {
        result += "  ";
    }

    string mark = "";
    if (node->markedLeft && node->markedRight) {
        mark = ", Mark: Left & Right";
    }
    else if (node->markedLeft) {
        mark = ", Mark: Left";
    }
    else if (node->markedRight) {
        mark = ", Mark: Right";
    }

    if (node->label == LEAF) {
        result += "Leaf: " + to_string(node->value) + mark + "\n";
    }
    else {
        result += "Label: " + to_string(node->label) + mark + "\n";
    }

    if (node->child) {
        result += generateTreeString(node->child, depth + 1);
    }
    if (depth > 0) {
        if (node->sibling) {
            result += generateTreeString(node->sibling, depth);
        }
    }

    return result;
}

/**
* A helper function to get a pre-Order of all leafs in the tree.
*
* @param node The current node.
* @param checkSiblings If the pre-Order of the nodes siblings should be returned.
* @return A pre-Order of the tree's leafs.
*/
vector<int> getPreOrderLeafsHelper(const TreeNode* node, bool checkSiblings) {
    vector<int> result;

    if (node == nullptr) {
        return result;
    }

    if (node->label == LEAF) {
        result.push_back(node->value);
    }

    auto childResult = getPreOrderLeafsHelper(node->child, true);
    result.insert(result.end(), childResult.begin(), childResult.end());

    if (checkSiblings) {
        auto siblingResult = getPreOrderLeafsHelper(node->sibling, true);
        result.insert(result.end(), siblingResult.begin(), siblingResult.end());
    }

    return result;
}

/**
* Return a pre-Order of the given tree's leafs.
*
* @param root The root node of the tree.
* @return A pre-Order of the tree's leafs.
*/
vector<int> getPreOrderLeafs(const TreeNode* root) {
    return getPreOrderLeafsHelper(root, false);
}

/**
* Resets all timestemps of the given tree to -1.
*
* @param root The root of the tree to reset the timestemps of.
*/
void resetTimestamps(TreeNode* node)
{
    node->timestamp = -1;
    if (node->child != nullptr) {
        resetTimestamps(node->child);
    }
    if (node->sibling != nullptr) {
        resetTimestamps(node->sibling);
    }
}

/**
* Returns the number of matching arguments of a vector and an unordered_set.
*
* @tparam T The type of elements in the vector and set.
* @param lhs The vector
* @param rhs The set
* @return The number of matching parameters.
*/
template <typename t>
int getNMatchingArguments(const vector<t>& lhs, const unordered_set<t>& rhs) {
    int count = 0;
    for (int element : lhs) {
        if (rhs.find(element) != rhs.end()) {
            count++;
        }
    }
    return count;
}

/**
* Calculates the maximal containing subtrees, based on a set X of leaf nodes. This means, that all nodes are returned,
* whose children can all be found in a set X. This cannot be true for the node's parent.
*
* @param treeNodes The set X of leaf nodes.
* @param currentTimestemp A timestemp that will be used for the calculation.
* @return The searched set of treeNodes
*/
unordered_set<TreeNode*> getMaxContSubTrees(const vector<TreeNode*>& treeNodes, int currentTimestamp) {

    unordered_set<TreeNode*> result;
    for (int i = 0; i < treeNodes.size(); i++) {
        TreeNode* currentNode = treeNodes[i];
        updateNodeInclusion(currentNode);
    }
    for (int i = 0; i < treeNodes.size(); i++) {
        TreeNode* currentNode = treeNodes[i];
        insertMaxContSubTrees(currentNode, result, currentTimestamp);
    }

    return result;
}

/**
* Updates a tree nodes inclusion information, based on the amount of its included children.
*
* @param node The tree node to update.
*/
void updateNodeInclusion(TreeNode* node)
{
    if (node->nChildNodes == node->nMarkedChildNodes) {
        node->includeNode = true;
        if (node->parent != nullptr) {
            node->parent->nMarkedChildNodes++;
            updateNodeInclusion(node->parent);
        }
    }
}

/**
* Inserts the correct node into a result - set, based on the requirements of the getMaxContSubTrees method.
* Starts at a leaf node and moves "upwards" in a tree.
*/
void insertMaxContSubTrees(TreeNode* node, unordered_set<TreeNode*>& result, int currentTimestamp)
{
    if (node->includeNode) {
        if (node->parent == nullptr) {
            result.insert(node);
        } else if (!node->parent->includeNode && node->parent->timestamp != currentTimestamp) {
            result.insert(node);
        }
        node->timestamp = currentTimestamp;
        node->includeNode = false;
        if (node->parent != nullptr && node->parent->nMarkedChildNodes > 0) {
            insertMaxContSubTrees(node->parent, result, currentTimestamp);
        }
    }

    node->nMarkedChildNodes = 0;
}

/**
 * Returns the maximum value of all leave descendents of the given node.
 *
 * @param node The given node
 * @return The maximum value of all leave descendants.
 */
int getMaxLeafValue(const TreeNode *node) {
    if (node->label == LEAF) {
        return node->value;
    }
    int maxLeafValue = -1;
    TreeNode* nextChild = node->child;
    while(nextChild != nullptr) {
        int currentValue = getMaxLeafValue(nextChild);
        if (currentValue > maxLeafValue) {
            maxLeafValue = currentValue;
        }
        nextChild = nextChild->sibling;
    }
    return maxLeafValue;
}
