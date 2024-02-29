#include "Util.h"

/**
 * Sorts the given tree alphabetically.
 *
 * @param tree The tree to sort.
 */
void Util::sortTree(MD_Tree& tree) {
    sortTreeHelper(tree.root);
}

/**
 * Rewrites the provided adjacency list, relabeling nodes consecutively and returning
 * a vector of differences between ASCII values of old and new node names.
 * The function modifies the input adjacency list to maintain the graph structure while
 * assigning consecutive labels to nodes.
 *
 * @param adjList A reference to the input adjacency list (modified in the method).
 * @return A vector of integers representing differences between ASCII values of old and new node names.
 */
vector<int> Util::rewriteAdjacencyList(string& adjList) {
    vector<int> indexMapping;
    string updatedAdjlist;

    istringstream iss(adjList);
    vector<string> lines;
    string line;

    while (getline(iss, line)) {
        lines.push_back(line);
    }

    sort(lines.begin(), lines.end(), [](const string& a, const string& b) {
        return stoi(a) < stoi(b);
    });

    for (int i = 0; i < lines.size(); i++) {
        int lineNum = stoi(lines[i]);
        indexMapping.push_back(lineNum);
        updatedAdjlist += to_string(i) + ": " + lines[i].substr(lines[i].find(':') + 2) + "\n";
    }

    adjList = updatedAdjlist;
    return indexMapping;
}

/**
 * Checks if the provided string of node names is consecutively ordered from 'a' onwards.
 *
 * @param input The string of node names.
 * @return True if consecutively ordered, false otherwise.
 */
bool Util::isConsecutivelyOrdered(const string& input) {
    istringstream iss(input);
    string line;
    int prevNumber = -1;

    while (getline(iss, line)) {
        std::istringstream lineStream(line);
        int currentNumber;
        char colon;
        if (lineStream >> currentNumber >> colon) {
            if (currentNumber != prevNumber + 1) {
                return false;
            }
            prevNumber = currentNumber;
        } else {
            return false;
        }
    }
    return true;
}

/**
 * Tests if the provided Modular Decomposition tree is valid for the given graph.
 *
 * @param graph The graph.
 * @param tree The Modular Decomposition tree.
 * @return True if the tree is valid, false otherwise.
 */
bool Util::testModularDecompositionTree(const Graph& graph, const MD_Tree& tree, const vector<TreeNode*>& nodeValueMapping)
{
    vector<vector<int>> adjlist = graph.getAdjlist();
    vector<vector<int>> isConnectedByPrimeNode(adjlist.size(), vector<int>(adjlist.size(), 0));

    for (int lhs = 0; lhs < adjlist.size() - 1; lhs++) {
        for (int rhs = lhs + 1; rhs < adjlist.size(); rhs++) {

            bool hasConnectionGraph = false;
            for (int adj : adjlist[lhs]) {
                if (adj == rhs) {
                    hasConnectionGraph = true;
                    break;
                }
            }

            bool hasConnectionTree = hasConnectionGraph;
            Ancestors ancestors = getCommonAncestor(nodeValueMapping, lhs, rhs);
            Label commonAncestorLabel = ancestors.commonAncestor->label;

            if (commonAncestorLabel == PARALLEL) {
                hasConnectionTree = false;
            }
            else if (commonAncestorLabel == SERIES) {
                hasConnectionTree = true;
            }
            else {
                if (isConnectedByPrimeNode[lhs][rhs] != 0) {
                    hasConnectionTree = (isConnectedByPrimeNode[lhs][rhs] == 1);
                }
                else {
                    vector<int> valuesLhs = getPreOrderLeafs(ancestors.ancestorChildLhs);
                    vector<int> valuesRhs = getPreOrderLeafs(ancestors.ancestorChildRhs);
                    int matrixValue = hasConnectionGraph ? 1 : -1;
                    for (int i : valuesLhs) {
                        for (int j : valuesRhs) {
                            isConnectedByPrimeNode[i][j] = matrixValue;
                            isConnectedByPrimeNode[j][i] = matrixValue;
                        }
                    }
                }
            }

            if (hasConnectionGraph != hasConnectionTree) {
                cout << "Error! Nodes " << lhs << " and " << rhs  << " have different connections in the graph and the MD_Tree!" << endl;
                return false;
            }
        }
    }
    return true;
}

/**
 * Checks the child node values of the provided Modular Decomposition tree.
 *
 * @param tree The Modular Decomposition tree.
 */
void Util::checkChildNodeValues(MD_Tree& tree)
{
    checkChildNodeValuesHelper(tree.root);
}

/**
 * @brief A helper function for the sortTree algorithm. It takes a node, sorts its children
 * alphabetically, and then recursively calls itself for the node's children and siblings.
 *
 * @param node The current node, of which the children should get sorted.
 */
void Util::sortTreeHelper(TreeNode* node) {
    vector<TreeNode*> children;
    TreeNode* next = node->child;

    while (next != nullptr) {
        children.push_back(next);
        next = next->sibling;
    }

    if (!children.empty()) {
        sort(children.begin(), children.end(), compareTreeNodePointers);
        for (int i = 0; i < children.size(); i++) {
            if (i == 0) {
                setChild(node, children[i]);
            }
            else {
                setSibling(children[i - 1], children[i]);
            }
        }
        children[children.size() - 1]->sibling = nullptr;
    }

    if (node->sibling != nullptr) {
        sortTreeHelper(node->sibling);
    }
    if (node->child != nullptr) {
        sortTreeHelper(node->child);
    }
}

/**
 * Gets the inverse of the provided vector.
 *
 * @param vec The input vector.
 * @return The inverse vector.
 */
vector<int> Util::getInverse(const vector<int>& vec)
{
    int maxElement = 0;
    for (int element : vec) {
        if (element > maxElement) {
            maxElement = element;
        }
    }
    vector<int> inverseVector(maxElement + 1, -1);
    for (int i = 0; i < vec.size(); i++) {
        inverseVector[vec[i]] = i;
    }
    return inverseVector;
}

/**
* Gets the common ancestor child node for the given nodes in the Modular Decomposition tree.
*
* @param tree The Modular Decomposition tree.
* @param lhs The left - hand side node.
* @param rhs The right - hand side node.
* @return The common ancestor child node.
*/
TreeNode* Util::getCommonAncestorChildLhs(const MD_Tree& tree, int lhs, int rhs)
{
    TreeNode* lhsNode = getCorrespondingTreeNode(tree.root, lhs);
    TreeNode* currentAncestor = lhsNode->parent;
    TreeNode* ancestorChildLhs = lhsNode;

    while (currentAncestor != nullptr) {
        vector<int> currentLeafNodes = getPreOrderLeafs(currentAncestor);
        if (find(currentLeafNodes.begin(), currentLeafNodes.end(), rhs) != currentLeafNodes.end()) {
            return ancestorChildLhs;
        }
        ancestorChildLhs = currentAncestor;
        currentAncestor = currentAncestor->parent;
    }
    return nullptr;
}

/**
 * Gets the corresponding tree node for the given node in the Modular Decomposition tree.
 *
 * @param currentNode The current node to search from.
 * @param node The target node.
 * @return The corresponding tree node.
 */
TreeNode* Util::getCorrespondingTreeNode(TreeNode* currentNode, int node)
{
    if (currentNode->value == node) {
        return currentNode;
    }
    TreeNode* out = nullptr;
    if (currentNode->child != nullptr) {
        out = getCorrespondingTreeNode(currentNode->child, node);
    }
    if (currentNode->sibling != nullptr && out == nullptr) {
        out = getCorrespondingTreeNode(currentNode->sibling, node);
    }
    return out;
}

/**
 * A helper function to check child node values of the Modular Decomposition tree.
 *
 * @param node The current node.
 */
void Util::checkChildNodeValuesHelper(TreeNode* node)
{
    TreeNode* currentChild = node->child;
    int childCount = 0;

    while (currentChild != nullptr) {
        childCount++;
        currentChild = currentChild->sibling;
    }

    if (childCount != node->nChildNodes) {
        cout << "Error at child count: Expected " << node->nChildNodes
             << " but found " << childCount << " nodes!";
    }

    if (node->sibling != nullptr) {
        checkChildNodeValuesHelper(node->sibling);
    }
    if (node->child != nullptr) {
        checkChildNodeValuesHelper(node->child);
    }
}

/**
 * Returns the number of vertices for a given graph.
 *
 * @param graph The given graph.
 * @return The number of vertices.
 */
int Util::getNumberVertices(const Graph &graph) {
    return static_cast<int>(graph.getAdjlist().size());
}

/**
 * Returns the number of edges for a given graph.
 *
 * @param graph The given graph.
 * @return The number of edges.
 */
int Util::getNumberEdges(const Graph &graph) {
    int count = 0;
    for (const auto & vertex : graph.getAdjlist()) {
        count += static_cast<int>(vertex.size());
    }
    return count/2;
}

/**
 * Generates a random modular decomposition tree.
 *
 * @param nVertices The number of vertices the random tree should have.
 * @param isCoGraph If the tree should be a cograph, i.e. have only PARALLEL and SERIES nodes
 * @return A randomly generated modular decomposition tree.
 */
MD_Tree Util::createRandomModularDecompositionTree(int nVertices, bool isCoGraph, vector<TreeNode*>& nodeValueMapping) {

    random_device rd;
    mt19937 gen(rd());
    int upperBoundNodeTypes = isCoGraph ? 1 : 2;
    int maxAmountChildren = max(5, nVertices/5);
    uniform_int_distribution<int> randomNodeType(0, upperBoundNodeTypes);
    uniform_int_distribution<int> randomRootNodeType(1, upperBoundNodeTypes);
    uniform_int_distribution<int> randomAmountChildNodes(2, maxAmountChildren);
    uniform_int_distribution<int> randomBool(0, 1);
    int vertexCount;
    TreeNode* rootNode;

    do {
        vertexCount = 0;
        vector<int> remainingVertices(nVertices);
        iota(remainingVertices.begin(), remainingVertices.end(), 0);
        queue<TreeNode*> remainingInnerNodes;
        rootNode = new TreeNode(getLabel(randomRootNodeType(gen)));
        remainingInnerNodes.push(rootNode);

        while (!remainingInnerNodes.empty()) {
            TreeNode* currentNode = remainingInnerNodes.front();
            remainingInnerNodes.pop();

            int nChildren = min(randomAmountChildNodes(gen), static_cast<int>(remainingVertices.size()));
            if (nChildren <= 3 && currentNode->label == PRIME) {
                currentNode->label = SERIES;
            }
            currentNode->nChildNodes = nChildren;
            TreeNode* currentChild;
            TreeNode* lastChild;
            for (int i = 0; i < nChildren; i++) {
                int isLeafNode = remainingVertices.size() >= 2 ? randomBool(gen) : 1;
                if (isLeafNode == 1) {
                    if (!remainingVertices.empty()) {
                        uniform_int_distribution<int> randomVertex(0, static_cast<int>(remainingVertices.size()) - 1);
                        int vertexIndex = randomVertex(gen);
                        int vertex = remainingVertices[vertexIndex];
                        remainingVertices.erase(remainingVertices.begin() + vertexIndex);
                        vertexCount++;
                        currentChild = new TreeNode(vertex);
                        nodeValueMapping[vertex] = currentChild;
                    }
                } else {
                    Label currentLabel = getLabel(randomNodeType(gen));
                    if (currentLabel == currentNode->label) {
                        currentLabel = generateDifferentLabel(currentNode->label, isCoGraph);
                    }
                    currentChild = new TreeNode(currentLabel);
                    remainingInnerNodes.push(currentChild);
                }
                if (i == 0) {
                    setChild(currentNode, currentChild);
                } else {
                    setSibling(lastChild, currentChild);
                }
                lastChild = currentChild;
            }
        }
        remainingVertices.clear();
    } while (vertexCount != nVertices);

    int res;
    do {
        res = removeFalseInnerNodes(rootNode, nullptr, false);
        updateParentPointersAndModuleTypes(rootNode, nullptr, false);
    } while (res == 1);

    if (res == -1) {
        return createRandomModularDecompositionTree(nVertices, isCoGraph, nodeValueMapping);
    }

    auto* result = new MD_Tree(rootNode);
    if (isValidModularDecompositionTree(*result)) {
        return *result;
    } else {
        return createRandomModularDecompositionTree(nVertices, isCoGraph, nodeValueMapping);
    }
}

/**
 * Creates a graph based on a modular decomposition tree. It is assumed, that in all PRIME modules the
 * children are connected via a simple path from the left-most to the right-most element.
 *
 * @param tree The tree to create the graph for.
 * @return The created graph.
 */
Graph Util::createGraphFromTree(const MD_Tree& tree, const vector<TreeNode*>& nodeValueMapping) {
    int nVertices = 0;
    getHighestVertexValue(tree.root, nVertices);
    nVertices++;
    vector<vector<int>> adjList(nVertices);

    for (int lhs = 0; lhs < nVertices - 1; lhs++) {
        for (int rhs = lhs + 1; rhs < nVertices; rhs++) {
            Ancestors ancestors = getCommonAncestor(nodeValueMapping, lhs, rhs);

            if (ancestors.commonAncestor->label == PRIME) {
                if (abs(ancestors.ancestorIndexLhs - ancestors.ancestorIndexRhs) == 1) {
                    adjList[lhs].push_back(rhs);
                    adjList[rhs].push_back(lhs);
                }
            } else if (ancestors.commonAncestor -> label == SERIES) {
                adjList[lhs].push_back(rhs);
                adjList[rhs].push_back(lhs);
            }
        }
    }
    return {adjList};
}

/**
 * Converts a number to the corresponding TreeNode-Label.
 *
 * @param n The number.
 * @return The label.
 */
Label Util::getLabel(int n) {
    switch (n) {
        case 0:
            return PARALLEL;
        case 1:
            return SERIES;
        case 2:
            return PRIME;
        default:
            return static_cast<Label>(NULL);
    }
}

/**
 * Removes inner nodes that have no children or only one child from a given modular decomposition tree.
 * Adjusts the label of an inner node if needed, as a PRIME node is not possible with only two children.
 *
 * @param currentNode The current root node of the tree.
 * @param callingNode The node that called this procedure recursively.
 * @param isParent If the calling node is the parent of the called node.
 * @return -1 if the tree is invalid, 0 if no node had to be removed and 1 if at least one node was removed.
 */
int Util::removeFalseInnerNodes(TreeNode* currentNode, TreeNode* callingNode, bool isParent) {
    if (currentNode == currentNode->sibling) {
        return -1;
    }

    bool removedNode = false;
    if (currentNode->sibling != nullptr) {
        removedNode |= removeFalseInnerNodes(currentNode->sibling, currentNode, false);
    }

    if (currentNode->label != LEAF) {
        int nChildren = 0;
        TreeNode* nextChild = currentNode->child;
        while(nextChild != nullptr) {
            nChildren++;
            nextChild = nextChild->sibling;
        }
        if (nChildren == 0) {
            if (isParent) {
                setChild(callingNode, currentNode->sibling);
            } else {
                setSibling(callingNode, currentNode->sibling);
            }
            return 1;
        } else if (nChildren == 1) {
            if (isParent) {
                setChild(callingNode, currentNode->child);
                setSibling(currentNode->child, currentNode->sibling);
            } else {
                setSibling(callingNode, currentNode->child);
                setSibling(currentNode->child, currentNode->sibling);
            }
            if (currentNode->parent != nullptr && currentNode->child->label == currentNode->parent->label) {
                currentNode->child->label = generateDifferentLabel(currentNode->parent->label, true);
            }
            return 1;
        } else if (nChildren < 4 && currentNode->label == PRIME) {
            currentNode->label = generateDifferentLabel(currentNode->label, false);
        }
    }

    if (currentNode->child != nullptr) {
        removedNode |= removeFalseInnerNodes(currentNode->child, currentNode, true);
    }
    return removedNode ? 1 : 0;
}

/**
 * Returns the highest leaf node value in the current subtree.
 *
 * @param node The node inducing the current subtree.
 * @param highestValue The currently highest value.
 */
void Util::getHighestVertexValue(TreeNode* node, int& highestValue) {
    if (node->label == LEAF) {
        if (node->value > highestValue) {
            highestValue = node->value;
        }
    }

    if (node->sibling != nullptr) {
        getHighestVertexValue(node->sibling, highestValue);
    }
    if (node->child != nullptr) {
        getHighestVertexValue(node->child, highestValue);
    }
}

/**
 * Generates a TreeNode Label that is different than the given one.
 *
 * @param label The given label.
 * @param isCoGraph If the current graph is a cograph.
 * @return A label different to the given one.
 */
Label Util::generateDifferentLabel(const Label &label, bool isCoGraph) {
    if (isCoGraph) {
        if (label == PARALLEL) {
            return SERIES;
        } else {
            return PARALLEL;
        }
    } else {
        if (label == PARALLEL) {
            return SERIES;
        } else if (label == SERIES) {
            return PRIME;
        } else {
            return PARALLEL;
        }
    }
}

/**
 * Computes if a given modular decomposition tree contains the correct amount of leave nodes and
 * has no parent and child node with the same label.
 *
 * @param tree The modular decomposition tree to check.
 * @return If the given tree is valid.
 */
bool Util::isValidModularDecompositionTree(const MD_Tree &tree) {
    unordered_set<int> leaves;
    return isValidTreeHelper(tree.root, leaves);
}

/**
 * A recursive helper function for the validTree calculation.
 *
 * @param node The current node in the tree.
 * @param leaves A set containing all leaves that have already been found.
 * @return If the current subtree is valid.
 */
bool Util::isValidTreeHelper(TreeNode* node, unordered_set<int>& leaves) {
    if (node->label == LEAF) {
        if (leaves.find(node->value) == leaves.end()) {
            leaves.insert(node->value);
            return true;
        } else {
            return false;
        }
    }
    if (node->parent != nullptr) {
        if (node->parent->label == node->label) {
            return false;
        }
    }

    if (node->sibling != nullptr) {
        if (!isValidTreeHelper(node->sibling, leaves)) {
            return false;
        }
    }
    if (node->child != nullptr) {
        if (!isValidTreeHelper(node->child, leaves)) {
            return false;
        }
    }
    return true;
}

/**
 * Updates the parent pointers of the given modular decomposition tree and adjusts the prime module types if there
 * are two or more PRIME modules that are children of the same node.
 *
 * @param currentNode The current node in the tree.
 * @param parentNode The parent node of the current node.
 * @param primeModule If there is already one PRIME module at the current level.
 */
void Util::updateParentPointersAndModuleTypes(TreeNode* currentNode, TreeNode* parentNode, bool primeModule) {
    if (parentNode != nullptr) {
        currentNode->parent = parentNode;
        if (currentNode->parent->label == currentNode->label) {
            currentNode->label = generateDifferentLabel(currentNode->label, true);
        }
    }

    if (currentNode->label == PRIME) {
        if (primeModule) {
            currentNode->label = generateDifferentLabel(currentNode->parent->label, true);
        }
        primeModule = true;
    }

    if (currentNode -> sibling != nullptr) {
        updateParentPointersAndModuleTypes(currentNode->sibling, parentNode, primeModule);
    }
    if (currentNode -> child != nullptr) {
        updateParentPointersAndModuleTypes(currentNode->child, currentNode, false);
    }
}

/**
 * Returns the common ancestor of two leaf-nodes in the given tree, represented by their values.
 *
 * @param nodeValueMapping The nodeValueMapping.
 * @param lhs The value representing the first node.
 * @param rhs The value representing the second node.
 * @return The common ancestor of the given node, as well as information about which of its children
 *      has the given nodes as descendants.
 */
Ancestors Util::getCommonAncestor(const vector<TreeNode*>& nodeValueMapping, int lhs, int rhs) {
    int timestampValueLhs = numeric_limits<int>::max() / 2;
    int timestampValueRhs = numeric_limits<int>::min() / 2;

    TreeNode* nodeLhs = nodeValueMapping[lhs];
    TreeNode* nodeRhs = nodeValueMapping[rhs];
    TreeNode* commonAncestor;
    TreeNode* ancestorChildLhs;
    TreeNode* ancestorChildRhs;
    int ancestorIndexLhs;
    int ancestorIndexRhs;

    TreeNode* currentLhs = nodeLhs;
    while (currentLhs != nullptr) {
        currentLhs->timestamp = timestampValueLhs;
        currentLhs = currentLhs->parent;
    }

    TreeNode* currentRhs = nodeRhs;
    while (currentRhs != nullptr) {
        if (currentRhs->timestamp == timestampValueLhs) {
            commonAncestor = currentRhs;
            currentRhs = nullptr;
        } else {
            currentRhs->timestamp = timestampValueRhs;
            currentRhs = currentRhs->parent;
        }
    }

    TreeNode* currentChild = commonAncestor->child;
    int currentIndex = 0;
    while (currentChild != nullptr) {
        if (currentChild->timestamp == timestampValueLhs) {
            ancestorIndexLhs = currentIndex;
            ancestorChildLhs = currentChild;
        }
        if (currentChild->timestamp == timestampValueRhs) {
            ancestorIndexRhs = currentIndex;
            ancestorChildRhs = currentChild;
        }
        currentChild = currentChild->sibling;
        currentIndex++;
    }
    resetTimestamps(commonAncestor);
    resetParentTimestamps(commonAncestor);

    return {commonAncestor, ancestorChildLhs, ancestorChildRhs, ancestorIndexLhs, ancestorIndexRhs};
}

/**
* Resets all timestamps of the given tree to -1.
*
* @param root The root of the tree to reset the timestamps of.
*/
void Util::resetTimestamps(TreeNode* node)
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
 * Resets all timestamps of all ancestors of the given node.
 *
 * @param node The given node.
 */
void Util::resetParentTimestamps(TreeNode *node) {
    node->timestamp = -1;
    if (node->parent != nullptr) {
        resetParentTimestamps(node->parent);
    }
}


