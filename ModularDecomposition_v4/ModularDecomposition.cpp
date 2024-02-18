#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <memory>
#include <queue>
#include <chrono>

#include "Graph.h"
#include "MD_Tree.h"
#include "Util.h"
#include "TreeList.h"


using namespace std;

MD_Tree getModularDecomposition(const Graph& graph);
MD_Tree getModularDecomposition(const Graph& graph, vector<TreeNode*>& nodeValueMapping);

/**
 * Returns the contents of a file as a string.
 *
 * @param filename The name of the file to be read.
 * @return A string containing the contents of the file, or an empty string if an error occurred.
 */
string readFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return "";
    }
    ostringstream oss;
    oss << file.rdbuf();
    file.close();

    return oss.str();
}


/**
* Prints a given modular - decomposition forest on the command line
*
* @param current The forest to be printed, passed by a vector of its trees
*/
void printForest(const vector<MD_Tree>& forest) {
    for (int i = 0; i < forest.size(); i++) {
        cout << "MD Tree " << (i + 1) << ": " << endl;
        printTree(forest[i].root);
        cout << endl;
    }
}


/**
* Marks a node and all its ancestor with either "left" or "right", based on the given parameter.
*
* @param node The node that should be marked (along with its ancestors).
* @param markLeft If the node should be marked with "left" ("right" otherwise).
*/
void markNodeAndAncestors(TreeNode* node, bool markLeft) {
    TreeNode* ancestor = node;
    if (markLeft) {
        while (ancestor != nullptr && !ancestor->markedLeft) {
            ancestor->markedLeft = true;
            ancestor = ancestor->parent;
        }
    }
    else {
        while (ancestor != nullptr && !ancestor->markedRight) {
            ancestor->markedRight = true;
            ancestor = ancestor->parent;
        }
    }
}

/**
* Marks all children of a given node with either "left" or "right", based on the given parameter.
*
* @param node The node which's children should be marked (along with its ancestors).
* @param markLeft If the node's children should be marked with "left" ("right" otherwise).
*/
void markChildren(TreeNode* node, bool markLeft) {
    TreeNode* next = node->child;
    while (next != nullptr) {
        if (markLeft) {
            next->markedLeft = true;
        }
        else {
            next->markedRight = true;
        }
        next = next->sibling;
    }
}

/**
* Constructs an MD_Tree, where the root has a given label and the given children.
*
* @param X A set containing all children that the root of the new tree should have.
* @param label The label of the new tree's root.
* @return The newly constructed MD_Tree
*/
MD_Tree constructTree(vector<TreeNode*> X, Label label) {
    TreeNode* T;
    if (X.size() == 1) {
        T = X[0];
        T->sibling = nullptr;
        T->parent = nullptr;
    }
    else {
        T = new TreeNode(label);
        for (int i = 0; i < X.size(); i++) {
            if (i == 0) {
                setChild(T, X[i]);
            }
            else {
                setSibling(X[i - 1], X[i]);
            }
        }
        X[X.size() - 1]->sibling = nullptr;
        T->nChildNodes = static_cast<int>(X.size());
    }
    return MD_Tree(T);
}

/**
 * Updates additional tree information like the placement in respect to the pivot or the markings of the root
 * of one tree based on this information on another tree.
 *
 * @param toUpdate The tree that has to be updated.
 * @param toCopy The tree to get the update information from.
 */
void updateTreeInformation(MD_Tree* toUpdate, MD_Tree* toCopy) {
    toUpdate->isLeftOfPivot = toCopy->isLeftOfPivot;
    toUpdate->root->markedLeft = toCopy->root->markedLeft;
    toUpdate->root->markedRight = toCopy->root->markedRight;
}

/**
* Refines a MD_Forest with a node that is not prime. For details on the refinement process,
* see the algorithm description.
*
* @param forest The forest that should be refined, passed as a list of MD_Trees.
* @param p The non-Prime tree node on base of which the forest should be refined.
* @param maxSubTrees A set containing all maximal Subtrees.
* @param isLeftSplit If a left-split should be used (right-split otherwise).
*/
void refineByNonPrimeNode(TreeList& forest, TreeNode* p,
                          const unordered_set<TreeNode*>& maxSubTrees, bool isLeftSplit) {

    vector<TreeNode*> A;
    vector<TreeNode*> B;

    TreeNode* next = p->child;
    while (next != nullptr) {
        if (maxSubTrees.find(next) != maxSubTrees.end()) {
            A.push_back(next);
        }
        else {
            B.push_back(next);
        }
        next = next->sibling;
    }

    if (A.size() > 0 && B.size() > 0) {
        MD_Tree* Ta = new MD_Tree(constructTree(A, p->label));
        MD_Tree* Tb = new MD_Tree(constructTree(B, p->label));

        if (p->parent == nullptr) {
            MD_Tree* currentTree = forest.getCorrespondingTree(p);
            updateTreeInformation(Ta, currentTree);
            updateTreeInformation(Tb, currentTree);

            if (isLeftSplit) {
                forest.replaceElement(currentTree, Ta, Tb);
            }
            else {
                forest.replaceElement(currentTree, Tb, Ta);
            }
        }
        else {
            setChild(p, Ta->root);
            setSibling(Ta->root, Tb->root);
            p->nChildNodes = 2;
        }
        markNodeAndAncestors(Ta->root, isLeftSplit);
        markNodeAndAncestors(Tb->root, isLeftSplit);
    }
}


/**
* Refines a MD_Forest with a node that is prime. For details on the refinement process,
* see the algorithm description.
*
* @param forest The forest that should be refined.
* @param isLeftSplit If a left-split should be used (right-split otherwise).
*/
void refineByPrimeNode(TreeNode* node, bool isLeftSplit, unordered_set<TreeNode*>& maxSubtrees) {
    markNodeAndAncestors(node, isLeftSplit);
    markChildren(node, isLeftSplit);
}

/**
* Refines a MD_Forest by a given set of active edges. For details on the refinement process see
* the algorithm description.
*
* @param forest The MD_Forest to be refined, passed by a list of MD_Trees.
* @param treeNodes The set of tree nodes that correspond to the active
    edges of the current refinement process.
* @param timestamp A timestamp that is used for the refinement process.
* @param nodeIsLeft If the node that was used for calculating the active edges can be found
*   on the left side of the pivot element.
*/
void refineBySet(TreeList& forest, vector<TreeNode*> treeNodes,
                 int timestamp, bool nodeIsLeft) {

    unordered_set<TreeNode*> maxSubtrees =
            getMaxContSubTrees(treeNodes, timestamp);

    unordered_set<TreeNode*> maxSubtreeParents;
    for (TreeNode* node : maxSubtrees) {
        if (node->parent != nullptr) {
            maxSubtreeParents.insert(node->parent);
        }
    }

    for (TreeNode* parentNode : maxSubtreeParents) {
        MD_Tree* correspondingTree = forest.getCorrespondingTree(parentNode);
        bool leftSplit = nodeIsLeft || correspondingTree->isLeftOfPivot;
        if (parentNode->label == PRIME) {
            refineByPrimeNode(parentNode, leftSplit, maxSubtrees);
        }
        else {
            refineByNonPrimeNode(forest, parentNode, maxSubtrees, leftSplit);
        }
    }
}

/**
* Executes the promotion algorithm for a specific tree, given by its root. For more details on promotion,
* see the algorithm description.
*
* @param root The root of the MD_Tree that should be promoted
* @return A list of all trees that are calculated in the promotion
*/
vector<MD_Tree> getPromotedTree(TreeNode* root) {
    vector<MD_Tree> forest;
    if (root->markedLeft) {
        TreeNode* markedChild = root->child;
        TreeNode* previous = root;

        while (markedChild != nullptr) {
            if (markedChild->markedLeft) {
                TreeNode* oldSibling = markedChild->sibling;
                if (previous == root) {
                    root->child = markedChild->sibling;
                }
                else {
                    previous->sibling = markedChild->sibling;
                }
                root->nChildNodes--;
                markedChild->parent = nullptr;
                markedChild->sibling = nullptr;
                vector<MD_Tree> left = getPromotedTree(markedChild);
                forest.insert(forest.end(), left.begin(), left.end());
                markedChild = oldSibling;
            } else {
                previous = markedChild;
                markedChild = markedChild->sibling;
            }
        }
    }
    forest.push_back(MD_Tree(root));
    if (root->markedRight) {
        TreeNode* markedChild = root->child;
        TreeNode* previous = root;

        while (markedChild != nullptr) {
            if (markedChild->markedRight) {
                TreeNode* oldSibling = markedChild->sibling;
                if (previous == root) {
                    root->child = markedChild->sibling;
                }
                else {
                    previous->sibling = markedChild->sibling;
                }
                root->nChildNodes--;
                markedChild->parent = nullptr;
                markedChild->sibling = nullptr;
                vector<MD_Tree> right = getPromotedTree(markedChild);
                forest.insert(forest.end(), right.begin(), right.end());
                markedChild = oldSibling;
            } else {
                previous = markedChild;
                markedChild = markedChild->sibling;
            }
        }
    }
    return forest;
}

/**
* Deletes all marks ("left" or "right"), starting at a given node.
*
* @param node The starting node for the mark-deletion process.
*/
void deleteMarks(TreeNode* node) {
    node->markedLeft = false;
    node->markedRight = false;

    if (node->sibling != nullptr) {
        deleteMarks(node->sibling);
    }
    if (node->child != nullptr) {
        deleteMarks(node->child);
    }
}

/**
* Cleans the forest after promotion is done. This includes:
*   - Every root with no child will be removed (Except the root is a leaf).
*   - Whenever a root as only one child, that child will take the place of the root.
*   - All markings will be delted.
*
* @param forest The forest to be cleaned up.
*/
void cleanUp(vector<MD_Tree>& forest) {
    vector<MD_Tree> newTreeList;
    for (MD_Tree& tree : forest) {
        if (tree.root->markedLeft || tree.root->markedRight) {
            if (tree.root->child != nullptr) {
                if (tree.root->child->sibling == nullptr) {
                    tree.root->child->parent = nullptr;
                    newTreeList.push_back(MD_Tree(tree.root->child));
                    delete tree.root;
                }
                else {
                    newTreeList.push_back(MD_Tree(tree.root));
                }
            }
            else if (tree.root->label == LEAF) {
                newTreeList.push_back(MD_Tree(tree.root));
            }
            else {
                delete tree.root;
            }
        }
        else {
            newTreeList.push_back(MD_Tree(tree.root));
        }
    }

    for (MD_Tree& tree : newTreeList) {
        deleteMarks(tree.root);
    }
    forest = newTreeList;
}

/**
* Updates the value of a given node based on a given list. Calls itself recursively.
*
* @param node The current node.
* @param updatedValues A vector containing the new values of the nodes.
*/
void updateTreeValues(TreeNode* node, const vector<int>& updatedValues) {
    if (node->label == LEAF && node->value < updatedValues.size()) {
        node->value = updatedValues[node->value];
    }
    if (node->sibling != nullptr) {
        updateTreeValues(node->sibling, updatedValues);
    }
    if (node->child != nullptr) {
        updateTreeValues(node->child, updatedValues);
    }
}

/**
* Updates the nodeValueMaping, based on a new indexMapping.
*
* @param newMapping This is used as the return of the function.
* @param oldMapping The current (now 'old') nodeValueMapping.
* @param indexMapping A new subgraphIndexMapping that is used to update the node values.
*/
void updateNodeValueMapping(vector<TreeNode*>& newMapping, vector<TreeNode*>& oldMapping,
                            const vector<int>& indexMapping) {

    for (int i = 0; i < oldMapping.size(); i++) {
        newMapping[indexMapping[i]] = oldMapping[i];
    }
}


/**
* Inserts the given index in the list of indices, if a given node has a value.
* Afterwards, this function recursively calls the nodes children and siblings.
* This is a helper function for getMaxModuleIndices.
*
* @param node The current node.
* @param indices The list to insert into.
* @param index The value to insert.
*/
void insertIndex(TreeNode* node, vector<int>& indices, int index) {
    if (node->label == LEAF) {
        indices[node->value] = index;
    }
    if (node->sibling != nullptr) {
        insertIndex(node->sibling, indices, index);
    }
    if (node->child != nullptr) {
        insertIndex(node->child, indices, index);
    }
}

/**
* Computes a list that provides information about every vertex in the graph:
* The index of the module it belongs to at the start of the assembly-process.
*
* @param forest The MD-forest
* @param graphSize The number of elements in the graph.
* @return A list that contains the searched information.
*/
vector<int> getMaxModuleIndices(const vector<MD_Tree>& forest, int graphSize) {
    vector<int> indices(graphSize, -1);
    for (int i = 0; i < forest.size(); i++) {
        insertIndex(forest[i].root, indices, i);
    }
    return indices;
}


/**
* Creates a large list of adajencies that contains all adjacencies of the elements in a given list.
*
* @param vertices The vertices to sum up the adjacencies of.
* @param graph The graph.
* @return The summed-up list of adjacencies.
*/
vector<int> getTotalAdjacencyList(const vector<int>& vertices, const Graph& graph) {
    vector<vector<int>> adjlist = graph.getAdjlist();
    vector<int> totalAdjacencies;
    for (int vertex : vertices) {
        totalAdjacencies.insert(totalAdjacencies.end(), adjlist[vertex].begin(), adjlist[vertex].end());
    }
    return totalAdjacencies;
}

/**
* Uses the connections given in the graph to insert left- and right pointers for every element in
* the forest. The left pointer of an element X is on the lowest index, so that all elements with lower
* indices are connected to X. The right pointer of this element is on the highest index, so that all
* elements with higher indices are disconnected to X.
*
* @param forest Contains the list of MD-trees.
* @param graph The graph.
*/
void insertLeftRightPointers(vector<MD_Tree>& forest, const Graph& graph) {
    int n = static_cast<int>(forest.size());
    vector<int> maxModuleIndices = getMaxModuleIndices(forest, static_cast<int>(graph.getAdjlist().size()));
    vector<bool> connections(n, false);

    for (int i = 0; i < n; i++) {
        vector<int> vertices = getPreOrderLeafs(forest[i].root);
        vector<int> totalAdjacencies = getTotalAdjacencyList(vertices, graph);
        int maxConnectionIndex = numeric_limits<int>::min();
        vector<int> connectionIndices;

        for (int j = 0; j < totalAdjacencies.size(); j++) {
            int connectedModule = maxModuleIndices[totalAdjacencies[j]];
            if (connectedModule != -1) {
                connections[connectedModule] = true;
                connectionIndices.push_back(connectedModule);
                if (connectedModule > maxConnectionIndex) {
                    maxConnectionIndex = connectedModule;
                }
            }
        }
        int leftPointer = 0;
        while (leftPointer < i && connections[leftPointer]) {
            leftPointer++;
        }

        // Left and Right pointers are always to the left of the element with the given index
        forest[i].leftIndex = leftPointer;
        forest[i].rightIndex = maxConnectionIndex + 1;

        // Reset the connections vector
        for (int connection : connectionIndices) {
            connections[connection] = false;
        }
    }
}

/**
* Checks, if the subtree induced by a given node is connected to the pivot element.
*
* @param node The node to check.
* @param pivotAdj The adjacency list of the pivot element.
* @return If there is a connection to the pivot.
*/
bool isConnectedToPivot(TreeNode* node, const vector<int>& pivotAdj) {
    if (node->label == LEAF) {
        return find(pivotAdj.begin(), pivotAdj.end(), node->value) != pivotAdj.end();
    }
    if (node->child != nullptr) {
        return isConnectedToPivot(node->child, pivotAdj);
    }
    return isConnectedToPivot(node->sibling, pivotAdj);
}

/**
* Checks, if the next module in a list of trees only extends to the right, meaning that it is parallel.
*
* @param forest Contains the list of trees.
* @param currentLeft The left index of the current module.
* @param currentRight The right index of the current module.
* @return if the next module is parallel.
*/
bool checkForParallel(vector<MD_Tree>& forest, int currentLeft, int currentRight) {
    if (currentRight >= forest.size()) {
        return false;
    }
    int i = currentRight;
    currentRight++;
    while (i < currentRight) {
        int leftPointer = forest[i].leftIndex;
        int rightPointer = forest[i].rightIndex;

        if (leftPointer < currentLeft) {
            return false;
        }
        if (rightPointer > currentRight) {
            currentRight = rightPointer;
        }
        i++;
    }
    return true;
}

/**
* A helper method for the construction of the modular decompositon
* tree in the assembly step of the algorithm.
*
* @param currentNode The node to be added to the MD tree.
* @param currentModuleType The type of the current module
* @param lastNode The last node that was added to the tree.
* @param moduleChildCounter A counter to keep track
*   of the amount of a modules children.
*/
void addToMDTree(TreeNode*& currentNode, Label& currentModuleType,
                 TreeNode*& lastNode, int& moduleChildCounter) {

    if (currentNode->label == currentModuleType && currentNode->child != nullptr) {
        setSibling(lastNode, currentNode->child);
        lastNode = currentNode->child;
        moduleChildCounter++;
        while (lastNode->sibling != nullptr) {
            setSibling(lastNode, lastNode->sibling);
            lastNode = lastNode->sibling;
            moduleChildCounter++;
        }
        delete currentNode;
    }
    else {
        setSibling(lastNode, currentNode);
        lastNode = currentNode;
        moduleChildCounter++;
    }
}

/**
* Performs the first step of the algorithm, the recursion. In this step, the vertices of
* the given graph are sorted by their distance to a given pivot element. The MD_Trees of each
* set are computed recursively and stored in a MD_Forest, that contains a list of theses MD_Trees.
* Furthermore, the active edges are computed as well.
*
* @param graph The graph on which the modular decomposition should be executed
* @param pivot The arbitrarly choosen vertex that works as pivot-element for the algorithm
* @param activeEdges This vector will contain the activeEdges of the given graph with the given pivot,
*   after the method has terminated
* @param leftNodes This vector will contain information about which nodes are in a tree to the
    left of the pivot, after the algorithm has terminated.
* @param nodeValueMapping This mapping stores the corresponding tree node
*   for every value of a graph's vertex.
* @return The resulting MD_Forest (list of MD_Trees)
*/
TreeList recursion(const Graph& graph, int pivot, vector<vector<int>>& activeEdges,
                   vector<bool>& leftNodes, vector<TreeNode*>& nodeValueMapping) {

    vector<vector<int>> adjlist = graph.getAdjlist();
    int adjlistSize = static_cast<int>(adjlist.size());
    vector<int> distancesToPivot(adjlistSize, -1);
    activeEdges.assign(adjlistSize, vector<int>());
    leftNodes.assign(adjlistSize, false);

    queue<int> bfsQueue;
    bfsQueue.push(pivot);
    distancesToPivot[pivot] = 0;
    int maxDistance = 0;

    while (!bfsQueue.empty()) {
        int current = bfsQueue.front();
        bfsQueue.pop();
        for (int neighbor : adjlist[current]) {
            if (distancesToPivot[neighbor] == -1) {
                distancesToPivot[neighbor] = distancesToPivot[current] + 1;
                if (distancesToPivot[current] + 1 > maxDistance) {
                    maxDistance = distancesToPivot[current] + 1;
                }
                bfsQueue.push(neighbor);
            }
            if (distancesToPivot[neighbor] != distancesToPivot[current]) {
                activeEdges[current].push_back(neighbor);
            }
            if (current == pivot) {
                leftNodes[neighbor] = true;
            }
        }
    }

    TreeList* output = new TreeList();
    vector<vector<int>> subgraphIndexMappings;
    vector<Graph> subgraphs = graph.getSubGraphs(distancesToPivot, subgraphIndexMappings, maxDistance);

    for (int i = 0; i < subgraphs.size(); i++) {
        Graph currentSubgraph = subgraphs[i];
        vector<int> currentIndexMapping = subgraphIndexMappings[i];
        vector<TreeNode*> nodeValueMap(currentSubgraph.getAdjlist().size());
        MD_Tree* tree = new MD_Tree(getModularDecomposition(currentSubgraph, nodeValueMap));
        if (i == 0) {
            tree->isLeftOfPivot = true;
        }
        updateNodeValueMapping(nodeValueMapping, nodeValueMap, currentIndexMapping);
        updateTreeValues(tree->root, currentIndexMapping);
        output->insert(tree);
    }

    return *output;
}


/**
* Performs the second step of the algorithm, the refinement. As step 1 (recursion) only
* used one node (the pivot) to * separate all vertices into modules, refinement takes care
* of using all other nodes to separate the vertices even further.
* For more details on the refinement process, see the algorithm description.
*
* @param graph The graph that should be modular decomposed.
* @param nodeValueMapping This mapping stores the corresponding tree node
*   for every value of a graph's vertex.
* @param previousPivot The pivot element that was selected in step 1 (recursion).
* @param forest The MD_Forest that was calculated in step 1 (recursion).
* @param activeEdges A set of all activeEdges.
* @param leftNode Contains information about which nodes are to the left of the pivot.
*/
void refinement(const Graph& graph, vector<TreeNode*>& nodeValueMapping, int previousPivot,
                TreeList& forest, const vector<vector<int>>& activeEdges, const vector<bool>& leftNodes) {

    for (int i = 0; i < graph.getAdjlist().size(); i++) {
        if (i != previousPivot) {
            vector<TreeNode*> treeNodes;
            for (int nodeValue : activeEdges[i]) {
                TreeNode* treeNode = nodeValueMapping[nodeValue];
                if (treeNode != nullptr) {
                    treeNodes.push_back(treeNode);
                }
            }
            refineBySet(forest, treeNodes, i, leftNodes[i]);
        }
    }
}

/**
* Executes the third step of the algorithm, the promotion. For more details on the promotion - process,
* see the algorithm description.
*
* @param forest The MD_Forest to execute promotion on.
* @return The promoted forest as a vector of MD_Trees
*/
vector<MD_Tree> promotion(TreeList& forest) {
    vector<MD_Tree> newTreeList;
    MD_Tree* current = forest.getStart();
    while (current != nullptr) {
        vector<MD_Tree> promotedList = getPromotedTree(current->root);
        newTreeList.insert(newTreeList.end(), promotedList.begin(), promotedList.end());
        current = current->right;
    }
    cleanUp(newTreeList);
    return newTreeList;
}

/**
* Executes the fourth (and final) step of the algorithm, the assembly. For more details on the
* assembly - process, see the algorithm description.
*
* @param forest The MD_Forest resulting from step 3: promotion.
* @param graph The graph.
* @param nodeValueMapping This mapping stores the corresponding tree node
*   for every value of a graph's vertex.
* @param pivot The previously chosen pivot element.
* @return The final assembled MD-tree.
*/
MD_Tree assembly(vector<MD_Tree>& forest, const Graph& graph, vector<TreeNode*>& nodeValueMapping, int pivot) {

    TreeNode* pivotNode = new TreeNode(pivot);
    MD_Tree pivotTree = MD_Tree(pivotNode);
    nodeValueMapping[pivot] = pivotNode;
    int pivotIndex = 1;
    while (pivotIndex < forest.size() && isConnectedToPivot(forest[pivotIndex].root,
                                                            graph.getAdjlist()[pivot])) {
        pivotIndex++;
    }
    forest.insert(forest.begin() + pivotIndex, pivotTree);
    insertLeftRightPointers(forest, graph);

    TreeNode* lastModule = pivotNode;

    int currentLeft = pivotIndex;
    int currentRight = pivotIndex + 1;
    int includedLeft = pivotIndex;
    int includedRight = pivotIndex + 1;

    do {
        bool addedRight = false;
        bool addedLeft = false;
        queue<int> maxModuleIndices;

        if (checkForParallel(forest, currentLeft, currentRight)) {
            maxModuleIndices.push(currentRight);
            addedRight = true;
            currentRight++;
        }
        else {
            maxModuleIndices.push(currentLeft - 1);
            addedLeft = true;
            currentLeft--;
        }

        do {
            int currentMaxModule = maxModuleIndices.front();
            maxModuleIndices.pop();

            int leftPointer = forest[currentMaxModule].leftIndex;
            int rightPointer = forest[currentMaxModule].rightIndex;

            if (leftPointer < currentLeft) {
                for (int i = currentLeft - 1; i >= leftPointer; i--) {
                    maxModuleIndices.push(i);
                }
                currentLeft = leftPointer;
                addedLeft = true;
            }
            if (rightPointer > currentRight) {
                for (int i = currentRight; i < rightPointer; i++) {
                    maxModuleIndices.push(i);
                }
                currentRight = rightPointer;
                addedRight = true;
            }
        } while (!maxModuleIndices.empty());

        Label moduleType = PARALLEL;
        if (addedLeft && addedRight) {
            moduleType = PRIME;
        }
        else if (addedLeft) {
            moduleType = SERIES;
        }
        TreeNode* moduleNode = new TreeNode(moduleType);
        setChild(moduleNode, lastModule);
        TreeNode* lastNode = lastModule;
        int moduleChildCounter = 1;

        for (int i = currentLeft; i < includedLeft; i++) {
            addToMDTree(forest[i].root, moduleType, lastNode, moduleChildCounter);
        }
        for (int i = includedRight; i < currentRight; i++) {
            addToMDTree(forest[i].root, moduleType, lastNode, moduleChildCounter);
        }

        includedLeft = currentLeft;
        includedRight = currentRight;
        lastModule = moduleNode;
        moduleNode->nChildNodes = moduleChildCounter;

    } while (currentLeft > 0 || currentRight < forest.size());

    return MD_Tree(lastModule);
}

/**
* Returns the Modular Decomposition for a disconnected Graph, by creating one PARALLEL node
* and setting the recursively computed MD_Trees of the graphs components as its children.
*
* @param graph The graph
* @param nodeValueMapping This mapping stores the corresponding tree node
*   for every value of a graph's vertex.
* @return The recursively computed MD_Tree
*/
MD_Tree getModularDecompositionDisconnectedGraph(const Graph& graph, vector<TreeNode*>& nodeValueMapping) {
    vector<vector<int>> components = graph.getConnectedComponents();
    TreeNode* rootNode = new TreeNode(PARALLEL);
    TreeNode* lastChild = nullptr;

    for (int i = 0; i < components.size(); i++) {
        vector<int> subgraphIndexMapping;
        Graph subgraph = graph.getSubGraph(components[i], subgraphIndexMapping);
        vector<TreeNode*> nodeValueMap(subgraph.getAdjlist().size());
        MD_Tree tree = getModularDecomposition(subgraph, nodeValueMap);
        updateNodeValueMapping(nodeValueMapping, nodeValueMap, subgraphIndexMapping);
        updateTreeValues(tree.root, subgraphIndexMapping);
        if (i == 0) {
            setChild(rootNode, tree.root);
            lastChild = tree.root;
        }
        else {
            setSibling(lastChild, tree.root);
            lastChild = tree.root;
        }
    }
    rootNode->nChildNodes = static_cast<int>(components.size());
    return MD_Tree(rootNode);
}

/**
* Returns the modular decomposition tree for a given graph. This is done by using the four steps:
*   - Recursion
*   - Refinement
*   - Promotion
*   - Assembly
* For more details on all of these steps and the algorithm itself, see the algorithm description.
*
* @param graph The graph to be modular decomposed.
* @param nodeValueMapping This mapping stores the corresponding tree node
*   for every value of a graph's vertex.
*/
MD_Tree getModularDecomposition(const Graph& graph, vector<TreeNode*>& nodeValueMapping) {

    if (graph.getAdjlist().size() == 0) {
        cerr << "You cannot get a modular decomposition of a graph with no vertices!" << endl;
    }
    else if (graph.getAdjlist().size() == 1) {
        TreeNode* root = new TreeNode(0);
        nodeValueMapping[0] = root;
        return MD_Tree(root);
    }

    if (graph.isConnected()) {
        int pivot = 0;
        vector<vector<int>> activeEdges;
        vector<bool> leftNodes;

        TreeList forest = recursion(graph, pivot, activeEdges, leftNodes, nodeValueMapping);

        if (false) {
            cout << "After Recursion: " << endl;
            forest.print();
            cout << endl;
        }

        refinement(graph, nodeValueMapping, pivot, forest, activeEdges, leftNodes);

        if (false) {
            cout << "After Refinement: " << endl;
            forest.print();
            cout << endl;
        }

        vector<MD_Tree> forestVec = promotion(forest);

        if (false) {
            cout << "After Promotion: " << endl;
            printForest(forestVec);
            cout << endl;
        }

        MD_Tree finalResult = assembly(forestVec, graph, nodeValueMapping, pivot);
        resetTimestamps(finalResult.root);
        //delete &forest;

        return finalResult;
    }
    else {
        MD_Tree finalTree = getModularDecompositionDisconnectedGraph(graph, nodeValueMapping);
        return finalTree;
    }
}

/**
* Returns the modular decomposition tree for a given graph. This is done by using the four steps:
*   - Recursion
*   - Refinement
*   - Promotion
*   - Assembly
* For more details on all of these steps and the algorithm itself, see the algorithm description.
*
* @param graph The graph to be modular decomposed.
*/
MD_Tree getModularDecomposition(const Graph& graph) {
    vector<TreeNode*> nodeValueMapping(graph.getAdjlist().size());
    return getModularDecomposition(graph, nodeValueMapping);
}


/**
* The main method.
*/
int main(int argc, char* argv[]) {
    string adjList = "";
    if (argc >= 2) {
        string filePath = argv[1];
        adjList = readFile(filePath);
    } else {
        cout << "Please provide valid arguments!" << endl;
        return 1;
    }
    vector<int> indexMapping;

    if (!Util::isConsecutivelyOrdered(adjList)) {
        vector<int> indexMapping = Util::rewriteAdjacencyList(adjList);
    }

    Graph graph = Graph(adjList);
    int nVertices = Util::getNumberVertices(graph);
    int nEdges = Util::getNumberEdges(graph);

    auto start = chrono::high_resolution_clock::now();
    MD_Tree mdTree = getModularDecomposition(graph);
    auto end = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start);

    /*
    Util::sortTree(mdTree);
    cout << "The final MD-tree: " << endl << endl;
    printTree(mdTree.root);

    cout << endl << "Time needed: " << duration.count() << " nanoseconds for a graph with " << Util::getNumberVertices(graph)
    << " vertices and " << Util::getNumberEdges(graph) << " edges" << endl << endl;
     */

    string output = "ModularDecomposition,";
    if (Util::testModularDecompositionTree(graph, mdTree)) {
        output += "0,";
    }
    else {
        output += "1,";
    }
    output += to_string(nVertices) + "," + to_string(nEdges) + "," + to_string(nVertices + nEdges) + ",";
    output += to_string(duration.count());
    cout << output << endl;
}

/**
 * Tests the correctness of the modular decomposition algorithm by randomly generating modular decomposition trees
 * of a specific size. These random trees get converted into adjacency-lists, which will be used by the modular decomposition
 * algorithm to compute a modular decomposition tree. This tree then gets compared to the initial MD-Tree. Multiple
 * repetitions of this process ensure the algorithms correctness.
 */
void testModularDecomposition () {
    int nRepetitions = 100;
    int nVertices = 25;
    bool useCoGraphs = false;

    int equalCounter = 0;
    for (int i = 0; i < nRepetitions; i++) {

        MD_Tree tree = Util::createRandomModularDecompositionTree(nVertices, useCoGraphs);
        Util::sortTree(tree);
        Graph graph = Util::createGraphFromTree(tree);
        MD_Tree mdTree = getModularDecomposition(graph);
        Util::sortTree(mdTree);

        string tree1Representation = generateTreeString(tree.root);
        string tree2Representation = generateTreeString(mdTree.root);
        if (tree1Representation == tree2Representation) {
            equalCounter++;
            if (equalCounter % 10 == 0) {
                cout << endl << "Equal: " << equalCounter << " of " << (i + 1) << endl;
            }
        } else {
            cout << "Not equal:" << endl;
            Util::testModularDecompositionTree(graph, mdTree);

            graph.print();
            cout << "Initial MD_Tree:" << endl;
            printTree(tree.root);
            cout << endl << endl;
            cout << "Calculated MD_Tree:" << endl;
            printTree(mdTree.root);
            cout << endl;
        }
    }
    cout << endl << "Equal: " << equalCounter << endl;
    cout << "Not Equal: " << (nRepetitions - equalCounter) << endl;
}
