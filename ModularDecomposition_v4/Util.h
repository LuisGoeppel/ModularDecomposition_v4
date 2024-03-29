#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <vector>
#include <random>

#include "MD_Tree.h"
#include "Graph.h"

using namespace std;

struct Ancestors{
    TreeNode* commonAncestor;
    TreeNode* ancestorChildLhs;
    TreeNode* ancestorChildRhs;
    int ancestorIndexLhs;
    int ancestorIndexRhs;
};

class Util {
public:
    static void sortTree(MD_Tree& tree);
    static vector<int> rewriteAdjacencyList(string& adjList);
    static bool isConsecutivelyOrdered(const string& input);
    static bool testModularDecompositionTree(const Graph& graph, const MD_Tree& tree, const vector<TreeNode*>& nodeValueMapping);
    static void checkChildNodeValues(MD_Tree& tree);
    static int getNumberVertices(const Graph& graph);
    static int getNumberEdges(const Graph& graph);
    static MD_Tree createRandomModularDecompositionTree(int nVertices, bool isCoGraph, vector<TreeNode*>& nodeValueMapping);
    static Graph createGraphFromTree(const MD_Tree& tree, const vector<TreeNode*>& nodeValueMapping);
    static bool isValidModularDecompositionTree(const MD_Tree& tree);
    static Ancestors getCommonAncestor(const vector<TreeNode*>& nodeValueMapping, int lhs, int rhs);

private:
    static void sortTreeHelper(TreeNode* node);
    static vector<int> getInverse(const vector<int>& vec);
    static TreeNode* getCommonAncestorChildLhs(const MD_Tree& tree, int lhs, int rhs);
    static TreeNode* getCorrespondingTreeNode(TreeNode* currentNode, int node);
    static void checkChildNodeValuesHelper(TreeNode* node);
    static int removeFalseInnerNodes(TreeNode* currentNode, TreeNode* callingNode, bool isParent);
    static void getHighestVertexValue(TreeNode* node, int& highestValue);
    static bool isValidTreeHelper(TreeNode* node, unordered_set<int>& leaves);
    static Label generateDifferentLabel(const Label& label, bool  isCoGraph);
    static void updateParentPointersAndModuleTypes(TreeNode* currentNode, TreeNode* parentNode, bool primeModule);
    static Label getLabel(int n);
    static void resetTimestamps(TreeNode* node);
    static void resetParentTimestamps(TreeNode* node);
};

