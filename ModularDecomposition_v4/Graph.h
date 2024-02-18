#pragma once
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

using namespace std;

class Graph
{
private:
    vector<vector<int>> adjlist;

public:
    Graph(vector<vector<int>>& adjlist);
    Graph(const string& graphString);
    const vector<vector<int>>& getAdjlist() const;
    Graph getSubGraph(const vector<int>& X, vector<int>& nodeMapping) const;
    vector<Graph> getSubGraphs(const vector<int>& distances,
                               vector<vector<int>>& nodeMapping, int nSubgraphs) const;
    bool isConnected() const;
    vector<vector<int>> getConnectedComponents() const;
    void print() const;

private:
    void dfs(int node, vector<bool>& visited, vector<int>& componentNodes) const;
    vector<string> splitString(const string& input) const;
};
