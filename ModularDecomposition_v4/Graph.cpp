#include "Graph.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <queue>

/**
* Constructor implementations
*/
Graph::Graph(vector<vector<int>>& adjlist)
{
    this->adjlist = adjlist;
}

Graph::Graph(const string& graphString)
{
    string graphStringCopy = graphString;
    istringstream iss(graphStringCopy);
    vector<string> lines;
    string line;

    while (getline(iss, line)) {
        lines.push_back(line);
    }

    for (int i = 0; i < lines.size(); i++) {
        vector<int> currentAdjacencies;
        vector<string> currentAdjacenciesString = splitString(lines[i]);
        for (int j = 1; j < currentAdjacenciesString.size(); j++) {
            currentAdjacencies.push_back(stoi(currentAdjacenciesString[j]));
        }
        adjlist.push_back(currentAdjacencies);
    }
}

/**
* Getter
*/
const vector<vector<int>>& Graph::getAdjlist() const
{
    return adjlist;
}

/**
* Returns a graph, that only contains the vertices in a given set X. All edges
* are calculated accordingly.
*
* @param X The set
* @param nodeMapping Mapping from old indices to new indices (updated in the method)
* @return A subgraph that only contains the vertices given in X.
*/
Graph Graph::getSubGraph(const vector<int>& X, vector<int>& nodeMapping) const
{
    vector<vector<int>> newAdjList;
    for (int i : X) {
        nodeMapping.push_back(i);
    }

    for (int i : X) {
        vector<int> neighbors;
        for (int neighbor : adjlist[i]) {
            auto it = find(X.begin(), X.end(), neighbor);
            if (it != X.end()) {
                neighbors.push_back(distance(X.begin(), it));
            }
        }
        newAdjList.push_back(neighbors);
    }

    return Graph(newAdjList);
}

/**
* Returns a list of subgraphs, where the vertices of G are distributed based on a vector with their
* distances to a pivot element.
*
* @param distances The elements of the graph will be distributed according to their values in this vector.
* @param nodeMapping Mapping from old indices to new indices (updated in the method).
* @param nSubgraphs The amount of subgraphs that have to be computed.
* @return A list of subgraphs based on the distances vector.
*/
vector<Graph> Graph::getSubGraphs(const vector<int>& distances, vector<vector<int>>& nodeMapping, int nSubgraphs) const
{
    vector<vector<vector<int>>> adjacencyLists(nSubgraphs);
    nodeMapping.assign(nSubgraphs, vector<int>());
    vector<int> inverseNodeMapping(distances.size());

    for (int i = 0; i < distances.size(); i++) {
        if (distances[i] > 0) {
            nodeMapping[distances[i] - 1].push_back(i);
            inverseNodeMapping[i] = nodeMapping[distances[i] - 1].size() - 1;
        }
    }

    for (int i = 0; i < distances.size(); i++) {
        if (distances[i] > 0) {
            vector<int> neighbors;
            for (int neighbor : adjlist[i]) {
                if (distances[i] == distances[neighbor]) {
                    neighbors.push_back(inverseNodeMapping[neighbor]);
                }
            }
            adjacencyLists[distances[i] - 1].push_back(neighbors);
        }
    }

    vector<Graph> output;
    for (vector<vector<int>> adjList : adjacencyLists) {
        output.push_back(Graph(adjList));
    }
    return output;
}


/**
* Uses BFS to check if the graph is connected
*/
bool Graph::isConnected() const {
    vector<bool> visited(adjlist.size(), false);
    queue<int> q;

    q.push(0);
    visited[0] = true;

    while (!q.empty()) {
        int current = q.front();
        q.pop();

        for (int neighbor : adjlist[current]) {
            if (!visited[neighbor]) {
                q.push(neighbor);
                visited[neighbor] = true;
            }
        }
    }

    return all_of(visited.begin(), visited.end(), [](bool v) { return v; });
}

/**
* Uses DFS to retrieve all components of the graph.
*/
vector<vector<int>> Graph::getConnectedComponents() const {
    vector<vector<int>> components;
    vector<bool> visited(adjlist.size(), false);

    for (int i = 0; i < adjlist.size(); ++i) {
        if (!visited[i]) {
            vector<int> componentNodes;
            dfs(i, visited, componentNodes);
            components.push_back(componentNodes);
        }
    }

    return components;
}

/**
* A simple recursive DFS algorithm.
*
* @param node The current node.
* @param visited Contains information about which nodes have already been visited.
* @param componentNodes This set is used to store the results.
*/
void Graph::dfs(int node, vector<bool>& visited, vector<int>& componentNodes) const {
    visited[node] = true;
    componentNodes.push_back(node);

    for (int neighbor : adjlist[node]) {
        if (!visited[neighbor]) {
            dfs(neighbor, visited, componentNodes);
        }
    }
}

vector<string> Graph::splitString(const string& input) const {
    vector<string> tokens;
    istringstream iss(input);

    string token;
    while (getline(iss, token, ' ')) {
        tokens.push_back(token);
    }
    return tokens;
}

void Graph::print() const {
    cout << "Graph:" << endl;
    for (int i = 0; i < adjlist.size(); i++) {
        cout << i << ": ";
        for (int j = 0; j < adjlist[i].size(); j++) {
            cout << adjlist[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

