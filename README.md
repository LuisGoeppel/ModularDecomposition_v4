## Simple, Linear-time Modular Decomposition Algorithm

### Introduction
This C++ project implements the *Simple, Linear-time Modular Decomposition* algorithm, as developed by *Marc Tedder*, *Derek Corneil*, *Michel Habib*, and *Christophe Paul*. The details of the algorithm are described in their paper titled *"Simple, Linear-time Modular Decomposition."* The authors are confident that they have achieved the first simple, linear-time solution for this problem.

### Modular Decomposition
*Modular Decomposition* is a fundamental concept in algorithmic graph theory. It aims to represent the modules of a graph through inner nodes in a Modular Decomposition (MD) tree. A (strong) module is characterized by the fact that every vertex outside of the module is either connected to all the vertices in the module or to none of them. This algorithm provides a method to efficiently compute the MD tree, which helps in understanding the modular structure of a graph.

### Algorithm Steps

#### Step 1: Recursion
The first step of the algorithm is the *Recursion* phase, which plays a crucial role in simplifying the complexity of modular decomposition. At the outset, a pivot element *x* is arbitrarily chosen. The remaining vertices are categorized into sets *N<sub>i</sub>* based on their distance from *x*. These sets include:

- *N<sub>0</sub>* = { *n* | starting at *x*, vertex *n* can be reached using only one edge},
- *N<sub>1</sub>* = { *n* | starting at *x*, the shortest path to *n* has exactly two edges},
- *N<sub>2</sub>* = { *n* | starting at *x*, the shortest path to *n* has exactly three edges},
- ...

The algorithm then recursively computes MD-trees *T(N<sub>i</sub>)* for these *N<sub>i</sub>*. Importantly, this separation of graph vertices aids in representing the modules of the graph through inner nodes in one of the recursively computed MD-trees.

*T(N<sub>0</sub> ), x, T(N<sub>1</sub> ), ... , T(N<sub>k</sub> )*

It is noteworthy that the base case for this recursion is trivial: if the input consists of only one vertex, a modular decomposition tree with that vertex as the only node is returned. Despite the assumption of a connected input graph, the algorithm gracefully handles disconnected graphs. In such cases, it determines all connected components separately and calculates their modular decomposition trees. These trees are then united under a common node labeled "PARALLEL," representing the root of the resulting modular decomposition tree.

#### Step 2: Refinement
Following the *Recursion* phase, the *Refinement* step refines the obtained modular decomposition trees, turning them into a forest of trees. Each module that does not contain the pivot element *x* is now represented by one node in this forest.

In essence, *Refinement* can be seen as the process of partitioning already partitioned sets further, aiming to include new information. The algorithm achieves this by refining the forest based on a set of active edges. An edge is considered active if and only if it is adjacent to *x* or connects two vertices from different *N<sub>i</sub>*.

The algorithm reexamines each vertex, except for *x*, and checks connections to nodes in *T(N<sub>k-1</sub>)* and *T(N<sub>k+1</sub>)* if they exist. This process ensures that all modules not containing *x* are correctly represented in the forest, as connections to all vertices (initially *x* in Step 1 and every other vertex in Step 2) have been thoroughly examined.

After refinement, strong modules not containing *x* appear consecutively in the forest. This ordered structure allows for a pre-order traversal of all trees in the forest, effectively listing all elements of modules not containing *x* one after another. Nodes in the forest without marked children correspond exactly to strong modules not containing *x*.

#### Step 3: Promotion
The *Promotion* step builds upon the work done in *Refinement*, aiming to split the trees of the forest so that each remaining tree corresponds to a component of the final MD tree. While much of the work has been completed in the previous step, some additional refinement is necessary.

*Refinement* introduced new nodes to separate siblings that do not belong to the same strong module and left markings on those nodes. *Promotion* utilizes these markings to determine where the different strong modules are placed in the forest and where the trees should be separated from each other.

The algorithm iterates through pairs of child- and parent nodes with the same label. As long as such pairs exist, the connection between these nodes is removed, resulting in two new trees. These trees replace the old one in the forest. The markings provide the information needed to determine the placement of the new trees – if the child- and parent nodes are marked with "left," the tree with the former child node as the root is placed on the left of the tree with the parent node. If marked with "right," it is placed on the right.

Upon completion of the promotion algorithm, the result of traversing all leaves of the forest in pre-order is a factorizing permutation.

#### Step 4: Assembly
In the final step of the algorithm, *Assembly*, the individual trees of the forest resulting from Step 3 are combined into a final modular decomposition tree. The process involves constructing the spine of the tree using the factorizing permutation and considering left- and right-pointers.

To construct the final tree, the algorithm starts at the position of the pivot element *x* and forms the MD-tree in multiple steps. In each step, a new inner node of the tree is created, having the roots of some trees from the forest as children. The resulting tree is built in such a way that *PARALLEL* nodes only add trees to the right of *x*, *SERIES* nodes only add trees to *x*'s left, and *PRIME* nodes have children on both sides.

The algorithm prioritizes the inclusion of parallel modules, attempting to include the element to the right of either *x* or the last-formed module. If a new tree is included in the current module, every tree up to its left- and up to its right-pointer is included as well. After no more trees can be added, the module is formed. If this process only added trees to the right of *x* (or to the right of the last-formed module), the module is successfully created. Otherwise, the algorithm starts again, this time by adding the tree to the left.

This process is repeated until every tree has been included. Before the finalization, small adjustments are made to unify nodes with the same label in the included trees. Nodes with the same label are replaced with a new node, preserving the label and incorporating both old nodes' children.

The result is the completed modular decomposition tree, ready to provide insights into the graph's modular structure.
### Usage
To use this C++ implementation of the *Simple, Linear-time Modular Decomposition* algorithm, follow these steps:

1. Clone the repository.
2. Provide the path to the Graph.txt-file as an argument in your configuration.
3. Insert your graph in the form of an adjacency list in the Graph.txt file.
4. Compile the source code using your preferred C++ compiler.
5. Execute the compiled program.

### Note
Ensure that the input graph is connected, simple and undirected. Otherwise, the algorithm might not produce any or faulty results.

### Correctness and Time-Complexity
This algorithm has been tested on more than 10 000 input graphs with sizes up to 10 000 vertices and edges. None of these tests showed signs of any wrong output. The method *testModularDecomposition*
can be used to test the algorithm manually.
Considering the time-complexity, the algorithm has shown to run within a quadratic time-bound. Unfortunately, this cannot keep up with the linear time-bounds promised by the authors of the algorithm.
The slow-down is due to a less complicated, but also less effective implementation of the recursion-step. Instead of performing a breath-first search over all recursive substeps, the implementation
performs one breadth-first search at every recursive substep. Hence, there is still room for improvement.

### Acknowledgments
Special thanks to *Marc Tedder*, *Derek Corneil*, *Michel Habib*, and *Christophe Paul* for their groundbreaking work on the *Simple, Linear-time Modular Decomposition* algorithm.

Feel free to explore and modify the source code to suit your specific needs. If you encounter any issues or have suggestions for improvements, please don't hesitate to open an issue or submit a pull request.

Happy graph analysis!
