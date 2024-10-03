
#ifndef Graph_h
#define Graph_h
/// @file Graph.h @brief Implements generic graph structure for fast neigbor search in arbitrary dimentions
/// @ingroup Topology


//  algorithm copied from
//    https://www.geeksforgeeks.org/bridge-in-a-graph/
// see also:
//      * https://cp-algorithms.com/graph/bridge-searching.html
//      * https://en.wikipedia.org/wiki/Bridge_(graph_theory)
//      * https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
//      * Find Bridges in a graph using Tarjans Algorithm https://www.youtube.com/watch?v=Rhxs4k6DyMM


// A C++ program to find bridges in a given undirected graph
//#include<iostream>
#include <list>
#define NIL -1
using namespace std;

// A class that represents an undirected graph
class Graph{ public:
	int V;                 // No. of vertices
	list<int> *adj;        // A dynamic array of adjacency lists
    vector<Vec2i> found;   // Array of lists    
    // ------ Functions
    void bridgeUtil(int v, bool visited[], int disc[], int low[], int parent[]);
	void addEdge(int v, int w); // to add an edge to graph
	void bridge(); // prints all bridges
	void fromBonds( int nb, const Vec2i* bonds );
    Graph(int V); // Constructor
};

Graph::Graph(int V){
	this->V = V;
	adj = new list<int>[V];
}

void Graph::addEdge(int v, int w){
	adj[v].push_back(w);
	adj[w].push_back(v); // Note: the graph is undirected
}

void Graph::fromBonds( int nb, const Vec2i* bonds ){
    for(int i=0; i<nb; i++){
        //if( (atypes[bonds[i].a]==ignore) || (atypes[bonds[i].b]==ignore)   ) continue;
        if( (bonds[i].a>=V) || (bonds[i].b>=V)  ) continue;
        addEdge(bonds[i].a, bonds[i].b);
    }
}

// A recursive function that finds and prints bridges using // DFS traversal
// u         --> The vertex to be visited next
// visited[] --> keeps track of visited vertices
// disc[]    --> Stores discovery times of visited vertices
// parent[]  --> Stores parent vertices in DFS tree
void Graph::bridgeUtil(int u, bool visited[], int disc[],  int low[], int parent[] ) {
	static int time = 0;          // A static variable is used for simplicity, we can avoid use of static variable by passing a pointer.
	visited[u]       = true;     // Mark the current node as visited
	disc[u] = low[u] = ++time;   // Initialize discovery time and low value
	list<int>::iterator i;       // Go through all vertices adjacent to this
	for (i = adj[u].begin(); i != adj[u].end(); ++i){
		int v = *i; // v is current adjacent of u
		if (!visited[v]){   // If v is not visited yet, then recur for it
			parent[v] = u;
			bridgeUtil( v, visited, disc, low, parent );
			low[u] = min(low[u], low[v]);                        // Check if the subtree rooted with v has a connection to one of the ancestors of u
			if (low[v] > disc[u]){   // If the lowest vertex reachable from subtree under v is below u in DFS tree, then u-v is a bridge
                printf( "%i %i \n", u, v );
                found.push_back( (Vec2i){u,v} );
            }
		} else if (v != parent[u])  low[u] = min(low[u], disc[v]);  // Update low value of u for parent function calls.
	}
}

// DFS based function to find all bridges. It uses recursive function bridgeUtil()
void Graph::bridge(){
    bool *visited = new bool[V];
    int *disc     = new int[V];
    int *low      = new int[V];
    int *parent   = new int[V];
    for (int i = 0; i < V; i++)	{   // Initialize parent and visited arrays
        parent[i] = NIL;
        visited[i] = false;          // Mark all the vertices as not visited
    }
    for (int i = 0; i < V; i++) if (visited[i] == false) bridgeUtil(i, visited, disc, low, parent);  // Call the recursive helper function to find Bridges in DFS tree rooted with vertex 'i'
}

int test_Graph_findBridge(){   // Driver program to test above function
	// Create graphs given in above diagrams
	printf( "\nBridges in first graph \n" );
	Graph g1(5);
	g1.addEdge(1, 0);
	g1.addEdge(0, 2);
	g1.addEdge(2, 1);
	g1.addEdge(0, 3);
	g1.addEdge(3, 4);
	g1.bridge();

	printf(  "\nBridges in second graph \n" );
	Graph g2(4);
	g2.addEdge(0, 1);
	g2.addEdge(1, 2);
	g2.addEdge(2, 3);
	g2.bridge();

	printf( "\nBridges in third graph \n" );
	Graph g3(7);
	g3.addEdge(0, 1);
	g3.addEdge(1, 2);
	g3.addEdge(2, 0);
	g3.addEdge(1, 3);
	g3.addEdge(1, 4);
	g3.addEdge(1, 6);
	g3.addEdge(3, 5);
	g3.addEdge(4, 5);
	g3.bridge();
	return 0;
}

#endif









