
#ifndef MolecularGraph_h
#define MolecularGraph_h

#include <math.h>
#include "fastmath.h"
#include "Vec3.h"

class MolecularGraph{ public:
    int natom=0;
    int nbond=0;
    //Vec3d* apos=0;
    Vec2i* bond2atom=0;
    int*   ngNs=0;
    int*   ngIs=0;
    int*   atom2bond=0;
    int*   atom2neigh=0;
    int*   amask=0;
    int*   bmask=0;

    int    if0=0,if1=0,if2=0;
    int*   front_=0;
    //int*   front2=0;
    int    endAtom=0;
    bool   endFound=0;

    void bindOrRealloc(int natom_, int nbond_, Vec2i* bond2atom_=0 ){
        natom=natom_;
        nbond=nbond_;
        if(bond2atom_){ bond2atom=bond2atom_; }else{ realloc(bond2atom,nbond); }
        _realloc( ngNs,   natom );
        _realloc( ngIs,   natom );
        _realloc( amask,  natom );
        _realloc( bmask,  nbond );
        _realloc( front_, natom );
        //_realloc( front2, natom );
        _realloc( atom2bond, nbond*2 );
        _realloc( atom2neigh, nbond*2 );
    }

    void makeNeighbors(){
        // NOTE : this algorithm is very similar to atoms->cells
        // count number of neighbors
        for(int i=0; i<natom;i++){ ngNs[i]=0; };
        //printf( "nbond %i \n", nbond );
        for(int i=0; i<nbond;i++){
            const Vec2i& b=bond2atom[i];
            ngNs[ b.a ]++;
            ngNs[ b.b ]++;
        }
        // set initial position of neighbors
        int n=0;
        //ngIs[0]=n;
        for(int i=0; i<natom; i++){
            //printf( "atom[%i] N %i i0 %i \n", i, ngNs[i], n );
            ngIs[i]=n;
            n+=ngNs[i];
            ngNs[i]=0;
        }
        //printf( "nTOT : %i \n", n );
        // bond2atom -> atom2bond
        for(int i=0; i<nbond; i++){
            const Vec2i& b=bond2atom[i];
            int j;
            j=ngIs[b.a]+ngNs[b.a]; atom2bond[j]=i; atom2neigh[j]=b.b; ngNs[b.a]++;
            j=ngIs[b.b]+ngNs[b.b]; atom2bond[j]=i; atom2neigh[j]=b.a; ngNs[b.b]++;
        }
    }

    inline int getOtherAtom( int ib, int ia  ){
        const Vec2i& b=bond2atom[ib];
        return (b.a==ia) ? b.b : b.a;
    }
    void cleanMask(){ for(int i=0; i<natom; i++){ amask[i]=-1; }; }

    inline void addToFront(int ja, int color){
        amask[ja]=color;
        front_[if2]=ja;
        if2++;
        if(endAtom)endFound=true;
    }

    void addNeigh2front( int ia, int color, int bg){
        int  ng =ngNs[ia];
        int* ngs=atom2neigh+ngIs[ia];
        for(int i=0; i<ng; i++){
            int ja = ngs[i];
            if( amask[ja] == bg ){  addToFront(ja, color); }
        }
    }
    inline void clearFront(){ if0=0;if1=0;if2=0; }
    inline void swapFronts(){
        if0=if1;if1=if2;
        //_swap(front1,front2);
        //_swap(nfront1,nfront2);
    }

    void fillSubGraph( int ia, int color ){
        int if1=1;
        front_[if0]=ia;
        int bg = amask[ia];
        do{
            for(int i=if0; i<if1; i++){
                addNeigh2front( front_[i], color, bg);
            }
            swapFronts();
        }while(if1>if0);
    }

    void splitByBond(int ib, int color){
        const Vec2i& b=bond2atom[ib];
        amask[ b.b ]=color;
        endAtom=b.b;
        fillSubGraph( b.a, color );
    }

    void maskCaps( int color ){
        for(int i=0; i<natom; i++){
            if(ngNs[i]<2){
                amask[ i     ]=color;
                bmask[ngIs[i]]=color;
            }
        }
    }

    void selectBondNeighs(){

    }

    void printNeighs(){
        for(int i=0; i<natom; i++){
            int  ng =ngNs[i];
            int i0=ngIs[i];
            printf( "atom[%i] nneigh %i \n", i, ng );
            for(int j=0; j<ng; j++){
                printf( "   ia %i ib %i \n", atom2neigh[i0+j], atom2bond[i0+j] );
            }
        }
    }

};




//  algorithm copied from 
//    https://www.geeksforgeeks.org/bridge-in-a-graph/


// A C++ program to find bridges in a given undirected graph
//#include<iostream>
#include <list>
#define NIL -1
using namespace std;

// A class that represents an undirected graph
class Graph{ public:
	int V; // No. of vertices
	list<int> *adj; // A dynamic array of adjacency lists
    // Array of lists // neighbors
	
    void bridgeUtil(int v, bool visited[], int disc[], int low[], int parent[]);
	Graph(int V); // Constructor
	void addEdge(int v, int w); // to add an edge to graph
	void bridge(); // prints all bridges
};

Graph::Graph(int V){
	this->V = V;
	adj = new list<int>[V];
}

void Graph::addEdge(int v, int w){
	adj[v].push_back(w);
	adj[w].push_back(v); // Note: the graph is undirected
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
			if (low[v] > disc[u]) printf( "%i %i \n", u, v );   // If the lowest vertex reachable from subtree under v is below u in DFS tree, then u-v is a bridge
		} else if (v != parent[u])  low[u] = min(low[u], disc[v]);  // Update low value of u for parent function calls.
	}
}

// DFS based function to find all bridges. It uses recursive function bridgeUtil()
void Graph::bridge(){
	bool *visited = new bool[V];
	int *disc = new int[V];
	int *low = new int[V];
	int *parent = new int[V];
	for (int i = 0; i < V; i++)	{   // Initialize parent and visited arrays
		parent[i] = NIL;
		visited[i] = false;          // Mark all the vertices as not visited
	}
	for (int i = 0; i < V; i++) if (visited[i] == false) bridgeUtil(i, visited, disc, low, parent);  // Call the recursive helper function to find Bridges in DFS tree rooted with vertex 'i'
}

// Driver program to test above function
int test(){
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


