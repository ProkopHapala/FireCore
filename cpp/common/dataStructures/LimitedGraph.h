
#ifndef LimitedGraph_h
#define LimitedGraph_h

#include <vector>
#include <array>

#include "macroUtils.h"


//  Modified from:   https://www.geeksforgeeks.org/bridge-in-a-graph/
// see also:
//      * https://cp-algorithms.com/graph/bridge-searching.html
//      * https://en.wikipedia.org/wiki/Bridge_(graph_theory)
//      * https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
//      * Find Bridges in a graph using Tarjans Algorithm https://www.youtube.com/watch?v=Rhxs4k6DyMM

//using namespace std;

template<unsigned int M>
class LimitedGraph{ public:
    static constexpr const int m = M; // max number of neighs
	using Vmi = std::array<int,M>;
    int  n =0;
    int* nneighs=0;
    Vmi* neighs=0;

    std::vector<Vec2i>  found;

    bool *visited=0;
    int  *disc=0;
    int  *low=0;
    int  *parent=0;

    bool bPrint=false;

// ============= Functions

void realloc(int n_){
    n=n_;
    _realloc( nneighs,n );  
    for(int i=0; i<n; i++){ 
        nneighs[i]=0;       
    }
    _realloc( neighs ,n );  
    for(int i=0; i<n; i++){ 
        neighs[i].fill(-1); 
    }

    _realloc( visited,n );
    _realloc( disc   ,n );
    _realloc( low    ,n );
    _realloc( parent ,n );
}

void dealloc(){
    _dealloc( nneighs);
    _dealloc( neighs);

    _dealloc( visited );
    _dealloc( disc    );
    _dealloc( low     );
    _dealloc( parent  );
}

LimitedGraph()=default;
LimitedGraph(int n_){
    realloc(n_);
}

bool addNeigh( int i, int j ){
    int ni=nneighs[i];
    if( ni<m ){ neighs[i][ni]=j; nneighs[i]=ni+1; return false; };
    return true;
}

bool addEdge( int i, int j ){ return addNeigh( i, j ) ||  addNeigh( j, i ); }

bool fromBonds( int nb, const Vec2i* bonds, bool bIgnoreAtoms, bool bExitOnErro=true ){
    bool err=0;
    for(int i=0; i<nb; i++){
        Vec2i b = bonds[i];
        if( (b.a>=n) || (b.b>=n) ) if(bIgnoreAtoms){ continue; }else[[unlikely]]{ printf( "ERROR in LimitedGraph::fromBonds() bond[%i](%i,%i) out of range(0..%i) \n => Exit(); \n", i, b.a,b.b, n ); exit(0); };
        bool err = addEdge(b.a, b.b);
        if(err && bExitOnErro)[[unlikely]]{ printf( "ERROR in LimitedGraph::fromBonds() cannot add bond[%i] neighs are filled nng[%i]=%i nng[%i]=%i  \n => Exit(); \n", i, b.a,nneighs[b.a], b.b,nneighs[b.b] ); exit(0); };
    }
    return err;
}

// A recursive function that finds and prints bridges using // DFS traversal
// u         --> The vertex to be visited next
// visited[] --> keeps track of visited vertices
// disc[]    --> Stores discovery times of visited vertices
// parent[]  --> Stores parent vertices in DFS tree
void bridgeUtil(int i, bool visited[], int disc[],  int low[], int parent[] ) {
	static int time  = 0;       
	visited[i] = true;    
	disc[i]    = low[i] = ++time;      
    int ni     = nneighs[i];
    const Vmi& ngsi = neighs[i]; 
	for(int ij=0; ij<ni; ij++){
		int j = ngsi[ij];
		if (!visited[j]){      
			parent[j] = i;
			bridgeUtil( j, visited, disc, low, parent );
			low[i] = _min(low[i], low[j]);         // Check if the subtree rooted with v has a connection to one of the ancestors of u
			if (low[j]>disc[i]){                   // If the lowest vertex reachable from subtree under v is below u in DFS tree, then u-v is a bridge
                if(bPrint)printf( "bridge %i %i \n", i,j );
                found.push_back( (Vec2i){i,j} );   // add bridge to found list
            }
		} else if ( j!= parent[i])  low[i] = _min(low[i], disc[j]);  // Update low value of u for parent function calls.
	}
}

// DFS based function to find all bridges. It uses recursive function bridgeUtil()
void bridge(){
    for (int i=0;i<n;i++){
        parent [i]  = -1;
        visited[i] = false;
    }
    for (int i=0;i<n;i++) if(visited[i]==false) bridgeUtil(i, visited, disc, low, parent);  // Call the recursive helper function to find Bridges in DFS tree rooted with vertex 'i'
}


void print(){
    printf( "LimitedGraph::print() n %i \n", n );
    for(int i=0; i<n; i++){
        printf( "node[%i] : ", i );
        for(int j=0; j<nneighs[i]; j++){ printf( " %i", neighs[i][j] ); }
        printf( "\n" );
    }
}

}; // ---- end class LimitedGraph


bool test_LimitedGraph_findBridge(){   // Driver program to test above function
	// Create graphs given in above diagrams
	printf( "\nBridges in first graph \n" );
	LimitedGraph<4> g1(5);
	g1.addEdge(1, 0);
	g1.addEdge(0, 2);
	g1.addEdge(2, 1);
	g1.addEdge(0, 3);
	g1.addEdge(3, 4);
	g1.bridge();

    /*
	printf(  "\nBridges in second graph \n" );
	LimitedGraph<4> g2(4);
	g2.addEdge(0, 1);
	g2.addEdge(1, 2);
	g2.addEdge(2, 3);
	g2.bridge();

	printf( "\nBridges in third graph \n" );
	LimitedGraph<4> g3(7);
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
    */

   return true;
}


#endif









