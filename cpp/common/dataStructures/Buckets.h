#ifndef  Buckets_h
#define  Buckets_h
/// @file Buckets.h @brief contains Buckets class, which is a class for storing indices of objects in buckets (groups) used typically to accelarate neighbourhood search
/// @defgroup Neighbours  Neighbours
/// @addtogroup Neighbours
/// @{

#include "macroUtils.h"

/// @brief Class for storing indices of objects in buckets (groups) used typically to accelarate neighbourhood search
/// @details contains bi-dierectional mapping between objects and buckets, 1 object is only in 1 bucket, 1 bucket may contain many objects, maping stored in simple array, with offset (start,end) of each bucket stored separately
class Buckets{ public:

    int  maxInBucket=0; /// maximum allowed number of objects in any bucket
    int  nobjSize=-1;

    int ncell,nobj;
    int* cellNs=0;   /// [ncell] number of objects contained in given cell, it may be not necessary as it can be computed as cellI0s[i+1]-cellI0s[i], but it is more convenient to have it here
    int* cellI0s=0;  /// [ncell] index of first object contained given cell in the array cell2obj
    int* cell2obj=0; /// [nobj]  indices of objects contained in given cell
    int* obj2cell=0;   // [nobj] index of cell containing given object, if object is not in any cell, then obj2cell[i]==-1
    //int  nobj_bind=-1;
    //int* obj2cell_bind=0; 

    // =========== Functions

    inline int addToCell( int icell, int iobj ){
        if(icell < 0) return -1;  // Don't add objects to negative cells
        int j = cellI0s[icell] + cellNs[icell];
        cell2obj[j] = iobj;
        cellNs[icell]++;
        return j;
    }

    inline void clean(){ for(int k=0; k<ncell; k++ ){ cellNs[k]=0; } }
    inline void cleanO2C( int icell=-1 ){ for(int i=0; i<nobj; i++ ){ obj2cell[i]=icell; } }

    /**
     * Counts the number of objects in each cell based on the given mapping. It stores the result in the array cellNs.
     *
     * @param nobj The number of objects.
     * @param obj2cell An array mapping each object to its corresponding cell.
     */
    inline void count( int nobj, int* obj2cell ){ 
        //printf( "Buckets::count() nobj=%i ncell=%i @obj2cell=%li \n", nobj, (long)obj2cell );
        for(int i=0; i<nobj; i++ ){ 
            //printf( "Buckets::count()[%i] \n", i );
            int ic = obj2cell[i]; 
            //if( (ic<0)||(ic>=ncell) ){ printf( "Buckets::count() ERROR i %i ic %i ncell %i \n", i, ic, ncell ); }
            //printf( "obj[%i ]-> cell %i \n", i, ic );
            if(ic>=0) cellNs[ ic ]++;  
        } 
    }

    /// @brief Updates the offsets of the buckets (cellI0s) based on the number of elements in each bucket (cellNs). It also determines the maximum number of elements in a bucket (maxInBucket).
    inline void updateOffsets(){   
        int ntot=0;
        maxInBucket=0;
        for(int k=0; k<ncell; k++ ){
            cellI0s[k] = ntot;
            int ni    = cellNs[k];
            ntot      += ni;
            if( ni>maxInBucket ) maxInBucket=ni;     // ? check most occupied cell
            cellNs[k]=0;
        }
    }

    /**
     * @brief Assigns objects to cells based on the given mapping. The result is stored in the array cell2obj.
     * 
     * @param nobj The number of objects.
     * @param obj2cell An array representing the mapping of objects to cells.
     */
    inline void objectsToCells( int nobj, int* obj2cell ){ 
        for(int i=0; i<nobj; i++){
            //printf( "[%i] ", i );
            int  ic =  obj2cell[i];
            if(ic<0) continue;
            if(ic>=ncell) {
                printf("ERROR: objectsToCells() ic=%i >= ncell=%i\n", ic, ncell);
                continue;
            }
            int& ni = cellNs[ic];
            int j   =  cellI0s[ic] + ni;
            if(j >= nobjSize) {
                printf("ERROR: objectsToCells() j=%i >= nobjSize=%i\n", j, nobjSize);
                continue;
            }
            //printf( " k %i j %i | nobj %i \n", k, j, nobj );
            cell2obj[j] = i;
            ni++;
        }
    }

    /// @brief Checks if every object from obj2cell can be found in cell2obj .
    int  checkObj2Cell( bool bPrint=false, bool bExit=true){
        int nbad=0;
        for(int i=0; i<nobj; i++){
            int ic = obj2cell[i];
            if(ic < 0) continue;  // Skip negative indices
            // Make sure this object appears in the corresponding cell
            bool found = false;
            int i0 = cellI0s[ic];
            int n  = cellNs [ic];
            for(int j=0; j<n; j++){
                if(cell2obj[i0+j] == i){
                    found = true;
                    break;
                }
            }
            if(!found){
                if(bPrint)printf("ERRORin Buckets::checkObj2Cell() object %i is assigned to cell %i but not found in cell list\n", i, ic);
                if(bExit) exit(1);
                nbad++;
            }
        }
        return nbad;
    }

    /// @brief Checks if all items stored in cells (cell2obj) can be found in obj2cell with the same cell index
    int  checkCell2Obj( int nobj_=-1, bool bPrint=false, bool bExit=true){
        if(nobj_<0) nobj_ = nobj;
        int nbad=0;
        for(int i=0; i<ncell; i++){
            int ni= cellNs[i];
            int i0= cellI0s[i];
            if(ni<=0) continue;
            
            for(int j=0; j<ni; j++){
                int iobj = cell2obj[i0+j];
                // Check if the object index is valid
                if( (iobj < 0) || (iobj >= nobj_) ){
                    if(bPrint)printf("ERRORin Buckets::checkCell2Obj() cell[%i|%i] has object index %i out of bounds [0..%i]\n", i, j, iobj, nobj_);
                    if(bExit) exit(1);
                    nbad++;
                    continue; // Skip to next object in cell
                }
                
                // Check if the object's assigned cell matches the current cell
                if(obj2cell[iobj] != i){
                    if(bPrint)printf("ERRORin Buckets::checkCell2Obj() cell[%i|%i] contains object %i but obj2cell[%i]=%i\n", 
                                     i, j, iobj, iobj, obj2cell[iobj]);
                    if(bExit) exit(1);
                    nbad++;
                }
            }
        }
        return nbad;
    }

    /**
     * @brief Updates the cells based on the given object-to-cell mapping.
     * This function performs the following steps:
     * 1. Cleans the cells.
     * 2. Counts the number of objects in each cell.
     * 3. Updates the offsets of the cells.
     * 4. Assigns objects to their respective cells.
     *
     * @param nobj_     number of objects.           If not specified, it uses the default number of objects.
     * @param obj2cell_ object-to-cell mapping.  If not specified, it uses the default mapping. (which have to be allocated or binded before).
     */
    inline void updateCells( int nobj_=-1, int* obj2cell_=0 ){
        if( obj2cell_==0 ) obj2cell_=obj2cell;
        if( nobj_<0      ) nobj_    =nobj;
        
        // Initialize cell2obj with -1 to ensure we don't have garbage values
        for(int i=0; i<nobjSize; i++){
            cell2obj[i] = -1;
        }
        
        clean         ();                      //printf("updateCells.clean() \n"  );
        count         ( nobj_, obj2cell_ );    //printf("updateCells.count() \n"  );
        updateOffsets ();                      //printf("updateCells.updateOffsets() \n"  );
        objectsToCells( nobj_, obj2cell_ );    //printf("updateCells.objectsToCells() \n"  );
    }

    inline bool resizeCells( int ncell_                 ){ bool b=(ncell_>ncell); if(b){ ncell=ncell_; _realloc(cellNs,ncell_); _realloc(cellI0s,ncell_); }; return b; };
    inline bool resizeObjs ( int nobj_, bool bO2C=false ){ 
        bool b=(nobj_>nobjSize);  // need to resize
        nobj=nobj_;  // true number of objects ( if nobj decreases we do not need to resize, but we want to itereate over just the used part of the array )
        if(b){       // if need to resize
            //printf( "Buckets::resizeObjs() nobjSize(%i) -> nobj_(%i) \n", nobjSize, nobj_ );
            nobjSize =nobj_; 
            // Initialize cell2obj with -1 to mark as unassigned
            _realloc0(cell2obj,nobjSize,-1); 
            if(bO2C)_realloc0(obj2cell,nobj_,-1);
        };
        //printf( "Buckets::resizeObjs() nobjSize(%i) nobj_(%i) \n", nobjSize, nobj_ ); 
        return b; 
    };
    inline void bindObjs   ( int nobj_, int* obj2cell_  ){ resizeObjs ( nobj_, false ); obj2cell=obj2cell_; }

    inline void realloc( int ncell_, int nobj_, bool bO2C=false ){ 
        ncell=-1; resizeCells(ncell_); 
        nobj=-1;  resizeObjs(nobj_,bO2C); 
    }

    inline int getInCell( int icell, int* out ){    
        if(icell < 0) return 0;  // Return 0 objects for negative cells
        // why we copy? we can just return pointer to cell2obj[ cellI0s[icell] ]
        const int i0 = cellI0s[icell];
        const int n  = cellNs [icell];
        const int i1 = i0+n;
        for(int i=0; i<n; i++){
            //printf("getInCell(%i) i %i(%i|%i) \n", icell, i, i0,i1);
            //printf("getInCell[%i] out[%i]=%i \n", icell, i, cell2obj[i] );
            out[i] = cell2obj[i+i0];
        }
        return n;
    }
    inline int* inCell( int icell ){ return cell2obj + cellI0s[icell]; }

    inline void printObjCellMaping(int i0=-1, int i1=-1, bool bBoundary=true){
        if(i0<0)i0=0; 
        if(i1<0)i1=nobj;
        printf( "Buckets::printObjCellMaping() ncell %i nobj %i \n", ncell, nobj );
        for(int i=i0; i<i1; i++){
            int ic = obj2cell[i];
            if(ic < 0) {
                // Skip negative cell indices, but optionally print them with a note
                // printf( "[%i] o2c %i (skipped - negative cell index)\n", i, ic );
                continue;
            }
            
            // Check if ic is within valid range
            if(ic >= ncell) {
                printf( "[%i] o2c %i (ERROR - cell index out of range)\n", i, ic );
                continue;
            }
            
            // Now find where this object is stored in the cell2obj array
            bool found = false;
            int location = -1;
            int i0_cell = cellI0s[ic];
            int n_cell = cellNs[ic];
            
            // Print basic information
            printf( "[%i] o2c %i (cell offset: %i, cell size: %i)\n", i, ic, i0_cell, n_cell );
            
            // Optionally debug the cell's content
            /*
            for(int j=0; j<n_cell; j++) {
                int obj_idx = cell2obj[i0_cell + j];
                if(obj_idx == i) {
                    found = true;
                    location = i0_cell + j;
                    printf( "  - Found at position %i in cell2obj\n", location );
                }
            }
            
            if(!found) {
                printf( "  - ERROR: Object not found in its assigned cell\n" );
            }
            */
        }
    }

    inline void printCells(int verb=1){
        printf( "Buckets::printCells() ncell %i nobj %i \n", ncell, nobj );
        for(int i=0; i<ncell; i++){
            printf( "cell[%i] n %i i0 %i \n", i, cellNs[i], cellI0s[i] );
            if(verb>0){
                int i0=cellI0s[i];
                int ni=cellNs [i];
                int i1=i0+ni;
                for(int j=0; j<ni; j++){
                    int obj_idx = cell2obj[i0+j];
                    if(obj2cell[obj_idx] < 0){
                        printf( "ERROR in printCells() i: %i obj_idx %i \n", i, obj_idx );
                        exit(0);
                    };
                    printf( "cell2obj[%i|%i,%i] = %i \n", i0+j, i,j, obj_idx );
                }
            }
        }
    }


};

/// @}
#endif
