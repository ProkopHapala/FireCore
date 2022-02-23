
#ifndef arrayAlloc_h
#define arrayAlloc_h

#include <vector>

class ArrayAllocator{ public:
    // this will allocate continuous array from several sub-arrays
    // it is usefull for example if we work with multiple degrees of freedom used from different force-fields and we want share them all in one DynamicOptimizer
    //
    // Example:
    // -----------
    //  mmff.ioff = allocator.demand( mmff.getByteSize() );
    //  piff.ioff = allocator.demand( piff.getByteSize() );
    //  allocator.alloc( 3 );
    //  mmff.apos    = (Vec3d*)allocator.get( 0, piff.ioff   );
    //  mmff.aforce  = (Vec3d*)allocator.get( 1, piff.ioff   );
    //  piff. plane0 = (Vec3d*)allocator.get( 0, piff.ioff_0 );
    //  piff. planeh = (Vec3d*)allocator.get( 0, piff.ioff_h );
    //  piff.fplane0 = (Vec3d*)allocator.get( 1, piff.ioff_0 );
    //  piff.fplaneh = (Vec3d*)allocator.get( 1, piff.ioff_h );
    //  opt.pos      = (Vec3d*)allocator.get( 0, 0 );
    //  opt.force    = (Vec3d*)allocator.get( 1, 0 );
    //  opt.vel      = (Vec3d*)allocator.get( 2, 0 );

    int nbtot=0;
    std::vector<char*> arrays;
    std::vector<int>   offsets;
    int demand( int nbytes ){
        int off=offsets.size();
        offsets.push_back(nbtot);
        nbtot+=nbytes;
    }
    char* get(int islot, int ioff){
        return arrays[islot] + offsets[ioff];
    }
    void alloc( int nslot ){
        arrays.resize( nslot );
        for(int i=0; i<arrays.size(); i++){ arrays[i]=new char[nbtot]; };        
    }
    void dealloc(){
        for( char* ar : arrays){ delete [] ar; };        
    }
}

#endif









