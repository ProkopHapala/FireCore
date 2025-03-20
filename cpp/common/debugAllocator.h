
#pragma once

#include <vector>
#include <unordered_map>
#include <string>

#include <cstdlib>
#include <stdio.h>

struct AllocItem{
    void*  ptr;     // pointer to allocated memory
    size_t n;       // number of elements
    size_t size_of; // size of each element
    std::string caption; // 
    AllocItem() = default;
    AllocItem( void* ptr_, size_t n_, size_t size_of_, const std::string& caption_ ): ptr(ptr_), n(n_), size_of(size_of_), caption{caption_}{};
};


class DebugAllocator { public:
    std::vector       <AllocItem>  items;
    std::unordered_map<void*, int> ptr2item;

    
    bool bPrintAlloc   = true;
    bool bPrintDealloc = true;
    bool bExitOnError = true;
    bool bLazyDealloc = true;

    inline int store( const AllocItem item ){
        if( bPrintAlloc ) printf( "DEBUG_ALLOCATOR: allocating @ %p n: %i at %s \n", item.ptr, item.n, item.caption.c_str() );
        items.push_back( item );
        int i = items.size()-1;
        ptr2item[item.ptr] = i;
        return i;
    }

    int add( void* ptr, size_t n, size_t size_of, const std::string& caption ){
        return store( AllocItem( ptr, n, size_of, caption ) );
    }

    void remove( void* ptr ){ 
        auto it = ptr2item.find( ptr );
        if( it == ptr2item.end() ){
            printf( "DebugAllocator::dealloc( %p ) ERROR : pointer not found \n", ptr );
            if(bExitOnError){ 
                print_allocated();
                exit(0);
            }
            return;
        };
        int i = it->second;
        if(bLazyDealloc){
            items[i].ptr = 0;
        }else{
            items.erase( items.begin() + i );
        }
        ptr2item.erase( ptr );
    }



    // template <typename T> T* alloc( size_t n, const std::string& caption ){
    //     T* ptr = new T[n];
    //     store( AllocItem( ptr, n, sizeof(T), caption ) );
    //     return ptr;
    // }

    // template <typename T> void dealloc( T* ptr ){ 
    //     auto it = ptr2item.find( ptr );
    //     if( it == ptr2item.end() ){
    //         printf( "DebugAllocator::dealloc( %p ) ERROR : pointer not found \n", ptr );
    //         if(bExitOnError) exit(0);
    //         return;
    //     };
    //     int i = it->second;
    //     if(bLazyDealloc){
    //         items[i].ptr = 0;
    //     }else{
    //         items.erase( items.begin() + i );
    //     }
    //     ptr2item.erase( ptr );
    //     delete [] ptr; 
    // }

    inline void clean_null_pointers(){
        // clean up the items storing null pointers
        int n=items.size();
        for(int i=0; i<n; i++){
            // if ptr==null we replace it with the last item and decrease the size of the vector
            while(items[i].ptr==0){
                items[i] = items[n-1];
                n--;
            }
        }
        items.resize(n);
    }

    // void clear(){
    //     for(int i=0; i<items.size(); i++){
    //         delete [] items[i].ptr;
    //     }
    //     items.clear();
    //     ptr2item.clear();
    // }

    inline void print_allocated(){
        printf( "DebugAllocator::print_allocated() \n" );
        for(int i=0; i<items.size(); i++){
            printf( "[%i] ptr @ %p n_items: %li item_size: %li allocated_at %s \n", i, items[i].ptr, items[i].n, items[i].size_of, items[i].caption.c_str() );
        }
    }

};

static DebugAllocator* debugAllocator = 0; 

inline void debugAllocator_init(){
    if(debugAllocator) return;
    debugAllocator = new DebugAllocator();
}

template <typename T> T* debug_alloc_store( T* ptr, size_t n, const std::string& caption ){
    if( ptr==0) return 0;
    if(debugAllocator){ 
        debugAllocator->add( ptr, n, sizeof(T), caption );
    }else{
        printf( "ERROR in debug_alloc() debugAllocator == null \n" );
        exit(0);
    }
    return ptr;
}

template <typename T> T* debug_alloc( DebugAllocator* debugAllocator,  size_t n, const std::string& caption ){
    T* ptr = new T[n];
    if( ptr==0) return 0;
    if(debugAllocator){ 
        debugAllocator->add( ptr, n, sizeof(T), caption );
    }else{
        printf( "ERROR in debug_alloc() debugAllocator == null \n" );
        exit(0);
    }
    return ptr;
}

template <typename T> void debug_dealloc( DebugAllocator* debugAllocator, T* ptr, const char* caption=0 ){ 
    if( ptr==0 ) return;
    if(debugAllocator){ 
        if(caption && debugAllocator->bPrintDealloc ) printf( "DEBUG_ALLOCATOR: deallocating @ %p at %s \n", ptr, caption );
        debugAllocator->remove( ptr );
    }else{
        printf( "ERROR in debug_dealloc() debugAllocator == null \n" );
        exit(0);
    }
}

inline void debugAllocator_print( ){
    if(debugAllocator){ 
        debugAllocator->print_allocated();
    }else{
        printf( "ERROR in debugAllocator_print() debugAllocator == null \n" );
        exit(0);
    }
} 



#define _STRINGIFY(x) #x
#define _TOSTRING(x) _STRINGIFY(x)
#define _CODE_LOCATION __FILE__ ":" _TOSTRING(__LINE__) 
//#define _CODE_LOCATION __FILE__ ":" _TOSTRING(__LINE__) __FUNCTION__

#ifdef DEBUG_ALLOCATOR
//#define _new(T,n ) debugAllocator->alloc<T>(n,_CODE_LOCATION )
// #define _delete(p) debugAllocator->dealloc(p)
//#define _new(T,n)  debug_alloc<T>(debugAllocator,n,_CODE_LOCATION ) 
#define _delete(p) debug_dealloc(debugAllocator,p)
#else
//#define _new(T,n)  new T[n]
#define _delete(p) delete[] p
#endif // DEBUG_ALLOCATOR;


