
#ifndef  macroUtils_h
#define  macroUtils_h

#include <vector>

//#undef  DEBUG_ALLOCATOR
//#define DEBUG_ALLOCATOR
#define _delete(p) delete[] p
#define _new(T,n)  new T[n]
#ifdef DEBUG_ALLOCATOR
    #include "debugAllocator.h"
#endif 


#define SWAP( a, b, TYPE ) { TYPE t = a; a = b; b = t; }
//#define _max(a,b)      ((a>b)?a:b)
//#define _min(a,b)      ((a<b)?a:b)
//#define _abs(a)        ((a>0)?a:-a)
//#define _clamp(x,a,b)  max(a, min(b, x))
//#define _clamp(n, lower, upper) if (n < lower) n= lower; else if (n > upper) n= upper
//#define _clamp(a, lower, upper) ((a>lower)?((a<upper)?a:upper):lower)

#define _minit( i, x, imin, xmin )  if( x<xmin ){ xmin=x; imin=i; }
#define _maxit( i, x, imax, xmax )  if( x>xmax ){ xmax=x; imax=i; }

#define _setmin( xmin, x )  if( x<xmin ){ xmin=x; }
#define _setmax( xmax, x )  if( x>xmax ){ xmax=x; }

#define _circ_inc( i, n )   i++; if(i>=n) i=0;
#define _circ_dec( i, n )   i--; if(i< 0) i=n-1;

//#define _realloc(TYPE,arr,n){ if(var) delete [] arr; arr=new TYPE[n]; }

// Get byte offset of a member in a struct ( call from outside the struct )
#define mapByteOffset(offset_map,instance,member) { offset_map[#member] = ((char*)&instance.member) - ((char*)&instance); }
// Get byte offset of a member in a struct ( call from inside the struct )
#define mapByteOffIn(offset_map,member) { offset_map[#member] = ((char*)&(this->member)) - ((char*)this); }


// ============= sorting
inline int selectMinHigher(int a0, int n, int* arr){
    int amin=0x7FFFFFFF; // 32 int max
    int imin=0;
    for(int i=0;i<n;i++){
        int a=arr[i];
        if((a>a0)&&(a<amin)){amin=a;imin=i;};
    }
    return imin;
};


template<typename  T>
struct Buf{
    T*   data=0;
    bool own=false;

    void bind(T* data_){ data=data_; };
    void alloc(int n  ){ own=true; data=new T[n]; };
    //bool realloc(int n){ bool b=false; if(data){ delete [] data; b=true; }; data=new T[n]; return b; };
    //Arr(T* data_){own=true; };
    Buf() = default;
    Buf(int n   ):own{true} {data=new T[n]; };
    Buf(T* data_):own{false},data{data_}{};
    ~Buf(){ if(own&&data)delete[]data; }

    T& operator[](int i){
        //if (index >= size) { printf( "Array index [%i] out of bound [%i] \n", i, n ); exit(0); }
        return data[i];
    }
};

template<typename  T>
struct Arr{
    T*   data;
    int  n;
    bool own;

    //Arr(T* data_){own=true; };
    Arr() = default;
    Arr(int n_           ){ own=true;  n=n_; data=new T[n]; };
    Arr(int n_ , T* data_){ own=false; n=n_; data=data_;    };
    ~Arr(){ if(own&&data)delete[]data; }

    T& operator[](int i){
        //if (index >= size) { printf( "Array index [%i] out of bound [%i] \n", i, n ); exit(0); }
        return data[i];
    }
};

//tempate<typename T> bool addFirstEmpty( T* arr, n, T what, T empty=-1 ){
inline bool addFirstEmpty( int* arr, int n, int what, int empty=-1 ){
    for(int i=0; i<n; i++){
        if(arr[i]==empty){ arr[i]=what; return true; }
    }
    return false;
};


#define WITH(x) auto& _=x;

#define _forN(i,n)         for(int i=0 ;i<n;i++)
#define _for0N(i,i0,n)     for(int i=i0;i<n;i++)
#define _for0Nd(i,i0,n,d)  for(int i=i0;i<n;i+=d)

//#define _list2array(T,n,lst,to) { T tmp[]=lst; for(int i=0 ;i<n;i++)to[i]=tmp[i]; }

#define _template_Func   template<typename Func>
#define _template_T      template<typename T>
#define _template_N      template<size_t N>
#define _template_TN     template<typename T,size_t N>
#define _inline_T        template<typename T> inline

#define _checkNull(var) if(var==nullptr){ printf("ERROR %s is NULL => Exit()\n", #var ); exit(0); }

template<typename T> void _vec2arr(T*out, const std::vector<T>& v){ for(int i=0;i<v.size();i++)out[i]=v[i];} 


_inline_T void _swap  (T& a, T& b) { T t=a; a=b; b=t;   }
_inline_T void _order (T& a, T& b) { if(a>b){T t=a; a=b; b=t; }; }
//_inline_T void _order (T& a, T& b) { if(a>b)_swap(a,b); }

_inline_T const T& _min  (const T& a, const T& b){ return (a>b)?b:a; }
_inline_T const T& _max  (const T& a, const T& b){ return (a<b)?b:a; }
_inline_T const T& _clamp(const T& a, const T& amin, const T& amax){ return _max(amin,_min(amax,a)); }

_inline_T T        _abs  (const T& a ){ return !(a<0)?a:-a; }
_inline_T int      signum(T val)      { return (T(0) < val) - (val < T(0)); }

// ======= allocation



// _inline_T bool _allocIfNull   (T*& arr, int n){ if(arr==0){ arr=_new(T,n); return true; } return false; }
// _inline_T void _alloc         (T*& arr, int n){ arr=_new(T,n); }
// _inline_T void __realloc       (T*& arr, int n){ if(arr){ _delete(arr); } arr=_new(T,n); }
// _inline_T void __realloc0      (T*& arr, int n, const T& v0){  _realloc(arr,n); for(int i=0;i<n;i++){ arr[i]=v0; } }
// _inline_T void _dealloc       (T*& arr       ){ if(arr){ _delete(arr); } arr=0;        }
// _inline_T bool _bindOrRealloc (int n, T* from, T*& arr ){ if(from){ arr=from; return false; }else{ _realloc(arr,n); return true; } }
// _inline_T T*   _allocPointer  (T**& pp, int n ){  if(pp){ if((*pp)==0)(*pp)=_new(T,n); return *pp; }; return 0; };
// _inline_T T*   _reallocPointer(T**& pp, int n){
//     if (pp){
//         if ((*pp))   _delete((*pp));
//         (*pp) = _new(T,n);
//         return (*pp);
//     }
//     return 0;
// }

// _inline_T bool _allocIfNull  (T*& arr, int n){ if(arr==0){ arr=new T[n]; return true; } return false; }
// _inline_T void _alloc        (T*& arr, int n){ arr=new T[n]; }
// _inline_T T* __realloc       (T*& arr, int n){ if(arr){ delete [] arr; } arr=new T[n]; return arr;  }
// _inline_T T* __realloc0      (T*& arr, int n, const T& v0){  __realloc(arr,n); for(int i=0;i<n;i++){ arr[i]=v0; } return arr; }
// _inline_T bool _dealloc      (T*& arr       ){ if(arr){ delete [] arr; arr=0; return true; } return false;  }
// _inline_T bool _bindOrRealloc(int n, T* from, T*& arr ){ if(from){ arr=from; return false; }else{ __realloc(arr,n); return true; } }
// _inline_T T* _allocPointer   (T**& pp, int n ){  if(pp){ if((*pp)==0)(*pp)=new T[n]; return *pp; }; return 0; };
// _inline_T T* _reallocPointer (T**& pp, int n){
//     if (pp) {
//         if ((*pp)) delete[] (*pp);
//         (*pp) = new T[n];
//         return (*pp);
//     }
//     return 0;
// }




_inline_T T*   __alloc       (T*& arr, int n){ arr=new T[n]; return arr;  }
_inline_T T*   __realloc     (T*& arr, int n){ if(arr){ _delete(arr); } arr=new T[n]; return arr;  }
_inline_T T*   __realloc0    (T*& arr, int n, const T& v0){  __realloc(arr,n); for(int i=0;i<n;i++){ arr[i]=v0; } return arr; }
_inline_T T*   __allocIfNull (T*& arr, int n){ if(arr==0){ arr=new T[n]; return arr; } return 0; }
_inline_T bool __dealloc     (T*& arr      ){ if(arr){ _delete(arr); arr=0; return true; } return false;  }
_inline_T bool _bindOrRealloc(int n, T* from, T*& arr ){ if(from){ arr=from; return false; }else{ __realloc(arr,n); return true; } }
_inline_T T* __allocPointer   (T**& pp, int n ){  if(pp){ if((*pp)==0)(*pp)=new T[n]; return *pp; }; return 0; };
_inline_T T* __reallocPointer (T**& pp, int n){
    if (pp) {
        if ((*pp)) delete[] (*pp);
        (*pp) = new T[n];
        return (*pp);
    }
    return 0;
}

#ifdef DEBUG_ALLOCATOR
#define _alloc(arr,n)          debug_alloc_store( __alloc(arr,n), n, _CODE_LOCATION )
#define _realloc(arr,n)        debug_alloc_store( __realloc(arr,n), n, _CODE_LOCATION )
#define _realloc0(arr,n,v0)    debug_alloc_store( __realloc0(arr,n,v0), n, _CODE_LOCATION )
#define _allocIfNull(arr,n)    debug_alloc_store( __allocIfNull(arr,n), n, _CODE_LOCATION )
#define _allocPointer(arr,n)   debug_alloc_store( __allocPointer(arr,n), n, _CODE_LOCATION )
#define _reallocPointer(arr,n) debug_alloc_store( __reallocPointer(arr,n), n, _CODE_LOCATION )
#define _dealloc(arr)          debug_dealloc( debugAllocator,arr, _CODE_LOCATION )
// fill in the missing
#else
#define _alloc(arr,n)          __alloc(arr,n)
#define _realloc(arr,n)        __realloc(arr,n)
#define _realloc0(arr,n,v0)    __realloc0(arr,n,v0)
#define _allocIfNull(arr,n)    __allocIfNull(arr,n)
#define _allocPointer(arr,n)   __allocPointer(arr,n)
#define _reallocPointer(arr,n) __reallocPointer(arr,n)
#define _dealloc(arr)          __dealloc(arr)
// fill in the missing
#endif // DEBUG_ALLOCATOR;







_inline_T  bool _clone( int i0, int imax, T* from, T*& arr, int n){
    bool ret = _allocIfNull(arr,n);
    for(int i=i0; i<imax; i++){ arr[i]=from[i-i0]; } // use mem copy instead ?
    return ret;
}
_inline_T bool _set( int i0, int imax, const T& from, T*& arr, int n){
    bool ret = _allocIfNull(arr,n);
    for(int i=i0; i<imax; i++){ arr[i]=from; }
    return ret;
}


_inline_T void copy( int n, int pitch1, int offset1, T* from, int pitch2, int offset2, T* out ){
    for(int i=0; i<n; i++ ){  out[pitch2*i+offset2] = from[pitch1*i+offset1]; }
}

#endif
