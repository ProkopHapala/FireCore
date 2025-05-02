
#ifndef  containers_h
#define  containers_h
/// @file containers.h @brief contains Buckets class, for accelarating neighbourhood search for particle-in-cell on rectangular 3D grids.
/// @ingroup Utilities

#include <string>
#include <vector>
#include <unordered_map>

typedef std::unordered_map<int,int>         IntMap;
typedef std::unordered_map<std::string,int> Dictionary;

// === Convenience functions for std::unordered_map

template<typename K,typename V>
V get( const std::unordered_map<K,V>& map, const K& key, const V& def ){
    auto it = map.find(key);
    if( it != map.end() ){ return it->second; }else{ return def; }
}

template<typename K,typename V>
bool set( std::unordered_map<K,V>& map, const K& key, const V& val, bool bReplace=false ){
    auto it = map.find(key);
    if( it != map.end() ){ if(bReplace){it->second = val;} return true; }else{ map.insert({key,val}); return false; }
}

template<typename K,typename V>
V* setp( std::unordered_map<K,V*>& map, const K& key, const V* val, bool bReplace=false ){
    auto it = map.find(key);
    if( it != map.end() ){ V* ret=it->second; if(bReplace){it->second=val;} return ret; }else{ map.insert({key,val}); return 0; }
}

// ================== Ditionary

/// @brief Dictionary with string keys and int values, generaly used as short-hand for python like dictionary
/// @tparam T type of the values stored in the dictionary (beside the index)
template<typename T>
class Dict{ public:
    std::unordered_map<std::string,int> map;
    std::vector<T>                     vec;

    int getId( const std::string& name, const int def = -1 )const{
        auto it = map.find(name);
        if( it != map.end() ){ return it->second; }else{ return def; }
    }

    T get( const std::string& name, const T def )const{
        auto it = map.find(name);
        if( it != map.end() ){ return vec[it->second]; }else{ return def; }
    }

    int getTo( const std::string& name, T& out )const{
        auto it = map.find(name);
        if( it != map.end() ){ out=vec[it->second]; return it->second; }else{ return -1; }
    }

    int add( const std::string& name, T t,  bool bDel=true, bool bReplace=true ){
        auto it = map.find( name );
        if( it == map.end() ){    // not found
            int i=vec.size(); vec.push_back(t); 
            map.insert({name,i}); 
            return i; 
        }else{    // found
            if constexpr (std::is_pointer<T>::value){
                if(bDel)delete vec[it->second]; 
            }
            if(bReplace){vec[it->second] = t; }
            return -it->second-1; 
        }
    }

    bool insert( const std::string& name, T t, int i, bool bDel=true ){
        if( vec.size(i)<=i ) vec.resize(i+1);
        vec[i] = t;
        auto it = map.find(name);
        if( it != map.end() ){  // found
            if constexpr (std::is_pointer<T>::value){
                if( i != it->second ){     
                    if(bDel) vec[it->second]=0; 
                }
            }
            map[name] = i; 
            vec[i]    = t; 
            return true;
        }else{ 
           map.insert({name,i}); 
           return false; 
        }
    }

    void rebind( const std::string& name, int i ){
        auto it = map.find(name);
        if( it != map.end() ){ 
            it->second = i; 
        }else{ 
            map.insert({name,i}); 
        }
    }

    bool clear(){
        map.clear();
        vec.clear();
        return true;
    }

};


#endif
