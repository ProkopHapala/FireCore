
#ifndef CircularArray_h
#define CircularArray_h
#include <vector>

template<typename T>
class CircularArray { public:
    std::vector<T> data;
    int head=0; // Index of the head element
    int n=0;
    // how to know if T is a pointer type ?
    bool isArray = false;
    bool bOwner  = true;

    CircularArray(int capacity_, bool isArray_, bool bOwner_ ):data(capacity_), isArray(isArray_), bOwner(bOwner_), n(0){
        for(int i=0; i<data.size(); i++){ data[i] = 0; }
        head = data.size()-1;
    }
    ~CircularArray(){
        for(int i=0; i<data.size(); i++){
            try_dealloc_field( i%data.size() );
        }
    }

    int size()const{ return n; } 

    bool try_dealloc_field(int i){
        if( data[i] == 0 ) return false;
        if( std::is_pointer<T>::value && bOwner ){
            if(isArray){ delete[] data[i]; }else{ delete data[i]; }
            return true;
        }
        return false;
    }

    void push(T value){
        try_dealloc_field(head);
        int j = head;
        data[j] = value;
        //head = (head+1)%data.size(); // Wrap around
        head = (head-1); 
        if( n<(data.size()-1) ) n++;
        //printf( "CircularArray.push[j=%i] n=%i sz=%i head=%i `%s`\n", j, n, data.size(), head, value );
    }
    
    T get(int i)const{
        if((i<0)||(i>=n)) return 0;
        int j = (head+i+1)%data.size();
        //printf( "CircularArray.get[j=%i] i=%i n=%i size=%i head=%i `%s`\n", j, i, n, data.size(), head, data[j] );
        return data[ j ];
    }


};

#endif









