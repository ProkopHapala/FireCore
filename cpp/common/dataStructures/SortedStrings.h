
#ifndef SortedStrings_h
#define SortedStrings_h


#include <vector>
//#include <iostream>
#include <string>

class SortedStrings{ public:

    std::vector<std::string> table;

    void sort() {
        // can we enforce isertion sort ?
        std::sort(table.begin(), table.end());
    }

    // Helper function to compare prefix with strings in the list
    static bool prefixMatches(const std::string& prefix, const std::string& str) {
        if (str.size() < prefix.size()) return false;
        return std::equal(prefix.begin(), prefix.end(), str.begin());
    }

    // Binary search to find the start position of matches
    int findPrefixStart( const std::string& s ){
        int lo = 0, hi = table.size() - 1;
        //printf( "findPrefixStart() lo(%i) hi(%i) table.size(%i)  s(%s)\n", lo, hi, table.size(), s.c_str() );
        while (lo <= hi) {
            int mid = lo + (hi-lo)/2;
            if  (table[mid].compare(0, s.size(), s) < 0) { lo = mid + 1; } 
            else                                         { hi = mid - 1; }
        }
        if( lo>=table.size() )lo=-1;
        return lo;
    }

    int findMatch( const std::string& s ){
        int i = findPrefixStart( s );
        if( i<0 ) return -1;
        //printf( "findMatch->i(%i) table.size(%i) s(%s)\n", i, table.size(), s.c_str() );
        if( prefixMatches(s, table[i]) ){
            return i;
        }
        return -1;
    }

    // Function to find all strings that start with the given prefix
    int findMatchEnd( const std::string& s, int& istart ){
        int n = 0;
        for (int i = istart; i < table.size() && prefixMatches(s, table[i]); ++i) {
            n++;
        }
        return n;
    }

};

#endif









