
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
    int findPrefixStart( const std::string& prefix ){
        int low = 0, high = table.size() - 1;
        while (low <= high) {
            int mid = low + (high - low) / 2;
            if (table[mid].compare(0, prefix.size(), prefix) < 0) {
                low  = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        return low;
    }

    int findMatch( const std::string& prefix ){
        int i = findPrefixStart( prefix );
        if( prefixMatches(prefix, table[i]) ){
            return i;
        }
        return -1;
    }

    // Function to find all strings that start with the given prefix
    int findMatchEnd( const std::string& prefix, int& istart ){
        int n = 0;
        for (int i = istart; i < table.size() && prefixMatches(prefix, table[i]); ++i) {
            n++;
        }
        return n;
    }

};

#endif









