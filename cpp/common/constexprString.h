#ifndef _CONSTEXPR_STRING_H_
#define _CONSTEXPR_STRING_H_

#include <cstddef>
#include <string>
#include <variant>

template<size_t N>
struct constexprString{
    char data[N] {};

    consteval constexprString(const char (&s)[N]){
        for (int i=0; i<N; i++) data[i] = s[i];
    }

    constexpr operator std::string() const {return std::string(data, N);}

    template<size_t M>
    consteval bool operator==(const constexprString<M>& other) const {
        if constexpr(N != M) return false;
        for (int i=0; i<N; i++) if (data[i] != other.data[i]) return false;
        return true;
    }

    consteval const char& operator[](size_t i) const { return data[i]; }
    consteval const char* c_str() const { return data; }
    consteval size_t size() const { return N; }
};

template<constexprString...list>
struct constexprStringList{
    static consteval size_t size() { return sizeof...(list); }
    using IdxSeq = std::make_index_sequence<size()>;

    template<constexprString name, size_t...Idx> static consteval bool _contains_impl(std::index_sequence<Idx...>){
        return ((name == list) || ...);
    }
    
    template<constexprString name, size_t...Idx> static consteval size_t _get_idx_impl(std::index_sequence<Idx...>){
        static_assert(contains<name>(), "ERROR: name not found. Check it is in this list with contains<name>() first.");
        return ((name == list ? Idx : 0)+...);
    }

    template<constexprString name> static consteval bool contains(){ return _contains_impl<name>(IdxSeq{}); }
    template<constexprString name> static consteval size_t get_idx(){ return _get_idx_impl<name>(IdxSeq{}); }
};

#endif // _CONSTEXPR_STRING_H_
