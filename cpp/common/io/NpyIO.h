/// @file NpyIO.h
/// @brief Parse NumPy .npy headers and payloads (little-endian, v1 header).
/// Decodes numeric/bool dtypes used by nanocrystal NPZ; `npy_try_parse_buffer` returns false for
/// unicode metadata (`<U*`) so NpzIO can skip non-geometry keys without aborting.
/// Parity reference: web/common_js/npzIO.js parseNpy / encodeNpy
#ifndef NpyIO_h
#define NpyIO_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

struct NpyArray { public:
    std::vector<uint8_t> storage;
    void* mmap_base = nullptr;
    size_t mmap_len = 0;
    const uint8_t* data = nullptr;
    size_t nbytes = 0;
    char descr[16] = {};
    int ndim = 0;
    int shape[8] = {};
    NpyArray() = default;
    NpyArray(const NpyArray&) = delete;
    NpyArray& operator=(const NpyArray&) = delete;
    NpyArray(NpyArray&& o) noexcept { *this = std::move(o); }
    NpyArray& operator=(NpyArray&& o) noexcept {
        if(mmap_base) munmap(mmap_base, mmap_len);
        storage = std::move(o.storage); mmap_base = o.mmap_base; mmap_len = o.mmap_len;
        data = o.data; nbytes = o.nbytes; strncpy(descr, o.descr, sizeof(descr));
        ndim = o.ndim; memcpy(shape, o.shape, sizeof(shape));
        o.mmap_base = nullptr; o.mmap_len = 0; o.data = nullptr;
        return *this;
    }
    ~NpyArray(){ if(mmap_base) munmap(mmap_base, mmap_len); }
    size_t count() const {
        if(ndim==0) return nbytes>0?1:0;
        size_t n = 1; for(int i=0;i<ndim;i++) n *= (size_t)shape[i];
        return n;
    }
};

inline void npy_fail(const char* msg){ printf("NpyIO ERROR: %s\n", msg); exit(1); }

inline bool npy_elem_size_opt(const char* descr, int& esz){
    if(!strcmp(descr,"<f8")) { esz=8; return true; }
    if(!strcmp(descr,"<f4")) { esz=4; return true; }
    if(!strcmp(descr,"<i8")) { esz=8; return true; }
    if(!strcmp(descr,"<i4")) { esz=4; return true; }
    if(!strcmp(descr,"|u1")) { esz=1; return true; }
    if(!strcmp(descr,"|b1")) { esz=1; return true; }
    return false;
}

inline int npy_elem_size(const char* descr){
    int esz = 0;
    if(!npy_elem_size_opt(descr, esz)){ printf("NpyIO ERROR: unsupported descr '%s'\n", descr); exit(1); }
    return esz;
}

inline void npy_parse_header(const uint8_t* buf, size_t len, int& hoff, int& hlen, char* descr_out, int descr_sz, int& ndim, int* shape, int shape_max){
    if(len < 10 || memcmp(buf, "\x93NUMPY", 6)) npy_fail("not npy magic");
    int major = buf[6];
    hoff = (major==1) ? 10 : 12;
    hlen = (major==1) ? (int)(buf[8] | (buf[9]<<8)) : (int)(buf[8] | (buf[9]<<8) | (buf[10]<<16) | (buf[11]<<24));
    if((size_t)(hoff+hlen) > len) npy_fail("header truncated");
    std::string header((const char*)(buf+hoff), hlen);
    size_t p1 = header.find("'descr':");
    size_t p2 = header.find("'shape':");
    if(p1==std::string::npos || p2==std::string::npos){ printf("NpyIO ERROR: bad header: %.120s\n", header.c_str()); exit(1); }
    size_t q1 = header.find('\'', p1+8); size_t q2 = header.find('\'', q1+1);
    if(q1==std::string::npos || q2==std::string::npos) npy_fail("bad descr in header");
    strncpy(descr_out, header.substr(q1+1, q2-q1-1).c_str(), descr_sz-1);
    size_t lp = header.find('(', p2); size_t rp = header.find(')', lp);
    if(lp==std::string::npos || rp==std::string::npos) npy_fail("bad shape in header");
    std::string inside = header.substr(lp+1, rp-lp-1);
    ndim = 0;
    size_t i0 = 0;
    while(i0 < inside.size()){
        while(i0<inside.size() && (inside[i0]==' '||inside[i0]==',')) i0++;
        if(i0>=inside.size()) break;
        size_t i1 = i0; while(i1<inside.size() && inside[i1]!=',') i1++;
        std::string tok = inside.substr(i0, i1-i0);
        while(!tok.empty() && tok.back()==' ') tok.pop_back();
        if(!tok.empty() && ndim < shape_max) shape[ndim++] = atoi(tok.c_str());
        i0 = i1;
    }
}

/// Parse npy buffer; return false for dtypes we do not decode (e.g. unicode metadata).
inline bool npy_try_parse_buffer(const uint8_t* buf, size_t len, bool copy_payload, NpyArray& out){
    int hoff=0, hlen=0, ndim=0, shape[8]={};
    char descr[16]={};
    npy_parse_header(buf, len, hoff, hlen, descr, sizeof(descr), ndim, shape, 8);
    int esz = 0;
    if(!npy_elem_size_opt(descr, esz)) return false;
    size_t data_off = (size_t)(hoff+hlen);
    if(data_off > len) npy_fail("data offset past end");
    size_t nbytes = len - data_off;
    if(nbytes % (size_t)esz) npy_fail("payload size not multiple of elem size");
    out = NpyArray();
    strncpy(out.descr, descr, sizeof(out.descr));
    out.ndim = ndim; memcpy(out.shape, shape, sizeof(out.shape));
    out.nbytes = nbytes;
    if(copy_payload){
        out.storage.assign(buf+data_off, buf+len);
        out.data = out.storage.data();
    }else{
        out.data = buf + data_off;
    }
    return true;
}

inline NpyArray npy_parse_buffer(const uint8_t* buf, size_t len, bool copy_payload=true){
    int hoff=0, hlen=0, ndim=0, shape[8]={};
    char descr[8]={};
    npy_parse_header(buf, len, hoff, hlen, descr, sizeof(descr), ndim, shape, 8);
    int esz = npy_elem_size(descr);
    size_t data_off = (size_t)(hoff+hlen);
    if(data_off > len) npy_fail("data offset past end");
    size_t nbytes = len - data_off;
    if(nbytes % (size_t)esz) npy_fail("payload size not multiple of elem size");
    NpyArray out;
    strncpy(out.descr, descr, sizeof(out.descr));
    out.ndim = ndim; memcpy(out.shape, shape, sizeof(out.shape));
    out.nbytes = nbytes;
    if(copy_payload){
        out.storage.assign(buf+data_off, buf+len);
        out.data = out.storage.data();
    }else{
        out.data = buf + data_off;
    }
    return out;
}

inline std::vector<uint8_t> npy_read_file_bytes(const char* path){
    FILE* f = fopen(path, "rb");
    if(!f){ printf("NpyIO ERROR: cannot open '%s'\n", path); exit(1); }
    fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
    if(sz < 0){ printf("NpyIO ERROR: ftell failed '%s'\n", path); exit(1); }
    std::vector<uint8_t> buf((size_t)sz);
    if(fread(buf.data(), 1, (size_t)sz, f) != (size_t)sz){ printf("NpyIO ERROR: fread '%s'\n", path); exit(1); }
    fclose(f);
    return buf;
}

inline NpyArray npy_load_file(const char* path, bool mmap_body=false){
    struct stat st;
    if(stat(path, &st) != 0){ printf("NpyIO ERROR: stat '%s'\n", path); exit(1); }
    if(S_ISDIR(st.st_mode)){ printf("NpyIO ERROR: '%s' is a directory, not .npy\n", path); exit(1); }
    if(!S_ISREG(st.st_mode)){ printf("NpyIO ERROR: '%s' is not a regular file\n", path); exit(1); }
    if(!mmap_body){
        std::vector<uint8_t> bytes = npy_read_file_bytes(path);
        return npy_parse_buffer(bytes.data(), bytes.size(), true);
    }
    int fd = open(path, O_RDONLY);
    if(fd < 0){ printf("NpyIO ERROR: open '%s'\n", path); exit(1); }
    size_t flen = (size_t)st.st_size;
    void* base = mmap(nullptr, flen, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if(base == MAP_FAILED){ printf("NpyIO ERROR: mmap '%s'\n", path); exit(1); }
    const uint8_t* buf = (const uint8_t*)base;
    int hoff=0, hlen=0, ndim=0, shape[8]={}; char descr[8]={};
    npy_parse_header(buf, flen, hoff, hlen, descr, sizeof(descr), ndim, shape, 8);
    size_t data_off = (size_t)(hoff+hlen);
    NpyArray out;
    strncpy(out.descr, descr, sizeof(out.descr));
    out.ndim = ndim; memcpy(out.shape, shape, sizeof(out.shape));
    out.nbytes = flen - data_off;
    out.mmap_base = base; out.mmap_len = flen;
    out.data = buf + data_off;
    return out;
}

inline const int32_t*  npy_as_i4(const NpyArray& a){ if(strcmp(a.descr,"<i4")) npy_fail("expected <i4"); return (const int32_t*)a.data; }
inline const int64_t*  npy_as_i8(const NpyArray& a){ if(strcmp(a.descr,"<i8")) npy_fail("expected <i8"); return (const int64_t*)a.data; }
inline const float*    npy_as_f4(const NpyArray& a){ if(strcmp(a.descr,"<f4")) npy_fail("expected <f4"); return (const float*)a.data; }
inline const double*   npy_as_f8(const NpyArray& a){ if(strcmp(a.descr,"<f8")) npy_fail("expected <f8"); return (const double*)a.data; }
inline const uint8_t*  npy_as_u1(const NpyArray& a){ if(strcmp(a.descr,"|u1")) npy_fail("expected |u1"); return (const uint8_t*)a.data; }
inline const uint8_t*  npy_as_b1(const NpyArray& a){ if(strcmp(a.descr,"|b1")) npy_fail("expected |b1"); return (const uint8_t*)a.data; }

#endif
