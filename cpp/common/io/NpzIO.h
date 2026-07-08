/// @file NpzIO.h
/// @brief Read NumPy .npz (zip deflate, method 8) into named NpyArray entries.
/// Loads all decodable arrays; warns and skips unsupported dtypes (pipeline unicode metadata).
/// Parity reference: web/common_js/npzIO.js readZipEntries / readNpzFile
#ifndef NpzIO_h
#define NpzIO_h

#include "NpyIO.h"
#include <map>
#include <string>
#include <zlib.h>

struct NpzFile { public:
    std::vector<uint8_t> file_bytes;
    std::map<std::string, NpyArray> arrays;
    const NpyArray* get(const char* key) const {
        auto it = arrays.find(key);
        return (it==arrays.end()) ? nullptr : &it->second;
    }
    NpyArray& require(const char* key){
        auto it = arrays.find(key);
        if(it==arrays.end()){ printf("NpzIO ERROR: missing key '%s'\n", key); exit(1); }
        return it->second;
    }
};

inline std::vector<uint8_t> npz_inflate_raw(const uint8_t* comp, size_t comp_sz, size_t raw_sz){
    std::vector<uint8_t> raw(raw_sz);
    z_stream strm = {};
    if(inflateInit2(&strm, -MAX_WBITS) != Z_OK){ printf("NpzIO ERROR: inflateInit2\n"); exit(1); }
    strm.next_in = (Bytef*)comp; strm.avail_in = (uInt)comp_sz;
    strm.next_out = raw.data(); strm.avail_out = (uInt)raw_sz;
    int ret = inflate(&strm, Z_FINISH);
    inflateEnd(&strm);
    if(ret != Z_STREAM_END || strm.total_out != raw_sz){
        printf("NpzIO ERROR: inflate failed ret=%d out=%lu expected=%zu\n", ret, (unsigned long)strm.total_out, raw_sz);
        exit(1);
    }
    return raw;
}

inline NpzFile npz_read_file(const char* path){
    std::vector<uint8_t> buf = npy_read_file_bytes(path);
    NpzFile out; out.file_bytes = std::move(buf);
    const uint8_t* b = out.file_bytes.data();
    size_t len = out.file_bytes.size();
    size_t off = 0;
    while(off + 30 <= len){
        uint32_t sig = b[off] | (b[off+1]<<8) | (b[off+2]<<16) | (b[off+3]<<24);
        if(sig==0x06054b50 || sig==0x02014b50) break;
        if(sig != 0x04034b50){ printf("NpzIO ERROR: bad zip sig 0x%08x at %zu\n", sig, off); exit(1); }
        uint16_t compM = b[off+8] | (b[off+9]<<8);
        uint32_t compSz = b[off+18] | (b[off+19]<<8) | (b[off+20]<<16) | (b[off+21]<<24);
        uint32_t rawSz  = b[off+22] | (b[off+23]<<8) | (b[off+24]<<16) | (b[off+25]<<24);
        uint16_t nameLen = b[off+26] | (b[off+27]<<8);
        uint16_t extraLen = b[off+28] | (b[off+29]<<8);
        std::string name((const char*)(b+off+30), nameLen);
        size_t dataStart = off + 30 + nameLen + extraLen;
        if(dataStart + compSz > len){ printf("NpzIO ERROR: entry '%s' truncated\n", name.c_str()); exit(1); }
        const uint8_t* comp = b + dataStart;
        std::vector<uint8_t> raw;
        if(compM == 0) raw.assign(comp, comp+compSz);
        else if(compM == 8) raw = npz_inflate_raw(comp, compSz, rawSz);
        else { printf("NpzIO ERROR: unsupported compression %u in %s\n", compM, name.c_str()); exit(1); }
        if(raw.size() != rawSz){ printf("NpzIO ERROR: size mismatch %s\n", name.c_str()); exit(1); }
        std::string key = name;
        if(key.size() > 4 && key.substr(key.size()-4)==".npy") key = key.substr(0, key.size()-4);
        NpyArray arr;
        if(npy_try_parse_buffer(raw.data(), raw.size(), true, arr)){
            out.arrays[key] = std::move(arr);
        }else{
            printf("NpzIO WARNING: skip unsupported array '%s' in '%s'\n", key.c_str(), path);
        }
        off = dataStart + compSz;
    }
    return out;
}

#endif
