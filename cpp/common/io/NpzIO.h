/// @file NpzIO.h
/// @brief Read NumPy .npz (zip stored or raw deflate) into named NpyArray entries.
/// Loads all decodable arrays; warns and skips unsupported dtypes (pipeline unicode metadata).
/// ZIP64: NumPy `np.savez` always opens entries with `force_zip64=True` (npyio.py, gh-10776),
/// so local 32-bit sizes are 0xFFFFFFFF and real sizes live in extra id 0x0001 (APPNOTE 4.5.3).
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

inline uint16_t npz_u16(const uint8_t* p){ return (uint16_t)(p[0] | (p[1]<<8)); }
inline uint32_t npz_u32(const uint8_t* p){ return (uint32_t)(p[0] | (p[1]<<8) | (p[2]<<16) | (p[3]<<24)); }
inline uint64_t npz_u64(const uint8_t* p){ return (uint64_t)npz_u32(p) | ((uint64_t)npz_u32(p+4)<<32); }

/// ZIP64 extra 0x0001: 8-byte fields appear in order only when the matching 32-bit local field is 0xFFFFFFFF (uncomp, then comp).
inline bool npz_zip64_sizes(const uint8_t* extra, size_t extraLen, uint32_t compSz32, uint32_t rawSz32, uint64_t& compSz, uint64_t& rawSz, const char* name){
    const bool needRaw  = (rawSz32  == 0xFFFFFFFFu);
    const bool needComp = (compSz32 == 0xFFFFFFFFu);
    if(!needRaw && !needComp){ compSz = compSz32; rawSz = rawSz32; return true; }
    size_t eoff = 0;
    while(eoff + 4 <= extraLen){
        uint16_t eid  = npz_u16(extra + eoff);
        uint16_t elen = npz_u16(extra + eoff + 2);
        if(eoff + 4 + elen > extraLen){ printf("NpzIO ERROR: extra field truncated in '%s' (eoff=%zu elen=%u extraLen=%zu)\n", name, eoff, elen, extraLen); return false; }
        if(eid == 0x0001){
            const uint8_t* p = extra + eoff + 4;
            size_t po = 0;
            rawSz = rawSz32; compSz = compSz32;
            if(needRaw){
                if(po + 8 > elen){ printf("NpzIO ERROR: ZIP64 extra missing uncompressed size in '%s'\n", name); return false; }
                rawSz = npz_u64(p + po); po += 8;
            }
            if(needComp){
                if(po + 8 > elen){ printf("NpzIO ERROR: ZIP64 extra missing compressed size in '%s'\n", name); return false; }
                compSz = npz_u64(p + po); po += 8;
            }
            return true;
        }
        eoff += 4 + elen;
    }
    printf("NpzIO ERROR: local sizes 0xFFFFFFFF but no ZIP64 extra 0x0001 in '%s' (extraLen=%zu)\n", name, extraLen);
    return false;
}

inline bool npz_inflate_raw(const uint8_t* comp, size_t comp_sz, size_t raw_sz, std::vector<uint8_t>& raw){
    raw.assign(raw_sz, 0);
    z_stream strm = {};
    if(inflateInit2(&strm, -MAX_WBITS) != Z_OK){ printf("NpzIO ERROR: inflateInit2\n"); return false; }
    strm.next_in = (Bytef*)comp; strm.avail_in = (uInt)comp_sz;
    strm.next_out = raw.data(); strm.avail_out = (uInt)raw_sz;
    int ret = inflate(&strm, Z_FINISH);
    inflateEnd(&strm);
    if(ret != Z_STREAM_END || strm.total_out != raw_sz){
        printf("NpzIO ERROR: inflate failed ret=%d out=%lu expected=%zu\n", ret, (unsigned long)strm.total_out, raw_sz);
        return false;
    }
    return true;
}

inline bool npz_try_read_file(const char* path, NpzFile& out){
    std::vector<uint8_t> buf = npy_read_file_bytes(path);
    out = NpzFile();
    out.file_bytes = std::move(buf);
    const uint8_t* b = out.file_bytes.data();
    size_t len = out.file_bytes.size();
    size_t off = 0;
    while(off + 30 <= len){
        uint32_t sig = npz_u32(b + off);
        if(sig==0x06054b50 || sig==0x02014b50) break; // EOCD / central directory
        if(sig != 0x04034b50){ printf("NpzIO ERROR: bad zip sig 0x%08x at %zu in '%s'\n", sig, off, path); return false; }
        uint16_t flags    = npz_u16(b + off + 6);
        uint16_t compM    = npz_u16(b + off + 8);
        uint32_t compSz32 = npz_u32(b + off + 18);
        uint32_t rawSz32  = npz_u32(b + off + 22);
        uint16_t nameLen  = npz_u16(b + off + 26);
        uint16_t extraLen = npz_u16(b + off + 28);
        if(off + 30 + nameLen + extraLen > len){
            printf("NpzIO ERROR: local header overflow at %zu in '%s' (nameLen=%u extraLen=%u fileLen=%zu)\n", off, path, nameLen, extraLen, len);
            return false;
        }
        std::string name((const char*)(b+off+30), nameLen);
        if(flags & 0x08){
            printf("NpzIO ERROR: entry '%s' in '%s' uses ZIP data-descriptor flag (0x08); local sizes unknown\n", name.c_str(), path);
            return false;
        }
        uint64_t compSz = 0, rawSz = 0;
        if(!npz_zip64_sizes(b + off + 30 + nameLen, extraLen, compSz32, rawSz32, compSz, rawSz, name.c_str())) return false;
        size_t dataStart = off + 30 + nameLen + extraLen;
        if(dataStart > len || compSz > (uint64_t)(len - dataStart)){
            printf("NpzIO ERROR: entry '%s' truncated in '%s' (dataStart=%zu compSz=%llu fileLen=%zu local32=%u zip64=%d)\n",
                name.c_str(), path, dataStart, (unsigned long long)compSz, len, compSz32, (int)(compSz32==0xFFFFFFFFu || rawSz32==0xFFFFFFFFu));
            return false;
        }
        const uint8_t* comp = b + dataStart;
        std::vector<uint8_t> raw;
        if(compM == 0) raw.assign(comp, comp + (size_t)compSz);
        else if(compM == 8){
            if(!npz_inflate_raw(comp, (size_t)compSz, (size_t)rawSz, raw)) return false;
        }else{
            printf("NpzIO ERROR: unsupported compression %u in '%s' of '%s'\n", compM, name.c_str(), path);
            return false;
        }
        if(raw.size() != (size_t)rawSz){ printf("NpzIO ERROR: size mismatch %s in '%s' (got %zu expected %llu)\n", name.c_str(), path, raw.size(), (unsigned long long)rawSz); return false; }
        std::string key = name;
        if(key.size() > 4 && key.substr(key.size()-4)==".npy") key = key.substr(0, key.size()-4);
        NpyArray arr;
        if(npy_try_parse_buffer(raw.data(), raw.size(), true, arr)){
            out.arrays[key] = std::move(arr);
        }else{
            printf("NpzIO WARNING: skip unsupported array '%s' in '%s'\n", key.c_str(), path);
        }
        off = dataStart + (size_t)compSz;
    }
    return true;
}

inline NpzFile npz_read_file(const char* path){
    NpzFile out;
    if(!npz_try_read_file(path, out)) exit(1);
    return out;
}

#endif
