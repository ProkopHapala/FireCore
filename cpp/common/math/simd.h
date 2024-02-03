
#ifndef  simd_h
#define  simd_h

#include <immintrin.h>

#include "datatypes.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

using Vec2sd = Vec2T <__m256d>;
using Vec3sd = Vec3T <__m256d>;
using Vec4sd = Quat4T<__m256d>;

void pack_Vec3d_simd( int n, const Vec3d* src, Vec3sd* dst ){
    int m = n/4;
    for(int i=0; i<m; i++){
        const int i4 = i*4;
        Vec3d v1 = src[i4  ];
        Vec3d v2 = src[i4+1];
        Vec3d v3 = src[i4+2];
        Vec3d v4 = src[i4+3];
        Vec3sd v{ 
            _mm256_set_pd( v1.x, v2.x, v3.x, v4.x ), 
            _mm256_set_pd( v1.y, v2.y, v3.y, v4.y ), 
            _mm256_set_pd( v1.z, v2.z, v3.z, v4.z ) 
        };
        dst[i] = v; 
    }
};

void pack_Vec4d_simd( int n, const Quat4d* src, Vec4sd* dst ){
    int m = n/4;
    for(int i=0; i<m; i++){
        const int i4 = i*4;
        Quat4d v1 = src[i4  ];
        Quat4d v2 = src[i4+1];
        Quat4d v3 = src[i4+2];
        Quat4d v4 = src[i4+3];
        Vec4sd v{ 
            _mm256_set_pd( v1.x, v2.x, v3.x, v4.x ), 
            _mm256_set_pd( v1.y, v2.y, v3.y, v4.y ), 
            _mm256_set_pd( v1.z, v2.z, v3.z, v4.z ),
            _mm256_set_pd( v1.w, v2.w, v3.w, v4.w ) 
        };
        dst[i] = v; 
    }
};

// https://stackoverflow.com/questions/49941645/get-sum-of-values-stored-in-m256d-with-sse-avx/49943540#49943540
inline double hsum_double_avx(__m256d v) {
    __m128d vlow   = _mm256_castpd256_pd128(v);
    __m128d vhigh  = _mm256_extractf128_pd(v, 1); // high 128
            vlow   = _mm_add_pd(vlow, vhigh);     // reduce down to 128
    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}


#endif



