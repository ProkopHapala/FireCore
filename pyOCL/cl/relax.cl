
#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

// https://www.khronos.org/registry/OpenCL/sdk/1.1/docs/man/xhtml/sampler_t.html

// see https://gist.github.com/likr/3735779
// https://www.khronos.org/registry/OpenCL/sdk/1.1/docs/man/xhtml/sampler_t.html
__constant sampler_t sampler_1 =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;
//__constant sampler_t sampler_1 =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;
//__constant sampler_t sampler_2 =  CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_NEAREST; // NOTE AMD-GPU seems to not accept CLK_NORMALIZED_COORDS_FALSE
__constant sampler_t sampler_2 =  CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_NEAREST;
//__constant sampler_t sampler_1 = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_REPEAT | CLK_FILTER_LINEAR;

/*
float3 tipForce( float3 dpos, float4 stiffness, float4 dpos0 ){
    float r = sqrt( dot( dpos,dpos) );
    return  (dpos-dpos0.xyz) * stiffness.xyz              // harmonic 3D
         + dpos * ( stiffness.w * (r-dpos0.w)/r );  // radial
}
*/

//#pragma OPENCL EXTENSION cl_amd_printf : enable
//#pragma OPENCL EXTENSION cl_intel_printf : enable
//#pragma OPENCL EXTENSION cl_amd_printf : enable

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.399644f  // [eV*Ang/e^2]

inline float3 rotMat( float3 v, float3 a, float3 b, float3 c ){ return (float3)(dot(v,a),dot(v,b),dot(v,c)); }
inline float3 rotMatT( float3 v,  float3 a, float3 b, float3 c  ){ return a*v.x + b*v.y + c*v.z; }

void print  (__constant char* s,float4 v){ printf("%s(%g,%g,%g|%g) " ,s,v.x,v.y,v.z,v.w); };
void println(__constant char* s,float4 v){ printf("%s(%g,%g,%g|%g)\n",s,v.x,v.y,v.z,v.w); };


float3 tipForce( float3 dpos, float4 stiffness, float4 dpos0 ){
    float r = sqrt( dot( dpos,dpos) );
    return  (dpos-dpos0.xyz) * stiffness.xyz        // harmonic 3D
         + dpos * ( stiffness.w * (r-dpos0.w)/r );  // radial
}

float4 read_imagef_trilin( __read_only image3d_t imgIn, float4 coord ){
    float4 d = (float4)(0.00666666666,0.00666666666,0.00666666666,1.0); 
    float4 icoord;
    float4 fc     =  fract( coord/d, &icoord );
    icoord*=d;
    float4 mc     = (float4)(1.0,1.0,1.0,1.0) - fc;
    // NOTE AMD-GPU seems to not accept CLK_NORMALIZED_COORDS_FALSE
    //return read_imagef( imgIn, sampler_2, icoord );
    //return read_imagef( imgIn, sampler_1, coord );
    return  
     (( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,0.0,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,0.0,0.0,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,d.y,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,d.y,0.0,0.0) ) * fc.x )*fc.y )*mc.z
    +(( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,0.0,d.z,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,0.0,d.z,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,d.y,d.z,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,d.y,d.z,0.0) ) * fc.x )*fc.y )*fc.z;
}; 


float4 interpFE( float3 pos, float3 dinvA, float3 dinvB, float3 dinvC, __read_only image3d_t imgIn ){
    const float4 coord = (float4)( dot(pos,dinvA),dot(pos,dinvB),dot(pos,dinvC), 0.0f );
    return read_imagef( imgIn, sampler_1, coord );
    //return coord;
}

float4 interpFE_prec( float3 pos, float3 dinvA, float3 dinvB, float3 dinvC, __read_only image3d_t imgIn ){
    const float4 coord = (float4)( dot(pos,dinvA),dot(pos,dinvB),dot(pos,dinvC), 0.0f );
    return read_imagef_trilin( imgIn, coord ); 
    // read_imagef( imgIn, sampler_1, coord );
    //return coord;
}

// this should be macro, to pass values by reference
void move_LeapFrog( float3 f, float3 p, float3 v, float2 RP ){
    v  =  f * RP.x + v*RP.y;
    p +=  v * RP.x;
}


//#define N_RELAX_STEP_MAX  64
#define N_RELAX_STEP_MAX  128
#define F2CONV  1e-8f

#define OPT_FIRE 1
#if OPT_FIRE 
#define FTDEC 0.5f
#define FTINC 1.1f
#define FDAMP 0.99f


//#define F2CONV  1e-6f
#define F2SAFE    1e-8f

float3 update_FIRE( float3 f, float3 v, float* dt, float* damp,    float dtmin, float dtmax, float damp0 ){
    // Bitzek, E., Koskinen, P., Gähler, F., Moseler, M., & Gumbsch, P. (2006). Structural Relaxation Made Simple. Physical Review Letters, 97(17), 170201. 
    // https://doi.org/10.1103/PhysRevLett.97.170201
    // http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf
    float ff = dot(f,f);
    float vv = dot(v,v);
    float vf = dot(v,f);
    if( vf < 0 ){ // if velocity along direction of force
        v      *= 0;
        (*dt)   = fmax( dtmin, (*dt) * FTDEC );
        (*damp) = damp0;
    }else{       // if velocity against direction of force
        // v = cV * v  + cF * F
        v       *= (1 - (*damp));
        v       +=  f * ( (*damp) * sqrt( vv / (ff + F2SAFE ) ) );
        (*dt)    = fmin( dtmax, (*dt) * FTINC );
        (*damp) *= FDAMP;
    }
    return v;
    //v  += f * dt;
    //p  += v * dt;
}

#endif  //  OPT_FIRE 1

__kernel void getFEinPoints(
    __read_only image3d_t  imgIn,
    __global  float4*      points,
    __global  float4*      FEs,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC
){
    //const float4 coord     = points[get_global_id(0)];
    //vals[get_global_id(0)] = read_imagef(imgIn, sampler_1, coord);
    FEs[get_global_id(0)]    = interpFE( points[get_global_id(0)].xyz, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
}

__kernel void getFEinPointsShifted(
    __read_only image3d_t  imgIn,
    __global  float4*      points,
    __global  float4*      FEs,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC,
    float4 dpos0
){
    FEs[get_global_id(0)] = interpFE( points[get_global_id(0)].xyz+dpos0.xyz, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
}

__kernel void getFEinStrokes(
    __read_only image3d_t  imgIn,    // 1
    __global  float4*      points,   // 2
    __global  float4*      FEs,      // 3
    float4 dinvA,                    // 4
    float4 dinvB,                    // 5
    float4 dinvC,                    // 6
    float4 dTip,                     // 7
    float4 dpos0,                    // 8
    int nz                           // 
){
    //if(get_global_id(0)==0){ printf( "GPU getFEinStrokes() nz %i dTip(%g,%g,%g) dpos0(%g,%g,%g)\n", nz, dTip.x,dTip.y,dTip.z,   dpos0.x,dpos0.y,dpos0.z ); }
    //if(get_global_id(0)==0){ printf( "GPU getFEinStrokes() dinvA(%g,%g,%g) dinvB(%g,%g,%g) dinvC(%g,%g,%g)\n", dinvA.x,dinvA.y,dinvA.z,  dinvB.x,dinvB.y,dinvB.z,  dinvC.x,dinvC.y,dinvC.z ); }
    float3 pos    =  points[get_global_id(0)].xyz + dpos0.xyz; 
    for(int iz=0; iz<nz; iz++){
        float4 fe  =  read_imagef( imgIn, sampler_1, (float4){pos.x,pos.y,pos.z,0} );
        //float4 fe  = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
        if(get_global_id(0)==100)printf( "GPU %li %i (%f,%f,%f) -> fe(%g,%g,%g,%g) \n", get_global_id(0), iz, pos.x, pos.y, pos.z, fe.x,fe.y,fe.z,fe.w );
        //if(get_global_id(0)==0)printf( "GPU iz %i (%f,%f,%f) -> fe(%g,%g,%g,%g) \n", iz, pos.x, pos.y, pos.z, fe.x,fe.y,fe.z,fe.w );
        FEs[get_global_id(0)*nz + iz] = fe;
        pos    += dTip.xyz;
    }
}

__kernel void getFEinStrokesTilted(
    __read_only image3d_t  imgIn,
    __global  float4*      points,
    __global  float4*      FEs,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC,
    float4 tipA,
    float4 tipB,
    float4 tipC,
    float4 dTip,
    float4 dpos0,
    int nz
){
    float3 pos    =  points[get_global_id(0)].xyz + dpos0.xyz; 
    for(int iz=0; iz<nz; iz++){
        //printf( " %li %i (%f,%f,%f) \n", get_global_id(0), iz, pos.x, pos.y, pos.z );
        float4 fe   = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
        float4 fe_  = fe;
        fe_.xyz = rotMat( fe.xyz, tipA.xyz, tipB.xyz, tipC.xyz );
        FEs[get_global_id(0)*nz + iz]    = fe_;
        pos    += dTip.xyz;
    }
}

__kernel void getZisoTilted(
    __read_only image3d_t  imgIn,
    __global  float4*      points,
    __global  float*       zMap,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC,
    float4 tipA,
    float4 tipB,
    float4 tipC,
    float4 dTip,
    float4 dpos0,
    int nz, float iso
){
    float3 pos     = points[get_global_id(0)].xyz + dpos0.xyz; 
    float4 ofe,fe;
    ofe     = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
    ofe.xyz = rotMat( ofe.xyz, tipA.xyz, tipB.xyz, tipC.xyz );
    for(int iz=1; iz<nz; iz++){
        pos    += dTip.xyz;
        fe     = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
        fe.xyz = rotMat( fe.xyz, tipA.xyz, tipB.xyz, tipC.xyz );
        //if( get_global_id(0) == 6050 ) printf( "iz %i fe %g iso %g \n", iz, fe.z, iso );
        if( fe.z/iso > 1.0 ){
            float t = (iso - ofe.z)/(fe.z - ofe.z);
            zMap[get_global_id(0)] = iz + t;
            return;
        }
        ofe      = fe;
    }
    zMap[get_global_id(0)] = -1;
}

__kernel void getZisoFETilted(
    __read_only image3d_t  imgIn,
    __read_only image3d_t  imgFE,
    __global  float4*      points,
    __global  float*       zMap,
    __global  float4*      feMap,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC,
    float4 tipA,
    float4 tipB,
    float4 tipC,
    float4 dTip,
    float4 dpos0,
    int nz, float iso
){
    float3 pos     = points[get_global_id(0)].xyz + dpos0.xyz; 
    float4 ofe,fe;
    ofe     = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
    ofe.xyz = rotMat( ofe.xyz, tipA.xyz, tipB.xyz, tipC.xyz );
    for(int iz=1; iz<nz; iz++){
        pos    += dTip.xyz;
        fe     = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
        fe.xyz = rotMat( fe.xyz, tipA.xyz, tipB.xyz, tipC.xyz );
        //if( get_global_id(0) == 6050 ) printf( "iz %i fe %g iso %g \n", iz, fe.z, iso );
        if( fe.z/iso > 1.0 ){
            float t = (iso - ofe.z)/(fe.z - ofe.z);
            zMap [get_global_id(0)] = iz + t;
            fe     = interpFE( pos+dTip.xyz*t, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgFE );
            fe.xyz = rotMat( fe.xyz, tipA.xyz, tipB.xyz, tipC.xyz );
            feMap[get_global_id(0)] = fe;
            return;
        }
        ofe      = fe;
    }
    zMap [get_global_id(0)] = -1;
    feMap[get_global_id(0)] = (float4)(0.0f,0.0f,0.0f,0.0f);
}

__kernel void relaxPoints(
    __read_only image3d_t  imgIn,
    __global  float4*      points,
    __global  float4*      FEs,
    float4 stiffness,
    float4 dpos0,
    float4 relax_params  // (dt,damp,tmin,tmax)
){

    float dt      = relax_params.x;
    float damp    = relax_params.y;

    float dtmax = dt;
    float dtmin = dtmax*0.1;
    float damp0 = damp;

    float3 tipPos = points[get_global_id(0)].xyz;
    float3 pos    = tipPos.xyz + dpos0.xyz; 
    float4 fe;
    float3 v    = 0.0f;
    for(int i=0; i<1000; i++){
        fe        = read_imagef( imgIn, sampler_1, (float4)(pos,0.0f) ); /// this would work only for unitary cell
        float3 f  = fe.xyz;
        f        += tipForce( pos-tipPos, stiffness, dpos0 );    

        #if OPT_FIRE
        v = update_FIRE( f, v, &dt, &damp, dtmin, dtmax, damp0 );
        #else
        v        *=    (1 - damp);
        #endif
        v        += f * dt;
        pos.xyz  += v * dt;

    }
    FEs[get_global_id(0)] = fe;
}

__kernel void relaxStrokes(
    __read_only image3d_t  imgIn,
    __global  float4*      points,
    __global  float4*      FEs,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC,
    float4 dTip,
    float4 stiffness,
    float4 dpos0,
    float4 relax_params,
    int nz
){
    float3 tipPos = points[get_global_id(0)].xyz;
    float3 pos    = tipPos.xyz + dpos0.xyz; 
    
    float dt      = relax_params.x;
    float damp    = relax_params.y;
    //printf( " %li (%f,%f,%f)  \n",  get_global_id(0), tipPos.x, tipPos.y, tipPos.z);
    
    float dtmax = dt;
    float dtmin = dtmax*0.1f;
    float damp0 = damp;

    for(int iz=0; iz<nz; iz++){
        float4 fe;
        float3 v   = 0.0f;
        for(int i=0; i<N_RELAX_STEP_MAX; i++){
            fe        = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
            float3 f  = fe.xyz;
            f        += tipForce( pos-tipPos, stiffness, dpos0 );
            
            #if OPT_FIRE
            v = update_FIRE( f, v, &dt, &damp, dtmin, dtmax, damp0 );
            #else
            v        *=    (1 - damp);
            #endif
            v        += f * dt;
            pos.xyz  += v * dt;

            if(dot(f,f)<F2CONV) break;
        }
        FEs[get_global_id(0)*nz + iz] = fe;
        //FEs[get_global_id(0)*nz + iz].xyz = pos;
        tipPos += dTip.xyz;
        pos    += dTip.xyz;
    }
}

__kernel void relaxStrokesTilted_debug(
    __read_only image3d_t  imgIn,
    __global  float4*      points,
    __global  float4*      FEs,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC,
    float4 tipA,
    float4 tipB,
    float4 tipC,
    float4 stiffness,
    float4 dpos0,
    float4 relax_params,
    float4 surfFF,
    int nz
){
    const float3 dTip   = tipC.xyz * tipC.w;
    const float4 dpos0_=dpos0; dpos0_.xyz= rotMatT( dpos0_.xyz , tipA.xyz, tipB.xyz, tipC.xyz );
    float3 tipPos = points[get_global_id(0)].xyz;
    float3 pos    = tipPos.xyz + dpos0_.xyz; 
    for(int iz=0; iz<nz; iz++){
        FEs[get_global_id(0)*nz + iz] = 1.0f;
    }
}


__kernel void relaxStrokesTilted(
    __read_only image3d_t  imgIn,   // 1
    __global  float4*      points,  // 2
    __global  float4*      FEs,     // 3
    float4 dinvA,                   // 4
    float4 dinvB,                   // 5
    float4 dinvC,                   // 6
    float4 tipA,                    // 7
    float4 tipB,                    // 8
    float4 tipC,                    // 9 
    float4 stiffness,               // 10
    float4 dpos0,                   // 11
    float4 relax_params,            // 12
    float4 surfFF,                  // 13
    int nz                          // 14
){

    const float3 dTip   = tipC.xyz * tipC.w;
    const float4 dpos0_=dpos0; dpos0_.xyz= rotMatT( dpos0_.xyz , tipA.xyz, tipB.xyz, tipC.xyz );

    float3 tipPos = points[get_global_id(0)].xyz;
    float3 pos    = tipPos.xyz + dpos0_.xyz; 

    float dt      = relax_params.x;
    float damp    = relax_params.y;

    float dtmax = dt;
    float dtmin = dtmax*0.1f;
    float damp0 = damp;

    // if( (get_global_id(0)==0) ){  
    //     printf( " dt %g damp %g \n", dt, damp );
    //     printf( " stiffness(%g,%g,%g|%g) dpos0(%g,%g,%g|%g) \n", stiffness.x,stiffness.y,stiffness.z,stiffness.w,  dpos0.x,dpos0.y,dpos0.z,dpos0.w  );
    //     printf( " relax_params(%g,%g,%g|%g) surfFF(%g,%g,%g|%g) \n", relax_params.x,relax_params.y,relax_params.z,relax_params.w,  surfFF.x,surfFF.y,surfFF.z,surfFF.w  );
    //     printf( " dinvA(%g,%g,%g|%g) tipA(%g,%g,%g|%g) \n", dinvA.x,dinvA.y,dinvA.z,dinvA.w,  tipA.x,tipA.y,tipA.z,tipA.w  );
    //     printf( " dinvB(%g,%g,%g|%g) tipB(%g,%g,%g|%g) \n", dinvB.x,dinvB.y,dinvB.z,dinvB.w,  tipB.x,tipB.y,tipB.z,tipB.w  );
    //     printf( " dinvc(%g,%g,%g|%g) tipC(%g,%g,%g|%g) \n", dinvC.x,dinvC.y,dinvC.z,dinvC.w,  tipC.x,tipC.y,tipC.z,tipC.w  );
    //     int i1=get_global_size(0)-1; printf( "pos0(%3.3f,%3.3f,%3.3f) pos1(%3.3f,%3.3f,%3.3f)\n", points[0].x,points[0].y,points[0].z, points[i1].x,points[i1].y,points[i1].z );
    //     //for(int i=0; i<get_global_size(0); i++ ){ printf( "pos[%i] (%3.3f,%3.3f,%3.3f)\n", i, points[i].x,points[i].y,points[i].z ); }
    // }
    //if( (get_global_id(0)==0) ){     float4 fe = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );  printf( " pos (%g,%g,%g) feImg(%g,%g,%g,%g) \n", pos.x, pos.y, pos.z, fe.x,fe.y,fe.z,fe.w );}

    for(int iz=0; iz<nz; iz++){
        float4 fe;
        float3 v   = (float3){0.f,0.f,0.f};
        
        for(int i=0; i<N_RELAX_STEP_MAX; i++){
            fe            = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
            float3 f      = fe.xyz;
            float3 dpos   = pos-tipPos;
            float3 dpos_  = rotMat  ( dpos, tipA.xyz, tipB.xyz, tipC.xyz );    // to tip-coordinates
            float3 ftip   = tipForce( dpos_, stiffness, dpos0 );

            f            += rotMatT ( ftip, tipA.xyz, tipB.xyz, tipC.xyz );      // from tip-coordinates
            f            +=  tipC.xyz * surfFF.x;                                // TODO: more sophisticated model of surface potential? Like Hamaker ?

            //f      +=  tipForce( dpos, stiffness, dpos0_ );  // Not rotated
            
            #if OPT_FIRE
            v = update_FIRE( f, v, &dt, &damp, dtmin, dtmax, damp0 );
            #else
            v        *=    (1 - damp);
            #endif
            v        += f * dt;
            pos.xyz  += v * dt;

            if(dot(f,f)<F2CONV) break;
        }
        fe            = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
        if(1){ // output tip-rotated force
            float4 fe_  = fe;
            fe_.xyz = rotMat( fe.xyz, tipA.xyz, tipB.xyz, tipC.xyz );
            fe_.w   = fe.w;
            FEs[get_global_id(0)*nz + iz] = fe_;
        }else{ // output molecule-rotated force 
            FEs[get_global_id(0)*nz + iz] = fe;
            //FEs[get_global_id(0)*nz + iz].xyz = pos;
        }
        tipPos += dTip.xyz;
        pos    += dTip.xyz;
        //if( (get_global_id(0)==0) ){ printf( "iz[%i] pos(%g,%g,%g) fe(%g,%g,%g|%g) \n", iz, pos.x, pos.y, pos.z, fe.x,fe.y,fe.z,fe.w ); }
        //if( (get_global_id(0)==0) ){ printf( "iz[%i] pos(%g,%g,%g) tipPos(%g,%g,%g) \n", iz, pos.x,pos.y,pos.z, tipPos.x,tipPos.y,tipPos.z ); }

    }
}



__kernel void relaxStrokesTilted_convZ(
    __read_only image3d_t  imgIn,
    __global  float4*      points,
    __constant  float*     weighs,
    __global  float4*      FEs,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC,
    float4 tipA,
    float4 tipB,
    float4 tipC,
    float4 stiffness,
    float4 dpos0,
    float4 relax_params,
    float4 surfFF,
    const int nz, const int nzout
){

    __local float  WEIGHTS[64];

    const float3 dTip   = tipC.xyz * tipC.w;
    const float4 dpos0_=dpos0; dpos0_.xyz= rotMatT( dpos0_.xyz , tipA.xyz, tipB.xyz, tipC.xyz );

    float3 tipPos = points[get_global_id(0)].xyz;
    float3 pos    = tipPos.xyz + dpos0_.xyz; 

    float dt      = relax_params.x;
    float damp    = relax_params.y;

    float dtmax = dt;
    float dtmin = dtmax*0.1f;
    float damp0 = damp;

    //if( (get_global_id(0)==0) ){     float4 fe = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );  printf( " pos (%g,%g,%g) feImg(%g,%g,%g,%g) \n", pos.x, pos.y, pos.z, fe.x,fe.y,fe.z,fe.w );}
    //if( (get_global_id(0)==0) ){ printf( "dt %g damp %g \n", dt, damp ); }; return;

    const int ioff = get_global_id(0)*nzout;
    const int nzw   = nz-nzout;
    const int iL=get_local_id(0);
    const int nL=get_local_size(0);
    for (int i=iL; i<nzw; i+=nL ){
        WEIGHTS[i] = weighs[i];
    }
    for (int iz=0; iz<nzout; iz++ ){
        FEs[ioff+iz] = 0.0f;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    int itr_tot = 0;

    for(int iz=0; iz<nz; iz++){
        float4 fe;
        float3 v   = 0.0f;
        for(int i=0; i<N_RELAX_STEP_MAX; i++){
        //for(int i=0; i<1; i++){ // DEBUG
            fe            = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
            //fe            = interpFE_prec( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
            float3 f      = fe.xyz;
            float3 dpos   = pos-tipPos;
            float3 dpos_  = rotMat  ( dpos, tipA.xyz, tipB.xyz, tipC.xyz );    // to tip-coordinates
            float3 ftip   = tipForce( dpos_, stiffness, dpos0 );

            f            += rotMatT ( ftip, tipA.xyz, tipB.xyz, tipC.xyz );      // from tip-coordinates
            f            +=  tipC.xyz * surfFF.x;                                // TODO: more sophisticated model of surface potential? Like Hamaker ?

            //f      +=  tipForce( dpos, stiffness, dpos0_ );  // Not rotated

            #if OPT_FIRE
            v = update_FIRE( f, v, &dt, &damp, dtmin, dtmax, damp0 );
            //if(get_global_id(0)==(64*128+64)){ printf( "itr,iz,i %i %i %i  |F| %g |v| %g <f,v> %g , (%g,%g,%g) (%g,%g,%g) damp %g dt %g \n", itr_tot, iz,i,  sqrt(dot(f,f)), sqrt(dot(v,v)),  dot(f,v),  fe.x,fe.y,fe.z, pos.x, pos.y, pos.z, damp, dt ); }
            #else
            v        *=    (1 - damp);
            //if(get_global_id(0)==(64*128+64)){ printf( "itr,iz,i %i %i %i  |F| %g |v| %g <f,v> %g , (%g,%g,%g) (%g,%g,%g) damp %g dt %g \n", itr_tot, iz,i,  sqrt(dot(f,f)), sqrt(dot(v,v)),  dot(f,v),  fe.x,fe.y,fe.z, pos.x, pos.y, pos.z, damp, dt ); }
            #endif
            v        += f * dt;
            pos.xyz  += v * dt;

            itr_tot++;
            if(dot(f,f)<F2CONV) break;
        }
        
        if(1){ // output tip-rotated force
            fe.xyz = rotMat( fe.xyz, tipA.xyz, tipB.xyz, tipC.xyz );
        }

        
        // do the convolution
        for(int izout=0;izout<nzout;izout++){
            int jzw = iz - izout;
            if((jzw<nzw)&&(jzw>0)){
                FEs[ ioff + izout] += fe * WEIGHTS[jzw];
            }
        }
        //if( iz<nzout ) FEs[ioff+iz] = fe;
        tipPos += dTip.xyz;
        pos    += dTip.xyz;
    }

}

__kernel void convolveZ(
    __global  float4* Fin,
    __global  float4* Fout,
    //__global  float*  weighs,
    __constant  float*  weighs,
    const int nzin, const int nzout
){
    const int ioffi = get_global_id(0)*nzin;
    const int ioffo = get_global_id(0)*nzout;
    const int nzw   = nzin-nzout;
    //if( get_global_id(0)==0 ) printf( "local size %i \n", get_local_size(0) );
    //if( get_global_id(0)==0 ) printf( "izo %i izi %i Fz %g W %g \n", nzin, nzout, nzw );

    __local float WEIGHTS[64];

    const int iL=get_local_id(0);
    const int nL=get_local_size(0);
    for (int i=iL; i<nzw; i+=nL ){
        if( i<nzw ) WEIGHTS[i] = weighs[i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    
    for(int izo=0; izo<nzout; izo++){
        float4 fe = 0.0f;
        for(int jz=0; jz<nzw; jz++){
            //fe += Fin[ ioffi + izo + jz ] * weighs[ jz ];
            fe += Fin[ ioffi + izo + jz ] * WEIGHTS[ jz ];
            //if( get_global_id(0)==0 ) printf( "izo %i izi %i Fz %g W %g \n", izo, jz, Fin[ ioffi + izo + jz ].z, weighs[ jz ] );
            //fe +=  tanh( Fin[ ioffi + izi ] ) * weighs[ izi - izo ];
        }
        //if( ioffi == 0 ){ printf( "izo %i w[i] %e \n", izo, weighs[ izo ] ); }
        //fe = (float)ioffo; // DEBUG
        Fout[ ioffo + izo ] = fe;
        //Fout[ ioffo + izo ] = weighs[ izo ];
        //Fout[ ioffo + izo ] = (float4) izo;
        //Fout[ ioffo + izo ] = Fin[ ioffi + izo ];
    }
}

__kernel void izoZ(
    __global  float4* Fin,
    __global  float*  zMap,
    int nz,   float iso
){
    int ioffi = get_global_id(0)*nz;
    float4 ofe = Fin[ ioffi ];
    for(int iz=1; iz<nz; iz++){
        float4 fe = Fin[ ioffi + iz ];
        // zMap[get_global_id(0)] = i;
        if( fe.z > iso ){
            float t = (iso - ofe.z)/(fe.z - ofe.z);
            zMap[get_global_id(0)] = iz + t;
            return;
        }
        ofe = fe;
    }
    zMap[get_global_id(0)] = -1;
}

// =========================================  
//           ForceField form FF.cl
// =========================================

float4 getCoulomb( float4 atom, float3 pos ){
     float3  dp  =  pos - atom.xyz;
     float   ir2 = 1.0f/( dot(dp,dp) +  R2SAFE );
     float   ir  = sqrt(ir2);
     float   E   = atom.w*sqrt(ir2);
     return (float4)(dp*(E*ir2), E );
}

float4 getLJ( float3 apos, float2 cLJ, float3 pos ){
     float3  dp  =  pos - apos;
     float   ir2 = 1.0f/( dot(dp,dp) + R2SAFE );
     float   ir6 = ir2*ir2*ir2;
     float   E   =  (    cLJ.y*ir6 -   cLJ.x )*ir6;
     float3  F   = (( 12.0f*cLJ.y*ir6 - 6.0f*cLJ.x )*ir6*ir2)*dp;
     return (float4)(F, E);
}

float4 getMorse( float3 dp, float3 REA ){
    //float3  dp  =  pos - apos;
    float   r     = sqrt( dot(dp,dp) + R2SAFE );
    float   expar = exp( REA.z*(r-REA.x) );
    float   E     = REA.y*expar*( expar - 2 );
    float   fr    = REA.y*expar*( expar - 1 )*2*REA.z;
    return (float4)(dp*(fr/r), E);
}

float8 getLJC( float4 atom, float2 cLJ, float3 pos ){
     float3  dp  =  pos - atom.xyz;
     float   ir2 = 1.0/( dot(dp,dp) +  R2SAFE );
     float   ir6 = ir2*ir2*ir2;
     float   ELJ =  (    cLJ.y*ir6 -   cLJ.x )*ir6;
     float3  FLJ = (( 12.0f*cLJ.y*ir6 - 6.0f*cLJ.x )*ir6*ir2)*dp;
     float   ir  = sqrt(ir2);
     float   Eel = atom.w*sqrt(ir2);
     return (float8)(FLJ, ELJ, dp*(Eel*ir2), Eel );
}

float getLorenz( float4 atom, float4 coefs, float3 pos ){
     float3  dp  =  pos - atom.xyz;
     return coefs.x/( dot(dp,dp) +  coefs.y*coefs.y );
     //return 1.0/( dot(dp,dp) +  0.000 );
}

__kernel void eval_GridFF(
    const int nAtoms, 
    __global float4*  atoms,
    __global float4*  REKQs,
    __global float4*  FE_Paul,
    __global float4*  FE_Lond,
    __global float4*  FE_Coul,
    int4 nGrid,
    float4 grid_p0,
    float4 grid_dA,
    float4 grid_dB,
    float4 grid_dC
){
    __local float4 LATOMS[32];
    __local float4 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;
    if(iG>nMax) return;
    float3 pos    = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;
    float4 fe_Paul = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    float4 fe_Lond = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    float4 fe_Coul = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = REKQs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                //fe += getLJC( LATOMS[j], LCLJS[j], pos );
                float4 REKQ = REKQs[j];
                float3 dp = pos - LATOMS[j].xyz;
                float  r2  = dot(dp,dp) +  R2SAFE;
                float ir2  = 1/r2; 
                float   r  = sqrt( r2 );
                // Electrostatic
                float   E  = REKQ.w*sqrt(ir2);
                fe_Coul   += (float4)(dp*(E*ir2), E );

                float   expar = exp( REKQ.z*(r-REKQ.x) );
                float   e     = REKQ.y*expar;
                float4  fe    = (float4)(dp*(e*REKQ.z), e );

                fe_Coul   += fe*2.0f;
                fe_Paul   += fe*expar;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE_Paul[iG] = fe_Paul;
    FE_Lond[iG] = fe_Lond;
    FE_Coul[iG] = fe_Coul;
}

__kernel void evalLJC_QZs(
    const int nAtoms,        // 1
    __global float4* atoms,  // 2
    __global float2*  cLJs,  // 3
    __global float4*    FE,  // 4
    int4 nGrid,              // 5
    float4 grid_p0,          // 6 
    float4 grid_dA,          // 7
    float4 grid_dB,          // 8
    float4 grid_dC,          // 9
    float4 Qs,               // 10
    float4 QZs               // 11
){
    __local float4 LATOMS[32];
    __local float2 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;

    //if (  get_global_id(0)==0 ) { printf("GPU evalLJC_QZs \n" ); }
    if(iG==0) printf( " Qs (%g,%g,%g,%g) QZs (%g,%g,%g,%g) \n", Qs.x,Qs.y,Qs.z,Qs.w,   QZs.x,QZs.y,QZs.z,QZs.w   );
    //if(iG==0) printf( " dA(%g,%g,%g) dB(%g,%g,%g) dC(%g,%g,%g) p0(%g,%g,%g)\n", grid_dA.x,grid_dA.y,grid_dA.z,   grid_dB.x,grid_dB.y,grid_dB.z,  grid_dC.x,grid_dC.y,grid_dC.z, grid_p0.x,grid_p0.y,grid_p0.z );
    if(iG>nMax) return;

    float3 pos    = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;

    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    
    Qs *= COULOMB_CONST;

    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs [i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                //fe += getLJC( LATOMS[j], LCLJS[j], pos );
                float4 xyzq = LATOMS[j];
                //if(iG==0) printf( "atom[%i](%g,%g,%g|%g) cLJ(%g,%g)\n", i, xyzq.x,xyzq.y,xyzq.z,  xyzq.w,   LCLJS[j].x, LCLJS[j].y );
                fe += getLJ     ( xyzq.xyz, LCLJS[j], pos );
                // ToDo : Electrostatics seems to be too strong in original forcefeidl
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.x) ) * Qs.x;
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.y) ) * Qs.y;
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.z) ) * Qs.z;
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.w) ) * Qs.w;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    float renorm = 100.0/fabs(fe.w);
    if( renorm<1.f ){ fe*=renorm; }
    //if ( (ia==nGrid.x/2)&&(ib==nGrid.y/2) ) { printf(" iz %i pos(%g,%g,%g) fe(%g,%g,%g|%g) \n", ic,  pos.x,pos.y,pos.z,  fe.x, fe.y, fe.z, fe.w ); }
    FE[iG] = fe;
}



__kernel void evalLJC_QZs_toImg(
    const int nAtoms,        // 1
    __global float4* atoms,  // 2
    __global float2*  cLJs,  // 3
    __write_only image3d_t  imgOut, // 4
    int4 nGrid,              // 5
    float4 grid_p0,          // 6 
    float4 grid_dA,          // 7
    float4 grid_dB,          // 8
    float4 grid_dC,          // 9
    float4 Qs,               // 10
    float4 QZs               // 11
){
    
    __local float4 LATOMS[32];
    __local float2 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;

    //if (  get_global_id(0)==0 ) { printf("GPU evalLJC_QZs \n" ); }
    // if(iG==0){
    //     printf( " nGrid(%i,%i,%i|%i)\n", nGrid.x,nGrid.y,nGrid.z,nGrid.w );
    //     printf( " grid_p0(%g,%g,%g|%g)\n", grid_p0.x,grid_p0.y,grid_p0.z,grid_p0.w );
    //     printf( " grid_dA(%g,%g,%g|%g)\n", grid_dA.x,grid_dA.y,grid_dA.z,grid_dA.w );
    //     printf( " grid_dB(%g,%g,%g|%g)\n", grid_dA.x,grid_dA.y,grid_dA.z,grid_dA.w );
    //     printf( " grid_dC(%g,%g,%g|%g)\n", grid_dB.x,grid_dB.y,grid_dB.z,grid_dB.w );
    //     printf( " dinvc(%g,%g,%g|%g)\n", grid_dC.x,grid_dC.y,grid_dC.z,grid_dC.w );
    //     printf( " Qs (%g,%g,%g|%g)\n", Qs.x,Qs.y,Qs.z,Qs.w );
    //     printf( " QZs(%g,%g,%g|%g)\n", QZs.x,QZs.y,QZs.z,QZs.w );
    //     for(int i=0; i<nAtoms; i++){
    //         printf( "atom(%g,%g,%g|%g) cLJ(%g,%g)\n", atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,  cLJs[i].x,cLJs[i].y );
    //     }
    // }
    if(iG>nMax) return;

    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    float3 pos    = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;
    
    Qs *= COULOMB_CONST;

    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs [i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                //fe += getLJC( LATOMS[j], LCLJS[j], pos );
                float4 xyzq = LATOMS[j];
                //if(iG==0) printf( "atom[%i](%g,%g,%g|%g) cLJ(%g,%g)\n", i, xyzq.x,xyzq.y,xyzq.z,  xyzq.w,   LCLJS[j].x, LCLJS[j].y );
                fe += getLJ     ( xyzq.xyz, LCLJS[j], pos );
                // ToDo : Electrostatics seems to be too strong in original forcefeidl
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.x) ) * Qs.x;
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.y) ) * Qs.y;
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.z) ) * Qs.z;
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.w) ) * Qs.w;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    float renorm = 100.0/fabs(fe.w);
    if( renorm<1.f ){ fe*=renorm; }
    //if ( (ia==nGrid.x/2)&&(ib==nGrid.y/2) ) { printf(" iz %i pos(%g,%g,%g) fe(%g,%g,%g|%g) \n", ic,  pos.x,pos.y,pos.z,  fe.x, fe.y, fe.z, fe.w ); }
    //imgOut[iG] = fe;
    //fe  = (float4){sin(ia*0.1), sin(ia*0.1), sin(ib*0.1), cos(ia*0.1)*cos(ib*0.1)*cos(ic*0.1) };
    write_imagef( imgOut, (int4){ia,ib,ic,0}, fe );
    //write_imagef( imgOut, (int4){0,0,0,0}, (float4){0.0f,0.0f,0.0f,0.0f} );
}