
// ===============================================================================================
// =========================  MINIMAL SIGMA-PI ALIGNMENT TEST KERNEL  ============================
// ===============================================================================================

inline float evalBond( float4 h, float dl, float k, __private float3* f ){
    float fr = dl*k;   // force magnitude
    *f = h.xyz * fr;   // force on atom a
    return fr*dl*0.5f;  // energy
}

// evaluate angular force and energy using cos(angle) formulation,    - faster, but not good for angles > 90 deg
inline float evalAngCos( const float4 hr1, const float4 hr2, float K, float c0, __private float3* f1, __private float3* f2 ){
    float  c = dot(hr1.xyz,hr2.xyz);
    float3 hf1,hf2;
    hf1 = hr2.xyz - hr1.xyz*c;
    hf2 = hr1.xyz - hr2.xyz*c;
    float c_   = c-c0;
    float E    = K*c_*c_;
    float fang = -K*c_*2;
    hf1 *= fang*hr1.w;
    hf2 *= fang*hr2.w;
    *f1=hf1;
    *f2=hf2;
    return E;
}

float2 sincosR2_taylor( float r2, const float2 cs ){
    const float c2 = -1.0/2;
    const float c3 = -1.0/6;
    const float c4 =  1.0/24;
    const float c5 =  1.0/120;
    const float c6 = -1.0/720;
    // r2  = w.x*w.x + w.y*w.y + w.z*w.z;
    return (float2){
      1  + r2*( c3 + c5*r2 ), // cos
      c2 + r2*( c4 + c6*r2 )  // sin
    };
}



float3 drotate_omega_csa( float3 p, float3 w, float2 cs ){
    float dx  =  w.y*w.z -  w.z*w.y;
    float dy  =  w.z*w.x -  w.x*w.z;
    float dz  =  w.x*w.y -  w.y*w.x;
    float ddx = dy*w.z - dz*w.y;
    float ddy = dz*w.x - dx*w.z;
    float ddz = dx*w.y - dy*w.x;
    p.x-=dx*cs.x + ddx*cs.y;
    p.y-=dy*cs.x + ddy*cs.y;
    p.z-=dz*cs.x + ddz*cs.y;
    return p;
}

float3 drotate_omega_taylor( float3 p, float3 w){
    float2 cs;
    sincosR2_taylor  ( dot(w,w), cs );
    return drotate_omega_csa( p,w,cs);
}

//================

// Pure Taylor: compute s = sin(r)/r and c = (1-cos(r))/r^2 as polynomials in r2
float2 sinc_div_r2_taylor(float r2){
    // series up to r^6 terms (i.e. up to r2^3)
    // s = sin(r)/r = 1 - r^2/6 + r^4/120 - r^6/5040
    // c = (1-cos r)/r^2 = 1/2 - r^2/24 + r^4/720 - r^6/40320
    const float r4 = r2*r2;
    const float r6 = r4*r2;
    float s = 1.0f + r2 * (-1.0f/6.0f) + r4 * (1.0f/120.0f) + r6 * (-1.0f/5040.0f);
    float c = 0.5f + r2 * (-1.0f/24.0f) + r4 * (1.0f/720.0f) + r6 * (-1.0f/40320.0f);
    return (float2){s, c};
}

// Pure-Taylor rotate (no sqrt). Good for small angles (e.g. |theta| < ~0.5...1.0)
float3 rotate_by_omega_taylor(float3 p, float3 w){
    float r2    = dot(w,w);
    float2 sc   = sinc_div_r2_taylor(r2); // sc.x = s, sc.y = c
    float3 wxp  = cross(w, p);
    float3 wwxp = cross(w, wxp);
    // p' = p + s*(w x p) + c*(w x (w x p))
    return (float3){ p.x + wxp.x*sc.x + wwxp.x*sc.y,
                     p.y + wxp.y*sc.x + wwxp.y*sc.y,
                     p.z + wxp.z*sc.x + wwxp.z*sc.y };
}

float3 rotate_by_omega(float3 p, float3 w){
    float r2    = dot(w,w);
    float  r    = sqrt(r2); 
    float2 sc   = (float2){cos(r), sin(r)}; // sc.x = s, sc.y = c
    float3 wxp  = cross(w, p);
    float3 wwxp = cross(w, wxp);
    // p' = p + s*(w x p) + c*(w x (w x p))
    return (float3){ p.x + wxp.x*sc.x + wwxp.x*sc.y,
                     p.y + wxp.y*sc.x + wwxp.y*sc.y,
                     p.z + wxp.z*sc.x + wwxp.z*sc.y };
}

// ====================


__kernel void simulate_sigma_pi(
    __global float4* pos,       // output: [nSteps,3] {node,cap,pi-rot}
    __global float4* vel,       // output: [nSteps,3] {node,cap,omega-pi}
    __global float4* force,     // output: [nSteps,3] {node,cap,torque-pi}
    const float       dt,       // time step
    const float       cdamp,    // linear damping factor
    const int         nSteps    // number of integration steps
){
    if(get_global_id(0) != 0) return;   // single-system kernel

    float3 pi   = pos [0].xyz;  // node position
    float3 pj   = pos [1].xyz;  // capping position
    float3 hpi  = pos [2].xyz;  // pi orientation

    float3 vi   = vel [0].xyz;
    float3 vj   = vel [1].xyz;
    float3 vpi  = vel [2].xyz;

    float bL0       = 1.5f; 
    float kBL       = 10.0f;
    float ksp       = 1.0f;
    float iM_i      = 1.0f;
    float iM_j      = 1.0f;
    float iI        = 1.0f;

    const float3 fallback_axis = (float3)(0.0f, 0.0f, 1.0f);

    for(int istep=0; istep<nSteps; istep++){

        const float3 dij    = pj - pi;
        const float  rij    = length(dij);
        const float  invl   = 1.f/rij;
        const float4 h      = (float4){dij*invl,invl};

        float3 fb;
        float Ebl = evalBond( h, rij-bL0, kBL, &fb );
        float3 fj = -fb;
        float3 fi =  fb; 

        float3 f1,f2;        
        float  Esp = evalAngCos( (float4){hpi,1.f}, h, ksp, 0.0f, &f1, &f2 );   //   pi-planarization (orthogonality), fpi is force on pi-orbital, fbs[i] is recoil force on i-th neighbor  
        fj        += f2;
        fi        -= f2; 
        float3 fpi = f1;

        fpi -= dot(fpi,hpi)*hpi; // remove component of force parallel to pi-vector

        // --- Integrate linear motion (semi-implicit Euler / damped)
        vi  += (fi  * iM_i ) * dt;  vi  *= cdamp;
        vj  += (fj  * iM_j ) * dt;  vj  *= cdamp;
        vpi += (fpi * iI   ) * dt;  vpi *= cdamp;

        vpi -= dot(vpi,hpi)*hpi; // remove component of velocity parallel to pi-vector
        pi  += vi  * dt;
        pj  += vj  * dt;
        hpi += vpi * dt;
        hpi  = normalize(hpi);

        // --- Store history for debugging (positions, orientation, forces)
        const int base  = istep*3;
        pos  [base+0] = (float4){ pi.x,  pi.y,  pi.z,  Ebl   };
        pos  [base+1] = (float4){ pj.x,  pj.y,  pj.z,  0.0f  };
        pos  [base+2] = (float4){ hpi.x, hpi.y, hpi.z, Esp   };

        vel  [base+0] = (float4){ vi.x,  vi.y,  vi.z,  0.0f  };
        vel  [base+1] = (float4){ vj.x,  vj.y,  vj.z,  0.0f  };
        vel  [base+2] = (float4){ vpi.x, vpi.y, vpi.z, 0.0f  };

        force[base+0] = (float4){ fi.x,  fi.y,  fi.z,  0.0f  };
        force[base+1] = (float4){ fj.x,  fj.y,  fj.z,  0.0f  };
        force[base+2] = (float4){ fpi.x, fpi.y, fpi.z, 0.0f  };
    }
}
    
__kernel void simulate_sigma_pi_rot(
    __global float4* pos,       // output: [nSteps,3] {node,cap,pi-rot}
    __global float4* vel,       // output: [nSteps,3] {node,cap,omega-pi}
    __global float4* force,     // output: [nSteps,3] {node,cap,torque-pi}
    const float       dt,       // time step
    const float       cdamp,    // linear damping factor
    const int         nSteps    // number of integration steps
){
    if(get_global_id(0) != 0) return;   // single-system kernel

    float3 pi   = pos [0].xyz;  // node position
    float3 pj   = pos [1].xyz;  // capping position
    float3 hpi  = pos [2].xyz;  // pi orientation

    float3 vi   = vel [0].xyz;
    float3 vj   = vel [1].xyz;
    float3 wi   = vel [2].xyz;

    float bL0       = 1.5f; 
    float kBL       = 10.0f;
    float ksp       = 1.0f;
    float iM_i      = 1.0f;
    float iM_j      = 1.0f;
    float iI        = 1.0f;

    const float3 fallback_axis = (float3)(0.0f, 0.0f, 1.0f);

    for(int istep=0; istep<nSteps; istep++){

        const float3 dij    = pj - pi;
        const float  rij    = length(dij);
        const float  invl   = 1.f/rij;
        const float4 h      = (float4){dij*invl,invl};

        float3 fb;
        float Ebl = evalBond( h, rij-bL0, kBL, &fb );
        float3 fj = -fb;
        float3 fi =  fb; 

        float3 f1,f2;        
        float  Esp = evalAngCos( (float4){hpi,1.f}, h, ksp, 0.0f, &f1, &f2 );   //   pi-planarization (orthogonality), fpi is force on pi-orbital, fbs[i] is recoil force on i-th neighbor  
        fj        += f2;
        fi        -= f2; 
        float3 tq  = cross(hpi, f1);

        // --- Integrate linear motion (semi-implicit Euler / damped)
        vi += (fi * iM_i ) * dt;  vi *= cdamp;
        vj += (fj * iM_j ) * dt;  vj *= cdamp;
        pi += vi * dt;
        pj += vj * dt;

        // --- Integrate rotational motion of pi vector
        wi   += (tq * iI) * dt;   wi  *= cdamp;
        //hpi  += cross(wi, hpi) * dt;

        float  w   = length(wi);
        float3 hwi = wi/w;
        w*=dt;
        drotate_omega_csa( hpi, hwi, (float2){cos(w), sin(w)});

        //hpi = rotate_by_omega(hpi,wi*dt);
        //hpi = drotate_omega_taylor( hpi, wi*dt );
        //hpi = rotate_by_omega_taylor( hpi, wi*dt );

        hpi = normalize(hpi);

        // --- Store history for debugging (positions, orientation, forces)
        const int base  = istep*3;
        pos  [base+0] = (float4){ pi.x, pi.y, pi.z, Ebl };
        pos  [base+1] = (float4){ pj.x, pj.y, pj.z, 0.0f };
        pos  [base+2] = (float4){ hpi.x, hpi.y, hpi.z, Esp   };

        vel  [base+0] = (float4){ vi.x, vi.y, vi.z, 0.0f   };
        vel  [base+1] = (float4){ vj.x, vj.y, vj.z, 0.0f   };
        vel  [base+2] = (float4){ wi.x, wi.y, wi.z, 0.0f   };

        force[base+0] = (float4){ fi.x, fi.y, fi.z, 0.0f   };
        force[base+1] = (float4){ fj.x, fj.y, fj.z, 0.0f   };
        force[base+2] = (float4){ tq.x, tq.y, tq.z, 0.0f   };
    }
}