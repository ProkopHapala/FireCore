
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


    inline float3 drotate_omega_csa( float3 p, float3 w, float ca, float sa){
        float dx  =  w.y*w.z -  w.z*w.y;
        float dy  =  w.z*w.x -  w.x*w.z;
        float dz  =  w.x*w.y -  w.y*w.x;
        float ddx = dy*w.z - dz*w.y;
        float ddy = dz*w.x - dx*w.z;
        float ddz = dx*w.y - dy*w.x;
        p.x-=dx*sa + ddx*ca;
        p.y-=dy*sa + ddy*ca;
        p.z-=dz*sa + ddz*ca;
        return p;
    }

    
__kernel void simulate_sigma_pi_single(
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
        float Ebl = evalBond( h, rij, 1.f, &fb );
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
        float  ca  = cos(w*dt);
        float  sa  = sin(w*dt);
        drotate_omega_csa( hpi, hwi, ca, sa);
        hpi   = normalize(hpi);

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