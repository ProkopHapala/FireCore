// ===============================================================================================
// ===========================  MINIMAL PI-PI ALIGNMENT TEST KERNEL  =============================
// ===============================================================================================

const float3 float3Zero = (float3)(0.0f, 0.0f, 0.0f);

inline float3 normalize_safe( const float3 v, const float3 fallback ){
    const float n2 = dot(v,v);
    if(n2 < 1.0e-12f) return fallback;
    const float inv = native_rsqrt(n2);
    return v * inv;
}

inline float3 project_onto_plane( const float3 v, const float3 n ){
    return v - n * dot(v,n);
}

inline float evalBond( const float4 h, const float dl, const float k, __private float3* f ){
    const float fr = dl * k;
    *f = h.xyz * fr;
    return fr * dl * 0.5f;
}

inline float evalAngCos( const float4 hr1, const float4 hr2, const float K, const float c0, __private float3* f1, __private float3* f2 ){
    const float c  = dot(hr1.xyz, hr2.xyz);
    const float3 hf1 = hr2.xyz - hr1.xyz * c;
    const float3 hf2 = hr1.xyz - hr2.xyz * c;
    const float cdiff = c - c0;
    const float E    = K * cdiff * cdiff;
    const float fang = -K * cdiff * 2.0f;
    *f1 = hf1 * (fang * hr1.w);
    *f2 = hf2 * (fang * hr2.w);
    return E;
}

inline float evalPiAlign( const float3 h1, const float3 h2, const float K, __private float3* f1, __private float3* f2 ){
    float  c   = dot(h1,h2);
    float3 hf1 = h2 - h1*c;
    float3 hf2 = h1 - h2*c;
    const float sgn  = (c < 0.0f) ? -1.0f : 1.0f;
    c = fabs(c);
    const float E    = -K * c;
    const float fang = K * sgn;
    hf1 *= fang;
    hf2 *= fang;
    *f1 = hf1;
    *f2 = hf2;
    return E;
}

float3 rotate_by_omega(float3 p, float3 w){
    float r2 = dot(w,w);
    float r = sqrt(r2);
    float s = sin(r) / r;
    float c = (1.0f - cos(r)) / r2;    // equals (1-cos r)/r^2
    float3 wxp  = cross(w, p);
    float3 wwxp = cross(w, wxp);
    return p + wxp * s + wwxp * c;
}

__kernel void simulate_pi_pi_rot(
    __global float4* pos,       // history: [nSteps,4] {nodeA,nodeB,piA,piB}
    __global float4* vel,       // history: [nSteps,4] {velA,velB,omegaA,omegaB}
    __global float4* force,     // history: [nSteps,4] {forceA,forceB,tqA,tqB}
    const float       dt,
    const float       cdamp,
    const int         nSteps
){
    if(get_global_id(0) != 0) return;   // single-system kernel for now

    float3 pa = pos[0].xyz;
    float3 pb = pos[1].xyz;
    float3 ha = normalize_safe(pos[2].xyz, (float3)(0.0f, 1.0f, 0.0f));
    float3 hb = normalize_safe(pos[3].xyz, (float3)(0.0f,-1.0f, 0.0f));

    float3 va = vel[0].xyz;
    float3 vb = vel[1].xyz;
    float3 wa = vel[2].xyz;
    float3 wb = vel[3].xyz;

    const float bL0  = 1.5f;
    const float kBL  = 10.0f;
    const float kPP  = 1.0f;
    const float kSP  = 1.0f;
    const float iM_A = 1.0f;
    const float iM_B = 1.0f;
    const float iI_A = 1.0f;
    const float iI_B = 1.0f;

    //const float3 fallback_axis = (float3)(0.0f, 0.0f, 1.0f);

    for(int istep=0; istep<nSteps; istep++){

        printf("OCL::simulate_pipi(): step %i dt %f ========= \n", istep, dt );

        float3 d = pb - pa;
        float  l = length(d);
        float inv_l = 1.0f/l;
        float3 h   = d * inv_l;
        //float3 h   = (rab > 1.0e-8f) ? dab * (1.0f / rab) : fallback_axis;
        float4 hbnd = (float4){ h.x, h.y, h.z, inv_l };

        // ====== Forcefield

        // ----- bond
        float3 fbond;
        float  Ebl = evalBond( hbnd, l - bL0, kBL, &fbond );
        float3 fA  =  fbond;
        float3 fB  = -fbond;

        float3 tqA = float3Zero, tqB = float3Zero;

        float Epp = 0.0f;
        //----- pi-pi align
        float3 fpa, fpb;
        Epp = evalPiAlign( ha, hb, kPP, &fpa, &fpb );
        tqA += cross(ha, fpa);
        tqB += cross(hb, fpb);
        //printf("OCL::simulate_pipi(): pi-pi:     iG %3i ing %3i epp %10.5f f(%10.5f,%10.5f,%10.5f) c %10.5f kpp %10.5f \n", 0, 4, Epp, fpa.x, fpa.y, fpa.z, dot(ha,h.xyz), kPP );
        printf("OCL::simu_pipi(): pi-pi:     iG %3i ing %3i epp %10.5f kpp %10.5f hpi(%10.5f,%10.5f,%10.5f) hpj(%10.5f,%10.5f,%10.5f)               f1(%10.5f,%10.5f,%10.5f) f2(%10.5f,%10.5f,%10.5f) \n", 0, 4, Epp, kPP,  ha.x,ha.y,ha.z,  hb.x,hb.y,hb.z,  fpa.x,fpa.y,fpa.z,  fpb.x,fpb.y,fpb.z );

        // ----- sigma-pi orthogonalization (each pi vs. bond axis)
        float3 f1, f2;
        float  EspA = evalAngCos( (float4){ha.x, ha.y, ha.z, 1.0f}, hbnd, kSP, 0.0f, &f1, &f2 );
        tqA += cross(ha, f1);
        fA  -= f2;
        fB  += f2;
        //printf("OCL::simulate_pipi(): pi-sigma:  iG %3i ing %3i esp %10.5f f1(%10.5f,%10.5f,%10.5f) c %10.5f ksp %10.5f \n", 0, 0, EspA, fpa_sig.x, fpa_sig.y, fpa_sig.z, dot(ha,h.xyz), kSP );
        printf("OCL::simu_pipi(): pi-sigma:  iG %3i ing %3i esp %10.5f ksp %10.5f hpi(%10.5f,%10.5f,%10.5f) h  (%10.5f,%10.5f,%10.5f|l:%10.5f|) f1(%10.5f,%10.5f,%10.5f) f2(%10.5f,%10.5f,%10.5f) \n", 0, 0, EspA, kSP,  ha.x,ha.y,ha.z,  h.x,h.y,h.z,l,  f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z );


        const float4 hbnd_rev = (float4){ -hbnd.x, -hbnd.y, -hbnd.z, hbnd.w };
        float  EspB = evalAngCos( (float4){hb.x, hb.y, hb.z, 1.0f}, hbnd_rev, kSP, 0.0f, &f1, &f2 );
        tqB += cross(hb, f1);
        fB  -= f2;
        fA  += f2;
        //printf("OCL::simulate_pipi(): pi-sigma:  iG %3i ing %3i esp %10.5f f1(%10.5f,%10.5f,%10.5f) c %10.5f ksp %10.5f \n", 1, 0, EspB, fpb_sig.x, fpb_sig.y, fpb_sig.z, dot(hb,h.xyz), kSP );
        printf("OCL::simu_pipi(): pi-sigma:  iG %3i ing %3i esp %10.5f ksp %10.5f hpi(%10.5f,%10.5f,%10.5f) h  (%10.5f,%10.5f,%10.5f|l:%10.5f|) f1(%10.5f,%10.5f,%10.5f) f2(%10.5f,%10.5f,%10.5f) \n", 1, 0, EspB, kSP,  hb.x,hb.y,hb.z,  h.x,h.y,h.z,l,   f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z );


        // ====== Propagation

        printf("OCL::simulate_pipi():         pi  iG %3i  fe(%10.5f,%10.5f,%10.5f)  pe(%10.5f,%10.5f,%10.5f)  ve(%10.5f,%10.5f,%10.5f) \n", 2, tqA.x,tqA.y,tqA.z, ha.x,ha.y,ha.z, wa.x,wa.y,wa.z );
        printf("OCL::simulate_pipi():         pi  iG %3i  fe(%10.5f,%10.5f,%10.5f)  pe(%10.5f,%10.5f,%10.5f)  ve(%10.5f,%10.5f,%10.5f) \n", 3, tqB.x,tqB.y,tqB.z, hb.x,hb.y,hb.z, wb.x,wb.y,wb.z );

        va += (fA * iM_A) * dt;
        vb += (fB * iM_B) * dt;
        va *= cdamp;
        vb *= cdamp;

        wa += (tqA * iI_A) * dt;
        wb += (tqB * iI_B) * dt;
        wa *= cdamp;
        wb *= cdamp;

        pa += va * dt;
        pb += vb * dt;

        // ha += cross(wa, ha) * dt;
        // hb += cross(wb, hb) * dt;
        // ha = normalize_safe(ha, project_onto_plane(fallback_axis, h));
        // hb = normalize_safe(hb, project_onto_plane(-fallback_axis, h));
        // wa -= ha * dot(wa, ha);
        // wb -= hb * dot(wb, hb);

        ha = rotate_by_omega( ha, wa*dt);  // rotate p by omega*dt
        hb = rotate_by_omega( hb, wb*dt);  // rotate p by omega*dt
        ha = normalize(ha);                // Normalize to correct for numerical drift
        hb = normalize(hb);                // Normalize to correct for numerical drift

        const int base = istep * 4;
        pos  [base+0] = (float4){ pa.x, pa.y, pa.z, Ebl + EspA };
        pos  [base+1] = (float4){ pb.x, pb.y, pb.z, Ebl + EspB };
        pos  [base+2] = (float4){ ha.x, ha.y, ha.z, Epp + EspA };
        pos  [base+3] = (float4){ hb.x, hb.y, hb.z, Epp + EspB };

        vel  [base+0] = (float4){ va.x, va.y, va.z, 0.0f };
        vel  [base+1] = (float4){ vb.x, vb.y, vb.z, 0.0f };
        vel  [base+2] = (float4){ wa.x, wa.y, wa.z, 0.0f };
        vel  [base+3] = (float4){ wb.x, wb.y, wb.z, 0.0f };

        force[base+0] = (float4){ fA.x, fA.y, fA.z, 0.0f };
        force[base+1] = (float4){ fB.x, fB.y, fB.z, 0.0f };
        force[base+2] = (float4){ tqA.x, tqA.y, tqA.z, 0.0f };
        force[base+3] = (float4){ tqB.x, tqB.y, tqB.z, 0.0f };
    }
}
