//>>>function  getLJQH (dp,REQ,ffpars)

inline float4 getLJQH( float3 dp, float4 REQ, float R2damp ){
    // ---- Electrostatic (damped Coulomb potential)
    float   r2    = dot(dp,dp);
    float   ir2_  = 1.f/(  r2 +  R2damp);              // inverse distance squared and damped
    float   Ec    =  COULOMB_CONST*REQ.z*sqrt( ir2_ ); // Ec = Q1*Q2/sqrt(r^2+R2damp)
    // --- Lennard-Jones and Hydrogen bond correction
    float  ir2 = 1.f/r2;          // inverse distance squared
    float  u2  = REQ.x*REQ.x*ir2; // u2 = (R0/r)^2
    float  u6  = u2*u2*u2;        // u6 = (R0/r)^6
    float vdW  = u6*REQ.y;        // vdW = E0*(R0/r)^6
    float E    =       (u6-2.f)*vdW     + Ec  ;     // E = E0*(R0/r)^6 - E0*(R0/r)^12 + Q1*Q2/sqrt(r^2+R2damp)
    float fr   = -12.f*(u6-1.f)*vdW*ir2 - Ec*ir2_;  // fr = -12*E0*( (R0/r)^8/r + 12*E0*(R0/r)^14) - Q1*Q2/(r^2+R2damp)^1.5
    return  (float4){ dp*fr, E };
}

//>>>function getMorseQH (dp,REQH,ffpars )

inline float4 getMorseQH( float3 dp,  float4 REQH, float K, float R2damp ){
    float r2    = dot(dp,dp);
    float ir2_  = 1/(r2+R2damp);
    float r     = sqrt( r2   );
    float ir_   = sqrt( ir2_ );     // ToDo: we can save some cost if we approximate r^2 = r^2 + R2damp;
    float e     = exp ( K*(r-REQH.x));
    //double e2    = e*e;
    //double fMors =  E0*  2*K*( e2 -   e ); // Morse
    //double EMors =  E0*      ( e2 - 2*e );
    float   Ae  = REQH.y*e;
    float fMors = Ae*  2*K*(e - 1); // Morse
    float EMors = Ae*      (e - 2);
    float Eel   = COULOMB_CONST*REQH.z*ir_;
    float fr    = fMors/r - Eel*ir2_ ;
    return  (float4){ dp*fr, EMors+Eel };
}

//>>>macro MODEL_LJQH2_PAIR
{
    // Distance safeguards
    float r_safe  = fmax(r, R2SAFE);
    float inv_r = 1.0f / r_safe;

    // Electrostatic
    float dE_dQ = inv_r * COULOMB_CONST;
    float Eel   = Q * dE_dQ;

    // Lennard-Jones 12-6 with H2 scaling of attraction
    float u   = R0 * inv_r;
    float u3  = u*u*u;
    float u6  = u3*u3;
    float u6p = (1.f + H) * u6;
    float dE_dE0 = u6 * (u6p - 2.f);
    float ELJ    =  E0 * dE_dE0;

    // Accumulate derivatives for atom i
    float dE_dR0 = 12.f * (E0/R0) * u6 * (u6p - 1.f);
    float dE_dH  = -E0 * u6 * u6;
    fREQi.x +=  -dE_dR0;
    fREQi.y +=  -dE_dE0 *        REQj.y;
    fREQi.z +=  -dE_dQ  *        REQj.z;
    fREQi.w +=  dE_dH  *        REQj.w * sH;

    // Accumulate energy
    Ei += ELJ + Eel;
}

//>>>macro MODEL_LJr8QH2_PAIR
{
    // Electrostatic
    float dE_dQ = inv_r * COULOMB_CONST;
    float Eel   = Q * dE_dQ;
    // r^-8 variant (8-6 like) with H2 scaling of r^-2 factor
    float u   = R0 * inv_r;
    float u2  = u*u;
    float u4  = u2*u2;
    float u6  = u4*u2;
    float u2p = (1.f + H) * u2;
    float dE_dE0 = u6 * (3.f * u2p - 4.f);
    float ELJ    =  E0 * dE_dE0;
    // Accumulate energy
    Ei += ELJ + Eel;

    // Accumulate derivatives for atom i
    float dE_dR0 = 24.f * (E0/R0) * u6 * (u2p - 1.f);
    float dE_dH  = -3.f * E0 * u6 * u2;
    fREQi.x +=  -dE_dR0;
    fREQi.y +=  -dE_dE0 *        REQj.y;
    fREQi.z +=  -dE_dQ  *        REQj.z;
    fREQi.w +=  dE_dH  *        REQj.w * sH;
}

//>>>macro MODEL_MorseQ_PAIR
{
    // Electrostatic
    float dE_dQ = inv_r * COULOMB_CONST;
    float Eel   = Q * dE_dQ;
    // Morse with alpha matching CPU kMorse = 1.8
    const float alpha = 1.8f;
    float e    = exp( -alpha * ( r - R0 ) );
    float e2   = e * e;
    float e2p  = (1.f + H) * e2;
    float dE_dE0 = e2p - 2.f * e;
    float ELJ    =  E0 * dE_dE0;
    Eij = ELJ + Eel;

    // Accumulate derivatives for atom i
    float dE_dR0 = 2.f * alpha * E0 * ( e2p - e );
    float dE_dH = - E0 * e2;
    // fREQi.x +=  dE_dR0;
    // fREQi.y +=  dE_dE0 *        REQj.y;
    // fREQi.z +=  dE_dQ  *        REQj.z;
    // fREQi.w +=  dE_dH *        REQj.w * sH;
    fij.x = -dE_dR0;                   // dEtot/dR0_i (match CPU)
    fij.y = -dE_dE0 * REQj.y;          // dEtot/dE0_i (match CPU)
    fij.z = -dE_dQ  * REQj.z;          // dEtot/dQi   (match CPU)
    fij.w =  dE_dH  * REQj.w * sH;     // dEtot/dH2i  (match CPU)
    // Accumulate energy
    
}

// Energy-only variants of the above models. These only accumulate Ei.
//>>>macro ENERGY_LJQH2_PAIR
{
    // Electrostatic
    float Eel   = Q * (inv_r * COULOMB_CONST);
    // Lennard-Jones 12-6 with H2 scaling on attraction
    float u   = R0 * inv_r;
    float u3  = u*u*u;
    float u6  = u3*u3;
    float u6p = (1.f + H)  *  u6;
    float ELJ =  E0 * ( u6 * (u6p - 2.f) );
    Ei += ELJ + Eel;
}

//>>>macro ENERGY_LJr8QH2_PAIR
{
    float Eel   = Q * (inv_r * COULOMB_CONST);
    float u   = R0 *   inv_r;
    float u2  = u*u;
    float u4  = u2*u2;
    float u6  = u4*u2;
    float u2p = (1.f + H) * u2;
    float ELJ =  E0 * ( u6 * (3.f * u2p - 4.f) );
    Ei += ELJ + Eel;
}

//>>>macro ENERGY_MorseQ_PAIR
{
    //float Eel   = Q*2.0 * (inv_r * COULOMB_CONST);
    float Eel   = Q * (inv_r * COULOMB_CONST);
    // Morse with alpha matching CPU kMorse = 1.8
    const float alpha = 1.8f;
    float e    = exp( -alpha * ( r - R0 ) );
    float e2   = e * e;
    float e2p  = (1.f + H) * e2;
    float ELJ  =  E0 * ( e2p - 2.f * e );
    Ei += ELJ + Eel;
}
