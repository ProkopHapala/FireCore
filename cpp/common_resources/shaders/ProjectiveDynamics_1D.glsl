precision highp float;

// --- Uniforms ---
uniform float dt;
uniform vec2 cmix; // .x = new_mix, .y = old_mix (heavy ball parameters)

// --- Textures ---
uniform sampler2D tPos;   // Current Position (p_k)
uniform sampler2D tRef;   // Predicted/Inertial Target (s_n)
uniform sampler2D tMom;   // Previous Momentum/Delta (d_k)
uniform sampler2D tBonds; // Bond Data (N x 16)

// --- Outputs (MRT) ---
layout(location = 0) out vec4 outPos; // p_{k+1}
layout(location = 1) out vec4 outMom; // d_{k+1} (stored for next step)

void main() {
    // 1. Get Atom Index
    // Since we use 1D textures (N x 1), Y is 0. X is the atom index.
    ivec2 coord = ivec2(gl_FragCoord.xy);
    int atomIdx = coord.x;

    // 2. Fetch Self Data
    vec4 posData = texelFetch(tPos, coord, 0);
    vec3 p_curr  = posData.rgb;
    float mass   = posData.a;
    
    // Fetch Inertial Reference (s_n)
    // This assumes tRef was pre-calculated: s_n = p_old + v*dt
    vec3 s_n = texelFetch(tRef, coord, 0).rgb;

    // 3. Handle Fixed Points (Infinite Mass)
    // Adjust threshold as needed for your specific "fixed" flag
    if (mass <= 0.0 || mass > 10000.0) {
        outPos = vec4(p_curr, mass);
        outMom = vec4(0.0);
        return;
    }

    // 4. Initialize Diagonal (Aii) and RHS (b)
    // The inertial term M/dt^2 corresponds to a spring pulling towards s_n
    float idt2 = 1.0 / (dt * dt);
    float w    = mass * idt2; 

    float Aii = w;       // Diagonal accumulation
    vec3  b   = w * s_n; // RHS accumulation

    // 5. Loop over Bonds (Constraint Projection)
    // tBonds is N wide, 16 high. 
    // Column = Atom Index, Row = Bond Index.
    for (int k = 0; k < 16; k++) {
        // Fetch bond data at column=atomIdx, row=k
        vec4 bond = texelFetch(tBonds, ivec2(atomIdx, k), 0);

        // Check if bond exists
        if (bond.a < 0.0) break; // -1 padding indicates end of list

        int neighIdx = int(bond.r);
        float l0     = bond.g;
        float kStiff = bond.b; // Stiffness

        // Fetch Neighbor Position (Linear Jacobi uses previous iteration's p)
        vec3 p_j = texelFetch(tPos, ivec2(neighIdx, 0), 0).rgb;

        // Projective Dynamics Constraint:
        // We want to satisfy distance constraint |p_i - p_j| = l0
        // The projection of p_i is: p_j + dir * l0
        
        vec3 dij   = p_curr - p_j;
        float dist = length(dij);
        
        // Safety for zero distance
        if (dist < 1e-9) {
            dij = vec3(0.0, 1.0, 0.0); // Arbitrary direction
            dist = 1.0;
        }

        // Projective Dynamics formulation:
        // Add stiffness to diagonal
        Aii += kStiff;
        
        // Add force/projection to RHS
        // b += k * projection_target
        vec3 projection = p_j + (dij / dist) * l0;
        b += kStiff * projection;
    }

    // 6. Jacobi Solve
    // p_jacobi = A^-1 * b
    vec3 p_jacobi = b / Aii;

    // 7. Heavy Ball / Momentum Acceleration
    // You requested: dpos = dpos_new * cmix.x + dpos_old * cmix.y
    
    // dpos_new is the change suggested by this Jacobi step
    vec3 dpos_new = p_jacobi - p_curr;
    
    // dpos_old is the change from the previous iteration
    vec3 dpos_old = texelFetch(tMom, coord, 0).rgb;

    // Apply Mixing
    vec3 dpos_final = dpos_new * cmix.x + dpos_old * cmix.y;

    // 8. Output
    vec3 p_final = p_curr + dpos_final;

    outPos = vec4(p_final, mass);
    outMom = vec4(dpos_final, 0.0); // Store momentum for next step
}