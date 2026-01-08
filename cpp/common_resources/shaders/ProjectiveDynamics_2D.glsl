precision highp float;

uniform sampler2D tPos;      // Current Iteration Positions (p_k) + mass in .a
uniform sampler2D tBonds;    // Bond Data (Wide texture: W*16 x H)
uniform sampler2D tMomentum; // Previous momentum (linsolve_yy)

uniform int texWidth;        // Width of atom grid
uniform float dt;
uniform float omega;         // Momentum parameter (mixer.get_bmix(i))

in vec2 vUV;

// We write to two textures simultaneously (MRT)
layout(location = 0) out vec4 outPos;      // p_{k+1}
layout(location = 1) out vec4 outMomentum; // d_{k+1}

ivec2 getUV(float index, int width) {
    int i = int(index + 0.5);
    return ivec2(i % width, i / width);
}

void main() {
    ivec2 centerUV = ivec2(gl_FragCoord.xy);
    
    // 1. Fetch Self Data (position xyz + mass in .a)
    vec4 pData = texelFetch(tPos, centerUV, 0);
    vec3 p_current = pData.rgb;
    float mass = pData.a;
    
    // Handle Fixed Points
    if (mass > 100000.0 || mass <= 0.0) {
        outPos = pData;
        outMomentum = vec4(0.0);
        return;
    }

    // 2. Setup PD Matrix Diagonal and RHS
    // Term 1: Inertia (M / dt^2)
    float idt2 = 1.0 / (dt * dt);
    float Aii = mass * idt2;
    
    // RHS starts with inertial target = current predicted position
    vec3 b_i = p_current * Aii; 

    // 3. Accumulate Constraints (The Loop)
    // Bond texture dimensions: width=maxBonds (columns per atom), height=N (one row per atom)
    int atomIndex = centerUV.x;
    
    for(int k = 0; k < 16; k++) {
        // Fetch bond info: R=NeighborIdx, G=L0, B=Stiffness, A=Active
        ivec2 bondUV = ivec2(k, atomIndex);
        vec4 bond = texelFetch(tBonds, bondUV, 0);

        if (bond.a < 0.0) break; // End of bond list

        float neighborIdx = bond.r;
        float l0          = bond.g;
        float stiffness   = bond.b;

        // Fetch Neighbor Position (from PREVIOUS iteration or current depending on sync)
        // In Jacobi, we use tPos which is the state at step k
        ivec2 neighborCoord = getUV(neighborIdx, texWidth);
        vec3 p_j = texelFetch(tPos, neighborCoord, 0).rgb;

        // Projective Dynamics Constraint Projection
        // We want to satisfy |p_i - p_j| = l0
        // The local step projection target for p_i relative to p_j is: p_j + dir * l0
        
        vec3 dij = p_current - p_j; // Vector from j to i
        float dist = length(dij);
        
        // Avoid division by zero
        if (dist < 1e-6) {
             dij = vec3(0.0, 1.0, 0.0); // Arbitrary safe dir
             dist = 1.0;
        }

        // Add to Diagonal: Aii += k
        Aii += stiffness;

        // Add to RHS: b_i += k * (p_j + (dij/dist) * l0)
        // This is the "Matrix Free" formulation of: b_i += k * projection
        vec3 projection = p_j + (dij / dist) * l0;
        b_i += stiffness * projection;
    }

    // 4. Jacobi Solve Step
    // p_jacobi = A^-1 * b
    vec3 p_jacobi = b_i / Aii;

    // 5. Heavy-Ball Momentum Acceleration
    // Corresponds to C++:
    // p = psb[j]; (p_jacobi)
    // p.add_mul( linsolve_yy[j], bmix ); 
    // linsolve_yy[j] = p - psa[j];
    
    vec3 momentum_prev = texelFetch(tMomentum, centerUV, 0).rgb;
    
    // Apply mixing
    vec3 p_accelerated = p_jacobi + momentum_prev * omega;
    
    // Calculate new momentum (delta)
    vec3 momentum_new = p_accelerated - p_current;

    // Output
    outPos = vec4(p_accelerated, mass);
    outMomentum = vec4(momentum_new, 0.0);
}