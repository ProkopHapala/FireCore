# USER

I made C++ library with ctypes python interface for fitting hydrogen bond correction they are of wo type H1,H2 are homogenous correction which decrease Pauli repulstion (repulsive part of Morse or Leanrd jones) and SR1,SR2,SR3 are applied to electron pairs/sigmaholes added as dummy atoms, which should introduce angular dependece. To better analyse / diagnose my fits and potential energy surfaces I thinking to decompose it to contribution of each atom pair and of different compomnent (consider Electrostatics+Morse as a base and  see how much is where energy decreased by H1,H2, resp Epairs/sigma hole). I'm thinking how to plot it best to be easy to navigate and understand for human scientest. 

I'm thinking to use stacked plots, these perhaps make sense only if components have same sign. And I think all the Hbond correction have same sign (negative - they are atractive, decrease energy), but the baseline is not zero but the Morse+Coulomb energy which is perhaps repulsive in relevant distance. So it would be like Positive baseline (repulsion) minus the stacked corrections, produce negative (bonding) potential profile. Is that fine? Can it be ploted easily?
https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.stackplot.html

Anay other idea how to plot it and analyse efficiently this problem to see what components contribute how much ?

```cpp
__attribute__((hot)) 
double evalEnergyDerivs ( int i0, int ni, int j0, int nj, int* types, Vec3d* ps, Quat4d* typeREQs, double* Qs, int* host, Quat4d* dEdREQs )const{
    double E_tot = 0.0;
//double E_Coul = 0.0, E_vdW = 0.0, E_Hcorr = 0.0, E_Epairs = 0.0;    
    //const double alpha   = kMorse;
//printf("JAMME     i     -dE_dR0_i       -dE_dRi         -dE_deps_i      -dE_dQ_i        -dE_dH_i        j     -dE_dR0_j       -dE_dRj         -dE_deps_j      -dE_dQ_j        -dE_dH_j        E_vdW           E_Hcorr         E_Coul          E_Epair\n");
    for(int ii=0; ii<ni; ii++){ // loop over all atoms[i] in system2
        const int     i     = i0+ii;
        const int     ih    = host[i];
        const bool    bEpi  = ih>=0;
        const Vec3d&  pi    = ps      [i ]; 
        const double  Qi    = Qs      [i ]; 
        const int     ti    = types   [i ];
        const Quat4d& REQi  = typeREQs[ti];
        Quat4d        fREQi = Quat4dZero;

        
        for(int jj=0; jj<nj; jj++){ // loop over all atoms[j] in system1
            const int     j    = j0+jj;
            const int     jh   = host[j];
            const bool    bEpj = jh>=0;
            const double  Qj   = Qs[j];
            const Vec3d   dij  = ps[j] - pi;
            const int     tj   = types[j];
            const Quat4d& REQj = typeREQs[tj];
            const double  R0   = REQi.x + REQj.x;
            const double  eps  = REQi.y * REQj.y; 
            const double  Q    = Qi     * Qj    ;
            double        H    = REQi.w * REQj.w;
            const double  sH   = (H<0.0) ? 1.0 : 0.0; 
            const double  r    = dij.norm();

            double Eij_Coul = 0.0, Eij_vdW = 0.0, Eij_Hcorr = 0.0, Eij_Epairs = 0.0;
            double dE_dR0 = 0.0, dE_deps = 0.0, dE_dQ = 0.0, dE_dH = 0.0;
            double fR = 0.0, fA = 0.0, fH1 = 0.0, fH2 = 0.0, fR0 = 0.0, fB = 1.0;
            double alpha = kMorse;
//double dE_dRi = 0.0, dE_dRj = 0.0; // JAMME

            if( bEpi ){  
                if(bEpj) continue; // dummy atoms should not interact with each other
/*
if(bEpj){
    double dE_dR = 0.0;
    if(iEpairs==1){ 
        Eij_Epairs = getSR_PN( r, H, REQj.x, dE_dH, dE_dR );
    }else if(iEpairs==2){ 
        Eij_Epairs = getSR2_PN( r, H, REQj.x, dE_dH, dE_dR );
    }
    dEdREQs[j].x -= dE_dR;
}
*/
                // --- Electron pair interaction
                double dE_dR = 0.0;
                if(iEpairs==1){ 
                    Eij_Epairs = getSR_PN( r, H, REQi.x, dE_dH, dE_dR );
                }else if(iEpairs==2){ 
                    Eij_Epairs = getSR2_PN( r, H, REQi.x, dE_dH, dE_dR );
                }else if(iEpairs==3){ 
                    Eij_Epairs = getSR3_PN( r, H, REQi.x, dE_dH, dE_dR );
                }
                dEdREQs[i].x -= dE_dR;
//dE_dRi = dE_dR; // JAMME
            }else if( bEpj ){
                // --- Electron pair interaction
                double dE_dR = 0.0;
                if(iEpairs==1){ 
                    Eij_Epairs = getSR_PN( r, H, REQj.x, dE_dH, dE_dR );
                }else if(iEpairs==2){ 
                    Eij_Epairs = getSR2_PN( r, H, REQj.x, dE_dH, dE_dR );
                }else if(iEpairs==3){ 
                    Eij_Epairs = getSR3_PN( r, H, REQj.x, dE_dH, dE_dR );
                }
                dEdREQs[j].x -= dE_dR;
//dE_dRj = dE_dR; // JAMME
            }else{
                // --- Electrostatic interaction
                if(iCoul==1){ // point charges
                    dE_dQ = COULOMB_CONST / r ;
                }else if(iCoul==2){ // point charges with softclamp
                    dE_dQ = dampCoulomb_SoftClamp(r, clamp_y1, clamp_y2) * COULOMB_CONST;
                }else if(iCoul>9){ // Boys clamping with different approximations
                    dE_dQ = dampCoulomb_Boys(r, boys_rmin, iCoul-10) * COULOMB_CONST;
                }
                Eij_Coul = Q * dE_dQ;                
                // --- Van der Waals interaction
                if(ivdW==1){ // Lennard-Jones 12-6
                    const double u  = R0 / r;
                    const double u3 = u * u * u;
                    fA              = u3 * u3;
                    fR              = fA * fA;
                    fH1             = 1.0;
                    fH2             = 2.0;
                    fR0             = 12.0 / R0;
                }else if(ivdW==2){ // Lennard-Jones 8-6
                    const double u  = R0 / r;
                    const double u2 = u * u;
                    fA              = u2 * u2 * u2;
                    fR              = fA * u2;
                    fH1             = 3.0;
                    fH2             = 4.0;
                    fR0             = 24.0 / R0;
                }else if(ivdW==3){ // Lennard-Jones 9-6
                    const double u  = R0 / r;
                    const double u3 = u * u * u;
                    fA              = u3 * u3;
                    fR              = fA * u3;
                    fH1             = 2.0;
                    fH2             = 3.0;
                    fR0             = 18.0 / R0;
                }else if(ivdW==4){ // Morse
                    if( kMorse < 0.0 ){ 
                        alpha = 6.0 / R0; 
                        fR0   = 2.0 * alpha * ( 1.0 + (r - R0) / R0);
                    }else{
                        fR0   = 2.0 * alpha;
                    }
                    fA                 = safe_exp2( -alpha * ( r - R0 ) );
                    fR                 = fA * fA;
                    fH1                = 1.0;
                    fH2                = 2.0;
                }else if(ivdW==5){ // Buckingham
                    if( kMorse < 0.0 ){ 
                        alpha = 6.0 / R0; 
                        fB    = 1.0 + (r - R0) / R0;
                    }
                    const double u     = R0 / r;
                    const double u3    = u * u * u;
                    const double e     = safe_exp2( -alpha * ( r - R0 ) );
                    fA                 = u3 * u3;
                    fR                 = e * e;
                    fH1                = 1.0;
                    fH2                = 2.0;
                    fR0                = 2.0 * alpha;
                }
                dE_deps = fH1 * fR - fH2 * fA;
                dE_dR0  = eps * fR0 * ( fB * fR - fA );
                Eij_vdW   = eps * dE_deps;
                // --- Hydrogen-bond corrections
                if(sH>0.0){
                    if(iHbond==1||iHbond==3){
                        const double f         = fH1 * fR;
                        const double dE_deps_H = H * f;
                        dE_deps       += dE_deps_H;
                        dE_dR0        += H * eps * fB * fR0 * fR;
                        dE_dH          = eps * f;
                        Eij_Hcorr      = eps * dE_deps_H;
                    }
                    if(iHbond==2||iHbond==3){
                        const double f         = fH2 * fA;
                        const double dE_deps_H = H * f;
                        dE_deps       += dE_deps_H;
                        dE_dR0        += H * eps * fR0 * fA;
                        dE_dH         += eps * f;
                        Eij_Hcorr     += eps * dE_deps_H;
                    }
                }
            }

            double Eij = Eij_Coul + Eij_vdW + Eij_Hcorr + Eij_Epairs;
//E_Coul += Eij_Coul; E_vdW += Eij_vdW; E_Hcorr += Eij_Hcorr; E_Epairs += Eij_Epairs;            
/*
printf("JAMME %5i %15.9f %15.9f %15.9f %15.9f %15.9f %5i %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n",
j+1, -dE_dR0,-dE_dRj,-dE_deps*REQi.y,-dE_dQ*Qi,-dE_dH*REQi.w, 
i+1, -dE_dR0,-dE_dRi,-dE_deps*REQj.y,-dE_dQ*Qj,-dE_dH*REQj.w, Eij_vdW, Eij_Hcorr, Eij_Coul, Eij_Epairs);
*/
            dEdREQs[j].sub( Quat4d{ 
                dE_dR0, 
                dE_deps * REQi.y, 
                dE_dQ   * Qi, 
                dE_dH   * REQi.w 
            } );

            // --- Energy and forces
            E_tot    +=  Eij;
            fREQi.x -= dE_dR0;              // dEtot/dR0_i
            fREQi.y -= dE_deps * REQj.y;    // dEtot/dE0_i
            fREQi.z -= dE_dQ   * Qj;        // dEtot/dQ_i
            fREQi.w -= dE_dH   * REQj.w;    // dEtot/dH_i
        }  // end loop over all atoms[j] in system1

        dEdREQs[i].add(fREQi);
    }  // end loop over all atoms[i] in system2
/*
printf("JAMMETOT     Etot\n");
printf("JAMMETOT %15.9f\n", E_tot );
printf("JAMMETOT     i     x_i             y_i             z_i             R0_i            eps_i           q_i             H_i             F_i[1]          F_i[2]          F_i[3]          F_i[4]\n");
for(int jj=0; jj<nj; jj++){ // loop over all atoms[j] in system1
    const int     j     = j0+jj;
    const int     tj    = types   [j ];
    const Quat4d& REQj  = typeREQs[tj];
    printf("JAMMETOT %5i %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
        j+1, ps[j].x,ps[j].y,ps[j].z, REQj.x,REQj.y,Qs[j],REQj.w, dEdREQs[j].x,dEdREQs[j].y,dEdREQs[j].z,dEdREQs[j].w );
}
for(int ii=0; ii<ni; ii++){ // loop over all atoms[i] in system2
    const int     i     = i0+ii;
    const int     ti    = types   [i ];
    const Quat4d& REQi  = typeREQs[ti];
    printf("JAMMETOT %5i %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
        i+1, ps[i].x,ps[i].y,ps[i].z, REQi.x,REQi.y,Qs[i],REQi.w, dEdREQs[i].x,dEdREQs[i].y,dEdREQs[i].z,dEdREQs[i].w );
}
*/
//exit(0);
    return E_tot;
}
```

---

# Gemini-3-pro

This is a classic problem in forcefield development: visualizing how "add-on" terms modify the base potential.

Your idea of a stacked plot is conceptually right (you want to see the sum of parts), but **standard `stackplot` is dangerous here**. `stackplot` usually assumes all components stack cumulatively from $y=0$. Since your Morse/LJ potential crosses zero (positive repulsion, negative attraction) and Coulomb can be either, a standard stack plot often results in a visual mess where "negative" layers cross "positive" layers.

Here is the recommended approach to visualize this effectively for human analysis, along with the Python code to do it.

### 1. The "Decomposition" Plot (Best for Analysis)

Instead of a raw stack, use a **Cumulative Line + Fill** approach.
You visualize the potential in "stages":
1.  **Baseline** (Morse + Coulomb)
2.  **+ H-bond Corrections** (The effect of H1, H2)
3.  **+ Angular/Dummy** (The effect of SR1, SR2, SR3) -> **Total**

This allows you to see exactly how the corrections deepen the well or shift the minimum.

#### Python Implementation using `fill_between`
This is superior to `stackplot` because it handles sign changes gracefully.

```python
import numpy as np
import matplotlib.pyplot as plt

def plot_potential_decomposition(r, E_coul, E_vdw, E_hcorr, E_epair):
    """
    r: array of distances
    E_*: arrays of energy components
    """
    
    # 1. Define the cumulative sums
    base_energy = E_coul + E_vdw
    base_plus_H = base_energy + E_hcorr
    total_energy = base_plus_H + E_epair
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 2. Plot the Total and Base as lines
    ax.plot(r, base_energy, 'k--', linewidth=1.5, label='Baseline (Morse+Coul)', alpha=0.6)
    ax.plot(r, total_energy, 'k-', linewidth=2.5, label='Total Energy')
    
    # 3. Fill the contributions
    # Contribution of H-bond corrections (Difference between Base and Base+H)
    ax.fill_between(r, base_energy, base_plus_H, 
                    color='cyan', alpha=0.4, label='H-bond Corr (H1/H2)')
    
    # Contribution of Epairs/SigmaHoles (Difference between Base+H and Total)
    ax.fill_between(r, base_plus_H, total_energy, 
                    color='orange', alpha=0.4, label='E-pair/Sigma (SR terms)')

    # Optional: Plot the components individually as thin lines to see shape
    # ax.plot(r, E_hcorr, 'c:', linewidth=1)
    # ax.plot(r, E_epair, 'orange', linestyle=':', linewidth=1)

    # Styling
    ax.set_ylim(-10.0, 10.0) # Set reasonable limits to ignore extreme repulsion
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.set_xlabel('Distance [Å]')
    ax.set_ylabel('Energy [kcal/mol]')
    ax.set_title('Decomposition of Interaction Energy')
    ax.legend()
    
    return fig, ax

# --- Dummy Data Generation for testing ---
r = np.linspace(1.5, 5.0, 100)
# Fake Morse
E_vdw = 5.0 * ((1 - np.exp(-1.5*(r-3.0)))**2 - 1) 
# Fake Coulomb (slightly attractive)
E_coul = -332.0 * 0.2 / r 
# Fake H-bond (short range attraction)
E_hcorr = -2.0 * np.exp(-2.0 * (r - 2.8)**2)
# Fake Epair (angular/short range)
E_epair = -1.5 * np.exp(-4.0 * (r - 2.5)**2)

plot_potential_decomposition(r, E_coul, E_vdw, E_hcorr, E_epair)
plt.show()
```

### 2. The "Interaction Matrix" (Best for Atom-Pair Diagnosis)

Since you mentioned decomposing by **atom pair**, a 1D plot gets messy if you have 20 atoms. Use a **Heatmap** to see which specific pairs are dominating the energy.

This helps you spot if a "dummy atom" (host < 0) is interacting too strongly with a real atom.

```python
def plot_interaction_matrix(n_atoms, interaction_data):
    """
    interaction_data: 2D array (N, N) where value is E_total for that pair
    """
    fig, ax = plt.subplots(figsize=(8,8))
    im = ax.imshow(interaction_data, cmap='coolwarm', vmin=-5, vmax=5)
    plt.colorbar(im, label='Pair Energy')
    
    ax.set_title("Atom-Atom Interaction Contributions")
    ax.set_xlabel("Atom Index System A")
    ax.set_ylabel("Atom Index System B")
    plt.show()
```

### 3. Modifying C++ for Analysis

Your current C++ function returns `E_tot`. To generate these plots, you can't just modify the return value. You have two efficient options to get the data out without slowing down your fitting loop:

#### Option A: The "Grid Mode" (Recommended)
Add a function to your library specifically for plotting. It takes a trajectory (a list of distances/positions) and returns the decomposed energy arrays (Coul, vdW, H, Epair) directly to Python.

#### Option B: The "Debug Pointer"
Modify your existing function to accept an optional debug pointer.

```cpp
// Inside C++
double evalEnergyDerivs(..., double* debug_out) {
    // ... inside loop ...
    double Eij_Coul = ..., Eij_vdW = ..., Eij_Hcorr = ..., Eij_Epairs = ...;
    
    // Accumulate total components if debug_out is not null
    if (debug_out != nullptr) {
        debug_out[0] += Eij_Coul;
        debug_out[1] += Eij_vdW;
        debug_out[2] += Eij_Hcorr;
        debug_out[3] += Eij_Epairs;
    }
    // ...
}
```

### 4. Handling the Angular Dependence (SR terms)

Since your SR terms (`E_Epair`) are angular (sigma holes), a standard Energy vs Distance ($R$) plot is misleading because it assumes a fixed angle.

**Recommendation:** Create 2 plots.
1.  **Radial Scan:** Fix the angle at the "ideal" hydrogen bond angle (e.g., 180° or 120°) and scan distance $R$. Use the "Fill Between" plot above.
2.  **Angular Scan:** Fix distance $R$ at the minimum ($R_{min}$) and scan the angle $\theta$.
    *   Plot `Total Energy` vs `Angle`.
    *   This will show how the `E_epair` term provides the directionality that isotropic Morse/LJ lacks.

### Summary
*   **Do not use standard stackplots**; the zero-crossing of the baseline makes them unintuitive.
*   **Use `fill_between`** (Base vs. Base+Corr) to visually highlight the "energy gain" provided by your H-bond corrections.
*   **Group your terms:**
    *   **Gray:** Physics Baseline (Coulomb + Morse/LJ)
    *   **Cyan:** Homogeneous Corrections (H1/H2)
    *   **Orange:** Directional Corrections (SR/Epairs)
*   The resulting plot will clearly show a "grey bucket" (vdW) getting deeper and shifting position due to the "colored layers" (corrections).

---

# USER

Great, now I realized what I actually wat to do visually. My potential maps are typically 2D poslar plots (angle,distance), I'm visualizing it using matplotlib imshow having 3 panels (DFT reference, my Model, and difference between the two, symatric 'bwr' colormap normalized by +/- energy of minimum - I cut out very repulsive region close to nuclei). 

Now I would like to make interactive version of this plot, where user can click certain point (pixel in the 2D map), which correspond to certain configuration, for this particular configuration it should plot 

1) 1D stacked plot we were talking above along angle and along distance coortinate (basically Data[iang,:] and Data[:,idist] slices ), it should show line without correction, final potential and in between stacked controbitions of H1/H2 and SR1/SR2/SR3 corections (they are exclusive) and eventually also decomposed to contribution of individual atoms or atom pairs (for small system like H2O it is possible, typically just few particles will me relevant as these are short-range exponentially decayiung contributions). Note that the 2D map imshow neglect true dimensions especially along distance (we have non uniform sampling, step is like dr=0.05A when we are close and potential is steep, while when we are far and potential is flat we have dr=5A ), this is fine for 2D plots, but tor the 1D stacked section we want to have true distance and angle on x-axis

2) I also want to visualize the geometry (also with added dummy atoms, electron pairs etc.) eventually with lines showing distaces between them (especially distance of electron pair to the atom from opposite fragment - we always have two fragments in each configuration and we only consider pairwise interaction between particles from different fragments in our model, as you can see in atached file there is n0 3 saying the 2nd fragment start from 3rd particle H2O-H2O=3+3 atoms, resp 5+5 after we add electron pairs on each oxygen or 7+7 when we add also sigma holes on each hydrogen ). 

I already have all configuration loaded in memory and wraped to numpy array over ctypes interface.

I want to make this visualization using matplotlib parhaps. But I'm also considering Vispy and pyQT as it is faster and nicer for plotting 3D geometry. 

Think it through how to do it best

this is example how looks my typical reference DFT data file in .xyz format:

### H2O-D1_H2O-A1-y.xyz

```
6
# n0 3 Etot .10510525787964335624 x0 01.50 z 00
O_3 -0.96870070846928291 5.5511151231257827E-017 0.0000000000000000 -0.67931
H_O -0.0000000000000000 0.0000000000000000 0.0000000000000000 0.339655
H_O -1.1823711352607018 0.94484177930196833 0.0000000000000000 0.339655
O_3 1.5000000000000000 0.0000000000000000 0.0000000000000000 -0.67931
H_O 2.1047306698843653 -0.75675747732684973 0.0000000000000000 0.339655
H_O 2.1047306623479072 0.75675746789575704 0.0000000000000000 0.339655
6
# n0 3 Etot .02838501404725735000 x0 01.55 z 00
O_3 -0.96870070846928291 5.5511151231257827E-017 0.0000000000000000 -0.67931
H_O -0.0000000000000000 0.0000000000000000 0.0000000000000000 0.339655
H_O -1.1823711352607018 0.94484177930196833 0.0000000000000000 0.339655
O_3 1.5500000000000000 0.0000000000000000 0.0000000000000000 -0.67931
H_O 2.1547306698843656 -0.75675747732684973 0.0000000000000000 0.339655
```

---

# Gemini-3.1-pro

This is an excellent idea for a diagnostic tool. Building an interactive dashboard is the best way to debug and understand parameterized forcefields. 

Here is a comprehensive plan and a complete, runnable **pure Matplotlib prototype** to achieve exactly what you want. 

### 1. Technology Choice: Matplotlib vs. PyQt/Vispy
*   **Start with Pure Matplotlib:** Since you already use it and have your arrays in memory, Matplotlib's event handling (`mpl_connect`) is perfectly capable of doing this layout (2D maps + 1D slices + 3D geometry). For small systems like water dimers (even with dummy atoms, ~14 particles max), Matplotlib's 3D engine is fast enough.
*   **When to upgrade to PyQt + PyVista:** If you eventually want GUI widgets (checkboxes to toggle specific atom pairs on/off, sliders for forcefield parameters to see live fitting updates), or if the 3D rotation in Matplotlib feels too laggy, you should move to **PyQt6 + pyqtgraph** (for lightning-fast 1D/2D) and **PyVista** (for incredible 3D rendering). 

### 2. Solving the `imshow` Non-Uniform Grid Problem
You mentioned `imshow` neglects true dimensions because your sampling is non-uniform ($dr=0.05$ close, $dr=0.5$ far). 
**Do not use `imshow` for this.** Use **`pcolormesh`**. 
`pcolormesh(angles, distances, Energy)` automatically stretches the cells according to your true coordinate arrays. This means your 2D maps will accurately reflect physical space, and clicking them will return actual angles and distances.

### 3. The Interactive Dashboard Code

Here is a complete, modular script. It creates mock data mimicking your setup, plots the 3 2D maps, the 2 1D slices, and a 3D geometry view, and ties them together with an interactive click event.

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

class PotentialDashboard:
    def __init__(self, angles, distances, E_dft, E_base, E_hcorr, E_sr, geometries, n0):
        self.angles = angles
        self.distances = distances
        
        # Energy arrays (Shape: len(angles), len(distances))
        self.E_dft = E_dft
        self.E_base = E_base
        self.E_hcorr = E_hcorr
        self.E_sr = E_sr
        self.E_model = E_base + E_hcorr + E_sr
        self.E_diff = self.E_model - self.E_dft
        
        self.geometries = geometries # Shape: (len(angles), len(distances), N_atoms, 3)
        self.n0 = n0 # Index where fragment 2 starts
        
        # Current selected indices
        self.curr_a = len(angles) // 2
        self.curr_d = len(distances) // 2
        
        self.setup_ui()
        self.update_plots()

    def setup_ui(self):
        self.fig = plt.figure(figsize=(16, 10))
        gs = GridSpec(2, 3, figure=self.fig)
        
        # --- Top Row: 2D Maps ---
        self.ax_dft = self.fig.add_subplot(gs[0, 0])
        self.ax_model = self.fig.add_subplot(gs[0, 1], sharex=self.ax_dft, sharey=self.ax_dft)
        self.ax_diff = self.fig.add_subplot(gs[0, 2], sharex=self.ax_dft, sharey=self.ax_dft)
        
        vmin, vmax = -10, 5 # Adjust based on your min energy
        
        # Using pcolormesh instead of imshow handles non-uniform grids natively!
        A, D = np.meshgrid(self.angles, self.distances, indexing='ij')
        self.ax_dft.pcolormesh(A, D, self.E_dft, cmap='coolwarm', vmin=vmin, vmax=vmax)
        self.ax_model.pcolormesh(A, D, self.E_model, cmap='coolwarm', vmin=vmin, vmax=vmax)
        self.ax_diff.pcolormesh(A, D, self.E_diff, cmap='bwr', vmin=-2, vmax=2)
        
        titles = ['DFT Reference', 'Model', 'Difference']
        for i, ax in enumerate([self.ax_dft, self.ax_model, self.ax_diff]):
            ax.set_title(titles[i])
            ax.set_xlabel('Angle')
            ax.set_ylabel('Distance')
            self.fig.canvas.mpl_connect('button_press_event', lambda event, a=ax: self.onclick(event, a))
            
        # Crosshairs for 2D maps
        self.crosshairs = []
        for ax in [self.ax_dft, self.ax_model, self.ax_diff]:
            ch, = ax.plot([], [], 'k+', markersize=15, markeredgewidth=2)
            self.crosshairs.append(ch)

        # --- Bottom Row: 1D Slices & 3D Geometry ---
        self.ax_slice_dist = self.fig.add_subplot(gs[1, 0])
        self.ax_slice_ang = self.fig.add_subplot(gs[1, 1])
        self.ax_geom = self.fig.add_subplot(gs[1, 2], projection='3d')
        
        plt.tight_layout()

    def onclick(self, event, ax):
        # Only respond if clicked inside one of the 2D maps
        if event.inaxes not in [self.ax_dft, self.ax_model, self.ax_diff]: return
        
        # Find nearest indices for the non-uniform grid
        self.curr_a = np.argmin(np.abs(self.angles - event.xdata))
        self.curr_d = np.argmin(np.abs(self.distances - event.ydata))
        
        self.update_plots()
        self.fig.canvas.draw_idle()

    def update_plots(self):
        ang_val = self.angles[self.curr_a]
        dist_val = self.distances[self.curr_d]
        
        # Update Crosshairs on 2D maps
        for ch in self.crosshairs:
            ch.set_data([ang_val], [dist_val])

        # --- Update 1D Distance Slice (Fixed Angle) ---
        self.ax_slice_dist.clear()
        self.ax_slice_dist.set_title(f"Distance Scan @ Angle = {ang_val:.1f}")
        self.ax_slice_dist.set_xlabel("Distance [Å]")
        self.ax_slice_dist.set_ylabel("Energy")
        
        base_d = self.E_base[self.curr_a, :]
        hcorr_d = self.E_hcorr[self.curr_a, :]
        sr_d = self.E_sr[self.curr_a, :]
        
        self.ax_slice_dist.plot(self.distances, base_d, 'k--', label='Baseline')
        self.ax_slice_dist.fill_between(self.distances, base_d, base_d + hcorr_d, alpha=0.5, label='+ Hcorr')
        self.ax_slice_dist.fill_between(self.distances, base_d + hcorr_d, base_d + hcorr_d + sr_d, alpha=0.5, label='+ SR (E-pair)')
        self.ax_slice_dist.plot(self.distances, self.E_dft[self.curr_a, :], 'r.', label='DFT')
        self.ax_slice_dist.axvline(dist_val, color='gray', linestyle=':')
        self.ax_slice_dist.set_ylim(-10, 5)
        self.ax_slice_dist.legend()

        # --- Update 1D Angle Slice (Fixed Distance) ---
        self.ax_slice_ang.clear()
        self.ax_slice_ang.set_title(f"Angular Scan @ Distance = {dist_val:.2f}")
        self.ax_slice_ang.set_xlabel("Angle")
        
        base_a = self.E_base[:, self.curr_d]
        hcorr_a = self.E_hcorr[:, self.curr_d]
        sr_a = self.E_sr[:, self.curr_d]
        
        self.ax_slice_ang.plot(self.angles, base_a, 'k--', label='Baseline')
        self.ax_slice_ang.fill_between(self.angles, base_a, base_a + hcorr_a, alpha=0.5)
        self.ax_slice_ang.fill_between(self.angles, base_a + hcorr_a, base_a + hcorr_a + sr_a, alpha=0.5)
        self.ax_slice_ang.plot(self.angles, self.E_dft[:, self.curr_d], 'r.')
        self.ax_slice_ang.axvline(ang_val, color='gray', linestyle=':')
        self.ax_slice_ang.set_ylim(-10, 5)

        # --- Update 3D Geometry ---
        self.ax_geom.clear()
        geom = self.geometries[self.curr_a, self.curr_d]
        
        # Plot Fragment 1 (e.g., first 3 atoms)
        f1 = geom[:self.n0]
        self.ax_geom.scatter(f1[:,0], f1[:,1], f1[:,2], c='blue', s=100, label='Frag 1')
        
        # Plot Fragment 2
        f2 = geom[self.n0:]
        self.ax_geom.scatter(f2[:,0], f2[:,1], f2[:,2], c='red', s=100, label='Frag 2')
        
        # Example: Draw interaction lines between specific particles (e.g., dummy to atom)
        # Assuming atom 0 is O in Frag 1, and atom self.n0+1 is H in Frag 2
        p1 = geom[0]
        p2 = geom[self.n0 + 1]
        self.ax_geom.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 'k:', linewidth=2)
        
        self.ax_geom.set_title("Configuration")
        self.ax_geom.set_box_aspect([1,1,1]) # Equal aspect ratio
        self.ax_geom.legend()

# ==========================================
# Mock Data Generation (To make script runnable)
# ==========================================
Na, Nd = 20, 30
angles = np.linspace(0, 90, Na)
# Non-uniform distances (dr is small at close range, large at long range)
distances = np.concatenate([np.linspace(1.5, 3.0, 20), np.linspace(3.2, 8.0, 10)]) 

# Create fake energy landscapes
A, D = np.meshgrid(angles, distances, indexing='ij')
E_base = 10.0 * np.exp(-1.5 * (D - 2.5)) - 2.0 / D
E_hcorr = -4.0 * np.exp(-3.0 * (D - 2.8)**2)
E_sr = -3.0 * np.exp(-2.0 * (D - 2.6)**2) * np.cos(np.radians(A))**2
E_dft = E_base + E_hcorr + E_sr + np.random.normal(0, 0.2, A.shape)

# Create fake geometries (Na, Nd, 6 atoms, 3 coords)
geometries = np.zeros((Na, Nd, 6, 3))
for i in range(Na):
    for j in range(Nd):
        # Fake coordinates moving apart based on distance and angle
        geometries[i, j, :3, 0] = -D[i,j] / 2
        geometries[i, j, 3:, 0] = D[i,j] / 2
        geometries[i, j, 3:, 1] = np.sin(np.radians(A[i,j]))

# Run Dashboard
dashboard = PotentialDashboard(angles, distances, E_dft, E_base, E_hcorr, E_sr, geometries, n0=3)
plt.show()
```

### Key Highlights of this Approach:
1.  **Handling the `n0 3` fragment logic:** Look at the `update_plots()` function under the 3D Geometry section. You pass `n0=3` to the class, and it naturally slices the coordinate array `geom[:self.n0]` and `geom[self.n0:]` to color Fragment 1 and Fragment 2 differently.
2.  **Interaction Lines:** I added a small piece of code in the 3D viewer that draws a dotted line (`k:`) between particle 0 (Oxygen) and particle 4 (Hydrogen of second fragment). You can loop through your interacting pairs here to draw lines visualizing the short-range forces.
3.  **Crosshairs:** When you click on *any* of the 3 maps, a black crosshair moves to your clicked location to orient you.
4.  **Instant Refresh:** Because the arrays are all loaded into memory via `ctypes`, indexing `self.E_hcorr[self.curr_a, :]` takes microseconds. The plots will update instantly when you click.

### Showing Individual Pair Contributions
You mentioned wanting to see decomposition by individual atoms/pairs in the 1D slice.
If you have an array shaped `(Na, Nd, N_pairs)` containing the energy of each pair, you can easily add this to the `ax_slice_dist` plot. 
Instead of a stacked fill, plot them as thin dotted lines under the main curve:
```python
# Assuming E_pairs is shaped (len(angles), len(distances), num_interacting_pairs)
for p in range(num_interacting_pairs):
    pair_energy = self.E_pairs[self.curr_a, :, p]
    # Only plot if it contributes significantly (e.g., < -0.1 kcal/mol)
    if np.min(pair_energy) < -0.1: 
        self.ax_slice_dist.plot(self.distances, pair_energy, ':', alpha=0.7)
```

---

# USER

NO, I misunderstood me. I like imshow, I prefer it over pcolormesh. The non uniformity is not problem for 2D map it is a feature, it allows to view the whole map but zooming the importaint region close to atoms, it is also much faster, I like it. But for 1D sections I want to see true distances and angles (true shape of potential not distorted by non uniform sampling), so there we use normal plot / stacked plot resp the fill_between()  function you proposed. Is that clear? Update your vision and desing. Also put more focs on how to render the connections between the particles for different fragments with clear labels of distances, we should have some memu where we can choose for which particle we want to plot conetions to other particles and also write distances on the lines as lables, we should swithc beteen modes - plot distance, or energy component ( Morse,Clulomb,Homogenous H-bond correction, SR component). Think this through deeply 

---

# Gemini-3.1-pro

I completely understand your vision now. The non-uniform grid visual stretching provided by `imshow` is a powerful tool for inspecting the dense regions (short-range) while still seeing the long-range tails in the same window. 

For the 1D plots, we will use the **true `distances` and `angles` arrays** for the x-axes, so the physical shape of the potential well is preserved without distortion.

For the 3D visualization, we will build a **diagnostic UI** directly into the Matplotlib window. It will include:
1.  **Radio buttons to switch the Label Mode** (Distance, Morse, Coulomb, H-bond, SR).
2.  **Buttons to cycle through the "Active Particle"**. It will draw lines connecting the active particle to all particles in the *other* fragment and write the value (energy or distance) directly on the midpoint of the line.

Here is the completely redesigned, runnable prototype.

### The Interactive Diagnostic Dashboard (Pure Python/Matplotlib)

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import RadioButtons, Button

class PotentialDashboard:
    def __init__(self, angles, distances, E_dft, E_base, E_hcorr, E_sr, geometries, pair_data, n0):
        self.angles = angles          # True physical angles
        self.distances = distances    # True physical distances
        
        # 2D Grids (Shape: len(angles), len(distances))
        self.E_dft = E_dft
        self.E_base = E_base
        self.E_hcorr = E_hcorr
        self.E_sr = E_sr
        self.E_model = E_base + E_hcorr + E_sr
        self.E_diff = self.E_model - self.E_dft
        
        # 3D geometries and pair interaction data
        self.geometries = geometries  # Shape: (Na, Nd, N_atoms, 3)
        self.pair_data = pair_data    # Shape: (Na, Nd, N_atoms, N_atoms, 4 components)
        self.n0 = n0                  # Index where Fragment 2 starts
        self.n_atoms = geometries.shape[2]
        
        # Application State
        self.curr_a = len(angles) // 2
        self.curr_d = len(distances) // 4
        self.focus_atom = 0           # The particle we are inspecting connections for
        self.label_mode = 'Distance'  # What to write on the 3D lines
        
        self.setup_ui()
        self.update_all()

    def setup_ui(self):
        self.fig = plt.figure(figsize=(18, 10))
        # Create a layout: 3 rows, 4 columns
        gs = GridSpec(3, 4, height_ratios=[1.2, 1, 0.4], width_ratios=[1, 1, 1.2, 0.4])
        
        # --- Top Row: 2D Imshow Maps (Distorted/Non-uniform) ---
        self.ax_dft = self.fig.add_subplot(gs[0, 0])
        self.ax_model = self.fig.add_subplot(gs[0, 1])
        self.ax_diff = self.fig.add_subplot(gs[0, 2])
        
        vmin, vmax = -10, 5
        # imshow treats rows as Y (angles) and cols as X (distances)
        self.im_dft = self.ax_dft.imshow(self.E_dft, origin='lower', aspect='auto', cmap='coolwarm', vmin=vmin, vmax=vmax)
        self.im_model = self.ax_model.imshow(self.E_model, origin='lower', aspect='auto', cmap='coolwarm', vmin=vmin, vmax=vmax)
        self.im_diff = self.ax_diff.imshow(self.E_diff, origin='lower', aspect='auto', cmap='bwr', vmin=-2, vmax=2)
        
        # Format imshow ticks to show TRUE physical values instead of pixel indices
        for ax, title in zip([self.ax_dft, self.ax_model, self.ax_diff], ['DFT', 'Model', 'Difference']):
            ax.set_title(title)
            # Distance ticks (X)
            xticks = np.linspace(0, len(self.distances)-1, 5, dtype=int)
            ax.set_xticks(xticks)
            ax.set_xticklabels([f"{self.distances[i]:.2f}" for i in xticks])
            # Angle ticks (Y)
            yticks = np.linspace(0, len(self.angles)-1, 5, dtype=int)
            ax.set_yticks(yticks)
            ax.set_yticklabels([f"{self.angles[i]:.1f}" for i in yticks])
            
            # Connect click event
            self.fig.canvas.mpl_connect('button_press_event', lambda event, a=ax: self.onclick_2d(event, a))
            
        self.crosshairs = [ax.plot([], [], 'k+', markersize=12, markeredgewidth=2)[0] for ax in [self.ax_dft, self.ax_model, self.ax_diff]]

        # --- Middle Row: 1D True-Scale Slices ---
        self.ax_slice_dist = self.fig.add_subplot(gs[1, 0])
        self.ax_slice_ang = self.fig.add_subplot(gs[1, 1])
        
        # --- Right Side: 3D Geometry ---
        self.ax_geom = self.fig.add_subplot(gs[1:, 2], projection='3d')
        
        # --- UI Panel ---
        # Label Mode Selection
        self.ax_radio = self.fig.add_subplot(gs[1, 3], facecolor='lightgoldenrodyellow')
        self.ax_radio.set_title("Line Labels")
        self.radio = RadioButtons(self.ax_radio, ('Distance', 'Morse', 'Coulomb', 'H-corr', 'SR (Epairs)'))
        self.radio.on_clicked(self.on_mode_change)
        
        # Atom Selection Buttons
        self.ax_btn_prev = self.fig.add_subplot(gs[2, 3])
        self.ax_btn_prev.set_position([0.86, 0.15, 0.05, 0.05]) # Manual tweak
        self.btn_prev = Button(self.ax_btn_prev, '< Prev')
        self.btn_prev.on_clicked(lambda e: self.change_atom(-1))
        
        self.ax_btn_next = self.fig.add_subplot(gs[2, 3])
        self.ax_btn_next.set_position([0.92, 0.15, 0.05, 0.05])
        self.btn_next = Button(self.ax_btn_next, 'Next >')
        self.btn_next.on_clicked(lambda e: self.change_atom(1))
        
        self.txt_atom = self.fig.text(0.89, 0.22, "", ha='center', fontsize=12, fontweight='bold')
        plt.tight_layout()

    def onclick_2d(self, event, ax):
        if event.inaxes not in [self.ax_dft, self.ax_model, self.ax_diff]: return
        # event.xdata/ydata are raw pixel indices in imshow
        x_idx = int(round(event.xdata))
        y_idx = int(round(event.ydata))
        
        # Clamp to bounds
        self.curr_d = max(0, min(len(self.distances) - 1, x_idx))
        self.curr_a = max(0, min(len(self.angles) - 1, y_idx))
        
        self.update_all()
        self.fig.canvas.draw_idle()

    def on_mode_change(self, label):
        self.label_mode = label
        self.update_3d()
        self.fig.canvas.draw_idle()
        
    def change_atom(self, step):
        self.focus_atom = (self.focus_atom + step) % self.n_atoms
        self.update_3d()
        self.fig.canvas.draw_idle()

    def update_all(self):
        # Update 2D Crosshairs
        for ch in self.crosshairs:
            ch.set_data([self.curr_d], [self.curr_a])
            
        self.update_1d_slices()
        self.update_3d()

    def update_1d_slices(self):
        ang_val = self.angles[self.curr_a]
        dist_val = self.distances[self.curr_d]
        
        # --- Distance Scan (Uses True Physical Distance Array) ---
        self.ax_slice_dist.clear()
        self.ax_slice_dist.set_title(f"Distance Scan @ Angle = {ang_val:.1f}°")
        self.ax_slice_dist.set_xlabel("True Distance [Å]")
        self.ax_slice_dist.set_ylabel("Energy [kcal/mol]")
        
        b_d = self.E_base[self.curr_a, :]
        h_d = self.E_hcorr[self.curr_a, :]
        s_d = self.E_sr[self.curr_a, :]
        
        # Stacked Fill
        self.ax_slice_dist.plot(self.distances, b_d, 'k--', label='Baseline')
        self.ax_slice_dist.fill_between(self.distances, b_d, b_d + h_d, color='cyan', alpha=0.5, label='+ Hcorr')
        self.ax_slice_dist.fill_between(self.distances, b_d + h_d, b_d + h_d + s_d, color='orange', alpha=0.5, label='+ SR (E-pair)')
        self.ax_slice_dist.plot(self.distances, self.E_dft[self.curr_a, :], 'r.', label='DFT')
        
        self.ax_slice_dist.axvline(dist_val, color='gray', linestyle=':')
        self.ax_slice_dist.set_ylim(-10, 5)
        self.ax_slice_dist.axhline(0, color='k', linewidth=0.5)
        self.ax_slice_dist.legend(loc='lower right', fontsize='small')

        # --- Angular Scan (Uses True Physical Angle Array) ---
        self.ax_slice_ang.clear()
        self.ax_slice_ang.set_title(f"Angular Scan @ Distance = {dist_val:.2f} Å")
        self.ax_slice_ang.set_xlabel("True Angle [deg]")
        
        b_a = self.E_base[:, self.curr_d]
        h_a = self.E_hcorr[:, self.curr_d]
        s_a = self.E_sr[:, self.curr_d]
        
        self.ax_slice_ang.plot(self.angles, b_a, 'k--')
        self.ax_slice_ang.fill_between(self.angles, b_a, b_a + h_a, color='cyan', alpha=0.5)
        self.ax_slice_ang.fill_between(self.angles, b_a + h_a, b_a + h_a + s_a, color='orange', alpha=0.5)
        self.ax_slice_ang.plot(self.angles, self.E_dft[:, self.curr_d], 'r.')
        
        self.ax_slice_ang.axvline(ang_val, color='gray', linestyle=':')
        self.ax_slice_ang.set_ylim(-10, 5)
        self.ax_slice_ang.axhline(0, color='k', linewidth=0.5)

    def update_3d(self):
        self.ax_geom.clear()
        geom = self.geometries[self.curr_a, self.curr_d]
        
        # Determine fragments
        frag1_idx = list(range(self.n0))
        frag2_idx = list(range(self.n0, self.n_atoms))
        
        # Plot Fragments
        f1 = geom[:self.n0]
        f2 = geom[self.n0:]
        self.ax_geom.scatter(f1[:,0], f1[:,1], f1[:,2], c='blue', s=150, alpha=0.6, label='Fragment 1')
        self.ax_geom.scatter(f2[:,0], f2[:,1], f2[:,2], c='red', s=150, alpha=0.6, label='Fragment 2')
        
        # Highlight Focus Atom
        focus_pos = geom[self.focus_atom]
        frag_name = "Frag 1" if self.focus_atom < self.n0 else "Frag 2"
        self.txt_atom.set_text(f"Focus: Atom {self.focus_atom}\n({frag_name})")
        self.ax_geom.scatter(*focus_pos, c='gold', s=250, edgecolors='black', linewidth=2, zorder=5)

        # Draw connections to opposite fragment
        target_atoms = frag2_idx if self.focus_atom < self.n0 else frag1_idx
        
        mode_idx_map = {'Morse': 0, 'Coulomb': 1, 'H-corr': 2, 'SR (Epairs)': 3}
        
        for tgt in target_atoms:
            tgt_pos = geom[tgt]
            
            # Plot line
            self.ax_geom.plot([focus_pos[0], tgt_pos[0]], 
                              [focus_pos[1], tgt_pos[1]], 
                              [focus_pos[2], tgt_pos[2]], 
                              'k--', alpha=0.7)
            
            # Calculate value for text label
            if self.label_mode == 'Distance':
                val = np.linalg.norm(focus_pos - tgt_pos)
                txt = f"{val:.2f} Å"
            else:
                comp_idx = mode_idx_map[self.label_mode]
                # Make sure we read symmetric pair matrix correctly:
                i, j = min(self.focus_atom, tgt), max(self.focus_atom, tgt)
                val = self.pair_data[self.curr_a, self.curr_d, i, j, comp_idx]
                # Only show text if interaction is non-zero
                if abs(val) < 1e-4: continue 
                txt = f"{val:+.2f}"
                
            # Place text at midpoint
            mid = (focus_pos + tgt_pos) / 2
            self.ax_geom.text(mid[0], mid[1], mid[2], txt, 
                              color='darkgreen', fontsize=10, fontweight='bold',
                              bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

        self.ax_geom.set_title(f"Geometry & Pair {self.label_mode}")
        self.ax_geom.set_box_aspect([1,1,1]) # Keep aspect ratio square
        self.ax_geom.set_axis_off()          # Hide 3D grid for cleaner view


# ==========================================
# Mock Data Generation (To make script runnable)
# ==========================================
Na, Nd = 20, 50
angles = np.linspace(0, 90, Na)
# Deliberately Non-Uniform Distance array
distances = np.concatenate([np.linspace(1.5, 3.0, 30), np.linspace(3.2, 8.0, 20)]) 

A, D = np.meshgrid(angles, distances, indexing='ij')

# Fake Energy Grids
E_base = 10.0 * np.exp(-1.5 * (D - 2.5)) - 2.0 / D
E_hcorr = -4.0 * np.exp(-3.0 * (D - 2.8)**2)
E_sr = -3.0 * np.exp(-2.0 * (D - 2.6)**2) * np.cos(np.radians(A))**2
E_dft = E_base + E_hcorr + E_sr + np.random.normal(0, 0.2, A.shape)

# Fake Geometries (Na, Nd, 6 atoms, 3 coords)
# n0 = 3 (First 3 atoms are Fragment 1, Next 3 are Fragment 2)
geometries = np.zeros((Na, Nd, 6, 3))
for i in range(Na):
    for j in range(Nd):
        geometries[i, j, :3, 0] = -D[i,j] / 2
        geometries[i, j, :3, 1] = [0, 0.5, -0.5] # Spread them out
        geometries[i, j, 3:, 0] = D[i,j] / 2
        # Angle changes Fragment 2's orientation
        geometries[i, j, 3:, 1] = [np.sin(np.radians(A[i,j])), 0.5, -0.5]

# Fake Pairwise Data: Shape (Na, Nd, Natoms, Natoms, 4 components)
# Components: [Morse, Coulomb, Hcorr, SR]
pair_data = np.zeros((Na, Nd, 6, 6, 4))
for i in range(Na):
    for j in range(Nd):
        # Fake an interaction between atom 0 (Frag 1) and atom 3 (Frag 2)
        pair_data[i, j, 0, 3, 0] = E_base[i, j] * 0.5    # Morse
        pair_data[i, j, 0, 3, 1] = E_base[i, j] * 0.5    # Coulomb
        pair_data[i, j, 0, 3, 2] = E_hcorr[i, j]         # Hcorr
        
        # Fake an SR interaction between atom 1 (Frag 1) and atom 4 (Frag 2)
        pair_data[i, j, 1, 4, 3] = E_sr[i, j]            # SR

# Run Dashboard
dashboard = PotentialDashboard(angles, distances, E_dft, E_base, E_hcorr, E_sr, geometries, pair_data, n0=3)
plt.show()
```

### Key Changes Addressing Your Feedback:

1. **`imshow` with Proper Ticks**: The 2D maps now use `imshow` (which gives you that non-uniform visual stretch mapping). However, I added code that calculates correct custom ticks. Instead of showing pixel indices `[0, 10, 20]`, the axes will show the physical values `[1.50, 3.20, 8.00]`.
2. **True Scale 1D Plots**: In `update_1d_slices()`, the `plot()` and `fill_between()` functions explicitly receive the physical `self.distances` and `self.angles` arrays. This guarantees the actual curve shape (the true Morse/LJ well profile) is plotted accurately on a linear X-axis.
3. **Interactive Menu for Particle Connections**:
    *   There is a large gold highlighted atom in the 3D plot indicating the **"Focus Atom"**.
    *   I added **`< Prev`** and **`Next >`** buttons in the bottom right corner. Clicking them cycles the focus atom from `0` up to `N_atoms - 1`.
    *   The code strictly separates fragments using your `n0` rule. It automatically figures out which fragment the Focus Atom is in, and loops over the particles in the *other* fragment to draw the connections.
4. **Interactive Component Switching**:
    *   There is a `RadioButtons` menu (highlighted yellow).
    *   Clicking **"Distance"** calculates `norm(pos_i - pos_j)` and draws it on the line.
    *   Clicking **"Morse", "Coulomb", "H-corr", or "SR (Epairs)"** looks up the specific `pair_data[iang, idist, atom_i, atom_j, mode_index]` array, pulling out the exact interaction strength for that specific pair at that specific point in the map, and writes it directly onto the middle of the 3D line in the plot.

---

# USER

now review the whole dicusion in @FitREQ_discussion.md 
and write detaied tutorial how to use these new visualization tools @visualize_static.py @visualizer_gui.py 
and also how tey works (what we added in @FitREQ_PN_lib.cpp @FitREQ_PN.h @FitREQ_PN.py as new functions? first review the files read the git status and then write notes in this chate

---

# GPT-5.1-Codex

Here’s a concise tutorial and review of the new visualization tools and underlying changes.

## What changed under the hood
- `cpp/common/molecular/FitREQ_PN.h`: Core FitREQ_PN logic; `addAndReorderEpairs` sets `atoms->n0` including epairs; `evalEnergyDerivs` already decomposes Coulomb, vdW/Morse, Hcorr, SR (epair).
- [cpp/libs/Molecular/FitREQ_PN_lib.cpp](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/cpp/libs/Molecular/FitREQ_PN_lib.cpp:0:0-0:0): Added [getSampleGeom](cci:1://file:///home/prokophapala/git/FireCore-fitREQH/cpp/libs/Molecular/FitREQ_PN_lib.cpp:99:0-112:1) returning geometry + n0, and [evalSamplePairs](cci:1://file:///home/prokophapala/git/FireCore-fitREQH/pyBall/FitREQ_PN.py:184:0-187:19) for per-pair components; exposed via C API.
- [pyBall/FitREQ_PN.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/pyBall/FitREQ_PN.py:0:0-0:0): ctypes bindings updated so [getSampleGeom](cci:1://file:///home/prokophapala/git/FireCore-fitREQH/cpp/libs/Molecular/FitREQ_PN_lib.cpp:99:0-112:1) returns n0; wrappers for pairwise evaluation and geometry retrieval.
- [tests/tFitREQ_PN/visualize_static.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/visualize_static.py:0:0-0:0): Non-interactive visualizer with label modes, matrix plots, diagnostics, argparse (input, DOF, types, label_mode, Lepairs, kMorse). Default input `wb97m-split/H2O-A1_H2O-D1-z.xyz`, `Lepairs=0.8`.
- [tests/tFitREQ_PN/visualizer_gui.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/visualizer_gui.py:0:0-0:0): Interactive GUI with label mode switcher, anchored interaction matrix + click-to-focus, equal 3D aspect, smaller atom markers, argparse (same options), default `Lepairs=0.8`. Uses per-pair data to label lines in 3D (distance/Morse/Coulomb/Hcorr/SR).

## How the visualization works
- Data prep: load types/DOF selection, load XYZ (with epairs), get energies/components ([getEs](cci:1://file:///home/prokophapala/git/FireCore-fitREQH/pyBall/FitREQ_PN.py:160:0-167:23), [getEs_components](cci:1://file:///home/prokophapala/git/FireCore-fitREQH/cpp/libs/Molecular/FitREQ_PN_lib.cpp:264:0-281:1)), map samples to angle–distance grid via [parse_xyz_mapping](cci:1://file:///home/prokophapala/git/FireCore-fitREQH/pyBall/FitREQ_PN.py:361:0-428:45).
- Geometry: [getSampleGeom](cci:1://file:///home/prokophapala/git/FireCore-fitREQH/cpp/libs/Molecular/FitREQ_PN_lib.cpp:99:0-112:1) supplies positions, types, charges, host mapping, and `n0` fragment split (including epairs).
- Per-pair: [evalSamplePairs](cci:1://file:///home/prokophapala/git/FireCore-fitREQH/pyBall/FitREQ_PN.py:184:0-187:19) returns a tensor [natoms, natoms, 4] for Morse(vdW), Coulomb, Hcorr, SR.
- Label modes: `numeric`, `chem`, `chem_no_unders`, sourced from AtomTypes.

## CLI usage (both scripts)
Common options:
- `--input` (default `wb97m-split/H2O-A1_H2O-D1-z.xyz`)
- `--dof_selection` (default [data/dofSelection_run.dat](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/data/dofSelection_run.dat:0:0-0:0))
- `--atom_types` (default [data/AtomTypes.dat](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/data/AtomTypes.dat:0:0-0:0))
- `--element_types` (default `data/ElementTypes.dat`)
- `--label_mode` (`chem_no_unders`|`chem`|`numeric`, default `chem_no_unders`)
- `--lepairs` (default `0.8`)
- `--kmorse` (default `1.8`)

Run examples (from repo root):
```bash
cd tests/tFitREQ_PN

# Static (saves PNGs in cwd)
python3 visualize_static.py \
  --input wb97m-split/H2O-A1_H2O-D1-z.xyz \
  --lepairs 0.8 --kmorse 1.8 --label_mode chem_no_unders

# GUI (interactive)
python3 visualizer_gui.py \
  --input wb97m-split/H2O-A1_H2O-D1-z.xyz \
  --lepairs 0.8 --kmorse 1.8 --label_mode chem_no_unders
```

## Static visualizer ([visualize_static.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/visualize_static.py:0:0-0:0))
- Generates 2D maps (DFT, model, diff) with physical ticks.
- 1D slices for distance/angle (baseline, +Hcorr, +SR, DFT).
- 3D geometry snapshot with fragment split and labels.
- Interaction matrices for components.
- Prints diagnostics (minima, component tables, etc.).
- Saves PNGs: [plot_2d_maps.png](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/plot_2d_maps.png:0:0-0:0), `plot_1d_slices.png`, `plot_3d_geom.png`, `plot_mat_*.png`.

## GUI ([visualizer_gui.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/visualizer_gui.py:0:0-0:0))
- Panels: DFT/model/diff maps; distance/angle slices; 3D geometry; interaction matrix.
- Click maps to change point; click matrix to focus atom and highlight corresponding 3D lines.
- Radio: interaction mode (Distance, Morse, Coulomb, H-corr, SR); label mode selector.
- Buttons: focus atom prev/next; anchored colorbar; equal 3D aspect; small atom markers.

## Key defaults and behaviors
- Default test XYZ: `H2O-A1_H2O-D1-z.xyz`.
- Epairs added and placed at `Lepairs=0.8`.
- Fragment split uses returned `n0` (includes epairs).
- Pair labels in 3D respect chosen interaction mode; distances shown when in Distance mode.

## Git status (snapshot)
Modified/added: [cpp/libs/Molecular/FitREQ_PN_lib.cpp](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/cpp/libs/Molecular/FitREQ_PN_lib.cpp:0:0-0:0), [pyBall/FitREQ_PN.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/pyBall/FitREQ_PN.py:0:0-0:0), [tests/tFitREQ_PN/visualize_static.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/visualize_static.py:0:0-0:0), [tests/tFitREQ_PN/visualizer_gui.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/visualizer_gui.py:0:0-0:0), [doc/DevNotes/FitREQ_discussion.md](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/doc/DevNotes/FitREQ_discussion.md:0:0-0:0), [doc/DevNotes/FitREQ_visualizer_discussion.md](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/doc/DevNotes/FitREQ_visualizer_discussion.md:0:0-0:0), plus build/temp files and generated plots.

## Quick start (recommended)
1) `cd tests/tFitREQ_PN`
2) `python3 visualize_static.py` (inspect generated PNGs)
3) `python3 visualizer_gui.py` (click around; matrix click highlights 3D lines)

If you want the tutorial inserted into a docs file (e.g., `doc/DevNotes/FitREQ_visualizer_tutorial.md`), I can add it.