import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import RadioButtons, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable

sys.path.append("/home/prokophapala/git/FireCore-fitREQH/")
from pyBall import FitREQ_PN as fit

def load_atom_types(filepath):
    type_names = {}
    idx = 0
    try:
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or not line: continue
                parts = line.split()
                if len(parts) > 1:
                    type_names[idx] = parts[0]
                    idx += 1
    except Exception as e:
        print(f"Warning: could not load atom types from {filepath}: {e}")
    return type_names

def get_atom_labels(types, host, type_names, mode="chem"):
    labels = []
    for i, (t, h) in enumerate(zip(types, host)):
        if mode == "numeric":
            lbl = f"{i}T{t}"
            if h >= 0:
                lbl += f"(E{h})"
            labels.append(lbl)
        else:
            chem_name = type_names.get(t, f"T{t}")
            if mode == "chem_no_unders":
                chem_name = chem_name.replace("_", "")
            labels.append(chem_name)
    return labels

class PotentialDashboard:
    def __init__(self, angles, distances, seq, E_dft, E_base, E_hcorr, E_sr, isamp_map, type_names):
        self.angles = angles
        self.distances = distances
        self.seq = seq
        self.type_names = type_names
        
        self.E_dft = E_dft            
        self.E_base = E_base          
        self.E_hcorr = E_hcorr        
        self.E_sr = E_sr              
        self.E_model = self.E_base + self.E_hcorr + self.E_sr
        self.E_diff = self.E_model - self.E_dft
        
        self.isamp_map = isamp_map

        valid = ~np.isnan(self.E_dft)
        min_idx = np.unravel_index(np.nanargmin(self.E_dft), self.E_dft.shape)
        self.curr_d, self.curr_a = min_idx
        self.Emin = self.E_dft[self.curr_d, self.curr_a]
        
        isamp0 = self.isamp_map.get((self.curr_d, self.curr_a), 0)
        pos0, types0, _, host0, n0_sample = fit.getSampleGeom(isamp0)
        self.n_atoms = len(pos0)
        self.types0 = types0
        self.host0 = host0
        self.n0 = n0_sample
        
        self.label_style = 'chem_no_unders'
        self.labels = get_atom_labels(self.types0, self.host0, self.type_names, self.label_style)
        
        self.focus_atom = 0
        self.label_mode = 'Distance'
        
        self.geom_cache = {}
        self.pair_cache = {}
        
        self.setup_ui()
        self.update_all()

    def get_current_data(self):
        isamp = self.isamp_map.get((self.curr_d, self.curr_a), -1)
        if isamp == -1: return None, None
        
        if isamp not in self.geom_cache:
            pos, _, _, _, _ = fit.getSampleGeom(isamp)
            self.geom_cache[isamp] = pos
            pair_out = fit.evalSamplePairs(isamp, self.n_atoms)
            self.pair_cache[isamp] = pair_out
            
        return self.geom_cache[isamp], self.pair_cache[isamp]

    def setup_ui(self):
        self.fig = plt.figure(figsize=(20, 12))
        gs = GridSpec(3, 4, height_ratios=[1.2, 1.2, 1.2], width_ratios=[1, 1, 1.2, 0.4])
        
        self.ax_dft = self.fig.add_subplot(gs[0, 0])
        self.ax_model = self.fig.add_subplot(gs[0, 1])
        self.ax_diff = self.fig.add_subplot(gs[0, 2])
        
        vmin, vmax = self.Emin, -self.Emin
        diff_scale = 0.25
        vmin_diff, vmax_diff = diff_scale * self.Emin, -diff_scale * self.Emin
        
        self.im_dft = self.ax_dft.imshow(self.E_dft, origin='lower', aspect='auto', cmap='coolwarm', vmin=vmin, vmax=vmax)
        self.im_model = self.ax_model.imshow(self.E_model, origin='lower', aspect='auto', cmap='coolwarm', vmin=vmin, vmax=vmax)
        self.im_diff = self.ax_diff.imshow(self.E_diff, origin='lower', aspect='auto', cmap='bwr', vmin=vmin_diff, vmax=vmax_diff)
        
        for ax, title in zip([self.ax_dft, self.ax_model, self.ax_diff], ['DFT', 'Model', f'Difference [scale={diff_scale}]']):
            ax.set_title(title)
            xticks = np.linspace(0, len(self.angles)-1, 5, dtype=int)
            ax.set_xticks(xticks)
            ax.set_xticklabels([f"{self.angles[i]:.1f}" for i in xticks])
            yticks = np.linspace(0, len(self.distances)-1, 5, dtype=int)
            ax.set_yticks(yticks)
            ax.set_yticklabels([f"{self.distances[i]:.2f}" for i in yticks])
            ax.set_xlabel("Angle [deg]")
            ax.set_ylabel("Distance [A]")
            
        self.crosshairs = [ax.plot([], [], 'k+', markersize=15, markeredgewidth=2)[0] for ax in [self.ax_dft, self.ax_model, self.ax_diff]]

        self.ax_slice_dist = self.fig.add_subplot(gs[1, 0])
        self.ax_slice_ang = self.fig.add_subplot(gs[1, 1])
        self.ax_matrix = self.fig.add_subplot(gs[2, 0:2])
        self.ax_geom = self.fig.add_subplot(gs[1:, 2], projection='3d')

        # Dedicated colorbar axes for the matrix, anchored to the matrix via divider (no layout drift)
        divider = make_axes_locatable(self.ax_matrix)
        self.ax_cbar_matrix = divider.append_axes("right", size="3%", pad=0.2)
        self.cbar_matrix = None
        
        plt.tight_layout()

        self.ax_radio = self.fig.add_axes([0.86, 0.65, 0.12, 0.15], facecolor='lightgoldenrodyellow')
        self.ax_radio.set_title("Interaction Mode")
        self.radio = RadioButtons(self.ax_radio, ('Distance', 'Morse', 'Coulomb', 'H-corr', 'SR (Epairs)'))
        self.radio.on_clicked(self.on_mode_change)

        self.ax_radio_label = self.fig.add_axes([0.86, 0.45, 0.12, 0.12], facecolor='lightcyan')
        self.ax_radio_label.set_title("Label Mode")
        self.radio_label = RadioButtons(self.ax_radio_label, ('chem_no_unders', 'chem', 'numeric'))
        self.radio_label.on_clicked(self.on_label_mode_change)
        
        self.txt_atom = self.fig.text(0.89, 0.35, "", ha='center', fontsize=12, fontweight='bold')
        
        self.ax_btn_prev = self.fig.add_axes([0.85, 0.25, 0.05, 0.05])
        self.btn_prev = Button(self.ax_btn_prev, '< Prev')
        self.btn_prev.on_clicked(lambda e: self.change_atom(-1))
        
        self.ax_btn_next = self.fig.add_axes([0.91, 0.25, 0.05, 0.05])
        self.btn_next = Button(self.ax_btn_next, 'Next >')
        self.btn_next.on_clicked(lambda e: self.change_atom(1))

        self.fig.canvas.mpl_connect('button_press_event', self.onclick)

    def onclick(self, event):
        if event.inaxes in [self.ax_dft, self.ax_model, self.ax_diff]:
            x_idx = int(round(event.xdata))
            y_idx = int(round(event.ydata))
            self.curr_a = max(0, min(len(self.angles) - 1, x_idx))
            self.curr_d = max(0, min(len(self.distances) - 1, y_idx))
            self.update_all()
            self.fig.canvas.draw_idle()
        elif event.inaxes == self.ax_matrix:
            if event.ydata is not None:
                self.focus_atom = max(0, min(self.n_atoms - 1, int(round(event.ydata))))
                self.update_3d()
                self.update_matrix()
                self.fig.canvas.draw_idle()

    def on_mode_change(self, label):
        self.label_mode = label
        self.update_3d()
        self.update_matrix()
        self.fig.canvas.draw_idle()

    def on_label_mode_change(self, label):
        self.label_style = label
        self.labels = get_atom_labels(self.types0, self.host0, self.type_names, self.label_style)
        self.update_3d()
        self.update_matrix()
        self.fig.canvas.draw_idle()
        
    def change_atom(self, step):
        self.focus_atom = (self.focus_atom + step) % self.n_atoms
        self.update_3d()
        self.update_matrix()
        self.fig.canvas.draw_idle()

    def update_all(self):
        for ch in self.crosshairs:
            ch.set_data([self.curr_a], [self.curr_d])
        self.update_1d_slices()
        self.update_3d()
        self.update_matrix()

    def update_1d_slices(self):
        ang_val = self.angles[self.curr_a]
        dist_val = self.distances[self.curr_d]
        
        y_lims = [1.2 * self.Emin, -1.2 * self.Emin]
        dist_xlim = [min(self.distances), 6.0]
        
        self.ax_slice_dist.clear()
        self.ax_slice_dist.set_title(f"Distance Scan @ Angle = {ang_val:.1f}")
        self.ax_slice_dist.set_xlabel("True Distance [A]")
        self.ax_slice_dist.set_ylabel("Energy [kcal/mol]")
        
        b_d = self.E_base[:, self.curr_a]
        h_d = self.E_hcorr[:, self.curr_a]
        s_d = self.E_sr[:, self.curr_a]
        dft_d = self.E_dft[:, self.curr_a]
        
        self.ax_slice_dist.plot(self.distances, dft_d, 'k-', linewidth=2, label='DFT')
        self.ax_slice_dist.plot(self.distances, b_d, color='gray', linestyle=':', linewidth=1, label='Baseline')
        self.ax_slice_dist.fill_between(self.distances, b_d, b_d + h_d, color='cyan', alpha=0.5, label='+ Hcorr')
        self.ax_slice_dist.fill_between(self.distances, b_d + h_d, b_d + h_d + s_d, color='orange', alpha=0.5, label='+ SR (E-pair)')
        
        self.ax_slice_dist.axvline(dist_val, color='gray', linestyle=':')
        self.ax_slice_dist.set_ylim(y_lims)
        self.ax_slice_dist.set_xlim(dist_xlim)
        self.ax_slice_dist.axhline(0, color='k', linewidth=0.5)
        self.ax_slice_dist.legend(loc='lower right', fontsize='small')

        self.ax_slice_ang.clear()
        self.ax_slice_ang.set_title(f"Angular Scan @ Distance = {dist_val:.2f} A")
        self.ax_slice_ang.set_xlabel("True Angle [deg]")
        
        b_a = self.E_base[self.curr_d, :]
        h_a = self.E_hcorr[self.curr_d, :]
        s_a = self.E_sr[self.curr_d, :]
        dft_a = self.E_dft[self.curr_d, :]
        
        self.ax_slice_ang.plot(self.angles, dft_a, 'k-', linewidth=2, label='DFT')
        self.ax_slice_ang.plot(self.angles, b_a, color='gray', linestyle=':', linewidth=1, label='Baseline')
        self.ax_slice_ang.fill_between(self.angles, b_a, b_a + h_a, color='cyan', alpha=0.5)
        self.ax_slice_ang.fill_between(self.angles, b_a + h_a, b_a + h_a + s_a, color='orange', alpha=0.5)
        
        self.ax_slice_ang.axvline(ang_val, color='gray', linestyle=':')
        self.ax_slice_ang.set_ylim(y_lims)
        self.ax_slice_ang.axhline(0, color='k', linewidth=0.5)

    def update_3d(self):
        self.ax_geom.clear()
        geom, pair_data = self.get_current_data()
        if geom is None:
            self.ax_geom.set_title("No sample configuration available at this point")
            return
            
        frag1_idx = list(range(self.n0))
        frag2_idx = list(range(self.n0, self.n_atoms))
        
        f1 = geom[:self.n0]
        f2 = geom[self.n0:]
        self.ax_geom.scatter(f1[:,0], f1[:,1], f1[:,2], c='blue', s=30, alpha=0.8, marker='o', label='Fragment 1')
        self.ax_geom.scatter(f2[:,0], f2[:,1], f2[:,2], c='red', s=30, alpha=0.8, marker='o', label='Fragment 2')
        
        for i, lbl in enumerate(self.labels):
            self.ax_geom.text(geom[i,0], geom[i,1], geom[i,2], lbl, color='black', fontsize=9, zorder=10)
        
        focus_pos = geom[self.focus_atom]
        frag_name = "Frag 1" if self.focus_atom < self.n0 else "Frag 2"
        self.txt_atom.set_text(f"Focus: {self.labels[self.focus_atom]}\n({frag_name})")
        self.ax_geom.scatter(*focus_pos, c='gold', s=250, edgecolors='black', linewidth=2, zorder=5)

        target_atoms = frag2_idx if self.focus_atom < self.n0 else frag1_idx
        mode_idx_map = {'Morse': 0, 'Coulomb': 1, 'H-corr': 2, 'SR (Epairs)': 3}
        
        for tgt in target_atoms:
            tgt_pos = geom[tgt]
            self.ax_geom.plot([focus_pos[0], tgt_pos[0]], [focus_pos[1], tgt_pos[1]], [focus_pos[2], tgt_pos[2]], 'k--', alpha=0.7)
            if self.label_mode == 'Distance':
                val = np.linalg.norm(focus_pos - tgt_pos)
                txt = f"{val:.2f} A"
            else:
                comp_idx = mode_idx_map[self.label_mode]
                i, j = min(self.focus_atom, tgt), max(self.focus_atom, tgt)
                val = pair_data[i, j, comp_idx]
                if abs(val) < 1e-4: continue 
                txt = f"{val:+.2f}"
                
            mid = (focus_pos + tgt_pos) / 2
            self.ax_geom.text(mid[0], mid[1], mid[2], txt,  color='darkgreen', fontsize=10, fontweight='bold', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

        self.ax_geom.set_title(f"Geometry & Pair {self.label_mode}")
        # Enforce equal scaling on x,y,z
        mins = geom.min(axis=0)
        maxs = geom.max(axis=0)
        spans = maxs - mins
        span = spans.max()
        if span <= 1e-9:
            span = 1.0
        centers = (maxs + mins) / 2.0
        self.ax_geom.set_xlim(centers[0]-span/2, centers[0]+span/2)
        self.ax_geom.set_ylim(centers[1]-span/2, centers[1]+span/2)
        self.ax_geom.set_zlim(centers[2]-span/2, centers[2]+span/2)
        self.ax_geom.set_box_aspect([1,1,1])

    def update_matrix(self):
        self.ax_matrix.clear()
        geom, pair_data = self.get_current_data()
        if geom is None: return

        if self.label_mode == 'Distance':
            mat = np.zeros((self.n_atoms, self.n_atoms))
            for i in range(self.n_atoms):
                for j in range(self.n_atoms):
                    mat[i,j] = np.linalg.norm(geom[i] - geom[j])
            cmap = 'viridis'
            vmin, vmax = 0, np.max(mat)
        else:
            mode_idx_map = {'Morse': 0, 'Coulomb': 1, 'H-corr': 2, 'SR (Epairs)': 3}
            comp_idx = mode_idx_map[self.label_mode]
            mat = pair_data[:, :, comp_idx]
            cmap = 'bwr'
            vmax = np.max(np.abs(mat))
            if vmax < 1e-6: vmax = 1.0
            vmin = -vmax

        im = self.ax_matrix.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)
        
        # Update or create the matrix colorbar on the anchored axes
        if self.cbar_matrix is not None:
            self.cbar_matrix.update_normal(im)
        else:
            self.cbar_matrix = self.fig.colorbar(im, cax=self.ax_cbar_matrix)

        self.ax_matrix.set_xticks(np.arange(len(self.labels)))
        self.ax_matrix.set_yticks(np.arange(len(self.labels)))
        self.ax_matrix.set_xticklabels(self.labels, rotation=90, fontsize=8)
        self.ax_matrix.set_yticklabels(self.labels, fontsize=8)
        self.ax_matrix.set_title(f"{self.label_mode} Matrix")

        self.ax_matrix.axhline(self.focus_atom, color='gold', linewidth=2, alpha=0.5)
        self.ax_matrix.axvline(self.focus_atom, color='gold', linewidth=2, alpha=0.5)

def main():
    parser = argparse.ArgumentParser(description="Interactive visualization for FitREQ_PN surfaces")
    parser.add_argument("--input", type=str, default="wb97m-split/H2O-A1_H2O-D1-z.xyz", help="Input XYZ path (relative to script dir or absolute)")
    parser.add_argument("--dof_selection", type=str, default="data/dofSelection_run.dat", help="DOF selection file")
    parser.add_argument("--atom_types", type=str, default="data/AtomTypes.dat", help="AtomTypes.dat path")
    parser.add_argument("--element_types", type=str, default="data/ElementTypes.dat", help="ElementTypes.dat path")
    parser.add_argument("--label_mode", type=str, default="chem_no_unders", choices=["chem_no_unders", "chem", "numeric"], help="Initial label style")
    parser.add_argument("--lepairs", type=float, default=0.8, help="Epair distance from host")
    parser.add_argument("--kmorse", type=float, default=1.8, help="Morse curvature")
    args = parser.parse_args()

    here = os.path.dirname(__file__)
    def resolve(p):
        return p if os.path.isabs(p) else os.path.join(here, p)

    fname = resolve(args.input)
    dof_selection = resolve(args.dof_selection)
    atom_types_file = resolve(args.atom_types)
    element_types_file = resolve(args.element_types)

    type_names = load_atom_types(atom_types_file)
    
    fit.setVerbosity(verbosity=1)
    fit.setModel(ivdW=4, iCoul=1, iHbond=3, Epairs=1, iEpairs=2, kMorse=args.kmorse, Lepairs=args.lepairs, bPN=True)
    fit.loadTypes(fEtypes=element_types_file, fAtypes=atom_types_file)
    fit.loadDOFSelection(fname=dof_selection)
    
    n_total = fit.loadXYZ(fname, bAddEpairs=True, bOutXYZ=False, bEvalOnlyCorrections=False, bAppend=False)
    fit.getBuffs()
    
    Erefs, x0s = fit.read_xyz_data(fname)
    Eerr, Es_model, _ = fit.getEs(bDOFtoTypes=True, bEs=True)
    _, Es_Coul, Es_vdW, Es_Hcorr, Es_Epairs = fit.getEs_components()
    
    Es_base = Es_Coul + Es_vdW
    
    Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(fname)
    
    shape = Gref.shape
    Gmodel = np.full(shape, np.nan)
    Gbase = np.full(shape, np.nan)
    Ghcorr = np.full(shape, np.nan)
    Gsr = np.full(shape, np.nan)
    
    isamp_map = {}
    nmap = min(len(Es_model), len(seq))
    for k in range(nmap):
        idist, iang = seq[k]
        Gmodel[idist, iang] = Es_model[k]
        Gbase[idist, iang] = Es_base[k]
        Ghcorr[idist, iang] = Es_Hcorr[k]
        Gsr[idist, iang] = Es_Epairs[k]
        isamp_map[(idist, iang)] = k
        
    dashboard = PotentialDashboard(angles, distances, seq, Gref, Gbase, Ghcorr, Gsr, isamp_map, type_names)
    dashboard.label_style = args.label_mode
    dashboard.labels = get_atom_labels(dashboard.types0, dashboard.host0, dashboard.type_names, dashboard.label_style)
    plt.show()

if __name__ == "__main__":
    main()
