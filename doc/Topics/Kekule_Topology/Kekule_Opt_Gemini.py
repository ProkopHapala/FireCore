import numpy as np
import matplotlib.pyplot as plt
import argparse

class ZigzagRibbonSolver:
    def __init__(self, nx=10, ny=4, t0=2.7, alpha=4.1, K=21.0, r0=1.42):
        self.nx = nx   
        self.ny = ny   
        self.t0, self.alpha, self.K, self.r0 = t0, alpha, K, r0
        self.generate_zigzag_lattice()
        
    def generate_zigzag_lattice(self):
        """Generates a proper honeycomb zigzag ribbon with PBC in X."""
        coords = []
        for i in range(self.nx):
            for j in range(self.ny):
                # Standard zigzag geometry unit cell (4 atoms)
                x0 = i * 3
                y0 = j * np.sqrt(3)
                coords.append([x0,       y0])
                coords.append([x0 + 1.0, y0])
                coords.append([x0 + 1.5, y0 + np.sqrt(3)/2])
                coords.append([x0 + 2.5, y0 + np.sqrt(3)/2])
        
        self.coords = np.array(coords)
        self.num_atoms = len(coords)
        self.adj = []
        
        # Proper Honeycomb Connectivity with PBC
        cell_width = self.nx * 3
        for i in range(self.num_atoms):
            for j in range(i + 1, self.num_atoms):
                dx = abs(self.coords[i][0] - self.coords[j][0])
                if dx > cell_width / 2: dx = cell_width - dx # PBC
                dy = abs(self.coords[i][1] - self.coords[j][1])
                dist = np.sqrt(dx**2 + dy**2)
                
                if 0.9 < dist < 1.1:
                    self.adj.append((i, j))
        
        y_coords = self.coords[:, 1]
        self.bottom_edge = np.where(np.isclose(y_coords, np.min(y_coords)))[0]
        self.top_edge = np.where(np.isclose(y_coords, np.max(y_coords)))[0]

    def solve(self, edge_configs, max_iter=200):
        num_bonds = len(self.adj)
        r_ij = np.full(num_bonds, self.r0)
        
        onsite = np.zeros(self.num_atoms)
        for idx, cfg in edge_configs.items():
            if cfg['type'] == 'N': onsite[idx] = -1.5
            elif cfg['type'] == 'O': onsite[idx] = -2.5

        for _ in range(max_iter):
            H = np.diag(onsite)
            for idx, (i, j) in enumerate(self.adj):
                t = self.t0 - self.alpha * (r_ij[idx] - self.r0)
                # Apply pinning for quinoid-like bonds
                if (i in edge_configs and edge_configs[i].get('pinned')) or \
                   (j in edge_configs and edge_configs[j].get('pinned')):
                    t = 3.5 
                H[i, j] = H[j, i] = -t
            
            evals, evecs = np.linalg.eigh(H)
            mid = self.num_atoms // 2
            # Bond Order P_ij
            P_new = np.zeros(num_bonds)
            for idx, (i, j) in enumerate(self.adj):
                P_new[idx] = np.sum(evecs[i, :mid] * evecs[j, :mid])
            
            r_next = self.r0 - (2 * self.alpha / self.K) * P_new
            if np.linalg.norm(r_ij - r_next) < 1e-6: break
            r_ij = 0.8 * r_ij + 0.2 * r_next

        self.final_r = r_ij
        # LDOS at Fermi Level (average of states near mid-gap)
        self.ldos_fermi = np.abs(evecs[:, mid-1])**2 + np.abs(evecs[:, mid])**2
        self.plot_ribbon(edge_configs)

    def plot_ribbon(self, edge_configs):
        plt.figure(figsize=(14, 6))
        
        # 1. Draw Bonds (width mapped to bond order)
        for idx, (i, j) in enumerate(self.adj):
            dr = self.r0 - self.final_r[idx]
            w = 1.0 + 15 * dr # Thicker for double bonds
            color = 'black' if dr > 0.02 else 'gray'
            # Adjust for PBC: don't draw lines across the whole plot
            if abs(self.coords[i][0] - self.coords[j][0]) < 2.0:
                plt.plot([self.coords[i][0], self.coords[j][0]], 
                         [self.coords[i][1], self.coords[j][1]], 
                         color=color, linewidth=max(0.2, w), alpha=0.5, zorder=1)
        
        # 2. Draw Bulk Atoms (colored by LDOS)
        bulk_indices = [i for i in range(self.num_atoms) if i not in edge_configs]
        sc = plt.scatter(self.coords[bulk_indices, 0], self.coords[bulk_indices, 1], 
                         c=self.ldos_fermi[bulk_indices], cmap='viridis', 
                         s=100, edgecolors='none', zorder=2, label='Carbon (LDOS)')
        plt.colorbar(sc, label='LDOS @ Fermi Level (Conductivity)')
        
        # 3. Draw Heteroatoms (Discrete colors)
        for idx, cfg in edge_configs.items():
            color = 'red' if cfg['type'] == 'O' else 'blue'
            plt.scatter(self.coords[idx, 0], self.coords[idx, 1], 
                        color=color, s=150, edgecolors='black', 
                        linewidth=2, zorder=3, label=cfg['type'] if cfg['type'] not in plt.gca().get_legend_handles_labels()[1] else "")

        plt.title(f"Zigzag Ribbon: Peierls/Kekule Optimization\n(Red=O, Blue=N, PBC in X)")
        plt.axis('equal')
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nx", type=int, default=12)
    parser.add_argument("--ny", type=int, default=4)
    parser.add_argument("--K", type=float, default=21.0)
    parser.add_argument("--alpha", type=float, default=4.1)
    args = parser.parse_args()

    solver = ZigzagRibbonSolver(nx=args.nx, ny=args.ny, K=args.K, alpha=args.alpha)
    
    # Define your edge atoms here. 
    # Example: Oxygen at start of bottom edge, Nitrogen at end of bottom edge.
    configs = {
        solver.bottom_edge[0]:  {'type': 'O', 'pinned': True},
        solver.bottom_edge[-1]: {'type': 'N', 'pinned': False}
    }
    
    solver.solve(configs)