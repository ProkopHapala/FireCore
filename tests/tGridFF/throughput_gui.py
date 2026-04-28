#!/usr/bin/env python3
"""
GUI interface for Throughput MD simulation parameters
Enhanced version with improved graphics and parameter loading
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import subprocess
import os
import threading
import re

class ThroughputGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("🧬 Throughput MD Simulation Suite")
        self.root.geometry("900x1000")
        self.root.configure(bg='#f0f0f0')

        # Set modern style
        style = ttk.Style()
        style.theme_use('clam')

        # Configure custom styles
        style.configure('Title.TLabel', font=('Arial', 14, 'bold'), foreground='#2c3e50')
        style.configure('Heading.TLabel', font=('Arial', 12, 'bold'), foreground='#34495e')
        style.configure('Info.TLabel', font=('Arial', 9), foreground='#7f8c8d')
        style.configure('Success.TLabel', font=('Arial', 10), foreground='#27ae60')
        style.configure('Error.TLabel', font=('Arial', 10), foreground='#e74c3c')

        # Control variables
        self.stop_requested = False
        
        # Create header
        header_frame = tk.Frame(root, bg='#2c3e50', height=60)
        header_frame.pack(fill=tk.X)
        header_frame.pack_propagate(False)

        title_label = tk.Label(header_frame, text="🧬 Molecular Dynamics Simulation Suite",
                              font=('Arial', 16, 'bold'), fg='white', bg='#2c3e50')
        title_label.pack(expand=True)

        # Create main frame with scrollbar
        main_frame = ttk.Frame(root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=10)

        # Create notebook for tabs
        notebook = ttk.Notebook(main_frame)
        notebook.pack(fill=tk.BOTH, expand=True)
        
        # Basic Parameters Tab
        basic_frame = ttk.Frame(notebook)
        notebook.add(basic_frame, text="Basic Parameters")
        self.create_basic_params(basic_frame)
        
        # Advanced Parameters Tab
        advanced_frame = ttk.Frame(notebook)
        notebook.add(advanced_frame, text="Advanced Parameters")
        self.create_advanced_params(advanced_frame)
        
        # Local Memory Parameters Tab
        memory_frame = ttk.Frame(notebook)
        notebook.add(memory_frame, text="Local Memory")
        self.create_memory_params(memory_frame)

        # Visualization Tab
        viz_frame = ttk.Frame(notebook)
        notebook.add(viz_frame, text="Visualization")
        self.create_visualization_tab(viz_frame)
        
        # Control buttons with improved styling
        control_frame = ttk.LabelFrame(main_frame, text="🎮 Control Panel")
        control_frame.pack(fill=tk.X, pady=(10, 0))

        # Main action buttons
        action_frame = ttk.Frame(control_frame)
        action_frame.pack(fill=tk.X, padx=10, pady=5)

        self.run_button = ttk.Button(action_frame, text="🚀 Run Optimization", command=self.run_simulation)
        self.run_button.pack(side=tk.LEFT, padx=(0, 10))

        self.stop_button = ttk.Button(action_frame, text="⏹️ Stop", command=self.stop_simulation, state=tk.DISABLED)
        self.stop_button.pack(side=tk.LEFT, padx=(0, 10))

        # Parameter management buttons
        param_frame = ttk.Frame(control_frame)
        param_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Button(param_frame, text="💾 Save Parameters", command=self.save_parameters).pack(side=tk.LEFT, padx=(0, 10))

        # Load parameters with dropdown
        load_frame = ttk.Frame(param_frame)
        load_frame.pack(side=tk.LEFT, padx=(0, 10))

        ttk.Button(load_frame, text="📂 Load Custom", command=self.load_parameters).pack(side=tk.LEFT)

        # Progress bar
        progress_frame = ttk.Frame(main_frame)
        progress_frame.pack(fill=tk.X, pady=(5, 0))

        ttk.Label(progress_frame, text="Progress:").pack(side=tk.LEFT)
        self.progress_var = tk.StringVar(value="Ready")
        ttk.Label(progress_frame, textvariable=self.progress_var).pack(side=tk.LEFT, padx=(10, 0))

        self.progress_bar = ttk.Progressbar(progress_frame, mode='determinate')
        self.progress_bar.pack(side=tk.RIGHT, fill=tk.X, expand=True, padx=(10, 0))
        
        # Output text area
        output_frame = ttk.LabelFrame(main_frame, text="Output")
        output_frame.pack(fill=tk.BOTH, expand=True, pady=(10, 0))
        
        self.output_text = tk.Text(output_frame, height=10)
        scrollbar = ttk.Scrollbar(output_frame, orient=tk.VERTICAL, command=self.output_text.yview)
        self.output_text.configure(yscrollcommand=scrollbar.set)
        
        self.output_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
    def create_basic_params(self, parent):
        # Basic simulation parameters with improved styling
        ttk.Label(parent, text="⚙️ Basic Simulation Parameters", style='Title.TLabel').pack(anchor=tk.W, pady=(10, 15))
        
        # Force field flags with improved layout
        flags_frame = ttk.LabelFrame(parent, text="🔬 Force Field Configuration")
        flags_frame.pack(fill=tk.X, pady=(0, 15), padx=5)

        # Create two columns for better organization
        flags_left = ttk.Frame(flags_frame)
        flags_left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=5)

        flags_right = ttk.Frame(flags_frame)
        flags_right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=5)

        # Left column
        self.dovdW = tk.IntVar(value=1)
        ttk.Checkbutton(flags_left, text="✓ Van der Waals (dovdW)", variable=self.dovdW).pack(anchor=tk.W, pady=2)

        self.doSurfAtoms = tk.IntVar(value=1)
        ttk.Checkbutton(flags_left, text="✓ Surface Atoms (doSurfAtoms)", variable=self.doSurfAtoms).pack(anchor=tk.W, pady=2)

        # Right column - Grid FF selection
        ttk.Label(flags_right, text="Grid Force Field:", style='Heading.TLabel').pack(anchor=tk.W)
        grid_frame = ttk.Frame(flags_right)
        grid_frame.pack(fill=tk.X, pady=5)

        self.bGridFF = tk.IntVar(value=6)
        ttk.Radiobutton(grid_frame, text="🔲 Disabled (0)", variable=self.bGridFF, value=0).pack(anchor=tk.W)
        ttk.Radiobutton(grid_frame, text="📈 Linear (1)", variable=self.bGridFF, value=1).pack(anchor=tk.W)
        ttk.Radiobutton(grid_frame, text="🎯 bSpline (6)", variable=self.bGridFF, value=6).pack(anchor=tk.W)

        self.bTex = tk.IntVar(value=0)
        ttk.Checkbutton(flags_right, text="🎨 Texture (bTex)", variable=self.bTex).pack(anchor=tk.W, pady=2)

        self.bSaveToDatabase = tk.IntVar(value=-1)
        ttk.Checkbutton(flags_right, text="💾 Save to Database (bSaveToDatabase)", variable=self.bSaveToDatabase).pack(anchor=tk.W, pady=2)
        
        # File paths
        files_frame = ttk.LabelFrame(parent, text="File Paths")
        files_frame.pack(fill=tk.X, pady=(0, 10))
        
        # XYZ file
        xyz_frame = ttk.Frame(files_frame)
        xyz_frame.pack(fill=tk.X, pady=2)
        ttk.Label(xyz_frame, text="XYZ File:").pack(side=tk.LEFT)
        self.xyz_name = tk.StringVar(value="data/xyz/xylitol_WO_gridFF")
        ttk.Entry(xyz_frame, textvariable=self.xyz_name, width=40).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(xyz_frame, text="Browse", command=self.browse_xyz).pack(side=tk.RIGHT)

        # Surface file
        surf_frame = ttk.Frame(files_frame)
        surf_frame.pack(fill=tk.X, pady=2)
        ttk.Label(surf_frame, text="Surface File:").pack(side=tk.LEFT)
        self.surf_name = tk.StringVar(value="data/xyz/surfaces_for_throughput/NaCl_${N}x${N}_L3")
        ttk.Entry(surf_frame, textvariable=self.surf_name, width=40).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(surf_frame, text="Browse", command=self.browse_surface).pack(side=tk.RIGHT)

        # Help text for surface file placeholders
        help_frame = ttk.Frame(files_frame)
        help_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Label(help_frame, text="💡 Use ${N} as placeholder for surface size (e.g., NaCl_${N}x${N}_L3)",
                 font=("Arial", 8), foreground="gray").pack(side=tk.LEFT)
        
        # Convergence criteria
        conv_frame = ttk.LabelFrame(parent, text="Convergence")
        conv_frame.pack(fill=tk.X, pady=(0, 10))
        
        fconv_frame = ttk.Frame(conv_frame)
        fconv_frame.pack(fill=tk.X)
        ttk.Label(fconv_frame, text="Force Convergence:").pack(side=tk.LEFT)
        self.Fconv = tk.StringVar(value="1e-4")
        ttk.Entry(fconv_frame, textvariable=self.Fconv, width=15).pack(side=tk.LEFT, padx=5)
        
    def create_advanced_params(self, parent):
        ttk.Label(parent, text="Advanced Parameters", font=("Arial", 12, "bold")).pack(anchor=tk.W, pady=(0, 10))
        
        # Simulation parameters arrays
        sim_frame = ttk.LabelFrame(parent, text="Simulation Parameters")
        sim_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Replicas
        rep_frame = ttk.Frame(sim_frame)
        rep_frame.pack(fill=tk.X, pady=2)
        ttk.Label(rep_frame, text="Replicas:").pack(side=tk.LEFT)
        self.replicas = tk.StringVar(value="1000,5000")
        ttk.Entry(rep_frame, textvariable=self.replicas, width=30).pack(side=tk.LEFT, padx=5)
        ttk.Label(rep_frame, text="(comma-separated)").pack(side=tk.LEFT)
        
        # Per frames
        pf_frame = ttk.Frame(sim_frame)
        pf_frame.pack(fill=tk.X, pady=2)
        ttk.Label(pf_frame, text="Per Frames:").pack(side=tk.LEFT)
        self.perframes = tk.StringVar(value="20,500")
        ttk.Entry(pf_frame, textvariable=self.perframes, width=30).pack(side=tk.LEFT, padx=5)
        ttk.Label(pf_frame, text="(comma-separated)").pack(side=tk.LEFT)
        
        # Per VF
        pvf_frame = ttk.Frame(sim_frame)
        pvf_frame.pack(fill=tk.X, pady=2)
        ttk.Label(pvf_frame, text="Per VF:").pack(side=tk.LEFT)
        self.perVF = tk.StringVar(value="20,50")
        ttk.Entry(pvf_frame, textvariable=self.perVF, width=30).pack(side=tk.LEFT, padx=5)
        ttk.Label(pvf_frame, text="(comma-separated)").pack(side=tk.LEFT)
        
        # PBC parameters
        pbc_frame = ttk.LabelFrame(parent, text="PBC Parameters")
        pbc_frame.pack(fill=tk.X, pady=(0, 10))
        
        npbc_frame = ttk.Frame(pbc_frame)
        npbc_frame.pack(fill=tk.X, pady=2)
        ttk.Label(npbc_frame, text="nPBC:").pack(side=tk.LEFT)
        self.nPBC = tk.StringVar(value="(1,1,0),(2,2,0),(3,3,0),(0,0,0)")
        ttk.Entry(npbc_frame, textvariable=self.nPBC, width=30).pack(side=tk.LEFT, padx=5)
        ttk.Label(npbc_frame, text="(comma-separated, e.g., (1,1,0),(2,2,0))").pack(side=tk.LEFT)
        
        # Surface parameters with range support
        surf_frame = ttk.LabelFrame(parent, text="🏔️ Surface Parameters")
        surf_frame.pack(fill=tk.X, pady=(0, 15), padx=5)

        ns_frame = ttk.Frame(surf_frame)
        ns_frame.pack(fill=tk.X, pady=5, padx=10)
        ttk.Label(ns_frame, text="Surface Sizes (N):", style='Heading.TLabel').pack(side=tk.LEFT)
        self.Ns = tk.StringVar(value="1-16")
        ttk.Entry(ns_frame, textvariable=self.Ns, width=25).pack(side=tk.LEFT, padx=5)

        # Help text for surface ranges
        help_frame = ttk.Frame(surf_frame)
        help_frame.pack(fill=tk.X, padx=10, pady=(0, 5))
        help_text = "Examples: '1-16' (range), '1,4,8,16' (specific), '16' (single)"
        ttk.Label(help_frame, text=help_text, style='Info.TLabel').pack(side=tk.LEFT)
        
    def create_memory_params(self, parent):
        ttk.Label(parent, text="Local Memory Parameters", font=("Arial", 12, "bold")).pack(anchor=tk.W, pady=(0, 10))
        
        mem_frame = ttk.LabelFrame(parent, text="Local Memory Settings")
        mem_frame.pack(fill=tk.X, pady=(0, 10))
        
        # nlocMMFFs
        mmff_frame = ttk.Frame(mem_frame)
        mmff_frame.pack(fill=tk.X, pady=2)
        ttk.Label(mmff_frame, text="nlocMMFFs:").pack(side=tk.LEFT)
        self.nlocMMFFs = tk.StringVar(value="32")
        ttk.Entry(mmff_frame, textvariable=self.nlocMMFFs, width=30).pack(side=tk.LEFT, padx=5)
        
        # nlocmoves
        moves_frame = ttk.Frame(mem_frame)
        moves_frame.pack(fill=tk.X, pady=2)
        ttk.Label(moves_frame, text="nlocmoves:").pack(side=tk.LEFT)
        self.nlocmoves = tk.StringVar(value="32")
        ttk.Entry(moves_frame, textvariable=self.nlocmoves, width=30).pack(side=tk.LEFT, padx=5)
        
        # nlocNBFFs
        nbff_frame = ttk.Frame(mem_frame)
        nbff_frame.pack(fill=tk.X, pady=2)
        ttk.Label(nbff_frame, text="nlocNBFFs:").pack(side=tk.LEFT)
        self.nlocNBFFs = tk.StringVar(value="--")
        ttk.Entry(nbff_frame, textvariable=self.nlocNBFFs, width=30).pack(side=tk.LEFT, padx=5)
        
        # nlocSurfs
        surfs_frame = ttk.Frame(mem_frame)
        surfs_frame.pack(fill=tk.X, pady=2)
        ttk.Label(surfs_frame, text="nlocSurfs:").pack(side=tk.LEFT)
        self.nlocSurfs = tk.StringVar(value="--")
        ttk.Entry(surfs_frame, textvariable=self.nlocSurfs, width=30).pack(side=tk.LEFT, padx=5)
        
        # nlocGridFFs
        gridff_frame = ttk.Frame(mem_frame)
        gridff_frame.pack(fill=tk.X, pady=2)
        ttk.Label(gridff_frame, text="nlocGridFFs:").pack(side=tk.LEFT)
        self.nlocGridFFs = tk.StringVar(value="--")
        ttk.Entry(gridff_frame, textvariable=self.nlocGridFFs, width=30).pack(side=tk.LEFT, padx=5)
        
        # nlocGridFFbSplines
        spline_frame = ttk.Frame(mem_frame)
        spline_frame.pack(fill=tk.X, pady=2)
        ttk.Label(spline_frame, text="nlocGridFFbSplines:").pack(side=tk.LEFT)
        self.nlocGridFFbSplines = tk.StringVar(value="--")
        ttk.Entry(spline_frame, textvariable=self.nlocGridFFbSplines, width=30).pack(side=tk.LEFT, padx=5)

    def create_visualization_tab(self, parent):
        ttk.Label(parent, text="Data Visualization", font=("Arial", 12, "bold")).pack(anchor=tk.W, pady=(0, 10))

        # Results file selection
        file_frame = ttk.LabelFrame(parent, text="Results File")
        file_frame.pack(fill=tk.X, pady=(0, 10))

        file_select_frame = ttk.Frame(file_frame)
        file_select_frame.pack(fill=tk.X, pady=5)
        ttk.Label(file_select_frame, text="Results file:").pack(side=tk.LEFT)
        self.results_file = tk.StringVar(value="results/all/results.dat")
        ttk.Entry(file_select_frame, textvariable=self.results_file, width=40).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(file_select_frame, text="Browse", command=self.browse_results_file).pack(side=tk.RIGHT)

        # Plot options
        options_frame = ttk.LabelFrame(parent, text="Plot Options")
        options_frame.pack(fill=tk.X, pady=(0, 10))

        # X-axis mode
        xaxis_frame = ttk.Frame(options_frame)
        xaxis_frame.pack(fill=tk.X, pady=5)
        ttk.Label(xaxis_frame, text="X-axis mode:").pack(side=tk.LEFT)
        self.use_nxn = tk.BooleanVar(value=False)
        ttk.Radiobutton(xaxis_frame, text="Atom Count", variable=self.use_nxn, value=False).pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(xaxis_frame, text="NxN Notation", variable=self.use_nxn, value=True).pack(side=tk.LEFT, padx=10)

        # Auto-open browser
        browser_frame = ttk.Frame(options_frame)
        browser_frame.pack(fill=tk.X, pady=5)
        self.auto_open_browser = tk.BooleanVar(value=True)
        ttk.Checkbutton(browser_frame, text="Auto-open plot in browser", variable=self.auto_open_browser).pack(side=tk.LEFT)

        # Plot generation buttons
        buttons_frame = ttk.Frame(parent)
        buttons_frame.pack(fill=tk.X, pady=(10, 0))

        ttk.Button(buttons_frame, text="Generate Interactive Plot", command=self.generate_plot).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(buttons_frame, text="View Results Table", command=self.view_results_table).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(buttons_frame, text="Open Results Folder", command=self.open_results_folder).pack(side=tk.LEFT)

        # Status display
        status_frame = ttk.LabelFrame(parent, text="Status")
        status_frame.pack(fill=tk.BOTH, expand=True, pady=(10, 0))

        self.viz_status = tk.Text(status_frame, height=8, wrap=tk.WORD)
        viz_scrollbar = ttk.Scrollbar(status_frame, orient=tk.VERTICAL, command=self.viz_status.yview)
        self.viz_status.configure(yscrollcommand=viz_scrollbar.set)

        self.viz_status.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        viz_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
    def browse_xyz(self):
        filename = filedialog.askopenfilename(
            title="Select XYZ file",
            filetypes=[("XYZ files", "*.xyz"), ("All files", "*.*")]
        )
        if filename:
            self.xyz_name.set(filename)

    def browse_surface(self):
        filename = filedialog.askopenfilename(
            title="Select Surface file",
            filetypes=[("Surface files", "*.xyz"), ("All files", "*.*")]
        )
        if filename:
            self.surf_name.set(filename)

    def browse_results_file(self):
        filename = filedialog.askopenfilename(
            title="Select results file",
            filetypes=[("Data files", "*.dat"), ("All files", "*.*")]
        )
        if filename:
            self.results_file.set(filename)

    def generate_plot(self):
        """Generate interactive plot using plot_data.py"""
        def plot_thread():
            try:
                self.viz_status.delete(1.0, tk.END)
                self.viz_status.insert(tk.END, "Generating plot...\n")
                self.viz_status.update()

                # Check if results file exists
                results_file = self.results_file.get()
                if not os.path.exists(results_file):
                    self.viz_status.insert(tk.END, f"Error: Results file '{results_file}' not found!\n")
                    self.viz_status.insert(tk.END, "Please run optimization first or select existing results file.\n")
                    return

                # Copy results file to expected location for plot_data.py
                import shutil
                shutil.copy(results_file, "results.dat")

                # Build command for plot_data.py
                cmd = ["python3", "plot_data.py"]
                if self.use_nxn.get():
                    cmd.append("--nxn")

                self.viz_status.insert(tk.END, f"Command: {' '.join(cmd)}\n")
                self.viz_status.update()

                # Run plot generation
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")

                if result.returncode == 0:
                    self.viz_status.insert(tk.END, "Plot generated successfully!\n")
                    self.viz_status.insert(tk.END, result.stdout)

                    # Auto-open browser if requested
                    if self.auto_open_browser.get():
                        mode_label = "NxN" if self.use_nxn.get() else "AtomCount"
                        html_file = f"interactive_performance_plot_{mode_label}.html"
                        if os.path.exists(html_file):
                            import webbrowser
                            webbrowser.open(f"file://{os.path.abspath(html_file)}")
                            self.viz_status.insert(tk.END, f"Opened {html_file} in browser\n")
                else:
                    self.viz_status.insert(tk.END, "Error generating plot:\n")
                    self.viz_status.insert(tk.END, result.stderr)

            except Exception as e:
                self.viz_status.insert(tk.END, f"Error: {str(e)}\n")
                import traceback
                self.viz_status.insert(tk.END, traceback.format_exc())

            self.viz_status.see(tk.END)

        threading.Thread(target=plot_thread, daemon=True).start()

    def view_results_table(self):
        """Open results table CSV if it exists"""
        try:
            csv_file = "results_table.csv"
            if os.path.exists(csv_file):
                import webbrowser
                webbrowser.open(f"file://{os.path.abspath(csv_file)}")
                self.viz_status.delete(1.0, tk.END)
                self.viz_status.insert(tk.END, f"Opened {csv_file} in default application\n")
            else:
                self.viz_status.delete(1.0, tk.END)
                self.viz_status.insert(tk.END, "Results table not found. Generate plot first.\n")
        except Exception as e:
            self.viz_status.delete(1.0, tk.END)
            self.viz_status.insert(tk.END, f"Error opening results table: {str(e)}\n")

    def open_results_folder(self):
        """Open results folder in file manager"""
        try:
            results_dir = "results"
            if os.path.exists(results_dir):
                import subprocess
                import platform

                if platform.system() == "Windows":
                    subprocess.run(["explorer", results_dir])
                elif platform.system() == "Darwin":  # macOS
                    subprocess.run(["open", results_dir])
                else:  # Linux
                    subprocess.run(["xdg-open", results_dir])

                self.viz_status.delete(1.0, tk.END)
                self.viz_status.insert(tk.END, f"Opened {results_dir} folder\n")
            else:
                self.viz_status.delete(1.0, tk.END)
                self.viz_status.insert(tk.END, "Results folder not found. Run optimization first.\n")
        except Exception as e:
            self.viz_status.delete(1.0, tk.END)
            self.viz_status.insert(tk.END, f"Error opening results folder: {str(e)}\n")
    
    def run_simulation(self):
        """Run the simulation in a separate thread"""
        def run_thread():
            self.stop_requested = False
            self.run_button.config(state=tk.DISABLED)
            self.stop_button.config(state=tk.NORMAL)
            self.progress_var.set("Initializing...")
            self.progress_bar['value'] = 0

            self.output_text.delete(1.0, tk.END)
            self.output_text.insert(tk.END, "Starting parameter optimization...\n")
            self.output_text.update()

            try:
                # Parse parameter arrays
                replicas = [int(x.strip()) for x in self.replicas.get().split(',')]
                perframes = [int(x.strip()) for x in self.perframes.get().split(',')]
                perVFs = [int(x.strip()) for x in self.perVF.get().split(',')]

                # Parse surface sizes with range support
                surface_names = self.parse_surface_range(self.Ns.get())
                Ns = []
                for surf_name in surface_names:
                    # Extract N from "NaCl_NxN" format
                    if 'x' in surf_name:
                        n_str = surf_name.split('_')[1].split('x')[0]
                        Ns.append(int(n_str))
                    else:
                        Ns.append(1)  # fallback

                # Parse local memory parameters
                nlocMMFFs = [x.strip() for x in self.nlocMMFFs.get().split(',')]
                nlocmoves = [x.strip() for x in self.nlocmoves.get().split(',')]
                nlocNBFFs = [x.strip() for x in self.nlocNBFFs.get().split(',')]
                nlocSurfs = [x.strip() for x in self.nlocSurfs.get().split(',')]
                nlocGridFFs = [x.strip() for x in self.nlocGridFFs.get().split(',')]
                nlocGridFFbSplines = [x.strip() for x in self.nlocGridFFbSplines.get().split(',')]

                # Parse nPBC values correctly - handle parentheses
                npbc_input = self.nPBC.get().strip()
                nPBCs = []
                if npbc_input:
                    # Split by "),(" to handle multiple nPBC values like "(1,1,0),(2,2,0)"
                    if '),(' in npbc_input:
                        parts = npbc_input.split('),(')
                        for i, part in enumerate(parts):
                            if i == 0:
                                part = part + ')'  # Add closing parenthesis to first part
                            elif i == len(parts) - 1:
                                part = '(' + part  # Add opening parenthesis to last part
                            else:
                                part = '(' + part + ')'  # Add both parentheses to middle parts
                            nPBCs.append(part.strip())
                    else:
                        # Single nPBC value
                        nPBCs.append(npbc_input)

                # Initialize best result tracking
                best_value = 0
                best_params = {}
                total_combinations = 0
                current_combination = 0

                # Calculate total combinations for progress tracking
                total_combinations = (len(nPBCs) * len(nlocMMFFs) * len(nlocmoves) *
                                    len(nlocNBFFs) * len(nlocSurfs) * len(nlocGridFFs) *
                                    len(nlocGridFFbSplines) * len(replicas) * len(perframes) *
                                    len(perVFs) * len(Ns))

                self.output_text.insert(tk.END, f"Total combinations to test: {total_combinations}\n\n")
                self.output_text.update()

                # Create results directories
                os.makedirs("results/best", exist_ok=True)
                os.makedirs("results/all", exist_ok=True)

                # Main optimization loops (following the bash script structure)
                for nPBC in nPBCs:
                    for nlocMMFF in nlocMMFFs:
                        self.set_nloc_param("nlocMMFF", nlocMMFF)
                        for nlocmove in nlocmoves:
                            self.set_nloc_param("nlocmove", nlocmove)
                            for nlocNBFF in nlocNBFFs:
                                if nlocNBFF != "--":
                                    self.set_nloc_param("nlocNBFF", nlocNBFF)
                                for nlocSurf in nlocSurfs:
                                    if nlocSurf != "--":
                                        self.set_nloc_param("nlocSurf", nlocSurf)
                                    for nlocGridFF in nlocGridFFs:
                                        if nlocGridFF != "--":
                                            self.set_nloc_param("nlocGridFF", nlocGridFF)
                                        for nlocGridFFbSpline in nlocGridFFbSplines:
                                            if nlocGridFFbSpline != "--":
                                                self.set_nloc_param("nlocGridFFbSpline", nlocGridFFbSpline)
                                            for nSys in replicas:
                                                for perframe in perframes:
                                                    for pvf in perVFs:
                                                        if pvf > perframe:
                                                            continue

                                                        if self.stop_requested:
                                                            self.output_text.insert(tk.END, "Optimization stopped by user.\n")
                                                            break

                                                        current_combination += 1
                                                        progress = (current_combination / total_combinations) * 100

                                                        self.progress_var.set(f"{current_combination}/{total_combinations} ({progress:.1f}%)")
                                                        self.progress_bar['value'] = progress

                                                        self.output_text.insert(tk.END,
                                                            f"Progress: {current_combination}/{total_combinations} ({progress:.1f}%)\n")
                                                        self.output_text.insert(tk.END,
                                                            f"Testing: nSys={nSys}, perframe={perframe}, perVF={pvf}\n")
                                                        self.output_text.update()

                                                        # Run simulations for each surface size
                                                        for N in Ns:
                                                            result_value = self.run_single_simulation(
                                                                nSys, perframe, pvf, N, nPBC,
                                                                nlocMMFF, nlocmove, nlocNBFF, nlocSurf,
                                                                nlocGridFF, nlocGridFFbSpline
                                                            )

                                                            if result_value and result_value > best_value:
                                                                best_value = result_value
                                                                best_params = {
                                                                    'nSys': nSys, 'perframe': perframe, 'perVF': pvf,
                                                                    'N': N, 'nPBC': nPBC, 'nlocMMFF': nlocMMFF,
                                                                    'nlocmove': nlocmove, 'nlocNBFF': nlocNBFF,
                                                                    'nlocSurf': nlocSurf, 'nlocGridFF': nlocGridFF,
                                                                    'nlocGridFFbSpline': nlocGridFFbSpline
                                                                }
                                                                self.output_text.insert(tk.END,
                                                                    f"NEW BEST: {best_value} with params: {best_params}\n")
                                                                self.output_text.update()

                                                        if self.stop_requested:
                                                            break
                                                    if self.stop_requested:
                                                        break
                                                if self.stop_requested:
                                                    break
                                            if self.stop_requested:
                                                break
                                        if self.stop_requested:
                                            break
                                    if self.stop_requested:
                                        break
                                if self.stop_requested:
                                    break
                            if self.stop_requested:
                                break
                        if self.stop_requested:
                            break
                    if self.stop_requested:
                        break

                # Final results
                if not self.stop_requested:
                    self.output_text.insert(tk.END, "\n" + "="*50 + "\n")
                    self.output_text.insert(tk.END, "OPTIMIZATION COMPLETE!\n")
                    self.output_text.insert(tk.END, f"Best value: {best_value}\n")
                    self.output_text.insert(tk.END, f"Best parameters: {best_params}\n")
                    self.output_text.insert(tk.END, "="*50 + "\n")
                    self.progress_var.set("Complete!")
                else:
                    self.progress_var.set("Stopped")

            except Exception as e:
                self.output_text.insert(tk.END, f"Error: {str(e)}\n")
                import traceback
                self.output_text.insert(tk.END, traceback.format_exc())
                self.progress_var.set("Error")

            finally:
                # Re-enable buttons
                self.run_button.config(state=tk.NORMAL)
                self.stop_button.config(state=tk.DISABLED)

            self.output_text.see(tk.END)

        threading.Thread(target=run_thread, daemon=True).start()

    def stop_simulation(self):
        """Stop the running simulation"""
        self.stop_requested = True
        self.progress_var.set("Stopping...")

    def set_nloc_param(self, param_name, param_value):
        """Set local memory parameter using the shell script"""
        try:
            if param_value != "--":
                subprocess.run(["./set_nloc_params.sh", param_name, param_value],
                             check=True, cwd=".")
        except subprocess.CalledProcessError as e:
            self.output_text.insert(tk.END, f"Warning: Failed to set {param_name}={param_value}: {e}\n")

    def run_single_simulation(self, nSys, perframe, pvf, N, nPBC,
                             nlocMMFF, nlocmove, nlocNBFF, nlocSurf,
                             nlocGridFF, nlocGridFFbSpline):
        """Run a single simulation with given parameters"""
        try:
            # Clear any existing minima.dat file before running simulation
            if os.path.exists("minima.dat"):
                os.remove("minima.dat")

            # Build command
            cmd = [
                "python3", "run_throughput_MD.py",
                "--nSys", str(nSys),
                "--xyz_name", self.xyz_name.get(),
                "--dovdW", str(self.dovdW.get()),
                "--doSurfAtoms", str(self.doSurfAtoms.get()),
                "--GridFF", str(self.bGridFF.get()) if self.doSurfAtoms.get() else "-1",
                "--Fconv", self.Fconv.get(),
                "--perframe", str(perframe),
                "--perVF", str(pvf),
                "--bSaveToDatabase", str(self.bSaveToDatabase.get())
            ]

            # Add surface file if using surface atoms
            if self.doSurfAtoms.get() and N > 0:
                # Use surface file template with placeholder substitution
                surf_template = self.surf_name.get()
                surf_name = surf_template.replace("${N}", str(N))
                cmd.extend(["--gridnPBC", nPBC])
                cmd.extend(["--surf_name", surf_name])

            # Run simulation
            self.output_text.insert(tk.END, f"Running simulation with command: {' '.join(cmd)}\n")
            self.output_text.update()
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")

            # if result.returncode != 0:
            #     self.output_text.insert(tk.END, f"Simulation failed: {result.stderr}\n")
            #     return None

            # Check if minima.dat was created
            if not os.path.exists("minima.dat"):
                self.output_text.insert(tk.END, f"Warning: minima.dat was not created by simulation\n")
                return None

            # Generate result file name (following bash script naming convention)
            dovdW_bit = 1 if self.dovdW.get() else 0
            doSurfAtoms_bit = 1 if self.doSurfAtoms.get() else 0
            bGridFF_bit = 1 if self.bGridFF.get() > 0 else 0
            bSpline_bit = 1 if self.bGridFF.get() in [5, 6] else 0
            bTex_bit = 1 if self.bTex.get() else 0
            flags_bitnum = (dovdW_bit << 4) | (doSurfAtoms_bit << 3) | (bGridFF_bit << 2) | (bSpline_bit << 1) | bTex_bit

            # Ensure nPBC keeps the full format with parentheses and commas
            # Replace problematic characters for filename safety
            nPBC_safe = nPBC.replace('(', '_').replace(')', '_').replace(',', '_')
            name = (f"minima__{flags_bitnum:04b}_surf:_NaCl_{N}x{N}_nPBC_{nPBC_safe}_"
                   f"nloc:_MMFF_{nlocMMFF}_move_{nlocmove}_NBFF_{nlocNBFF}_surf_{nlocSurf}_"
                   f"gridFF_{nlocGridFF}_gridFFbSpline_{nlocGridFFbSpline}___"
                   f"replica:_{nSys}_perframe:_{perframe}_perVF:_{pvf}")

            print
            # Process minima.dat file
            # Read the result data from minima.dat
            with open("minima.dat", 'r') as f:
                lines = f.readlines()
                if not lines:
                    self.output_text.insert(tk.END, f"Warning: minima.dat is empty\n")
                    return None

                # Find the last data line (skip comments)
                last_line = None
                for line in reversed(lines):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        last_line = line
                        break

                if not last_line:
                    self.output_text.insert(tk.END, f"Warning: No data found in minima.dat\n")
                    return None

            # Move minima.dat to results/all with proper name
            target_file = f"results/all/{name}.dat"
            os.rename("minima.dat", target_file)

            # Extract throughput value from column 6 (0-based indexing)
            parts = last_line.split()
            if len(parts) > 5:
                try:
                    value = float(parts[5])

                    # Append entry to results.dat file
                    with open("results/all/results.dat", 'a') as rf:
                        rf.write(f"{name} {last_line}\n")

                    return value
                except ValueError:
                    self.output_text.insert(tk.END, f"Warning: Could not parse throughput value from: {parts[5]}\n")
                    return None
            else:
                self.output_text.insert(tk.END, f"Warning: Insufficient data columns in result line\n")
                return None

        except Exception as e:
            self.output_text.insert(tk.END, f"Error in simulation: {str(e)}\n")
            return None

    def save_parameters(self):
        """Save current parameters to a file"""
        filename = filedialog.asksaveasfilename(
            title="Save parameters",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            try:
                with open(filename, 'w') as f:
                    f.write(f"dovdW={self.dovdW.get()}\n")
                    f.write(f"doSurfAtoms={self.doSurfAtoms.get()}\n")
                    f.write(f"bGridFF={self.bGridFF.get()}\n")
                    f.write(f"bTex={self.bTex.get()}\n")
                    f.write(f"xyz_name={self.xyz_name.get()}\n")
                    f.write(f"surf_name={self.surf_name.get()}\n")
                    f.write(f"Fconv={self.Fconv.get()}\n")
                    f.write(f"replicas={self.replicas.get()}\n")
                    f.write(f"perframes={self.perframes.get()}\n")
                    f.write(f"perVF={self.perVF.get()}\n")
                    f.write(f"nPBC={self.nPBC.get()}\n")
                    f.write(f"Ns={self.Ns.get()}\n")
                    f.write(f"nlocMMFFs={self.nlocMMFFs.get()}\n")
                    f.write(f"nlocmoves={self.nlocmoves.get()}\n")
                    f.write(f"nlocNBFFs={self.nlocNBFFs.get()}\n")
                    f.write(f"nlocSurfs={self.nlocSurfs.get()}\n")
                    f.write(f"nlocGridFFs={self.nlocGridFFs.get()}\n")
                    f.write(f"nlocGridFFbSplines={self.nlocGridFFbSplines.get()}\n")
                messagebox.showinfo("Success", "Parameters saved successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save parameters: {str(e)}")
    
    def parse_surface_range(self, surface_str):
        """Parse surface range string like '1-16' or '1,4,8,16' into list of surface names"""
        surface_str = surface_str.strip()
        surfaces = []

        if '-' in surface_str and ',' not in surface_str:
            # Range format like "1-16"
            try:
                start, end = map(int, surface_str.split('-'))
                for i in range(start, end + 1):
                    surfaces.append(f"NaCl_{i}x{i}")
            except ValueError:
                # Fallback to single value
                surfaces.append(f"NaCl_{surface_str}x{surface_str}")
        else:
            # Comma-separated or single value
            parts = surface_str.split(',')
            for part in parts:
                part = part.strip()
                if part:
                    surfaces.append(f"NaCl_{part}x{part}")

        return surfaces

    def load_preset(self, filename):
        """Load parameters from a preset file"""
        try:
            if os.path.exists(filename):
                self._load_parameters_from_file(filename)
                messagebox.showinfo("Success", f"Loaded preset: {filename}")
            else:
                messagebox.showerror("Error", f"Preset file not found: {filename}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load preset: {str(e)}")

    def load_parameters(self):
        """Load parameters from a file"""
        filename = filedialog.askopenfilename(
            title="Load parameters",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            self._load_parameters_from_file(filename)
            messagebox.showinfo("Success", "Parameters loaded successfully!")

    def _load_parameters_from_file(self, filename):
        """Internal method to load parameters from file"""
        try:
            with open(filename, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        if '=' in line:
                            key, value = line.split('=', 1)
                            key = key.strip()
                            value = value.strip()

                            if hasattr(self, key):
                                var = getattr(self, key)
                                if isinstance(var, tk.IntVar):
                                    var.set(int(value))
                                elif isinstance(var, tk.StringVar):
                                    var.set(value)
        except Exception as e:
            raise Exception(f"Failed to load parameters: {str(e)}")

def main():
    root = tk.Tk()
    app = ThroughputGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
