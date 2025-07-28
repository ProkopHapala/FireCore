import sys
import json
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QWidget, QFileDialog
from PyQt5 import QtWidgets
from PyQt5.QtCore import QTimer
from PyQt5.QtCore import Qt

from ..OGL.BaseGUI import BaseGUI
from .GLCLGUI_old import GLCLWidget
from .OCLsystem import OCLSystem
from .OGLsystem_old import OGLSystem
import numpy as np
import pyopencl as cl
import os

class GLCLBrowser(BaseGUI):
    def __init__(self, json_filepath=None):
        print("GLCLBrowser() json_filepath: ", json_filepath)
        super().__init__()
        self.setWindowTitle("GLCL Browser")
        self.setGeometry(100, 100, 1600, 900) # Increased window size

        self.ogl_system = OGLSystem()
        self.ocl_system = OCLSystem()

        self.glcl_widget = GLCLWidget(self)
        self.glcl_widget.set_systems(self.ogl_system, self.ocl_system)

        main_layout = QHBoxLayout()
        self.main_widget.setLayout(main_layout)

        self.param_widgets = {} # Initialize param_widgets

        control_panel  = QtWidgets.QWidget()
        control_layout = QtWidgets.QVBoxLayout()
        control_panel.setLayout(control_layout)
        control_layout.setAlignment(Qt.AlignTop)
        control_panel.setMaximumWidth(300) # Narrow right column

        # Add a button to load JSON config
        self.button("Load Simulation JSON", self.load_simulation_json, layout=control_layout)

        # Parameters section (will be populated dynamically)
        self.params_layout = self.add_form_layout(control_layout) # BaseGUI method

        # Auto-load nbody_simulation.json if path is provided or default exists
        if json_filepath is not None:
            with open(json_filepath, 'r') as f:
                config = json.load(f)
            self.apply_simulation_config(config, filepath=json_filepath)
            print(f"Auto-loaded simulation config from: {json_filepath}")

        main_layout.addWidget(self.glcl_widget, 1) # GLCLWidget takes most space
        main_layout.addWidget(control_panel) # Control panel on the right

    def load_simulation_json(self):
        print("load_simulation_json()")
        file_dialog = QFileDialog()
        filepath, _ = file_dialog.getOpenFileName(self, "Load Simulation JSON", "", "JSON Files (*.json)")
        if filepath:
            with open(filepath, 'r') as f:
                config = json.load(f)
            self.apply_simulation_config(config, filepath=filepath)
            print(f"Loaded simulation config from: {filepath}")

    def apply_simulation_config(self, config, filepath=None):
        print("!!!!! apply_simulation_config()")
        # This method will parse the JSON and configure OGLSystem, OCLSystem, and GLCLWidget
        # For now, let's just populate parameters in the GUI
        if "parameters" in config:
            # The structure of 'parameters' in JSON should match what populate_params_from_json expects
            # e.g., {"param_name": ["type", [default_value], step]}
            # For now, let's assume a simple structure for testing:
            # {"param_name": [default_value, step]}
            # We need to convert it to the format expected by BaseGUI's populate_params_from_json
            # BaseGUI expects: {"name": (typ, defaults, step)}
            # Let's assume for now that 'typ' is always 'scalar' and 'defaults' is a list of one item for single values
            
            # Example expected JSON 'parameters' format:
            # "parameters": {
            #     "gravity": {"value": 9.81, "step": 0.1, "type": "float"},
            #     "num_particles": {"value": 1000, "step": 1, "type": "int"}
            # }

            # Convert to BaseGUI format: {"name": ("type", [default_value], step)}
            params_for_gui = {}
            for name, p_config in config["parameters"].items():
                value = p_config.get("value")
                step = p_config.get("step", 0.1) # Default step if not provided
                param_type = p_config.get("type", "float") # Default type
                
                # BaseGUI's populate_params_from_json expects 'defaults' as a list
                if isinstance(value, list):
                    defaults = value
                else:
                    defaults = [value]
                
                params_for_gui[name] = (param_type, defaults, step)
            
            self.populate_params_from_json(params_for_gui, callback_func=self.update_sim_uniforms)

        # path to this python source file
        root_path = os.path.dirname(__file__)

        # Configure OCLSystem and OGLSystem based on JSON
        if "opencl_kernels" in config:
            for kernel_info in config["opencl_kernels"]:
                name = kernel_info["name"]
                filepath = os.path.join(root_path, kernel_info["filepath"]) # Adjust path
                self.ocl_system.load_program(name, filepath)

        if "opengl_shaders" in config:
            for shader_info in config["opengl_shaders"]:
                name = shader_info["name"]
                # make absolute path for relative paths
                vertex_filepath   = os.path.join( root_path, shader_info["vertex_filepath"]) ; 
                fragment_filepath = os.path.join( root_path, shader_info["fragment_filepath"]) ; print("fragment_filepath ", fragment_filepath)
                
                with open(vertex_filepath, 'r') as f: vertex_src = f.read()
                with open(fragment_filepath, 'r') as f: fragment_src = f.read()
                
                # Store shader sources in the widget for deferred loading
                self.glcl_widget.set_shader_sources(name, vertex_src, fragment_src)

        # Initialize particle data and buffers
        particle_count = int(self.param_widgets["particle_count"].value()) if "particle_count" in self.param_widgets else 2048
        dt = float(self.param_widgets["dt"].value()) if "dt" in self.param_widgets else 0.001

        # Initial particle positions and velocities
        positions = (np.random.rand(particle_count, 4) * 2 - 1).astype(np.float32)
        positions[:, 2] = 0.0; positions[:, 3] = 1.0 # Set z to 0 and w to 1
        velocities = (np.random.rand(particle_count, 4) * 0.1).astype(np.float32)

        # Create OpenCL buffers
        self.ocl_system.create_buffer("positions", positions.nbytes, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR)
        self.ocl_system.toGPU("positions", positions)
        self.ocl_system.create_buffer("velocities", velocities.nbytes, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR)
        self.ocl_system.toGPU("velocities", velocities)

        # Set kernel parameters
        self.ocl_system.set_kernel_param("dt", np.float32(dt))
        self.ocl_system.set_kernel_param("particle_count", np.int32(particle_count))

        # Pass initial data to GLCLWidget for OpenGL setup
        self.glcl_widget.set_particle_data(positions, particle_count)

        # Start simulation timer
        if not hasattr(self, 'sim_timer'):
            self.sim_timer = QTimer(self)
            self.sim_timer.timeout.connect(self.update_simulation)
        self.sim_timer.start(16) # ~60 FPS
        print("apply_simulation_config() DONE")

    def update_simulation(self):
        try:
            # Acquire GL buffer if sharing (not yet implemented, will use copy for now)
            # Execute OpenCL kernel
            self.ocl_system.execute_kernel("nbody_sim", "nbody_sim", (self.glcl_widget.particle_count,), local_size=None)

            # Copy results back to host or GL buffer
            # For now, copy to host and then update GL buffer
            updated_positions = np.empty_like(self.glcl_widget.positions)
            self.ocl_system.fromGPU("positions", updated_positions)
            self.glcl_widget.update_particle_vbo(updated_positions)
            self.glcl_widget.update()
        except Exception as e:
            print("GLCLBrowser: Exception in update_simulation(), stopping timer.")
            self.sim_timer.stop()
            raise e

    def update_sim_uniforms(self):
        # This method is called by BaseGUI when parameters are changed
        # Update OpenCL kernel arguments based on current GUI parameter values
        if "dt" in self.param_widgets:
            dt = self.param_widgets["dt"].value()
            self.ocl_system.set_kernel_param("dt", np.float32(dt))
        if "particle_count" in self.param_widgets:
            particle_count = self.param_widgets["particle_count"].value()
            self.ocl_system.set_kernel_param("particle_count", np.int32(particle_count))
            # If particle count changes, re-initialize positions and velocities
            # This is a simplified approach; a more robust solution might handle resizing buffers.
            # For now, we'll just re-create them if the count changes significantly.
            if particle_count != self.glcl_widget.particle_count:
                # Ensure particle_count is an integer for numpy array creation
                particle_count_int = int(particle_count)
                positions = (np.random.rand(particle_count_int, 4) * 2 - 1).astype(np.float32)
                positions[:, 2] = 0.0; positions[:, 3] = 1.0
                velocities = (np.random.rand(particle_count_int, 4) * 0.1).astype(np.float32)
                self.ocl_system.create_buffer("positions", positions.nbytes, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=positions)
                self.ocl_system.create_buffer("velocities", velocities.nbytes, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=velocities)
                self.glcl_widget.set_particle_data(positions, particle_count)

    # Main entry point for the application
if __name__ == '__main__':

    # run like this:
    #  python -u -m pyBall.GLCLBrowser
    #  __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia python -u -m pyBall.GLCLBrowser
    
    app = QApplication(sys.argv)
    json_path = os.path.join(os.path.dirname(__file__), "scripts/nbody.json")
    viewer = GLCLBrowser(json_filepath=json_path)
    viewer.show()
    sys.exit(app.exec_())
