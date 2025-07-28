import sys
import os
import importlib.util
import numpy as np
import pyopencl as cl

from PyQt5.QtWidgets import QApplication, QHBoxLayout, QWidget, QFileDialog, QVBoxLayout, QGroupBox
from PyQt5 import QtWidgets
from PyQt5.QtCore import QTimer, Qt

from ..OGL.BaseGUI import BaseGUI
from .GLCLGUI import GLCLWidget
from .OCLsystem import OCLSystem
from .OGLsystem import OGLSystem

class BakedKernelCall:
    """Encapsulates a pre-baked OpenCL kernel call with all parameters resolved."""
    
    def __init__(self, ocl_system, kernel_name, kernel_obj, global_size_resolver, local_size, arg_infos):
        self.ocl_system = ocl_system
        self.kernel_name = kernel_name
        self.kernel_obj = kernel_obj
        self.global_size_resolver = global_size_resolver
        self.local_size = local_size
        self.arg_infos = arg_infos  # List of (param_name, is_buffer, type_hint)

    def execute(self):
        """Execute the kernel with current parameters."""
        try:
            global_size = self.global_size_resolver()
            
            # Prepare arguments for kernel execution
            kernel_args = []
            for param_name, is_buffer, type_hint in self.arg_infos:
                if is_buffer:
                    kernel_args.append(self.ocl_system.get_buffer(param_name))
                else:
                    param_value = self.ocl_system.get_kernel_param(param_name)
                    # Type conversion based on hint
                    if type_hint == 'int':
                        kernel_args.append(np.int32(param_value))
                    elif type_hint == 'float':
                        kernel_args.append(np.float32(param_value))
                    else:
                        kernel_args.append(param_value)

            self.ocl_system.execute_kernel_with_args(
                self.kernel_obj, global_size, self.local_size, kernel_args
            )
        except Exception as e:
            print(f"Error executing kernel {self.kernel_name}: {e}")
            raise


class GLCLBrowser(BaseGUI):
    """Main browser class for loading and executing scientific simulation scripts."""
    
    def __init__(self, python_script_path=None):
        super().__init__()
        self.setWindowTitle("GLCL Browser - Scientific Simulation Framework")
        self.setGeometry(100, 100, 1600, 900)

        # Core systems
        self.ogl_system = OGLSystem()
        self.ocl_system = OCLSystem()
        
        # Runtime state
        self.current_script = None
        self.current_config = None
        self.baked_kernel_calls = []
        self.sim_params = {}
        self.param_widgets = {}
        
        # Create main widget and layout
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        
        # Setup GLCL widget
        self.glcl_widget = GLCLWidget(self)
        self.glcl_widget.set_systems(self.ogl_system, self.ocl_system)

        # Build UI
        self._build_ui()
        
        # Load initial script if provided
        if python_script_path is not None:
            self.load_and_apply_script(python_script_path)
            
    def _build_ui(self):
        """Build the main user interface."""
        main_layout = QHBoxLayout()
        self.main_widget.setLayout(main_layout)

        # Left panel - Controls
        control_panel = QWidget()
        control_layout = QVBoxLayout()
        control_panel.setLayout(control_layout)
        control_panel.setMaximumWidth(300)
        
        # File operations group
        file_group = QGroupBox("Simulation Script")
        file_layout = QVBoxLayout()
        file_group.setLayout(file_layout)
        
        self.button("Load Script", self.load_simulation_script, layout=file_layout)
        self.button("Reload Current", self.reload_current_script, layout=file_layout)
        
        if hasattr(self, 'script_label'):
            file_layout.addWidget(self.script_label)
        
        control_layout.addWidget(file_group)
        
        # Parameters group
        params_group = QGroupBox("Simulation Parameters")
        self.params_layout = QVBoxLayout()
        params_group.setLayout(self.params_layout)
        control_layout.addWidget(params_group)
        
        # Simulation controls group
        sim_group = QGroupBox("Simulation Control")
        sim_layout = QVBoxLayout()
        sim_group.setLayout(sim_layout)
        
        self.button("Start/Pause", self.toggle_simulation, layout=sim_layout)
        self.button("Reset", self.reset_simulation, layout=sim_layout)
        
        control_layout.addWidget(sim_group)
        control_layout.addStretch()

        # Add widgets to main layout
        main_layout.addWidget(self.glcl_widget, 1)
        main_layout.addWidget(control_panel)
        
        # Timer for simulation updates
        self.sim_timer = QTimer()
        self.sim_timer.timeout.connect(self.update_simulation)
        self.simulation_running = False

    def load_and_apply_script(self, script_path):
        """Load and apply a simulation script."""
        try:
            # Clear previous state
            self.baked_kernel_calls.clear()
            self.sim_params.clear()
            self.param_widgets.clear()
            
            # Load the module
            spec = importlib.util.spec_from_file_location("simulation_script", script_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            # Extract configuration
            config = getattr(module, 'config', None)
            init_func = getattr(module, 'init', None)

            if config is None:
                raise ValueError(f"Python script '{script_path}' must define a 'config' dictionary.")

            self.current_script = script_path
            self.current_config = config
            
            self.apply_simulation_config(config, init_func=init_func, script_path=script_path)
            print(f"Successfully loaded simulation script: {script_path}")
            
        except Exception as e:
            print(f"Error loading script {script_path}: {e}")
            import traceback
            traceback.print_exc()
            print(f"Error loading simulation script {script_path}: {e}")

    def load_simulation_script(self):
        """Load a simulation script via file dialog."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Simulation Script", "", "Python Files (*.py)"
        )
        if file_path:
            self.load_and_apply_script(file_path)

    def reload_current_script(self):
        """Reload the current simulation script."""
        if self.current_script:
            self.load_and_apply_script(self.current_script)

    def apply_simulation_config(self, config, init_func=None, script_path=None):
        """Apply the simulation configuration to both OpenCL and OpenGL systems."""
        print("Applying simulation configuration...")
        self.current_config = config
        
        # Clear previous state
        self.baked_kernel_calls.clear()
        self.ocl_system.clear_buffers()
        
        # Setup parameter controls
        self.create_parameter_controls(config)
        
        # Determine script directory for relative paths
        script_dir = os.path.dirname(script_path) if script_path else os.getcwd()
        
        # Initialize simulation data
        self.init_simulation_data(config, init_func, script_dir)
        
        # Setup OpenCL system
        self.setup_opencl_system(config, script_dir)
        
        # Setup OpenGL system
        self.setup_opengl_system(config, script_dir)
        
        # Bake kernels for execution
        self.bake_kernels(config)
        
        # Configure GLCL widget
        self.glcl_widget.set_config(config)
        
        # Start simulation if not already running
        if not self.simulation_running:
            self.toggle_simulation()

    def init_simulation_data(self, config, init_func=None, script_dir="."):
        """Initialize simulation data using init function or defaults."""
        parameters = config.get("parameters", {})
        buffers_config = config.get("buffers", {})
        
        # Initialize buffers
        for buffer_name, buffer_info in buffers_config.items():
            size_expr, stride, dtype_str = buffer_info
            
            # Resolve size expression
            if isinstance(size_expr, str) and size_expr in parameters:
                size = parameters[size_expr][0]
            else:
                size = int(size_expr)
            
            # Create numpy array
            if dtype_str == "f4":
                data = np.zeros((size, stride), dtype=np.float32)
            elif dtype_str == "f8":
                data = np.zeros((size, stride), dtype=np.float64)
            elif dtype_str == "i4":
                data = np.zeros((size, stride), dtype=np.int32)
            else:
                data = np.zeros((size, stride), dtype=np.float32)
            
            # Use init function data if available
            if init_func:
                init_data = init_func()
                if buffer_name in init_data:
                    data = init_data[buffer_name]
            
            # Create OpenCL buffer
            self.ocl_system.create_buffer(buffer_name, data)

    def setup_opencl_system(self, config, script_dir="."):
        """Setup OpenCL system based on configuration."""
        # Load OpenCL programs
        opencl_sources = config.get("opencl_source", [])
        for source_file in opencl_sources:
            filepath = os.path.join(script_dir, source_file)
            if os.path.exists(filepath):
                program_name = os.path.splitext(os.path.basename(source_file))[0]
                self.ocl_system.load_program(program_name, filepath)
            else:
                print(f"Warning: OpenCL source file not found: {filepath}")

    def setup_opengl_system(self, config, script_dir="."):
        """Setup OpenGL system based on configuration."""
        shaders_config = config.get("opengl_shaders", {})
        for shader_name, shader_info in shaders_config.items():
            vertex_file, fragment_file, uniforms = shader_info
            
            vertex_path = os.path.join(script_dir, vertex_file)
            fragment_path = os.path.join(script_dir, fragment_file)
            
            if os.path.exists(vertex_path) and os.path.exists(fragment_path):
                with open(vertex_path, 'r') as f:
                    vertex_src = f.read()
                with open(fragment_path, 'r') as f:
                    fragment_src = f.read()
                
                self.ogl_system.load_shader_program(
                    shader_name, vertex_path, fragment_path, uniforms
                )
            else:
                print(f"Warning: Shader files not found: {vertex_path}, {fragment_path}")

    def bake_kernels(self, config):
        """Bake OpenCL kernels for efficient execution."""
        kernels_config = config.get("kernels", {})
        pipeline = config.get("kernel_pipeline", [])
        parameters = config.get("parameters", {})
        
        self.baked_kernel_calls.clear()
        
        for kernel_name in pipeline:
            if kernel_name in kernels_config:
                kernel_info = kernels_config[kernel_name]
                local_size, global_size_expr, buffer_names, param_names = kernel_info
                
                # Get kernel object
                kernel_obj = self.ocl_system.get_kernel(kernel_name)
                if kernel_obj is None:
                    print(f"Warning: Kernel '{kernel_name}' not found")
                    continue
                
                # Build argument information
                arg_infos = []
                for buf_name in buffer_names:
                    arg_infos.append((buf_name, True, None))
                for param_name in param_names:
                    param_info = parameters.get(param_name, [None, "float", 0])
                    arg_infos.append((param_name, False, param_info[1]))
                
                # Create global size resolver
                def make_global_resolver(expr):
                    return lambda: self._resolve_global_size(expr, parameters)
                
                baked_call = BakedKernelCall(
                    self.ocl_system, kernel_name, kernel_obj,
                    make_global_resolver(global_size_expr), local_size, arg_infos
                )
                
                self.baked_kernel_calls.append(baked_call)

    def _resolve_global_size(self, size_expr, parameters):
        """Resolve global size expression to actual integer."""
        if isinstance(size_expr, (list, tuple)):
            resolved = []
            for item in size_expr:
                if isinstance(item, str) and item in parameters:
                    resolved.append(parameters[item][0])
                else:
                    resolved.append(int(item))
            return tuple(resolved)
        elif isinstance(size_expr, str) and size_expr in parameters:
            return parameters[size_expr][0]
        else:
            return int(size_expr)

    def toggle_simulation(self):
        """Toggle simulation running state."""
        self.simulation_running = not self.simulation_running
        if self.simulation_running:
            self.sim_timer.start(16)  # ~60 FPS
            print("Simulation started")
        else:
            self.sim_timer.stop()
            print("Simulation paused")

    def reset_simulation(self):
        """Reset simulation to initial state."""
        if self.current_script and self.current_config:
            print("Resetting simulation...")
            self.load_and_apply_script(self.current_script)

    def update_simulation(self):
        """Main simulation update loop."""
        try:
            # Update parameters from GUI
            self.update_sim_uniforms()
            
            # Execute baked kernel calls
            for kernel_call in self.baked_kernel_calls:
                kernel_call.execute()
                
            # Update OpenGL visualization
            self.glcl_widget.update()
            
        except Exception as e:
            print(f"Simulation error: {e}")
            import traceback
            traceback.print_exc()
            self.simulation_running = False

    def update_sim_uniforms(self):
        """Update simulation uniforms based on current parameters."""
        if not self.current_config:
            return
            
        parameters = self.current_config.get("parameters", {})
        
        # Update kernel parameters from GUI
        for param_name, widget in self.param_widgets.items():
            if param_name in parameters:
                param_info = parameters[param_name]
                type_str = param_info[1]
                value = widget.value()
                
                if type_str == "int":
                    self.ocl_system.set_kernel_param(param_name, np.int32(value))
                elif type_str == "float":
                    self.ocl_system.set_kernel_param(param_name, np.float32(value))



    
    # Main entry point for the application
if __name__ == '__main__':
    app = QApplication(sys.argv)
    
    # Determine default script path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    default_script = os.path.join(current_dir, "scripts", "nbody.py")
    
    # Create and show browser
    browser = GLCLBrowser(python_script_path=default_script)
    browser.show()
    
    sys.exit(app.exec_())
