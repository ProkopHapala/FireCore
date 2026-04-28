import sys
import os
import importlib.util
import numpy as np
import pyopencl as cl

from PyQt5.QtWidgets import QApplication, QHBoxLayout, QWidget, QFileDialog, QVBoxLayout, QGroupBox
from PyQt5 import QtWidgets
from PyQt5.QtCore import QTimer, Qt

import OpenGL.GL as GL

from ..OGL.BaseGUI import BaseGUI
from .GLCLGUI import GLCLWidget
from .OCLsystem import OCLSystem
from .OGLsystem import OGLSystem

class BakedKernelCall:
    """Simplified pre-baked kernel execution class."""
    
    def __init__(self, ocl_system, program_name, kernel_name, global_size, local_size, args_tuple):
        self.ocl_system = ocl_system
        self.program_name = program_name
        self.kernel_name = kernel_name
        self.global_size = global_size
        self.local_size = local_size
        self.args_tuple = args_tuple
    
    def execute(self):
        """Execute the kernel with pre-computed arguments."""
        try:
            self.ocl_system.execute_kernel(
                self.program_name, self.kernel_name, self.global_size, self.local_size, self.args_tuple
            )
        except Exception as e:
            print(f"Error executing kernel {self.kernel_name}: {e}")
            import traceback
            traceback.print_exc()
            raise


class GLCLBrowser(BaseGUI):
    """Main browser class for loading and executing scientific simulation scripts."""
    
    def __init__(self, python_script_path=None, bDebugCL=False, bDebugGL=False, nDebugFrames=5):
        super().__init__()
        self.setWindowTitle("GLCL Browser - Scientific Simulation Framework")
        self.setGeometry(100, 100, 1600, 900)

        # Core systems
        self.ogl_system = OGLSystem()
        self.ocl_system = OCLSystem()
        
        # Debug flags
        self.bDebugCL = bDebugCL
        self.bDebugGL = bDebugGL
        self.nDebugFrames = nDebugFrames
        self.debug_frame_counter = 0

        # Runtime state
        self.current_script = None
        self.current_config = None
        self.baked_kernel_calls = []
        self.sim_params = {}
        self.param_widgets = {}
        self.buffer_shapes = {}
        
        # Dynamic buffer management
        self.gl_objects = {}  # Dictionary mapping buffer names to GLobject instances
        self.render_pipeline_info = []  # List of (shader_name, element_count_name, vertex_buffer_name, index_buffer_name)
        self.buffer_data = {}  # Dictionary storing buffer data for deferred GLobject creation
        
        # Pre-baked buffer synchronization
        self.buffers_to_sync = []  # Pre-computed list of buffer names for GPU→CPU→GPU transfer
        
        # Create main widget and layout
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        
        # Setup GLCL widget
        self.glcl_widget = GLCLWidget(self, enable_opengl_debug=self.bDebugGL)
        self.glcl_widget.set_systems(self.ogl_system, self.ocl_system, self)

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

            # Debug: Print config parameters right after loading
            print("GLCLBrowser::load_and_apply_script() Config loaded from script:")
            params = config.get("parameters", {})
            for name, info in params.items():
                print(f"  {name}: {info}")

            self.current_script = script_path
            self.current_config = config
            
            self.apply_simulation_config(config, init_func=init_func, script_path=script_path)
            print(f"GLCLBrowser::load_and_apply_script() Successfully loaded simulation script: {script_path}")
            
        except Exception as e:
            print(f"GLCLBrowser::load_and_apply_script() Error loading script {script_path}: {e}")
            self.on_exception(e)

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
        print("GLCLBrowser::apply_simulation_config() Original config parameters:")
        params = config.get("parameters", {})
        for name, info in params.items():
            print(f"  {name}: {info}")
        
        # Create a deep copy of the config to prevent in-place modification
        import copy
        self.current_config = copy.deepcopy(config)
        
        self.baked_kernel_calls.clear()
        self.ocl_system.clear_buffers()
        self.param_widgets.clear()
        
        # Create parameter controls
        self.create_parameter_controls(self.current_config)
        
        # Initialize simulation data
        self.init_simulation_data(self.current_config, init_func=init_func, script_dir=os.path.dirname(script_path) if script_path else ".")
        
        # Setup OpenCL and OpenGL systems
        self.setup_opencl_system(self.current_config, script_dir=os.path.dirname(script_path) if script_path else ".")
        self.setup_opengl_system(self.current_config, script_dir=os.path.dirname(script_path) if script_path else ".")
        
        # Bake kernels for execution - use original config values
        self.bake_kernels(self.current_config)
        
        # Set initial kernel parameters from config
        self.set_initial_kernel_params()
        
        print("Simulation configuration applied successfully.")
        
        # Print summary of baked kernels and buffers before starting runtime loop
        print("GLCLBrowser::apply_simulation_config() Buffers:")
        for b, shp in self.buffer_shapes.items():
            print(f"  {b}: shape {shp}")
        print("GLCLBrowser::apply_simulation_config() Baked Kernels:")
        for kc in self.baked_kernel_calls:
            print(f"  {kc.kernel_name}: global_size={kc.global_size}, local_size={kc.local_size}")
        
        print("=====================")

        # Shaders are now handled in setup_opengl_system method
        pass
        
        # Create dynamic buffer registry from user script data
        self._create_dynamic_buffer_registry()
        
        # Now that dynamic buffer registry is created, precompute buffer sync list
        # This must happen after render_pipeline_info is populated
        self._precompute_buffer_sync_list(self.current_config)
        
        # Don't compile shaders here - wait until OpenGL context is ready
        # Shaders will be compiled in GLCLWidget.initializeGL()
        
        # Start simulation if not already running
        if not self.simulation_running:
            self.toggle_simulation()

    def init_simulation_data(self, config, init_func=None, script_dir="."):
        """Initialize simulation data using init function or defaults."""
        print("=== init_simulation_data called ===")
        print(f"config type: {type(config)}")
        print(f"config keys: {list(config.keys())}")
        print(f"full config: {config}")
        
        # Get buffer configuration
        buffers_config = config.get("buffers", {})
        parameters = config.get("parameters", {})
        
        print(f"GLCLBrowser::init_simulation_data() buffers_config:{buffers_config} parameters:{parameters} of type:{type(parameters)}")
        
        # Initialize data dictionary
        data = {}
        self.init_data = {}  # Store for later use
        
        # Call init function if provided
        if init_func is not None:
            try:
                print("Calling init function...")
                init_data = init_func()
                print(f"init_data returned: {type(init_data)} {init_data}")
                if init_data:
                    data.update(init_data)
                    self.init_data.update(init_data)  # Store for GLCLWidget
            except Exception as e:
                self.on_exception(e)
        
        # Create buffers based on configuration
        for buffer_name, buffer_info in buffers_config.items():
            print(f"GLCLBrowser::init_simulation_data() Processing buffer {buffer_name}: {buffer_info}")
            if len(buffer_info) == 3:
                size_expr, stride, dtype = buffer_info
                
                print(f"GLCLBrowser::init_simulation_data() size_expr: {size_expr} (type: {type(size_expr)})")
                print(f"GLCLBrowser::init_simulation_data() parameters accessible: {list(parameters.keys())}")
                
                # Resolve size expression
                if isinstance(size_expr, str):
                    if size_expr in parameters:
                        param_value = parameters[size_expr]
                        print(f"Found parameter {size_expr}: {param_value} (type: {type(param_value)})")
                        if isinstance(param_value, tuple) and len(param_value) >= 1:
                            size = param_value[0]
                        else:
                            size = param_value
                    else:
                        print(f"GLCLBrowser::init_simulation_data() ERROR: Parameter {size_expr} not found in parameters")
                        raise KeyError(f"Parameter {size_expr} not found")
                else:
                    size = int(size_expr)
                
                print(f"GLCLBrowser::init_simulation_data() Resolved size: {size} (type: {type(size)})")
                
                # Create buffer data
                if buffer_name not in data:
                    if dtype == "f4":
                        data[buffer_name] = np.zeros((size, stride), dtype=np.float32)
                    elif dtype == "f8":
                        data[buffer_name] = np.zeros((size, stride), dtype=np.float64)
                    else:
                        data[buffer_name] = np.zeros((size, stride), dtype=np.float32)
                
                print(f"GLCLBrowser::init_simulation_data() Creating buffer {buffer_name} with shape {data[buffer_name].shape}")
                
                # Create buffer and upload data to OpenCL
                buffer_size = data[buffer_name].nbytes
                print(f"GLCLBrowser::init_simulation_data() Buffer size in bytes: {buffer_size}")
                self.ocl_system.create_buffer(buffer_name, buffer_size)
                self.ocl_system.toGPU(buffer_name, data[buffer_name])
                
                print(f"GLCLBrowser::init_simulation_data() Successfully created buffer {buffer_name}: {data[buffer_name].shape}")
                self.buffer_shapes[buffer_name] = data[buffer_name].shape
                
                # Store data for potential GLCLWidget use
                if buffer_name not in self.init_data:
                    self.init_data[buffer_name] = data[buffer_name]

    def setup_opencl_system(self, config, script_dir="."):
        """Setup OpenCL system based on configuration."""
        # Load OpenCL programs
        opencl_sources = config.get("opencl_source", [])
        print(f"GLCLBrowser::setup_opencl_system() OpenCL sources: {opencl_sources}")
        print(f"GLCLBrowser::setup_opencl_system() Script directory: {script_dir}")
        
        for source_file in opencl_sources:
            filepath = os.path.join(script_dir, source_file)
            print(f"GLCLBrowser::setup_opencl_system() Loading OpenCL program from: {filepath}")
            
            if os.path.exists(filepath):
                program_name = os.path.splitext(os.path.basename(source_file))[0]
                print(f"GLCLBrowser::setup_opencl_system() Loading program '{program_name}' from {filepath}")
                self.ocl_system.load_program(program_name, filepath)
                print(f"GLCLBrowser::setup_opencl_system() Programs now loaded: {list(self.ocl_system.programs.keys())}")
            else:
                print(f"Warning: OpenCL source file not found: {filepath}")
                print(f"GLCLBrowser::setup_opencl_system() Current directory contents: {os.listdir(script_dir)}")

    def setup_opengl_system(self, config, script_dir="."):
        """Setup OpenGL system with all shaders from config."""
        # Store shader information for later compilation when OpenGL context is available
        shaders_config = config.get("opengl_shaders", {})
        self.shader_configs = shaders_config
        self.script_dir = script_dir
        print(f"GLCLBrowser::setup_opengl_system() Stored shader configs for later compilation: {list(shaders_config.keys())}")
        
    def compile_shaders(self):
        """Compile stored shaders now that OpenGL context is available."""
        if not hasattr(self, 'shader_configs'):
            print("GLCLBrowser::compile_shaders() No shader configs stored")
            return
            
        # Dictionary to store compiled shader programs
        self.shader_programs = {}
        
        # Compile all stored shaders
        for shader_name, shader_info in self.shader_configs.items():
            vertex_file, fragment_file, uniforms = shader_info
            
            vertex_path = os.path.join(self.script_dir, vertex_file)
            fragment_path = os.path.join(self.script_dir, fragment_file)
            
            try:
                if os.path.exists(vertex_path) and os.path.exists(fragment_path):
                    with open(vertex_path, 'r') as f:
                        vertex_src = f.read()
                    with open(fragment_path, 'r') as f:
                        fragment_src = f.read()
                    
                    # Load shader program into OGLSystem now that context is available
                    success = self.ogl_system.load_shader_program(shader_name, vertex_src, fragment_src)
                    shader_program = self.ogl_system.get_shader_program(shader_name)
                    
                    if success and shader_program:
                        self.shader_programs[shader_name] = shader_program
                        print(f"GLCLBrowser::compile_shaders() Loaded and compiled shader: {shader_name} (program ID: {shader_program})")
                        
                        # Set default uniforms if specified
                        if uniforms:
                            print(f"GLCLBrowser::compile_shaders() Setting default uniforms for shader {shader_name}: {uniforms}")
                    else:
                        print(f"GLCLBrowser::compile_shaders() Warning: Failed to compile shader program: {shader_name}")
                else:
                    print(f"GLCLBrowser::compile_shaders() Warning: Shader files not found: {vertex_path}, {fragment_path}")
            except Exception as e:
                print(f"GLCLBrowser::compile_shaders() Error loading shader {shader_name}: {e}")
                import traceback
                traceback.print_exc()
        
        # Pass compiled shaders to GLCLWidget
        if hasattr(self, 'glcl_widget'):
            self.glcl_widget.set_shader_programs(self.shader_programs)

    def bake_kernels(self, config):
        """Bake OpenCL kernels for efficient execution."""
        print(">>>> GLCLBrowser::bake_kernels() Baking kernels with config:", config)
        kernels_config = config.get("kernels", {})
        pipeline = config.get("kernel_pipeline", [])
        parameters = config.get("parameters", {})
        print("GLCLBrowser::bake_kernels() Parameters:")
        for name, info in parameters.items():
            print(f"  {name}: {info}")
        
        self.baked_kernel_calls.clear()
        
        # Note: Buffer sync precomputation moved to after dynamic buffer registry creation
        # Use parameter values directly from config, not from widgets
        for kernel_name in pipeline:
            if kernel_name in kernels_config:
                kernel_info = kernels_config[kernel_name]
                print(f"GLCLBrowser::bake_kernels() Processing kernel {kernel_name} with info: {kernel_info}")
                local_size, global_size_expr, buffer_names, param_names = kernel_info
                
                # Resolve global size once during baking
                global_size = self._resolve_global_size(global_size_expr, parameters)
                
                # Build pre-computed arguments tuple
                args_list = []
                
                # Add buffer arguments
                for buf_name in buffer_names:
                    buffer_obj = self.ocl_system.get_buffer(buf_name)
                    args_list.append(buffer_obj)
                
                # Add scalar arguments with pre-computed type conversion
                for param_name in param_names:
                    if param_name in parameters:
                        value, type_str, step = parameters[param_name]
                        print(f"  Using parameter '{param_name}' = {value} ({type_str})")
                        if type_str == "int":
                            args_list.append(np.int32(value))
                        elif type_str == "float":
                            args_list.append(np.float32(value))
                        else:
                            args_list.append(value)
                    else:
                        print(f"Warning: Parameter '{param_name}' not found, using 0.0f")
                        args_list.append(np.float32(0.0))
                
                args_tuple = tuple(args_list)
                
                # Create baked kernel call with pre-computed arguments
                baked_call = BakedKernelCall(
                    self.ocl_system, "nbody", kernel_name, global_size, local_size, args_tuple
                )
                self.baked_kernel_calls.append(baked_call)
                
                print(f"GLCLBrowser::bake_kernels() Successfully baked kernel {kernel_name} with {len(args_tuple)} arguments")
                for i, arg in enumerate(args_tuple):
                    if hasattr(arg, 'dtype') and arg.dtype.kind == 'f':
                        print(f"    arg[{i}] = {arg} (float)")
                    elif hasattr(arg, 'dtype') and arg.dtype.kind == 'i':
                        print(f"    arg[{i}] = {arg} (int)")
                    else:
                        print(f"    arg[{i}] = {arg}")
        print("<<<< GLCLBrowser::bake_kernels() DONE.")

    def setup_kernel_params_from_config(self, parameters):
        """Initialize all kernel parameters from the config dictionary."""
        print("GLCLBrowser::setup_kernel_params_from_config() Initializing kernel parameters...")
        for param_name, param_info in parameters.items():
            value, type_str, step = param_info
            print(f"  Setting parameter '{param_name}' = {value} ({type_str})")
            if type_str == "int":
                self.ocl_system.set_kernel_param(param_name, np.int32(value))
            elif type_str == "float":
                self.ocl_system.set_kernel_param(param_name, np.float32(value))
            else:
                self.ocl_system.set_kernel_param(param_name, value)

    def setup_kernel_params(self, parameters):
        """Initialize all kernel parameters from configuration."""
        print("GLCLBrowser::setup_kernel_params() Initializing kernel parameters...")
        
        for param_name, param_info in parameters.items():
            value, type_str, step = param_info
            print(f"  Setting parameter '{param_name}' = {value} ({type_str})")
            
            if type_str == "int":
                self.ocl_system.set_kernel_param(param_name, np.int32(value))
            elif type_str == "float":
                self.ocl_system.set_kernel_param(param_name, np.float32(value))
            else:
                self.ocl_system.set_kernel_param(param_name, value)

    def _resolve_global_size(self, size_expr, parameters):
        """Resolve global size expression to actual integer."""
        print(f"GLCLBrowser::_resolve_global_size() Resolving global size '{size_expr}' with parameters: {parameters}")
        if isinstance(size_expr, (list, tuple)):
            resolved = []
            for item in size_expr:
                if isinstance(item, str) and item in parameters:
                    resolved.append(parameters[item][0])
                else:
                    resolved.append(int(item))
            result = tuple(resolved)
            print(f"GLCLBrowser::_resolve_global_size() Resolved global size {size_expr} -> {result}")
            return result
        elif isinstance(size_expr, str) and size_expr in parameters:
            result_val = parameters[size_expr][0]
            result = (result_val,) if isinstance(result_val, int) else tuple(result_val)
            print(f"GLCLBrowser::_resolve_global_size() Resolved global size {size_expr} -> {result}")
            return result
        else:
            result_val = int(size_expr)
            result = (result_val,)
            print(f"GLCLBrowser::_resolve_global_size() Resolved global size {size_expr} -> {result}")
            return result

    def create_parameter_controls(self, config):
        """Create GUI controls for simulation parameters."""
        # Debug: Print parameters when creating controls
        print("GLCLBrowser::create_parameter_controls() Creating controls with parameters:")
        params = config.get("parameters", {})
        for name, info in params.items():
            print(f"  {name}: {info}")
            
        # Clear existing controls
        for i in reversed(range(self.params_layout.count())):
            item = self.params_layout.itemAt(i)
            if item.widget():
                item.widget().setParent(None)
        self.param_widgets.clear()
        
        parameters = config.get("parameters", {})
        for param_name, param_info in parameters.items():
            if len(param_info) == 3:
                value, type_str, step = param_info
            else:
                continue
                
            print(f"  Creating control for {param_name}: {param_info}")
                
            if type_str == "int":
                # Create integer input
                from PyQt5.QtWidgets import QSpinBox, QLabel, QHBoxLayout, QWidget
                container = QWidget()
                layout = QHBoxLayout(container)
                layout.setContentsMargins(0, 0, 0, 0)
                
                label = QLabel(param_name)
                widget = QSpinBox()
                widget.setRange(-1000000, 1000000)
                
                # Block signals BEFORE setting value to prevent immediate updates
                widget.blockSignals(True)
                widget.setValue(int(value))  # Initialize with config value
                widget.setSingleStep(int(step))
                
                layout.addWidget(label)
                layout.addWidget(widget)
                self.params_layout.addWidget(container)
                
            elif type_str == "float":
                # Create float input
                from PyQt5.QtWidgets import QDoubleSpinBox, QLabel, QHBoxLayout, QWidget
                container = QWidget()
                layout = QHBoxLayout(container)
                layout.setContentsMargins(0, 0, 0, 0)
                
                label = QLabel(param_name)
                widget = QDoubleSpinBox()
                widget.setRange(-1000000.0, 1000000.0)
                widget.setDecimals(6)  # Increased precision for small values
                
                # Block signals BEFORE setting value to prevent immediate updates
                widget.blockSignals(True)
                widget.setValue(float(value))  # Initialize with config value
                widget.setSingleStep(float(step))
                widget.setDecimals(4)
                
                layout.addWidget(label)
                layout.addWidget(widget)
                self.params_layout.addWidget(container)
            else:
                continue
                
            self.param_widgets[param_name] = widget
            
            # Store the original value to prevent reset
            widget.setProperty("original_value", float(value))
            
            # Connect value changes to update config and re-bake kernels
            if hasattr(widget, 'valueChanged'):
                if type_str == "int":
                    widget.valueChanged.connect(lambda v, p=param_name: self.update_parameter_from_widget(p, v, "int"))
                elif type_str == "float":
                    widget.valueChanged.connect(lambda v, p=param_name: self.update_parameter_from_widget(p, v, "float"))
            
            # Unblock signals after setup
            widget.blockSignals(False)

    def toggle_simulation(self):
        """Toggle simulation running state."""
        self.simulation_running = not self.simulation_running
        if self.simulation_running:
            self.sim_timer.start(16)  # ~60 FPS
            print("GLCLBrowser::toggle_simulation() Simulation started")
        else:
            self.sim_timer.stop()
            print("GLCLBrowser::toggle_simulation() Simulation paused")

    def reset_simulation(self):
        """Reset simulation to initial state."""
        if self.current_script and self.current_config:
            print("Resetting simulation...")
            self.load_and_apply_script(self.current_script)

    def update_simulation(self):
        """Main simulation update loop - simplified with pre-baked buffer management."""
        self.debug_frame_counter += 1
        try:
            # Execute baked kernel calls unless GL debugging is on
            if not self.bDebugGL:
                # Execute all pre-baked kernel calls
                for kernel_call in self.baked_kernel_calls:
                    kernel_call.execute()
                
                # Simple loop over pre-computed buffer sync list - no runtime checks needed
                for buf_name in self.buffers_to_sync:
                    # Download from OpenCL
                    host_data = np.empty(self.buffer_shapes[buf_name], dtype=np.float32)
                    self.ocl_system.fromGPU(buf_name, host_data)
                    
                    # Upload to OpenGL
                    gl_obj = self.glcl_widget.gl_objects[buf_name]
                    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, gl_obj.vbo)
                    GL.glBufferSubData(GL.GL_ARRAY_BUFFER, 0, host_data.nbytes, host_data)
                    
                    if self.bDebugCL:
                        print(f"[Frame {self.debug_frame_counter}] Updated OpenGL buffer '{buf_name}' with OpenCL data")
                
                # Print buffer snapshots for debugging
                if self.bDebugCL:
                    for buf_name in self.buffer_shapes.keys():
                        host = np.empty(self.buffer_shapes[buf_name], dtype=np.float32)
                        self.ocl_system.fromGPU(buf_name, host)
                        print(f"[Frame {self.debug_frame_counter}] Buffer '{buf_name}': 1st {host[:1]} last {host[-1:]}")
                
            self.glcl_widget.update()
            
            if self.debug_frame_counter >= self.nDebugFrames:
                self.simulation_running = False
                self.sim_timer.stop()
                print(f"GLCLBrowser::update_simulation() STOP! debug_frame_counter reached nDebugFrames: {self.nDebugFrames}")
                #exit()

        except Exception as e:
            self.on_exception(e)

    def on_exception(self, e):
        print(f"Simulation error: {e}")
        import traceback
        traceback.print_exc()
        self.simulation_running = False
        self.sim_timer.stop()
        os._exit(1)  # Forcefully terminate to avoid PyQt swallowing exceptions

    def update_parameter_from_widget(self, param_name, new_value, type_str):
        """Update parameter value from GUI widget."""
        print(f"GLCLBrowser::update_parameter_from_widget() {param_name} = {new_value} ({type_str})")
        
        if self.current_config is None:
            return
            
        parameters = self.current_config.get("parameters", {})
        if param_name in parameters:
            old_value, old_type, step = parameters[param_name]
            print(f"  Old value: {old_value}, New value: {new_value}")
            parameters[param_name] = (new_value, old_type, step)
            
            # Update OpenCL kernel parameter
            if type_str == "int":
                self.ocl_system.set_kernel_param(param_name, np.int32(new_value))
            elif type_str == "float":
                self.ocl_system.set_kernel_param(param_name, np.float32(new_value))
            
            # Re-bake kernels with updated parameters
            self.bake_kernels(self.current_config)
            
            # Update simulation uniforms from controls
            self.update_uniforms_from_controls()
            
            print(f"Parameter '{param_name}' updated to {new_value} ({type_str})")

    def set_initial_kernel_params(self):
        """Set initial kernel parameters from config without reading widgets."""
        if not self.current_config:
            return
            
        parameters = self.current_config.get("parameters", {})
        for param_name, param_info in parameters.items():
            value, type_str, step = param_info
            
            if type_str == "int":
                self.ocl_system.set_kernel_param(param_name, np.int32(value))
            elif type_str == "float":
                self.ocl_system.set_kernel_param(param_name, np.float32(value))
        
        # Re-bake kernels with current parameters
        if self.current_config:
            self.bake_kernels(self.current_config)

    def update_uniforms_from_controls(self):
        """Update parameters from widget controls (original purpose)."""
        if not self.current_config:
            return
            
        parameters = self.current_config.get("parameters", {})
        for param_name, widget in self.param_widgets.items():
            if param_name in parameters:
                param_info = parameters[param_name]
                value, type_str, step = param_info
                new_value = widget.value()
                parameters[param_name] = (new_value, type_str, step)
                
                if type_str == "int":
                    self.ocl_system.set_kernel_param(param_name, np.int32(new_value))
                elif type_str == "float":
                    self.ocl_system.set_kernel_param(param_name, np.float32(new_value))
        
        # Re-bake kernels with updated parameters
        if self.current_config:
            self.bake_kernels(self.current_config)

    def update_sim_uniforms(self):
        """Update simulation uniforms based on current parameters."""
        if not self.current_config:
            return
            
        parameters = self.current_config.get("parameters", {})
        for param_name, param_info in parameters.items():
            value, type_str, step = param_info
            
            if type_str == "int":
                self.ocl_system.set_kernel_param(param_name, np.int32(value))
            elif type_str == "float":
                self.ocl_system.set_kernel_param(param_name, np.float32(value))
        
        # Re-bake kernels with current parameters
        if self.current_config:
            self.bake_kernels(self.current_config)

    def _create_dynamic_buffer_registry(self):
        """Create dynamic buffer registry from user script data according to the pipeline approach."""
        #if not self.current_config:
        #    return
            
        print("Creating dynamic buffer registry using pipeline approach...")
        
        # Clear existing registry
        self.gl_objects.clear()
        self.render_pipeline_info.clear()
        self.buffer_data.clear()
        
        # 1. Get all buffer names used in render pipeline (both index buffers and vertex buffers)
        # to make sure we create only one buffer per unique name
        buffer_names = set()
        render_pipeline = self.current_config.get('render_pipeline', [])
        
        for render_pass in render_pipeline:
            # render_pass format: (shader_name, element_count_name, vertex_buffer_name, index_buffer_name)
            if len(render_pass) >= 3:
                shader_name, element_count_name, vertex_buffer_name, index_buffer_name = render_pass
                if vertex_buffer_name:
                    buffer_names.add(vertex_buffer_name)
                if index_buffer_name:
                    buffer_names.add(index_buffer_name)
        
        print(f"Found buffer names used in render pipeline: {list(buffer_names)}")
        
        # 2. Create all buffers in this set using the parameters from "buffers" 
        buffers_config = self.current_config.get('buffers', {})
        from .OGLsystem import GLobject
        
        for buffer_name in buffer_names:
            # Only create buffers that are actually used in render pipeline
            if buffer_name in buffers_config and hasattr(self, 'init_data') and self.init_data:
                data = self.init_data.get(buffer_name)
                if data is not None:
                    print(f"Creating buffer '{buffer_name}' with shape {data.shape}")
                    self.buffer_data[buffer_name] = {
                        'data': data.astype(np.float32),
                        'nelements': len(data),
                        'components': 1 if data.ndim == 1 else data.shape[1]
                    }
                else:
                    self.buffer_data[buffer_name] = {'data': None, 'nelements': 0, 'components': 3}
        
        # 3. Store render pipeline information for later use
        self.render_pipeline_info = render_pipeline
        
        # 4. Pass buffer data and render pipeline to GLCLWidget for deferred GLobject creation
        self.glcl_widget.set_buffer_data(self.buffer_data, self.render_pipeline_info)
        
        print(f"Dynamic buffer registry created with {len(self.buffer_data)} buffers")

    def _precompute_buffer_sync_list(self, config):
        """Pre-compute the list of buffer names that need GPU→CPU→GPU transfer.
        
        This method identifies the common subset of buffer names that appear in both
        the render pipeline and the kernel arguments. The complex logic is performed
        once during initialization, not every frame.
        """
        print(">>>> GLCLBrowser::_precompute_buffer_sync_list() Pre-computing buffer sync list...")
        
        # Reset the pre-computed list
        self.buffers_to_sync = []
        
        # Get buffers from render pipeline
        render_buffers = set()
        for render_pass in self.render_pipeline_info:
            if len(render_pass) >= 3 and render_pass[2]:  # vertex_buffer_name
                render_buffers.add(render_pass[2])
        
        # Get buffers from kernel arguments
        kernel_buffers = set()
        kernels_config = config.get("kernels", {})
        for kernel_name, kernel_info in kernels_config.items():
            if len(kernel_info) >= 3:  # local_size, global_size_expr, buffer_names, param_names
                buffer_names = kernel_info[2]
                kernel_buffers.update(buffer_names)
        
        # Find intersection - buffers that need sync
        buffers_to_sync = render_buffers.intersection(kernel_buffers)
        
        # Store as pre-computed list
        self.buffers_to_sync = list(buffers_to_sync)
        
        print(f"  Render pipeline buffers: {list(render_buffers)}")
        print(f"  Kernel buffers: {list(kernel_buffers)}")
        print(f"  Buffers to sync: {self.buffers_to_sync}")
        print("<<<< GLCLBrowser::_precompute_buffer_sync_list() DONE.")

    
    # Main entry point for the application
if __name__ == '__main__':
    import argparse
    
    # run like this:
    # __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia python -u -m pyBall.GLCL.GLCLBrowser
    # Enable OpenGL debug:
    # __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia python -u -m pyBall.GLCL.GLCLBrowser --debug-gl

    parser = argparse.ArgumentParser(description='GLCL Browser - Scientific Simulation Framework')
    parser.add_argument('--script', type=str, help='Path to simulation script')
    parser.add_argument('--debug-cl', action='store_true', help='Enable OpenCL debugging')
    parser.add_argument('--debug-gl', action='store_true', help='Enable OpenGL debugging')
    parser.add_argument('--debug-frames', type=int, default=1000000000, help='Number of debug frames to capture')
    
    args = parser.parse_args()
    
    app = QApplication(sys.argv)
    
    # Determine default script path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    default_script = os.path.join(current_dir, "scripts", "nbody.py")
    
    script_path = args.script or default_script

    ## Temporary hack for easier debugging
    #args.debug_cl = True
    #args.debug_frames = 10

    
    print(f"GLCLBrowser::__main__() Script path: {script_path}")
    print(f"GLCLBrowser::__main__() Current directory: {current_dir}")
    print(f"GLCLBrowser::__main__() Script exists: {os.path.exists(script_path)}")
    print(f"GLCLBrowser::__main__() OpenCL Debug: {args.debug_cl}")
    print(f"GLCLBrowser::__main__() OpenGL Debug: {args.debug_gl}")
    print(f"GLCLBrowser::__main__() Debug Frames: {args.debug_frames}")
    
    # Check if script exists and print its parameters
    if os.path.exists(script_path):
        import importlib.util
        spec = importlib.util.spec_from_file_location("simulation_script", script_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        
        config = getattr(module, 'config', None)
        if config:
            print("GLCLBrowser::__main__() Script config parameters:")
            params = config.get("parameters", {})
            for name, info in params.items():
                print(f"  {name}: {info}")
    
    # Create and show browser
    browser = GLCLBrowser(
        python_script_path=script_path, 
        bDebugCL=args.debug_cl, 
        bDebugGL=args.debug_gl, 
        nDebugFrames=args.debug_frames
    )
    browser.show()
    
    sys.exit(app.exec_())
