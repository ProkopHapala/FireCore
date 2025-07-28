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
    
    def __init__(self, ocl_system, kernel_name, kernel_obj, global_size, local_size, arg_infos):
        self.ocl_system = ocl_system
        self.kernel_name = kernel_name
        self.kernel_obj = kernel_obj
        self.global_size = global_size  # Already-resolved tuple
        self.local_size = local_size
        self.arg_infos = arg_infos  # List of (param_name, is_buffer, type_hint)

    def execute(self):
        """Execute the kernel with current parameters."""
        try:
            global_size = self.global_size
            
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

            # Use OCLSystem's execute_kernel method with correct program name
            # Program name is derived from source file name (without extension)
            program_name = "nbody"  # This matches the actual .cl file name
            self.ocl_system.execute_kernel(
                program_name, self.kernel_name, global_size, self.local_size, kernel_args
            )
        except Exception as e:
            print(f"Error executing kernel {self.kernel_name}: {e}")
            self.on_exception(e)


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
        
        # Create main widget and layout
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        
        # Setup GLCL widget
        self.glcl_widget = GLCLWidget(self)
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
        self.current_config = config
        
        # Clear previous state
        self.baked_kernel_calls.clear()
        self.ocl_system.clear_buffers()
        self.param_widgets.clear()
        
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
        
        # Bake OpenCL kernels for efficient execution
        self.bake_kernels(config)
        
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
        kernels_config = config.get("kernels", {})
        pipeline = config.get("kernel_pipeline", [])
        parameters = config.get("parameters", {})
        
        print(f"GLCLBrowser::bake_kernels() Baking kernels with config: {config}")
        print(f"GLCLBrowser::bake_kernels() Parameters: {parameters}")
        print(f"GLCLBrowser::bake_kernels() Kernels: {kernels_config}")
        print(f"GLCLBrowser::bake_kernels() Pipeline: {pipeline}")
        
        self.baked_kernel_calls.clear()
        
        for kernel_name in pipeline:
            if kernel_name in kernels_config:
                kernel_info = kernels_config[kernel_name]
                print(f"GLCLBrowser::bake_kernels() Processing kernel {kernel_name} with info: {kernel_info}")
                local_size, global_size_expr, buffer_names, param_names = kernel_info
                
                # Build argument information
                arg_infos = []
                for buf_name in buffer_names:
                    arg_infos.append((buf_name, True, None))
                for param_name in param_names:
                    if param_name in parameters:
                        param_info = parameters[param_name]
                        arg_infos.append((param_name, False, param_info[1]))
                    else:
                        print(f"Warning: Parameter '{param_name}' not found, using float")
                        arg_infos.append((param_name, False, "float"))
                
                # Create global size resolver
                # Resolve global size once during baking
                resolved_gs = self._resolve_global_size(global_size_expr, parameters)
                baked_call = BakedKernelCall(
                    self.ocl_system, kernel_name, kernel_name,
                    resolved_gs, local_size, arg_infos
                )
                
                self.baked_kernel_calls.append(baked_call)

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
                
            if type_str == "int":
                # Create integer input
                from PyQt5.QtWidgets import QSpinBox, QLabel, QHBoxLayout, QWidget
                container = QWidget()
                layout = QHBoxLayout(container)
                layout.setContentsMargins(0, 0, 0, 0)
                
                label = QLabel(param_name)
                widget = QSpinBox()
                widget.setRange(-1000000, 1000000)
                widget.setValue(int(value))
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
                widget.setValue(float(value))
                widget.setSingleStep(float(step))
                widget.setDecimals(4)
                
                layout.addWidget(label)
                layout.addWidget(widget)
                self.params_layout.addWidget(container)
            else:
                continue
                
            self.param_widgets[param_name] = widget
            
            # Connect value changes
            if hasattr(widget, 'valueChanged'):
                if type_str == "int":
                    widget.valueChanged.connect(lambda v, p=param_name: self.ocl_system.set_kernel_param(p, np.int32(v)))
                elif type_str == "float":
                    widget.valueChanged.connect(lambda v, p=param_name: self.ocl_system.set_kernel_param(p, np.float32(v)))

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
        """Main simulation update loop."""
        self.debug_frame_counter += 1
        try:
            # Update parameters from GUI
            self.update_sim_uniforms()
            # Execute baked kernel calls unless GL debugging is on
            if not self.bDebugGL:
                for kernel_call in self.baked_kernel_calls:
                    kernel_call.execute()
                    # If CL debugging, print buffer snapshots
                    if self.bDebugCL:
                        for buf_name in self.buffer_shapes.keys():
                            host = np.empty(self.buffer_shapes[buf_name], dtype=np.float32)
                            self.ocl_system.fromGPU(buf_name, host)
                            print(f"[Frame {self.debug_frame_counter}] Buffer '{buf_name}':\n{host[:5] if host.size>20 else host}")
                
            # Update OpenGL visualization if not CL-only debug
            if not self.bDebugCL:
                self.glcl_widget.update()
            
            if self.debug_frame_counter >= self.nDebugFrames:
                self.simulation_running = False
                self.sim_timer.stop()
                print("GLCLBrowser::update_simulation() STOP! debug_frame_counter reached nDebugFrames: {self.nDebugFrames}")
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

    def update_sim_uniforms(self):
        """Update simulation uniforms based on current parameters."""
        #if not self.current_config:
        #    return
            
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

    
    # Main entry point for the application
if __name__ == '__main__':
    app = QApplication(sys.argv)
    
    # Determine default script path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    default_script = os.path.join(current_dir, "scripts", "nbody.py")
    
    # Create and show browser
    #browser = GLCLBrowser(python_script_path=default_script, bDebugCL=False, bDebugGL=True, nDebugFrames=5)
    browser = GLCLBrowser(python_script_path=default_script, nDebugFrames=1000 )
    browser.show()
    
    sys.exit(app.exec_())
