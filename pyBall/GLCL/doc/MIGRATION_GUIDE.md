# Migration Guide: From Hardcoded to Modular OpenGL System

## Overview

This guide shows how to migrate from the hardcoded OpenGL system to the new modular system that provides maximum flexibility with minimal code.

## Current Problems

### Before (Hardcoded)
```python
# GLCLGUI_old.py - hardcoded buffer names
class GLCLWidget(QOpenGLWidget):
    def setup_particle_vbo(self):  # Hardcoded function name
        # Hardcoded buffer name "positions"
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glBufferData(GL_ARRAY_BUFFER, self.positions.nbytes, self.positions, GL_DYNAMIC_DRAW)
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, None)
        
    def update_particle_vbo(self, new_positions):  # Hardcoded again
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glBufferSubData(GL_ARRAY_BUFFER, 0, new_positions.nbytes, new_positions)
```

### After (Modular)
```python
# Generic buffer management - works with any name
render_obj.add_buffer("positions", positions, location=0, components=4)
render_obj.update_buffer("positions", new_positions)
```

## Key Improvements

1. **Generic Buffer Management**: No hardcoded buffer names
2. **Reusable Components**: Same code works for particles, atoms, bonds, etc.
3. **Factory Pattern**: Common patterns available as one-liners
4. **Resource Management**: Automatic cleanup
5. **Performance**: Efficient instanced rendering

## Migration Steps

### Step 1: Replace GLCLWidget

**Before:**
```python
class GLCLWidget(QOpenGLWidget):
    def setup_particle_vbo(self):
        # Hardcoded setup
        
    def update_particle_vbo(self, data):
        # Hardcoded update
```

**After:**
```python
from .ModularGL import RenderSystem, RenderFactory

class ModularGLCLWidget(QOpenGLWidget):
    def __init__(self):
        super().__init__()
        self.render_system = RenderSystem()
        
    def setup_render_objects(self, config):
        """Generic setup for any type of data"""
        positions = config.get('positions')
        if positions is not None:
            # Create particle system with generic name
            particles = RenderFactory.create_particle_system(
                config.get('name', 'particles'), positions
            )
            self.render_system.add_render_object('particles', particles)
```

### Step 2: Replace Buffer Updates

**Before:**
```python
# Hardcoded update
self.update_particle_vbo(new_positions)
```

**After:**
```python
# Generic update - works with any buffer name
self.render_objects['particles'].update_buffer('positions', new_positions)
```

### Step 3: Replace Shader Management

**Before:**
```python
# Hardcoded shader loading
self.shader_program = compile_shader_program(vertex_src, fragment_src)
```

**After:**
```python
# Modular shader management
shader = ShaderProgram("nbody", vertex_src, fragment_src)
shader.compile()
self.render_system.add_shader("nbody", shader)
```

## Usage Examples

### Example 1: N-Body Simulation
```python
# Before: 50+ lines of hardcoded setup
# After: 3 lines
positions = np.random.randn(2048, 4).astype(np.float32)
colors = np.random.rand(2048, 4).astype(np.float32)

widget.setup_particle_system(positions, colors)
```

### Example 2: Molecular Visualization
```python
# Create sphere instancer for atoms
atom_positions = mol.get_atom_positions()
widget.setup_sphere_instancer(atom_positions, radius=0.5)

# Create bond renderer
bond_vertices = mol.get_bond_vertices()
widget.setup_mesh_renderer(bond_vertices)
```

### Example 3: Mixed Rendering
```python
# Multiple systems can coexist
widget.setup_particle_system(positions1, colors1, name="electrons")
widget.setup_sphere_instancer(positions2, name="atoms")
widget.setup_mesh_renderer(vertices, name="bonds")
```

## API Reference

### Core Classes

- **RenderObject**: Generic render object for any data
- **InstancedRenderObject**: Efficient instanced rendering
- **ShaderProgram**: Shader management
- **RenderSystem**: High-level render management
- **RenderFactory**: Factory for common patterns

### Key Methods

```python
# Create objects
render_obj = RenderFactory.create_particle_system(name, positions, colors)
render_obj = RenderFactory.create_sphere_instancer(name, positions, radius)

# Update data
render_obj.update_buffer("positions", new_positions)
render_obj.update_buffer("colors", new_colors)

# Render
render_system.render_all()
```

## Migration Checklist

- [ ] Replace hardcoded buffer functions with generic ones
- [ ] Use RenderFactory for common patterns
- [ ] Replace hardcoded names with configurable ones
- [ ] Add proper resource cleanup
- [ ] Update all existing simulations to use new system
- [ ] Test with existing nbody.py simulation
- [ ] Verify performance is maintained or improved

## Backward Compatibility

The new system is designed to be backward compatible. You can:

1. **Gradual Migration**: Replace components one at a time
2. **Keep Existing Names**: Use same buffer names as before
3. **Mix Systems**: Use both old and new approaches temporarily

## Performance Benefits

- **Instanced Rendering**: Efficient for large particle systems
- **Resource Management**: Automatic cleanup prevents leaks
- **Memory Efficiency**: Shared buffers when possible
- **Minimal Overhead**: Thin abstraction layer

## Testing

Test the new system with:

```bash
# Test with existing simulations
python -m pyBall.GLCL.GLCLBrowser --script nbody.py

# Verify all buffer names work correctly
python -m pyBall.GLCL.GLCLBrowser --script your_simulation.py
```
