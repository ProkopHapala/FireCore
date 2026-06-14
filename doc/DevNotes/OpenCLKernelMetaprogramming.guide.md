# OpenCL Kernel Metaprogramming Guide

## 1. Introduction
This document describes a general approach for dynamic OpenCL kernel generation using Python-based textual substitution. The system enables:
- Runtime code substitution in kernel templates
- Generation of specialized kernel variants
- Debuggable kernel source output

(Example implementation shown using FireCore's molecular dynamics system)

## 2. Core Workflow

The general metaprogramming workflow:

1. **Template Preparation**: Create OpenCL kernel with substitution markers
   ```opencl
   // Example marker in template:
   result = //<<<CALCULATION>>>;
   ```

2. **Substitution Definition**: Define code variants in Python
   ```python
   substitutions = { 'macros': { 'CALCULATION': 'your_custom_code_here()' } }

3. **Preprocessing**: Generate specialized kernels
   ```python
   ocl.preprocess_opencl_source( 'template.cl', substitutions, 'specialized.cl' )
   ```

4. **Execution**: Compile and run generated kernels

## 3. Key Components (General)

### 3.1 Host Manager Class
Basic requirements for the host class:
- Inherits from opencl base class `pyBall/OCL/OpenCLBase.py`
- Provides preprocessing method
- Handles program compilation
- Manages kernel execution

(FireCore example: `MolecularDynamics` extends `OpenCLBase`)

### 3.2 Kernel Template
Template requirements:
- Valid OpenCL code with substitution markers
- Clear variable scoping
- Unique marker identifiers

(FireCore example: `relax_multi_mini.cl` uses `//<<<GET_FORCE_NONBOND`)

### 3.3 Test Harness
Typical harness functions:
- Defines substitution variants
- Manages test cases
- Benchmarks performance

(FireCore example: `run_scanNonBond.py`)

## 4. Implementation Patterns

### 4.1 Basic Substitution
```python
# Define substitution
expression = 'a*b + c'
subs = {'macros': {'CALCULATION': expression}}

# Apply to template
ocl.preprocess_opencl_source( 'math_kernel.cl', subs, 'math_kernel_specialized.cl' )
```

### 4.2 Multiple Substitutions
```python
subs = {
    'macros': {
        'INIT': 'float total = 0;',
        'LOOP': 'total += input[i];',
        'OUTPUT': 'result = total;'
    }
}
```

## 5. Best Practices

1. **Template Design**:
   - Use descriptive marker names
   - Document required variables
   - Keep templates modular

2. **Substitution Scope**:
   - Verify variable availability
   - Maintain consistent types
   - Check boundary conditions

3. **Debugging**:
   - Save generated kernels
   - Check compilation logs
   - Test incrementally

## 6. Case Study: FireCore Implementation

### Force Calculation Substitution
```opencl
// In relax_multi_mini.cl:
force = //<<<GET_FORCE_NONBOND;
```

### Python Substitution Definition
```python
# In run_scanNonBond.py:
potentials = {
    "Lennard-Jones": "compute_lj_force(r);",
    "Morse":         "compute_morse_force(r);"
}
```

### Execution Flow
```python
# For each potential:
md.preprocess_opencl_source(...)
md.load_program(...)
results = md.scanNonBond2(...)
```

## 7. Usage Examples

### 7.1 Adding a New Force Function
1. Add implementation to `Forces.cl`
2. Define substitution in `run_scanNonBond.py`
3. Add test case to potentials dictionary

### 7.2 Debugging Substitutions
1. Check generated `.cl` files in `tests/tmp/cl/`
2. Verify variable scope matches kernel requirements
3. Examine compilation logs

## 8. Debugging Tips

1. **Common Issues**:
   - Undefined variables in substitutions
   - Mismatched function signatures
   - Incorrect buffer sizes

2. **Debugging Tools**:
   - Examine generated `.cl` files
   - Check OpenCL compilation logs
   - Use `print` debugging in kernels

3. **Verification Steps**:
   - Compare against reference implementations
   - Test with simplified cases first
   - Validate energy conservation
