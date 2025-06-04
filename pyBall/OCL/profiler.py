import time
import numpy as np
import pyopencl as cl

class OpenCLProfiler:
    def __init__(self):
        # Timers for Python-level operations
        self.python_timers = {}
        self.python_active = {}
        
        # Timers for OpenCL operations
        self.cl_timers = {}
        self.cl_timer_counts = {}
        
        # Memory transfer statistics
        self.memory_transfers = {
            'host_to_device': 0,
            'device_to_host': 0
        }
        
    def start_python(self, name):
        """Start timing a Python operation"""
        self.python_active[name] = time.time()
        
    def stop_python(self, name):
        """Stop timing a Python operation"""
        if name in self.python_active:
            elapsed = time.time() - self.python_active[name]
            if name not in self.python_timers:
                self.python_timers[name] = 0
            self.python_timers[name] += elapsed
            del self.python_active[name]
            
    def record_cl_event(self, name, event):
        """Record time for an OpenCL event"""
        # Wait for the event to complete
        event.wait()
        
        # Calculate event execution time (nanoseconds -> seconds)
        start = event.get_profiling_info(cl.profiling_info.START)
        end = event.get_profiling_info(cl.profiling_info.END)
        elapsed = (end - start) * 1e-9  # Convert from nanoseconds to seconds
        
        # Store timing info
        if name not in self.cl_timers:
            self.cl_timers[name] = 0
            self.cl_timer_counts[name] = 0
        self.cl_timers[name] += elapsed
        self.cl_timer_counts[name] += 1
        
    def record_memory_transfer(self, direction, size_bytes):
        """Record memory transfer between host and device"""
        self.memory_transfers[direction] += size_bytes
    
    def reset(self):
        """Reset all timers"""
        self.python_timers = {}
        self.python_active = {}
        self.cl_timers = {}
        self.cl_timer_counts = {}
        self.memory_transfers = {
            'host_to_device': 0,
            'device_to_host': 0
        }
    
    def print_stats(self):
        """Print timing statistics"""
        print("\n===== OPENCL PROFILING STATISTICS =====\n")
        
        # Print Python operation timings
        if self.python_timers:
            print("Python Operations:")
            python_total = sum(self.python_timers.values())
            sorted_python = sorted(self.python_timers.items(), key=lambda x: x[1], reverse=True)
            for name, elapsed in sorted_python:
                print(f"  {name:30s}: {elapsed:.4f}s ({elapsed/python_total*100:.1f}%)")
            print(f"  {'Total Python Time':30s}: {python_total:.4f}s\n")
        
        # Print OpenCL kernel timings
        if self.cl_timers:
            print("OpenCL Kernels:")
            cl_total = sum(self.cl_timers.values())
            sorted_cl = sorted(self.cl_timers.items(), key=lambda x: x[1], reverse=True)
            for name, elapsed in sorted_cl:
                count = self.cl_timer_counts[name]
                avg = elapsed / count
                print(f"  {name:30s}: {elapsed:.4f}s ({elapsed/cl_total*100:.1f}%) - {count} calls, avg: {avg:.6f}s")
            print(f"  {'Total OpenCL Time':30s}: {cl_total:.4f}s\n")
        
        # Print memory transfer statistics
        print("Memory Transfers:")
        h2d = self.memory_transfers['host_to_device'] / (1024 * 1024)  # Convert to MB
        d2h = self.memory_transfers['device_to_host'] / (1024 * 1024)  # Convert to MB
        print(f"  {'Host to Device':30s}: {h2d:.2f} MB")
        print(f"  {'Device to Host':30s}: {d2h:.2f} MB")
        print(f"  {'Total Transfer':30s}: {h2d + d2h:.2f} MB\n")
        
        # Overall timing statistics
        if self.python_timers and self.cl_timers:
            overall_total = python_total + cl_total
            print("Overall Time Distribution:")
            print(f"  {'Python Overhead':30s}: {python_total:.4f}s ({python_total/overall_total*100:.1f}%)")
            print(f"  {'OpenCL Kernels':30s}: {cl_total:.4f}s ({cl_total/overall_total*100:.1f}%)")
            print(f"  {'Total Time':30s}: {overall_total:.4f}s")
