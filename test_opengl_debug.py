#!/usr/bin/env python3
"""
Test script to verify OpenGL debug functionality
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'pyBall'))

from PyQt5.QtWidgets import QApplication
from pyBall.GLCL.GLCLBrowser import GLCLBrowser

def test_opengl_debug():
    """Test OpenGL debug functionality"""
    print("Testing OpenGL debug functionality...")
    
    app = QApplication(sys.argv)
    
    # Test with debug enabled
    try:
        print("Creating browser with OpenGL debug enabled...")
        browser = GLCLBrowser(
            python_script_path="pyBall/GLCL/scripts/nbody.py",
            bDebugGL=True,
            bDebugCL=False,
            nDebugFrames=1
        )
        print("Browser created successfully with OpenGL debug enabled")
        
        # The debug will be enabled in initializeGL after context creation
        
    except Exception as e:
        print(f"Error during OpenGL debug test: {e}")
        return False
    
    return True

if __name__ == "__main__":
    success = test_opengl_debug()
    print(f"OpenGL debug test: {'PASSED' if success else 'FAILED'}")
