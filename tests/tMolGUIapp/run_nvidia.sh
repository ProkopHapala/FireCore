#!/bin/bash

# Force NVIDIA GPU with render offload
export __NV_PRIME_RENDER_OFFLOAD=1
export __GLX_VENDOR_LIBRARY_NAME=nvidia
export VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/nvidia_icd.json

# Disable multisampling to avoid GLX visual errors
# This overrides SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8) at runtime
export SDL_VIDEO_X11_VISUALID=

# Force X11 video driver
export SDL_VIDEODRIVER=x11

# Run the application with all arguments passed through
exec ./MolGUIapp "$@"
