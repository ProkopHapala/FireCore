#!/bin/bash
# Compile and install OpenMPI 4.1.8 from source
# This matches the version ORCA 6.1.1 was built against

set -e

INSTALL_DIR="$HOME/sw/openmpi-418"
BUILD_DIR="$HOME/tmp/openmpi-418-build"
VERSION="4.1.8"

echo "========================================"
echo "Installing OpenMPI $VERSION"
echo "Install dir: $INSTALL_DIR"
echo "========================================"

# Create directories
mkdir -p "$INSTALL_DIR"
mkdir -p "$BUILD_DIR"

# Download if not already present
cd "$BUILD_DIR"
if [ ! -f "openmpi-$VERSION.tar.gz" ]; then
    echo "Downloading OpenMPI $VERSION..."
    wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-$VERSION.tar.gz
fi

# Extract
if [ ! -d "openmpi-$VERSION" ]; then
    echo "Extracting..."
    tar -xzf openmpi-$VERSION.tar.gz
fi

# Configure
cd openmpi-$VERSION
echo "Configuring..."
./configure --prefix="$INSTALL_DIR" \
    --enable-mpi-thread-multiple \
    --enable-shared \
    --disable-static

# Build
echo "Building (this may take 10-20 minutes)..."
make -j$(nproc)

# Install
echo "Installing..."
make install

# Add to .bashrc if not already there
BASHRC="$HOME/.bashrc"
if ! grep -q "$INSTALL_DIR/bin" "$BASHRC"; then
    echo ""
    echo "Adding OpenMPI to PATH in ~/.bashrc..."
    echo "" >> "$BASHRC"
    echo "# OpenMPI 4.1.8 for ORCA" >> "$BASHRC"
    echo "export PATH=\"$INSTALL_DIR/bin:\$PATH\"" >> "$BASHRC"
    echo "export LD_LIBRARY_PATH=\"$INSTALL_DIR/lib:\$LD_LIBRARY_PATH\"" >> "$BASHRC"
    echo "export PKG_CONFIG_PATH=\"$INSTALL_DIR/lib/pkgconfig:\$PKG_CONFIG_PATH\"" >> "$BASHRC"
fi

echo ""
echo "========================================"
echo "OpenMPI $VERSION installed successfully!"
echo "========================================"
echo ""
echo "Source your ~/.bashrc or run:"
echo "  export PATH=\"$INSTALL_DIR/bin:\$PATH\""
echo "  export LD_LIBRARY_PATH=\"$INSTALL_DIR/lib:\$LD_LIBRARY_PATH\""
echo ""
echo "Then verify with:"
echo "  mpirun --version"
echo ""
