#!/bin/bash

# Exit on error
set -e

echo "Cleaning up C++ extension..."

# Close VS Code if it's running
echo "Please close VS Code before continuing."
read -p "Press Enter when VS Code is closed..."

# Remove the C++ extensions
echo "Removing C++ extensions..."
rm -rf ~/.vscode/extensions/ms-vscode.cpptools-*

# Clear the extension cache
echo "Clearing extension cache..."
rm -rf ~/.config/Code/CachedExtensionVSIXs/*cpptools*

# Create a backup of the C++ configuration
echo "Backing up C++ configuration..."
if [ -f ~/.vscode/c_cpp_properties.json ]; then
    cp ~/.vscode/c_cpp_properties.json ~/.vscode/c_cpp_properties.json.bak
fi

echo "Cleanup complete. Please restart VS Code and reinstall the C++ extension from the marketplace."
echo "After reinstalling, reload the window (Ctrl+Shift+P, then 'Reload Window')."
