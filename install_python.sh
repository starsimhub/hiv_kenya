#!/bin/bash
# Install Python dependencies for the HIV Kenya model.
# Run from the repo root:  bash install_python.sh

set -e

echo "Installing Python dependencies..."
pip install -e ".[test]"
echo "Done. You can now run: python hiv_model.py"
