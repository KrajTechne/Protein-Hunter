#!/bin/bash
set -e

echo "🚀 Setting up Boltz Design Environment..."

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Miniconda/Anaconda first."
    exit 1
fi

# Create and activate environment
echo "📦 Creating conda environment..."
conda create -n boltz_ph python=3.10 -y
source $(conda info --base)/etc/profile.d/conda.sh
conda activate boltz_ph

# Install boltz
if [ -d "boltz" ]; then
    echo "📂 Installing Boltz..."
    cd boltz
    pip install -e .
    cd ..
else
    echo "❌ boltz directory not found. Please run this script from the project root."
    exit 1
fi
# Install conda dependencies
echo "🔧 Installing conda dependencies..."
conda install -c anaconda ipykernel -y

# Install Python dependencies
echo "🔧 Installing Python dependencies..."
pip install matplotlib seaborn prody tqdm PyYAML requests pypdb py3Dmol logmd==0.1.45
pip install ml_collections

# Install PyRosetta
echo "⏳ Installing PyRosetta (this may take a while)..."
pip install pyrosettacolabsetup pyrosetta-installer
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'


# ✅ Fix NumPy + Numba compatibility (PyRosetta downgrades NumPy to 1.23)
echo "🩹 Fixing NumPy/Numba version for Boltz and diffusion..."
pip install --upgrade "numpy>=1.24,<1.27" numba


# Download Boltz weights and dependencies
echo "⬇️  Downloading Boltz weights and dependencies..."
python -c "
from boltz.main import download_boltz2
from pathlib import Path
cache = Path('~/.boltz').expanduser()
cache.mkdir(parents=True, exist_ok=True)
download_boltz2(cache)
print('✅ Boltz weights downloaded successfully!')
"

# Setup LigandMPNN if directory exists
if [ -d "LigandMPNN" ]; then
    echo "🧬 Setting up LigandMPNN..."
    cd LigandMPNN
    bash get_model_params.sh "./model_params"
    cd ..
fi

# Make DAlphaBall.gcc executable
chmod +x "boltz/utils/DAlphaBall.gcc" || { echo -e "Error: Failed to chmod DAlphaBall.gcc"; exit 1; }

# Setup Jupyter kernel for the environment
echo "📓 Setting up Jupyter kernel..."
python -m ipykernel install --user --name=boltz_ph --display-name="Boltz Protein Hunter"

echo "🎉 Installation complete! Activate environment with: conda activate boltz_ph"