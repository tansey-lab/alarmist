# Alarmist Installation Guide

## Quick Install (Recommended)

```bash
# 1. Activate bptf conda environment
conda activate bptf

# 2. Navigate to alarmist directory
cd /Users/jiayifan/tansey_lab/alarmist

# 3. Install in development mode
pip install -e .

# 4. Test installation
python -c "import alarmist; print(f'✅ Alarmist v{alarmist.__version__} installed')"

# 5. Test CLI commands
alarmist-patch-lri --help
alarmist-bptf --help
```

## Verify Installation

```bash
# Test Python API
python << 'EOF'
import alarmist as al

# Check version
print(f"Alarmist version: {al.__version__}")

# Check available functions
print(f"Available: {al.__all__}")

# Check BPTF
from alarmist.factorization.bptf import BPTF_AVAILABLE
print(f"BPTF available: {BPTF_AVAILABLE}")
EOF

# Test CLI commands are registered
which alarmist-patch-lri
which alarmist-bptf
which alarmist-bptf-viz
```

## Environment Requirements

The package is designed to work with the **bptf** conda environment which has:
- Python 3.9+
- numpy, pandas, scipy
- anndata, scanpy
- matplotlib, seaborn
- liana (for LRI databases)
- sparse, bptf (for matrix factorization)

## Alternative: Create New Environment

If you need a fresh environment:

```bash
# Create new environment
conda create -n alarmist python=3.9

# Activate
conda activate alarmist

# Install dependencies
pip install numpy pandas scipy scikit-learn anndata scanpy matplotlib seaborn liana h5py pyarrow tqdm networkx pillow

# Install BPTF (optional)
pip install git+https://github.com/aschein/bptf.git

# Install alarmist
cd /Users/jiayifan/tansey_lab/alarmist
pip install -e .
```

## Troubleshooting

### Import Error: NumPy version conflict
If you see "Numba needs NumPy 1.24 or less":
```bash
conda activate bptf
pip install "numpy<1.25"
```

### CLI Commands Not Found
If CLI commands are not found after installation:
```bash
# Reinstall with pip
pip install -e . --force-reinstall

# Or check pip installation location
pip show alarmist
```

### BPTF Not Available
If `BPTF_AVAILABLE = False`:
```bash
pip install git+https://github.com/aschein/bptf.git
```

## Next Steps

After installation, see:
- **TUTORIAL.md** - Usage examples
- **PACKAGE_README.md** - Package structure and features
