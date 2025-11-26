# 🚀 Quick Start Guide

## 📦 Installation

```bash
# 1. Clone the repository
git clone https://github.com/yourusername/LPS_NMA_Analysis.git
cd LPS_NMA_Analysis

# 2. Create conda environment
conda create -n lps_nma python=3.9
conda activate lps_nma

# 3. Install dependencies
pip install -r requirements.txt
```

## ⚡ 5-Minute Analysis

### Option 1: Using the Main Analyzer

```bash
# Basic analysis with your PDB files
python src/nma_analyzer.py \
    --apo your_apo_structure.pdb \
    --complex your_complex_structure.pdb \
    --all \
    --output results/
```

### Option 2: Using Example Data

```bash
# Test with provided examples
python src/nma_analyzer.py \
    --apo examples/no-peptide.pdb \
    --complex examples/complex.pdb \
    --all \
    --output test_results/
```

### Option 3: Advanced Analysis

```bash
# Copy your files (must be named exactly as shown)
cp your_apo_structure.pdb src/no-peptide.pdb
cp your_complex_structure.pdb src/complex.pdb

# Run comprehensive analysis
cd src/
python create_free_energy_landscape_nma.py
```

## 📊 What You'll Get

### Key Output Files:
- `eigenvalues_comparison.pdf` - Stiffness comparison
- `free_energy_landscape_nma.pdf` - 2D energy landscapes
- `autoinhibition_analysis.txt` - Quantitative metrics
- `gate_openness_advanced.pdf` - Gate dynamics
- `structure_trajectory_3d.gif` - 3D animations

### Analysis Results:
- **Eigenvalue ratios > 2**: Strong stiffening effect
- **Reduced conformational space**: Autoinhibition evidence
- **LPS position marking**: Spatial relationship analysis
- **Unified coordinate axes**: Accurate comparison

## 🎯 Command Line Options

```bash
# Basic eigenvalue comparison only
python src/nma_analyzer.py --apo apo.pdb --complex complex.pdb --eigenvalues

# Free energy landscape only
python src/nma_analyzer.py --apo apo.pdb --complex complex.pdb --fel --method nma

# All analyses with verbose output
python src/nma_analyzer.py --apo apo.pdb --complex complex.pdb --all --verbose

# Custom output directory
python src/nma_analyzer.py --apo apo.pdb --complex complex.pdb --all --output my_results/
```

## 🔧 Individual Analysis Scripts

```bash
# Comprehensive free energy landscape (requires specific filenames)
cd src/
python create_free_energy_landscape_nma.py

# Gate dynamics visualization
python create_advanced_gate_visualization.py

# 3D trajectory animations
python create_3d_trajectory_animation.py

# Complete NMA analysis suite
python nma_analysis.py
```

## 📁 Project Structure

```
LPS_NMA_Analysis/
├── src/                          # Main analysis scripts
│   ├── nma_analyzer.py          # Main CLI tool
│   ├── create_free_energy_landscape_nma.py
│   ├── create_advanced_gate_visualization.py
│   ├── create_3d_trajectory_animation.py
│   └── nma_analysis.py
├── examples/                     # Example data and results
│   ├── complex.pdb
│   ├── no-peptide.pdb
│   └── nma_results/
├── docs/                        # Documentation
│   └── USAGE_GUIDE.md
├── README.md                    # Main documentation
├── requirements.txt             # Dependencies
└── setup.py                    # Installation script
```

## 🎨 Customization

### Modify Analysis Parameters:
```python
# Edit in the scripts
n_samples = 1000        # Increase for better resolution
amplitude_range = 3.0   # Adjust sampling range
dpi = 300              # Change plot resolution
```

### Select Specific Protein Regions:
```python
# Focus on specific chains or residues
protein_selection = structure.select('protein and chain A and resnum 50 to 300')
```

## 🔍 Troubleshooting

### Common Issues:
1. **ProDy not found**: `pip install prody`
2. **Memory error**: Reduce `n_samples` or use smaller protein regions
3. **No hetero atoms**: Check PDB file contains LPS/sugar residues
4. **File not found**: Ensure correct file paths and names

### Performance Tips:
- Use `--verbose` flag to monitor progress
- Reduce sampling for large proteins (>5000 atoms)
- Run on systems with >8GB RAM for best performance

## 📚 Next Steps

1. **Read the full documentation**: `README.md`
2. **Check detailed usage guide**: `docs/USAGE_GUIDE.md`
3. **Understand the theory**: `examples/nma_results/PROJECTION_SPACE_EXPLANATION.md`
4. **Explore examples**: `examples/README.md`

## 🆘 Getting Help

- **Issues**: Check troubleshooting section in `README.md`
- **Examples**: Use provided example data to test installation
- **Theory**: Read explanation documents in `examples/nma_results/`

---

**Ready to analyze your protein-LPS system? Start with the examples and then use your own data!** 🧬
