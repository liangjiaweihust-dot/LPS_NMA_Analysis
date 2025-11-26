# LPS-Protein Normal Mode Analysis (NMA) Tool

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive computational tool for analyzing protein-lipopolysaccharide (LPS) interactions using Normal Mode Analysis (NMA). This tool provides insights into protein dynamics, conformational changes, and autoinhibition mechanisms upon LPS binding.

## 🚀 Features

- **Normal Mode Analysis (NMA)**: Calculate protein flexibility and dominant motion modes
- **Free Energy Landscape (FEL)**: Generate 2D energy landscapes using NMA or PCA projections
- **Gate Openness Analysis**: Quantify protein gate dynamics and opening mechanisms
- **Autoinhibition Studies**: Analyze peptide-induced conformational restrictions
- **3D Visualizations**: Create animated trajectories and motion representations
- **Unified Coordinate Systems**: Ensure accurate comparative analysis
- **Publication-Ready Plots**: Generate Nature-style figures with proper formatting

## 📋 Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Examples](#usage-examples)
- [Analysis Workflows](#analysis-workflows)
- [Output Files](#output-files)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)

## 🛠️ Installation

### Prerequisites

- Python 3.8 or higher
- conda or pip package manager

### Option 1: Clone and Install

```bash
# Clone the repository
git clone https://github.com/yourusername/LPS_NMA_Analysis.git
cd LPS_NMA_Analysis

# Create conda environment (recommended)
conda create -n lps_nma python=3.9
conda activate lps_nma

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

### Option 2: Direct Installation

```bash
pip install lps-nma-analysis
```

### Verify Installation

```bash
python -c "import prody; print('ProDy version:', prody.__version__)"
python src/nma_analyzer.py --help
```

## 🚀 Quick Start

### Basic Analysis

```bash
# Run complete NMA analysis
python src/nma_analyzer.py --apo apo_structure.pdb --complex complex_structure.pdb --all --output results/

# Generate only free energy landscape
python src/nma_analyzer.py --apo apo.pdb --complex complex.pdb --fel --method nma --output fel_results/

# Calculate eigenvalues comparison
python src/nma_analyzer.py --apo apo.pdb --complex complex.pdb --eigenvalues --output eigenvals/
```

### Individual Analysis Scripts

```bash
# Free energy landscape with LPS marking
python src/create_free_energy_landscape_nma.py

# Gate openness visualization
python src/create_advanced_gate_visualization.py

# 3D trajectory animation
python src/create_3d_trajectory_animation.py

# Complete NMA analysis with all plots
python src/nma_analysis.py
```

## 📊 Usage Examples

### 1. Basic NMA Comparison

```bash
python src/nma_analyzer.py \
    --apo examples/apo_structure.pdb \
    --complex examples/complex_structure.pdb \
    --eigenvalues \
    --output basic_analysis/
```

**Output**: Eigenvalue comparison table and plots

### 2. Free Energy Landscape Analysis

```bash
python src/nma_analyzer.py \
    --apo examples/apo_structure.pdb \
    --complex examples/complex_structure.pdb \
    --fel \
    --method nma \
    --output fel_analysis/
```

**Output**: 2D free energy landscapes with unified axes

### 3. Complete Analysis Pipeline

```bash
python src/nma_analyzer.py \
    --apo examples/apo_structure.pdb \
    --complex examples/complex_structure.pdb \
    --all \
    --verbose \
    --output complete_analysis/
```

**Output**: All analyses including FEL, gate dynamics, and autoinhibition

### 4. Advanced Free Energy Landscape

```bash
python src/create_free_energy_landscape_nma.py
```

**Requirements**: 
- `complex.pdb` and `no-peptide.pdb` in current directory
- Generates comprehensive FEL with LPS position marking

## 🔬 Analysis Workflows

### Workflow 1: Basic Protein Dynamics Comparison

1. **Prepare structures**: Ensure clean PDB files with proper atom naming
2. **Run NMA analysis**: Calculate normal modes for both systems
3. **Compare eigenvalues**: Identify stiffness changes
4. **Generate plots**: Create publication-ready figures

```bash
# Step-by-step workflow
python src/nma_analyzer.py --apo apo.pdb --complex complex.pdb --eigenvalues
python src/nma_analyzer.py --apo apo.pdb --complex complex.pdb --fel --method nma
```

### Workflow 2: Autoinhibition Analysis

1. **NMA calculation**: Perform normal mode analysis
2. **FEL generation**: Create free energy landscapes
3. **Comparative analysis**: Quantify conformational restrictions
4. **Statistical analysis**: Generate autoinhibition metrics

```bash
python src/create_free_energy_landscape_nma.py
# Generates: FEL plots + autoinhibition analysis report
```

### Workflow 3: Gate Dynamics Study

1. **Identify gate regions**: Automatically detect beta-barrel gates
2. **Calculate openness**: Quantify gate opening dynamics
3. **Visualize motion**: Create 3D animations
4. **Compare systems**: Analyze opening/closing trends

```bash
python src/create_advanced_gate_visualization.py
python src/create_3d_trajectory_animation.py
```

## 📁 Output Files

### Standard Outputs

| File | Description |
|------|-------------|
| `eigenvalues.txt` | Eigenvalue comparison table |
| `eigenvalues_comparison.pdf` | Eigenvalue bar plots |
| `free_energy_landscape_nma.pdf` | NMA-based FEL with unified axes |
| `autoinhibition_analysis.txt` | Quantitative autoinhibition metrics |
| `gate_openness_advanced.pdf` | Gate dynamics visualization |

### Advanced Outputs

| File | Description |
|------|-------------|
| `free_energy_landscape_autoinhibition.pdf` | 4-panel autoinhibition comparison |
| `structure_trajectory_3d.gif` | 3D animated protein motion |
| `gate_openness_3d_animation.gif` | Gate opening animation |
| `rms_displacement_comparison.pdf` | RMS displacement analysis |

### Documentation Files

| File | Description |
|------|-------------|
| `PROJECTION_SPACE_EXPLANATION.md` | Theory behind projection coordinates |
| `UNIFIED_AXIS_EXPLANATION.md` | Importance of unified coordinate axes |
| `COORDINATE_AXIS_EXPLANATION.md` | Axis range differences explanation |

## ⚙️ Advanced Usage

### Custom Analysis Parameters

```python
from src.nma_analyzer import NMAAnalyzer

# Initialize with custom settings
analyzer = NMAAnalyzer(
    apo_pdb="apo.pdb",
    complex_pdb="complex.pdb",
    output_dir="custom_results"
)

# Load and analyze
analyzer.load_structures()
analyzer.perform_nma()

# Custom FEL generation
analyzer.generate_free_energy_landscape(method='nma')
```

### Batch Processing

```bash
# Process multiple structure pairs
for apo in apo_*.pdb; do
    complex=${apo/apo_/complex_}
    output_dir="results_$(basename $apo .pdb)"
    python src/nma_analyzer.py --apo $apo --complex $complex --all --output $output_dir
done
```

### Configuration File

Create `config.yaml`:

```yaml
analysis:
  n_modes: 20
  n_samples: 1000
  amplitude_range: 3.0

visualization:
  dpi: 300
  figure_format: 'pdf'
  colormap: 'nature'

output:
  save_frames: true
  generate_animations: true
```

## 🔧 Troubleshooting

### Common Issues

1. **ProDy Import Error**
   ```bash
   conda install -c conda-forge prody
   # or
   pip install prody
   ```

2. **Memory Issues with Large Proteins**
   ```python
   # Reduce sampling in large systems
   analyzer._sample_conformational_space_nma(..., n_samples=500)
   ```

3. **Missing Hetero Atoms**
   ```bash
   # Check PDB file for LPS/sugar residues
   grep "HETATM" your_structure.pdb | head -10
   ```

4. **Coordinate System Issues**
   - Ensure both structures use the same coordinate system
   - Check for missing atoms or chains
   - Verify proper PDB formatting

### Performance Optimization

```bash
# For large proteins (>5000 atoms)
export OMP_NUM_THREADS=4  # Limit CPU cores
ulimit -m 8000000         # Set memory limit (8GB)

# Run with reduced precision
python src/nma_analyzer.py --apo apo.pdb --complex complex.pdb --fast
```

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
git clone https://github.com/yourusername/LPS_NMA_Analysis.git
cd LPS_NMA_Analysis
pip install -e .[dev]
pre-commit install
```

### Running Tests

```bash
pytest tests/
pytest --cov=src tests/  # With coverage
```

## 📚 Scientific Background

### Normal Mode Analysis (NMA)

NMA is a computational method to study protein flexibility by:
1. Building a Hessian matrix of second derivatives
2. Calculating eigenvalues (λ) and eigenvectors (modes)
3. Analyzing low-frequency modes representing functional motions

### Free Energy Landscape (FEL)

The FEL is calculated using:
- **Projection coordinates**: q₁ = ΔR·u₁, q₂ = ΔR·u₂
- **Energy formula**: E = ½(λ₁q₁² + λ₂q₂²)
- **Conformational sampling**: Monte Carlo sampling along normal modes

### Autoinhibition Analysis

Quantifies peptide-induced conformational restrictions through:
- Eigenvalue increases (higher stiffness)
- Reduced accessible conformational space
- Energy barrier changes

## 📖 Citation

If you use this tool in your research, please cite:

```bibtex
@software{lps_nma_analysis,
  title={LPS-Protein Normal Mode Analysis Tool},
  author={Your Name and Collaborators},
  year={2024},
  url={https://github.com/yourusername/LPS_NMA_Analysis},
  version={1.0.0}
}
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🆘 Support

- **Email**: 929458347@qq.com

## 🙏 Acknowledgments

- ProDy development team for the excellent structural analysis library
- Scientific community for NMA methodology development
- Beta testers and contributors

---

**Note**: This tool is designed for research purposes. Please validate results with experimental data and consult relevant literature for interpretation guidelines.
