# LPS-NMA Analysis Usage Guide

## 📋 Complete Analysis Workflow

### Step 1: Prepare Your Data

1. **Structure Files**: Ensure you have clean PDB files
   - `apo_structure.pdb`: Protein without LPS/peptide
   - `complex_structure.pdb`: Protein with LPS/peptide bound

2. **File Requirements**:
   - Proper atom naming (PDB standard)
   - Complete protein chains
   - Hetero atoms (LPS/sugars) properly labeled
   - No missing residues in critical regions

### Step 2: Basic Analysis

```bash
# Navigate to your working directory
cd /path/to/your/project

# Run basic eigenvalue comparison
python /path/to/LPS_NMA_Analysis/src/nma_analyzer.py \
    --apo apo_structure.pdb \
    --complex complex_structure.pdb \
    --eigenvalues \
    --output basic_results/
```

**Expected Output**:
- `basic_results/eigenvalues.txt`: Numerical comparison
- `basic_results/eigenvalues_comparison.pdf`: Bar plot visualization

### Step 3: Free Energy Landscape

```bash
# Generate NMA-based free energy landscape
python /path/to/LPS_NMA_Analysis/src/nma_analyzer.py \
    --apo apo_structure.pdb \
    --complex complex_structure.pdb \
    --fel \
    --method nma \
    --output fel_results/
```

**Expected Output**:
- `fel_results/free_energy_landscape_nma.pdf`: 2D energy landscape with unified axes

### Step 4: Complete Analysis

```bash
# Run all analyses
python /path/to/LPS_NMA_Analysis/src/nma_analyzer.py \
    --apo apo_structure.pdb \
    --complex complex_structure.pdb \
    --all \
    --verbose \
    --output complete_analysis/
```

**Expected Output**:
- All eigenvalue and FEL files
- Autoinhibition analysis report
- Publication-ready figures

## 🔬 Advanced Analysis Scripts

### 1. Comprehensive Free Energy Landscape

```bash
# Copy your PDB files to the script directory
cp apo_structure.pdb /path/to/LPS_NMA_Analysis/src/no-peptide.pdb
cp complex_structure.pdb /path/to/LPS_NMA_Analysis/src/complex.pdb

# Run advanced FEL analysis
cd /path/to/LPS_NMA_Analysis/src/
python create_free_energy_landscape_nma.py
```

**Features**:
- LPS position marking on FEL
- Autoinhibition quantification
- Multiple comparison plots
- Unified coordinate axes

**Output Files**:
- `nma_results/free_energy_landscape_nma.pdf`
- `nma_results/free_energy_landscape_autoinhibition.pdf`
- `nma_results/autoinhibition_analysis.txt`

### 2. Gate Dynamics Analysis

```bash
cd /path/to/LPS_NMA_Analysis/src/
python create_advanced_gate_visualization.py
```

**Features**:
- Automatic gate region detection
- Gate openness quantification
- Multi-panel visualization
- 3D motion representation

### 3. 3D Trajectory Animation

```bash
cd /path/to/LPS_NMA_Analysis/src/
python create_3d_trajectory_animation.py
```

**Features**:
- Animated protein motion
- Gate tracking
- Dual-system comparison
- High-quality GIF output

### 4. Complete NMA Analysis

```bash
cd /path/to/LPS_NMA_Analysis/src/
python nma_analysis.py
```

**Features**:
- All standard NMA analyses
- RMS displacement comparison
- Key residue motion analysis
- Summary tables and plots

## 📊 Understanding the Output

### 1. Eigenvalues Table (`eigenvalues.txt`)

```
Mode    Apo         Complex     Difference  Ratio
1       0.172283    0.445373    0.273090    2.586
2       0.280410    0.617239    0.336829    2.201
3       0.341562    0.723891    0.382329    2.120
...
```

**Interpretation**:
- Higher eigenvalues = higher stiffness
- Ratio > 1 indicates complex is stiffer
- Large differences suggest significant conformational changes

### 2. Free Energy Landscape

**Coordinate System**:
- `q₁`: Projection along Mode 1 (easiest motion)
- `q₂`: Projection along Mode 2 (second easiest motion)
- Origin (0,0): Reference structure conformation

**Energy Contours**:
- Blue regions (Apo): Lower energy, more accessible
- Red regions (Complex): Higher energy, less accessible
- Elliptical shapes: Anisotropic motion (different stiffness in different directions)

### 3. LPS Position Markers

**Gold Diamonds**: LPS/sugar position in projection space
- Near origin: LPS inside/close to beta-barrel
- Far from origin: LPS at beta-barrel-NTD interface

### 4. Autoinhibition Analysis

**Key Metrics**:
- Eigenvalue increase: % increase in stiffness
- Accessible area: Conformational space reduction
- Energy basin width: Flexibility restriction

**Interpretation**:
- >50% eigenvalue increase: Strong autoinhibition
- Reduced accessible area: Conformational restriction
- Narrower energy basins: Reduced flexibility

## 🛠️ Customization Options

### 1. Modify Sampling Parameters

Edit the script to change sampling:

```python
# In create_free_energy_landscape_nma.py
n_samples = 2000  # Increase for better resolution
amplitude_range = 4.0  # Increase for wider sampling
```

### 2. Adjust Beta-Barrel Region

```python
# Modify get_beta_barrel_region function
start_fraction = 0.15  # Start at 15% instead of 20%
end_fraction = 0.85    # End at 85% instead of 80%
```

### 3. Change Visualization Settings

```python
# Modify plotting parameters
figsize = (20, 10)  # Larger figures
dpi = 600          # Higher resolution
levels = 30        # More contour levels
```

## 🔍 Troubleshooting Common Issues

### Issue 1: "No hetero atoms found"

**Solution**:
```bash
# Check for hetero atoms in your PDB
grep "HETATM" your_structure.pdb

# If missing, ensure LPS/sugars are properly included
# Check residue names: EA2, L01, LIG, PCJ, etc.
```

### Issue 2: "Eigenvalue calculation failed"

**Possible Causes**:
- Incomplete protein structure
- Missing atoms
- Coordinate system issues

**Solution**:
```bash
# Check structure integrity
python -c "
import prody
structure = prody.parsePDB('your_structure.pdb')
protein = structure.select('protein')
print(f'Protein atoms: {protein.numAtoms()}')
print(f'Residues: {protein.numResidues()}')
"
```

### Issue 3: "Memory error during NMA"

**Solution**:
```bash
# Reduce system size or use more memory
export OMP_NUM_THREADS=2  # Use fewer CPU cores
ulimit -m 16000000        # Increase memory limit

# Or select specific chains/regions
python -c "
import prody
structure = prody.parsePDB('large_structure.pdb')
protein = structure.select('protein and chain A')  # Select specific chain
prody.writePDB('reduced_structure.pdb', protein)
"
```

### Issue 4: "Coordinate axes not unified"

This has been fixed in the latest version. Ensure you're using the updated scripts.

## 📈 Performance Tips

### 1. For Large Proteins (>5000 atoms)

```bash
# Use reduced sampling
n_samples = 500

# Select important regions only
protein_core = structure.select('protein and resnum 50 to 300')
```

### 2. Batch Processing

```bash
#!/bin/bash
# batch_analysis.sh

for apo_file in apo_*.pdb; do
    complex_file=${apo_file/apo_/complex_}
    output_dir="results_$(basename $apo_file .pdb)"
    
    echo "Processing: $apo_file vs $complex_file"
    python src/nma_analyzer.py \
        --apo $apo_file \
        --complex $complex_file \
        --all \
        --output $output_dir
done
```

### 3. Parallel Processing

```python
# For multiple structure pairs
from multiprocessing import Pool
import subprocess

def run_analysis(file_pair):
    apo, complex_file = file_pair
    cmd = f"python src/nma_analyzer.py --apo {apo} --complex {complex_file} --all"
    subprocess.run(cmd, shell=True)

file_pairs = [("apo1.pdb", "complex1.pdb"), ("apo2.pdb", "complex2.pdb")]
with Pool(4) as p:
    p.map(run_analysis, file_pairs)
```

## 📝 Best Practices

### 1. Structure Preparation

- Always check PDB file quality
- Ensure consistent atom naming
- Remove water molecules if not needed
- Check for missing residues

### 2. Analysis Workflow

- Start with basic eigenvalue comparison
- Generate FEL for detailed analysis
- Use 3D visualizations for presentation
- Validate results with experimental data

### 3. Result Interpretation

- Compare eigenvalue ratios between systems
- Look for consistent trends across multiple modes
- Consider biological relevance of observed changes
- Validate with literature and experimental data

### 4. Publication Guidelines

- Always use unified coordinate axes for comparisons
- Include eigenvalue tables in supplementary material
- Provide clear figure captions explaining coordinate systems
- Cite relevant NMA methodology papers

## 🎯 Example Analysis Pipeline

```bash
#!/bin/bash
# complete_pipeline.sh

# 1. Basic setup
mkdir -p analysis_results
cd analysis_results

# 2. Copy structures
cp ../apo_structure.pdb ./
cp ../complex_structure.pdb ./

# 3. Basic analysis
python ../LPS_NMA_Analysis/src/nma_analyzer.py \
    --apo apo_structure.pdb \
    --complex complex_structure.pdb \
    --eigenvalues \
    --output eigenvalue_analysis/

# 4. Free energy landscape
python ../LPS_NMA_Analysis/src/nma_analyzer.py \
    --apo apo_structure.pdb \
    --complex complex_structure.pdb \
    --fel \
    --method nma \
    --output fel_analysis/

# 5. Advanced analysis (requires specific filenames)
cp apo_structure.pdb no-peptide.pdb
cp complex_structure.pdb complex.pdb
python ../LPS_NMA_Analysis/src/create_free_energy_landscape_nma.py

# 6. Gate dynamics
python ../LPS_NMA_Analysis/src/create_advanced_gate_visualization.py

# 7. 3D animations
python ../LPS_NMA_Analysis/src/create_3d_trajectory_animation.py

echo "Analysis complete! Check the following directories:"
echo "- eigenvalue_analysis/"
echo "- fel_analysis/"
echo "- nma_results/"
```

This pipeline will generate all major analyses and visualizations for your protein-LPS system.
