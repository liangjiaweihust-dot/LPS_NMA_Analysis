# Example Data and Usage

This directory contains example PDB structures and analysis results to help you get started with the LPS-NMA Analysis tool.

## 📁 Files Included

### Input Structures
- `complex.pdb`: Protein-LPS complex structure
- `no-peptide.pdb`: Apo protein structure (without LPS/peptide)

### Example Results
- `nma_results/`: Complete analysis results directory
  - `free_energy_landscape_nma.pdf`: Free energy landscape with unified axes
  - `autoinhibition_analysis.txt`: Quantitative autoinhibition metrics
  - `eigenvalues_comparison.pdf`: Eigenvalue comparison plots
  - Various other analysis outputs

## 🚀 Quick Start Examples

### 1. Basic Analysis

```bash
# Navigate to the LPS_NMA_Analysis directory
cd /path/to/LPS_NMA_Analysis

# Run basic eigenvalue comparison
python src/nma_analyzer.py \
    --apo examples/no-peptide.pdb \
    --complex examples/complex.pdb \
    --eigenvalues \
    --output test_results/
```

### 2. Free Energy Landscape

```bash
# Generate free energy landscape
python src/nma_analyzer.py \
    --apo examples/no-peptide.pdb \
    --complex examples/complex.pdb \
    --fel \
    --method nma \
    --output fel_test/
```

### 3. Complete Analysis

```bash
# Run all analyses
python src/nma_analyzer.py \
    --apo examples/no-peptide.pdb \
    --complex examples/complex.pdb \
    --all \
    --output complete_test/
```

### 4. Advanced Analysis (Using Original Scripts)

```bash
# Copy example files to src directory
cp examples/no-peptide.pdb src/
cp examples/complex.pdb src/

# Run advanced analysis
cd src/
python create_free_energy_landscape_nma.py
```

## 📊 Expected Results

After running the examples, you should see:

1. **Eigenvalue Analysis**:
   - Complex shows higher eigenvalues (increased stiffness)
   - Ratios > 2 for first few modes indicate significant stiffening

2. **Free Energy Landscape**:
   - Apo: Broader energy distribution (higher flexibility)
   - Complex: More restricted conformational space (autoinhibition)
   - LPS position marked with gold diamonds

3. **Autoinhibition Metrics**:
   - Eigenvalue increase: >100% (strong autoinhibition)
   - Reduced accessible conformational area
   - Narrower energy basins in complex

## 🔬 Understanding the Example System

### Biological Context
- **Apo structure**: Protein in its native, flexible state
- **Complex structure**: Protein with bound LPS/peptide showing autoinhibition
- **Key finding**: LPS binding increases protein stiffness and reduces gate opening

### Structural Features
- **Beta-barrel region**: Main structural domain analyzed
- **Gate regions**: Specific residues controlling protein opening/closing
- **LPS binding site**: Located at beta-barrel-NTD interface

### Analysis Insights
- **Mode 1 & 2**: Dominant motion modes showing stiffening upon LPS binding
- **Projection space**: (q₁, q₂) coordinates representing conformational changes
- **Energy landscape**: 2D representation of conformational accessibility

## 🛠️ Customizing for Your Data

### 1. Replace Example Structures

```bash
# Copy your structures to examples directory
cp your_apo_structure.pdb examples/apo_structure.pdb
cp your_complex_structure.pdb examples/complex_structure.pdb

# Run analysis with your data
python src/nma_analyzer.py \
    --apo examples/apo_structure.pdb \
    --complex examples/complex_structure.pdb \
    --all \
    --output your_results/
```

### 2. Modify Analysis Parameters

Edit the scripts to adjust:
- Sampling parameters (`n_samples`, `amplitude_range`)
- Beta-barrel region definition
- Visualization settings
- Output formats

### 3. Batch Processing Multiple Systems

```bash
# Create batch script
for i in {1..5}; do
    python src/nma_analyzer.py \
        --apo examples/apo_system_${i}.pdb \
        --complex examples/complex_system_${i}.pdb \
        --all \
        --output results_system_${i}/
done
```

## 📈 Performance Notes

- **Runtime**: ~5-15 minutes for typical protein systems (2000-5000 atoms)
- **Memory**: ~2-8 GB RAM depending on protein size
- **Output size**: ~50-200 MB including all plots and animations

## 🔍 Troubleshooting

If you encounter issues with the examples:

1. **Check dependencies**: Ensure all required packages are installed
2. **Verify file paths**: Make sure PDB files are accessible
3. **Check ProDy installation**: Test with `python -c "import prody"`
4. **Memory issues**: Reduce sampling parameters for large proteins

## 📚 Further Reading

- See `docs/USAGE_GUIDE.md` for detailed usage instructions
- Check `nma_results/PROJECTION_SPACE_EXPLANATION.md` for theoretical background
- Review `nma_results/UNIFIED_AXIS_EXPLANATION.md` for coordinate system details

## 🤝 Contributing Examples

To contribute additional example systems:

1. Ensure PDB files are clean and properly formatted
2. Include both apo and complex structures
3. Provide biological context and expected results
4. Test with the analysis pipeline
5. Submit via pull request with documentation
