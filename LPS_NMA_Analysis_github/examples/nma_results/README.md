# Example Results

This directory would contain analysis results. Due to GitHub file size limits, large files (>25MB) are not included in this repository.

## Missing Files (Available in Full Release)

- `structure_trajectory_3d.gif` (13MB) - 3D protein trajectory animation
- `mode_trajectory_comparison.pdf` (13MB) - Mode trajectory comparison plots  
- `gate_openness_3d_animation.gif` (6.9MB) - Gate opening animation
- `gate_openness_3d_clear.gif` (8.7MB) - Clear gate dynamics visualization

## Available Files

- `free_energy_landscape_nma.pdf` (212KB) - Free energy landscape
- `eigenvalues_comparison.pdf` (643KB) - Eigenvalue comparison
- `autoinhibition_analysis.txt` - Quantitative analysis results
- Various other analysis outputs

## How to Generate Results

Run the analysis tools to generate all results:

```bash
# Generate all results
python src/nma_analyzer.py --apo examples/no-peptide.pdb --complex examples/complex.pdb --all --output results/

# Generate specific analyses
python src/create_free_energy_landscape_nma.py
python src/create_3d_trajectory_animation.py
```

## Download Full Results

For the complete results including large visualization files, download the full release package from the releases section.

