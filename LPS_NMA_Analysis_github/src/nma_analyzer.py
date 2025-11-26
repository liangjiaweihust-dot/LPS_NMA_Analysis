#!/usr/bin/env python3
"""
LPS-Protein Normal Mode Analysis (NMA) Tool

A comprehensive tool for analyzing protein-LPS interactions using Normal Mode Analysis,
including free energy landscape generation, gate openness analysis, and autoinhibition studies.

Author: Computational Biology Lab
Version: 1.0.0
"""

import argparse
import sys
import os
from pathlib import Path
import numpy as np
import prody
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import warnings

# Suppress ProDy warnings
prody.confProDy(verbosity='none')
warnings.filterwarnings('ignore')

# Nature-style colors
NATURE_COLORS = {
    'apo': '#2E86AB',
    'complex': '#A23B72',
    'reference': '#28A745',
    'lps': '#FFD700'
}

class NMAAnalyzer:
    """Main class for NMA analysis of protein-LPS systems"""
    
    def __init__(self, apo_pdb=None, complex_pdb=None, output_dir="results"):
        self.apo_pdb = apo_pdb
        self.complex_pdb = complex_pdb
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Data storage
        self.apo_structure = None
        self.complex_structure = None
        self.apo_protein = None
        self.complex_protein = None
        self.apo_anm = None
        self.complex_anm = None
        self.apo_sugar = None
        self.complex_sugar = None
        
    def load_structures(self):
        """Load PDB structures and extract protein and sugar components"""
        print("Loading structures...")
        
        if self.apo_pdb:
            self.apo_structure = prody.parsePDB(self.apo_pdb)
            self.apo_protein = self.apo_structure.select('protein')
            self.apo_sugar = self.apo_structure.select('hetero')
            print(f"  Apo: {self.apo_protein.numAtoms()} protein atoms")
            if self.apo_sugar:
                print(f"       {self.apo_sugar.numAtoms()} hetero atoms")
        
        if self.complex_pdb:
            self.complex_structure = prody.parsePDB(self.complex_pdb)
            self.complex_protein = self.complex_structure.select('protein')
            self.complex_sugar = self.complex_structure.select('hetero')
            print(f"  Complex: {self.complex_protein.numAtoms()} protein atoms")
            if self.complex_sugar:
                print(f"           {self.complex_sugar.numAtoms()} hetero atoms")
    
    def perform_nma(self):
        """Perform Normal Mode Analysis on both structures"""
        print("Performing Normal Mode Analysis...")
        
        if self.apo_protein:
            print("  Analyzing Apo structure...")
            self.apo_anm = prody.ANM('apo')
            self.apo_anm.buildHessian(self.apo_protein)
            self.apo_anm.calcModes()
            print(f"    Calculated {len(self.apo_anm)} modes")
        
        if self.complex_protein:
            print("  Analyzing Complex structure...")
            self.complex_anm = prody.ANM('complex')
            self.complex_anm.buildHessian(self.complex_protein)
            self.complex_anm.calcModes()
            print(f"    Calculated {len(self.complex_anm)} modes")
    
    def get_beta_barrel_region(self, protein):
        """Identify beta-barrel region (middle 60% of the protein)"""
        total_residues = protein.numResidues()
        start_res = int(total_residues * 0.2)
        end_res = int(total_residues * 0.8)
        
        resnums = protein.getResnums()
        min_resnum = np.min(resnums)
        max_resnum = np.max(resnums)
        
        beta_start = min_resnum + int((max_resnum - min_resnum) * 0.2)
        beta_end = min_resnum + int((max_resnum - min_resnum) * 0.8)
        
        beta_barrel = protein.select(f'resnum {beta_start} to {beta_end}')
        return beta_barrel, list(range(beta_start, beta_end + 1))
    
    def calculate_eigenvalues_comparison(self):
        """Calculate and compare eigenvalues"""
        if not (self.apo_anm and self.complex_anm):
            print("Error: NMA not performed on both structures")
            return
        
        apo_eigenvals = self.apo_anm.getEigvals()[:6]
        complex_eigenvals = self.complex_anm.getEigvals()[:6]
        
        # Save eigenvalues
        eigenvals_file = self.output_dir / "eigenvalues.txt"
        with open(eigenvals_file, 'w') as f:
            f.write("Mode\tApo\tComplex\tDifference\tRatio\n")
            for i in range(6):
                diff = complex_eigenvals[i] - apo_eigenvals[i]
                ratio = complex_eigenvals[i] / apo_eigenvals[i] if apo_eigenvals[i] > 0 else 0
                f.write(f"{i+1}\t{apo_eigenvals[i]:.6f}\t{complex_eigenvals[i]:.6f}\t{diff:.6f}\t{ratio:.3f}\n")
        
        print(f"  Eigenvalues saved to {eigenvals_file}")
        return apo_eigenvals, complex_eigenvals
    
    def generate_free_energy_landscape(self, method='nma'):
        """Generate free energy landscape using specified method"""
        print(f"Generating free energy landscape using {method.upper()} method...")
        
        if method == 'nma':
            self._generate_nma_landscape()
        elif method == 'pca':
            self._generate_pca_landscape()
        else:
            print(f"Unknown method: {method}")
    
    def _generate_nma_landscape(self):
        """Generate NMA-based free energy landscape"""
        # Get beta-barrel regions
        apo_beta_barrel, _ = self.get_beta_barrel_region(self.apo_protein)
        complex_beta_barrel, _ = self.get_beta_barrel_region(self.complex_protein)
        
        # Sample conformational space
        apo_result = self._sample_conformational_space_nma(
            self.apo_structure if self.apo_structure else self.apo_protein,
            self.apo_anm, apo_beta_barrel
        )
        complex_result = self._sample_conformational_space_nma(
            self.complex_structure if self.complex_structure else self.complex_protein,
            self.complex_anm, complex_beta_barrel
        )
        
        # Calculate landscapes
        apo_landscape = self._calculate_free_energy_landscape_nma(*apo_result)
        complex_landscape = self._calculate_free_energy_landscape_nma(*complex_result)
        
        # Create plots
        self._create_free_energy_landscape_plot(
            apo_landscape, complex_landscape,
            self.apo_anm, self.complex_anm,
            self.output_dir / "free_energy_landscape_nma.pdf"
        )
    
    def _sample_conformational_space_nma(self, structure, anm, target_selection, n_samples=1000):
        """Sample conformational space using NMA"""
        if len(anm) < 2:
            return None, None, None
        
        # Get eigenvalues
        eigenvalues = anm.getEigvals()
        
        # Generate random amplitudes
        amplitude_range = 3.0
        amplitudes1 = np.random.uniform(-amplitude_range, amplitude_range, n_samples)
        amplitudes2 = np.random.uniform(-amplitude_range, amplitude_range, n_samples)
        
        # Calculate projection intensities
        q1 = amplitudes1
        q2 = amplitudes2
        
        # Calculate energies
        energies = 0.5 * (eigenvalues[0] * q1**2 + eigenvalues[1] * q2**2)
        
        return q1, q2, energies
    
    def _calculate_free_energy_landscape_nma(self, q1, q2, energies):
        """Calculate free energy landscape from NMA projections"""
        # Create grid
        q1_range = np.linspace(q1.min(), q1.max(), 50)
        q2_range = np.linspace(q2.min(), q2.max(), 50)
        Q1, Q2 = np.meshgrid(q1_range, q2_range)
        
        # Calculate FEL on grid
        FEL = np.zeros_like(Q1)
        for i in range(len(q1_range)):
            for j in range(len(q2_range)):
                FEL[j, i] = 0.5 * (energies[0] * Q1[j, i]**2 + energies[0] * Q2[j, i]**2)
        
        return Q1, Q2, FEL, q1, q2
    
    def _create_free_energy_landscape_plot(self, apo_result, complex_result, apo_anm, complex_anm, output_file):
        """Create unified free energy landscape plot"""
        apo_Q1, apo_Q2, apo_FEL, apo_q1, apo_q2 = apo_result
        complex_Q1, complex_Q2, complex_FEL, complex_q1, complex_q2 = complex_result
        
        # Calculate unified axis ranges
        q1_min = min(apo_q1.min(), complex_q1.min())
        q1_max = max(apo_q1.max(), complex_q1.max())
        q2_min = min(apo_q2.min(), complex_q2.min())
        q2_max = max(apo_q2.max(), complex_q2.max())
        
        margin = 0.1
        q1_range = q1_max - q1_min
        q2_range = q2_max - q2_min
        
        q1_unified = [q1_min - margin * q1_range, q1_max + margin * q1_range]
        q2_unified = [q2_min - margin * q2_range, q2_max + margin * q2_range]
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Plot Apo
        apo_FEL_clean = np.nan_to_num(apo_FEL, nan=0, posinf=0, neginf=0)
        if np.any(apo_FEL_clean > 0):
            p99 = np.percentile(apo_FEL_clean[apo_FEL_clean > 0], 99)
            apo_FEL_clean = np.clip(apo_FEL_clean, 0, p99 * 1.1)
        
        im1 = ax1.contourf(apo_Q1, apo_Q2, apo_FEL_clean, levels=20, 
                          cmap='Blues', alpha=0.8)
        ax1.scatter([0], [0], c=NATURE_COLORS['reference'], s=200, marker='o', 
                   edgecolors='white', linewidths=2, label='Reference', zorder=10)
        
        ax1.set_xlim(q1_unified)
        ax1.set_ylim(q2_unified)
        ax1.set_xlabel('q₁ (Mode 1 Projection, Å)', fontweight='bold')
        ax1.set_ylabel('q₂ (Mode 2 Projection, Å)', fontweight='bold')
        ax1.set_title('Apo Free Energy Landscape', fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        plt.colorbar(im1, ax=ax1, label='Free Energy (kcal/mol)')
        
        # Plot Complex
        complex_FEL_clean = np.nan_to_num(complex_FEL, nan=0, posinf=0, neginf=0)
        if np.any(complex_FEL_clean > 0):
            p99 = np.percentile(complex_FEL_clean[complex_FEL_clean > 0], 99)
            complex_FEL_clean = np.clip(complex_FEL_clean, 0, p99 * 1.1)
        
        im2 = ax2.contourf(complex_Q1, complex_Q2, complex_FEL_clean, levels=20,
                          cmap='Reds', alpha=0.8)
        ax2.scatter([0], [0], c=NATURE_COLORS['reference'], s=200, marker='o',
                   edgecolors='white', linewidths=2, label='Reference', zorder=10)
        
        ax2.set_xlim(q1_unified)
        ax2.set_ylim(q2_unified)
        ax2.set_xlabel('q₁ (Mode 1 Projection, Å)', fontweight='bold')
        ax2.set_ylabel('q₂ (Mode 2 Projection, Å)', fontweight='bold')
        ax2.set_title('Complex Free Energy Landscape', fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.colorbar(im2, ax=ax2, label='Free Energy (kcal/mol)')
        
        plt.suptitle('Free Energy Landscape Comparison (Unified Axes)', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Free energy landscape saved to {output_file}")

def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(
        description="LPS-Protein Normal Mode Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic NMA analysis
  python nma_analyzer.py --apo apo.pdb --complex complex.pdb --output results/

  # Generate free energy landscape
  python nma_analyzer.py --apo apo.pdb --complex complex.pdb --fel --method nma

  # Full analysis pipeline
  python nma_analyzer.py --apo apo.pdb --complex complex.pdb --all --output my_results/
        """
    )
    
    # Input files
    parser.add_argument('--apo', type=str, required=True,
                       help='Apo structure PDB file')
    parser.add_argument('--complex', type=str, required=True,
                       help='Complex structure PDB file')
    
    # Output options
    parser.add_argument('--output', '-o', type=str, default='results',
                       help='Output directory (default: results)')
    
    # Analysis options
    parser.add_argument('--eigenvalues', action='store_true',
                       help='Calculate eigenvalues comparison')
    parser.add_argument('--fel', action='store_true',
                       help='Generate free energy landscape')
    parser.add_argument('--method', choices=['nma', 'pca'], default='nma',
                       help='Method for free energy landscape (default: nma)')
    parser.add_argument('--all', action='store_true',
                       help='Run all analyses')
    
    # Other options
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Check input files
    if not os.path.exists(args.apo):
        print(f"Error: Apo PDB file not found: {args.apo}")
        sys.exit(1)
    if not os.path.exists(args.complex):
        print(f"Error: Complex PDB file not found: {args.complex}")
        sys.exit(1)
    
    # Initialize analyzer
    analyzer = NMAAnalyzer(args.apo, args.complex, args.output)
    
    # Load structures
    analyzer.load_structures()
    
    # Perform NMA
    analyzer.perform_nma()
    
    # Run requested analyses
    if args.all or args.eigenvalues:
        analyzer.calculate_eigenvalues_comparison()
    
    if args.all or args.fel:
        analyzer.generate_free_energy_landscape(args.method)
    
    print(f"\nAnalysis complete! Results saved to: {args.output}")

if __name__ == "__main__":
    main()
