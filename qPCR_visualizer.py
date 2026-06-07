"""
qPCR Visualization Generator
Generate plots with different colormaps (viridis, magma, plasma, inferno, cividis)
Useful for neurodivergent-friendly colormap selection
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def load_results_csv(csv_file):
    """Load previously calculated results from CSV."""
    return pd.read_csv(csv_file)


def plot_delta_ct_with_colormap(delta_ct_df, colormap='viridis', output_filename=None):
    """
    Create bar plot for ΔCt values with specified colormap.
    
    Parameters:
    -----------
    delta_ct_df : pandas.DataFrame
        DataFrame with Delta Ct values
    colormap : str
        Colormap to use ('viridis', 'magma', 'plasma', 'inferno', 'cividis')
    output_filename : str, optional
        Filename to save the plot. If None, uses colormap name.
    """
    if output_filename is None:
        output_filename = f"qPCR_results_delta_ct_{colormap}.png"
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create color palette
    colors = plt.colormaps.get_cmap(colormap)(np.linspace(0, 1, len(delta_ct_df)))
    
    # Create bar plot with error bars
    x_pos = np.arange(len(delta_ct_df))
    bars = ax.bar(x_pos, delta_ct_df['Delta_Ct'], 
                   yerr=delta_ct_df['Delta_Ct_STD'],
                   capsize=5, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
    
    # Formatting
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('ΔCt (Target - Reference)', fontsize=12, fontweight='bold')
    ax.set_title(f'qPCR Delta Ct Values ({colormap.capitalize()} Colormap)', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    if 'Sample_Name' in delta_ct_df.columns:
        ax.set_xticklabels(delta_ct_df['Sample_Name'], rotation=45, ha='right')
    else:
        ax.set_xticklabels([f"S{int(n)}" for n in delta_ct_df['Sample_Number']], rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_filename}")
    plt.close()


def plot_fold_change_with_colormap(fold_change_df, colormap='viridis', output_filename=None):
    """
    Create bar plot for Fold Change values with specified colormap.
    
    Parameters:
    -----------
    fold_change_df : pandas.DataFrame
        DataFrame with Fold Change values
    colormap : str
        Colormap to use ('viridis', 'magma', 'plasma', 'inferno', 'cividis')
    output_filename : str, optional
        Filename to save the plot. If None, uses colormap name.
    """
    if output_filename is None:
        output_filename = f"qPCR_results_fold_change_{colormap}.png"
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create color palette
    colors = plt.colormaps.get_cmap(colormap)(np.linspace(0, 1, len(fold_change_df)))
    
    # Create bar plot with error bars
    x_pos = np.arange(len(fold_change_df))
    bars = ax.bar(x_pos, fold_change_df['Fold_Change'],
                   yerr=fold_change_df['Fold_Change_STD'],
                   capsize=5, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
    
    # Formatting
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Fold Change (2^-ΔΔCt)', fontsize=12, fontweight='bold')
    ax.set_title(f'qPCR Fold Change Values ({colormap.capitalize()} Colormap)', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    if 'Sample_Name' in fold_change_df.columns:
        ax.set_xticklabels(fold_change_df['Sample_Name'], rotation=45, ha='right')
    else:
        ax.set_xticklabels([f"S{int(n)}" for n in fold_change_df['Sample_Number']], rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=1, color='red', linestyle='--', linewidth=2, alpha=0.5, label='No change')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_filename}")
    plt.close()


def generate_all_colormaps():
    """
    Generate plots with all neurodivergent-friendly colormaps.
    """
    # Available neurodivergent-friendly colormaps
    colormaps = ['viridis', 'magma', 'plasma', 'inferno', 'cividis']
    
    print("\n" + "="*70)
    print("GENERATING PLOTS WITH ALL COLORMAPS")
    print("="*70)
    
    # Load results
    delta_ct_df = load_results_csv("qPCR_results_delta_ct.csv")
    fold_change_df = load_results_csv("qPCR_results_fold_change.csv")
    
    print(f"\nLoaded {len(delta_ct_df)} samples")
    print(f"Available colormaps: {', '.join(colormaps)}")
    print("\nGenerating visualizations...\n")
    
    # Generate plots for each colormap
    for colormap in colormaps:
        print(f"Generating {colormap}...")
        plot_delta_ct_with_colormap(delta_ct_df, colormap)
        plot_fold_change_with_colormap(fold_change_df, colormap)
    
    print("\n" + "="*70)
    print("✓ ALL VISUALIZATIONS COMPLETE")
    print("="*70)
    print("\nColormaps explained:")
    print("  • viridis  - Perceptually uniform, works for colorblind")
    print("  • magma    - Dark to bright, good contrast")
    print("  • plasma   - Vibrant, high contrast")
    print("  • inferno  - Dark to bright with warm tones")
    print("  • cividis - Optimized for colorblind vision (deuteranopia)")


if __name__ == "__main__":
    # Generate plots with all colormaps
    generate_all_colormaps()
