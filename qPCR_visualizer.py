"""
qPCR Visualization Generator
Generate plots with different colormaps (viridis, magma, plasma, inferno, cividis)
Useful for neurodivergent-friendly colormap selection
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def get_output_directory(output_dir=None):
    """Return the desktop output directory used for qPCR results and plots."""
    if output_dir is None:
        output_dir = Path.home() / "Desktop" / "qPCR_Ergebnisse"
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def save_figure_to_output(fig, output_filename, output_dir=None):
    """Save a matplotlib figure to the qPCR output directory and return the path."""
    target_path = Path(output_filename)
    if not target_path.is_absolute():
        target_path = get_output_directory(output_dir) / target_path

    target_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {target_path}")
    plt.close(fig)
    return target_path


def load_results_csv(csv_file):
    """Load previously calculated results from CSV."""
    return pd.read_csv(csv_file)


def plot_delta_ct_with_colormap(delta_ct_df, colormap='viridis', output_filename=None, sample_name=None):
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
    title = f'qPCR Delta Ct Values ({colormap.capitalize()} Colormap)'
    if sample_name:
        title = f'{sample_name} - {title}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    if 'Sample_Name' in delta_ct_df.columns:
        ax.set_xticklabels(delta_ct_df['Sample_Name'], rotation=45, ha='right')
    else:
        ax.set_xticklabels([f"S{int(n)}" for n in delta_ct_df['Sample_Number']], rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.5)
    
    plt.tight_layout()
    save_figure_to_output(fig, output_filename)


def plot_fold_change_with_colormap(fold_change_df, colormap='viridis', output_filename=None, sample_name=None):
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
    title = f'qPCR Fold Change Values ({colormap.capitalize()} Colormap)'
    if sample_name:
        title = f'{sample_name} - {title}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    if 'Sample_Name' in fold_change_df.columns:
        ax.set_xticklabels(fold_change_df['Sample_Name'], rotation=45, ha='right')
    else:
        ax.set_xticklabels([f"S{int(n)}" for n in fold_change_df['Sample_Number']], rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=1, color='red', linestyle='--', linewidth=2, alpha=0.5, label='No change')
    ax.legend()
    
    plt.tight_layout()
    save_figure_to_output(fig, output_filename)


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


def save_combined_publication_quality_plots(delta_ct_df, fold_change_df, output_filename, sample_name=None):
    """Save a combined publication-quality figure containing delta Ct and fold change plots."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 14), dpi=300, constrained_layout=True)

    x_pos = np.arange(len(delta_ct_df))
    sample_labels = delta_ct_df['Sample_Name'].tolist() if 'Sample_Name' in delta_ct_df.columns else [f"S{int(n)}" for n in delta_ct_df['Sample_Number']]

    # Delta Ct plot
    colors = plt.colormaps.get_cmap('viridis')(np.linspace(0, 1, len(delta_ct_df)))
    axes[0].bar(x_pos, delta_ct_df['Delta_Ct'],
                yerr=delta_ct_df['Delta_Ct_STD'],
                capsize=6, color=colors, edgecolor='black', linewidth=1.25, alpha=0.9)
    axes[0].set_title('qPCR ΔCt Analysis', fontsize=18, fontweight='bold')
    axes[0].set_ylabel('ΔCt (Target − Reference)', fontsize=14, fontweight='bold')
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(sample_labels, rotation=45, ha='right', fontsize=10)
    axes[0].grid(axis='y', alpha=0.25, linestyle='--')
    axes[0].axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.6)

    # Fold Change plot
    colors_fc = plt.colormaps.get_cmap('plasma')(np.linspace(0, 1, len(fold_change_df)))
    axes[1].bar(x_pos, fold_change_df['Fold_Change'],
                yerr=fold_change_df['Fold_Change_STD'],
                capsize=6, color=colors_fc, edgecolor='black', linewidth=1.25, alpha=0.9)
    axes[1].set_title('qPCR Fold Change (2^-ΔΔCt)', fontsize=18, fontweight='bold')
    axes[1].set_ylabel('Fold Change', fontsize=14, fontweight='bold')
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels([f"S{int(n)}" for n in fold_change_df['Sample_Number']], rotation=45, ha='right', fontsize=10)
    axes[1].grid(axis='y', alpha=0.25, linestyle='--')
    axes[1].axhline(y=1, color='red', linestyle='--', linewidth=2, alpha=0.6, label='No change')
    axes[1].legend(fontsize=11)

    for ax in axes:
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    title = 'qPCR Publication-Quality Summary'
    if sample_name:
        title = f'{sample_name} - {title}'
    fig.suptitle(title, fontsize=22, fontweight='bold')
    save_figure_to_output(fig, output_filename)


if __name__ == "__main__":
    # Generate plots with all colormaps
    generate_all_colormaps()
