"""
qPCR Delta CT Calculation and Visualization
Calculates ΔCt (Delta Ct) = Target Ct - Reference Ct
and creates accessible visualizations with neurodivergent-friendly colormaps
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


# Import the data loader module
import sys
sys.path.insert(0, str(Path(__file__).parent))
from qPCR_data_loader import load_qpcr_data, separate_genes, get_sample_groups


def calculate_delta_ct(target_df, reference_df):
    """
    Calculate ΔCt = Target Ct - Reference Ct
    
    Parameters:
    -----------
    target_df : pandas.DataFrame
        Target gene data (e.g., gene of interest)
    reference_df : pandas.DataFrame
        Reference gene data (e.g., housekeeping gene)
    
    Returns:
    --------
    delta_ct_df : pandas.DataFrame
        DataFrame with calculated ΔCt values
    """
    target_df = target_df.copy()
    reference_df = reference_df.copy()
    
    # Use parsed sample numbers from the loader if available
    if 'Sample_Number' not in target_df:
        target_df['Sample_Number'] = target_df['Well1'].str.extract(r'(\d+)$')[0].astype(int)
    if 'Sample_Number' not in reference_df:
        reference_df['Sample_Number'] = reference_df['Well1'].str.extract(r'(\d+)$')[0].astype(int)
    
    target_df = target_df.sort_values('Sample_Number').reset_index(drop=True)
    reference_df = reference_df.sort_values('Sample_Number').reset_index(drop=True)

    n = min(len(target_df), len(reference_df))
    if n == 0:
        return pd.DataFrame(columns=['Sample_Number', 'Sample_Name', 'Day', 'Target_Ct', 'Reference_Ct', 'Target_STD', 'Reference_STD', 'Delta_Ct', 'Delta_Ct_STD'])

    target_df = target_df.iloc[:n].copy()
    reference_df = reference_df.iloc[:n].copy()
    
    delta_ct_df = pd.DataFrame({
        'Sample_Number': target_df['Sample_Number'],
        'Sample_Name': target_df.get('Sample_Name', target_df['Well1'] + '_' + target_df['Well2']),
        'Day': target_df.get('Day', target_df['Sample_Name'].str.extract(r'^(D[12])', expand=False).fillna('Unknown')),
        'Target_Ct': target_df['MeanCp'].values,
        'Reference_Ct': reference_df['MeanCp'].values,
        'Target_STD': target_df['STD_Cp'].values,
        'Reference_STD': reference_df['STD_Cp'].values,
    })
    
    valid = (delta_ct_df['Target_Ct'] > 0) & (delta_ct_df['Reference_Ct'] > 0)
    delta_ct_df.loc[valid, 'Delta_Ct'] = delta_ct_df.loc[valid, 'Target_Ct'] - delta_ct_df.loc[valid, 'Reference_Ct']
    delta_ct_df.loc[~valid, 'Delta_Ct'] = np.nan
    
    delta_ct_df['Delta_Ct_STD'] = np.sqrt(
        delta_ct_df['Target_STD']**2 + delta_ct_df['Reference_STD']**2
    )
    delta_ct_df.loc[~valid, 'Delta_Ct_STD'] = np.nan
    
    return delta_ct_df


def calculate_fold_change(delta_ct_df, reference_sample=None):
    """
    Calculate Fold Change = 2^(-ΔΔCt)
    
    Parameters:
    -----------
    delta_ct_df : pandas.DataFrame
        DataFrame with ΔCt values
    reference_sample : int or str, optional
        Sample number or sample name to use as control.
        If None, uses DMSO control per day when available.
    
    Returns:
    --------
    fold_change_df : pandas.DataFrame
        DataFrame with calculated fold change values
    """
    fold_change_df = delta_ct_df.copy()
    if 'Day' not in fold_change_df:
        fold_change_df['Day'] = fold_change_df['Sample_Name'].str.extract(r'^(D[12])', expand=False).fillna('Unknown')
    
    if reference_sample is None:
        dmso_mask = fold_change_df['Sample_Name'].str.contains('DMSO', case=False, na=False)
        dmso_baseline = fold_change_df.loc[dmso_mask].groupby('Day')['Delta_Ct'].first()
        fold_change_df['Baseline_Delta_Ct'] = fold_change_df['Day'].map(dmso_baseline)
        if fold_change_df['Baseline_Delta_Ct'].isna().all():
            reference_sample = fold_change_df['Sample_Number'].iloc[0]
        else:
            fold_change_df['Delta_Delta_Ct'] = fold_change_df['Delta_Ct'] - fold_change_df['Baseline_Delta_Ct']
    
    if reference_sample is not None:
        if isinstance(reference_sample, str):
            if reference_sample.lower() == 'dmso':
                dmso_mask = fold_change_df['Sample_Name'].str.contains('DMSO', case=False, na=False)
                dmso_baseline = fold_change_df.loc[dmso_mask].groupby('Day')['Delta_Ct'].first()
                fold_change_df['Baseline_Delta_Ct'] = fold_change_df['Day'].map(dmso_baseline)
                fold_change_df['Delta_Delta_Ct'] = fold_change_df['Delta_Ct'] - fold_change_df['Baseline_Delta_Ct']
            else:
                ref_delta_ct = fold_change_df.loc[fold_change_df['Sample_Name'] == reference_sample, 'Delta_Ct'].iloc[0]
                fold_change_df['Delta_Delta_Ct'] = fold_change_df['Delta_Ct'] - ref_delta_ct
        else:
            ref_delta_ct = fold_change_df.loc[fold_change_df['Sample_Number'] == reference_sample, 'Delta_Ct'].iloc[0]
            fold_change_df['Delta_Delta_Ct'] = fold_change_df['Delta_Ct'] - ref_delta_ct
    
    if 'Delta_Delta_Ct' not in fold_change_df:
        fold_change_df['Delta_Delta_Ct'] = np.nan
    
    fold_change_df['Fold_Change'] = 2 ** (-fold_change_df['Delta_Delta_Ct'])
    fold_change_df.loc[fold_change_df['Delta_Delta_Ct'].isna(), 'Fold_Change'] = np.nan
    fold_change_df['Fold_Change_STD'] = fold_change_df['Fold_Change'] * np.log(2) * fold_change_df['Delta_Ct_STD']
    fold_change_df.loc[fold_change_df['Fold_Change'].isna(), 'Fold_Change_STD'] = np.nan
    
    return fold_change_df


def save_results(delta_ct_df, fold_change_df=None, output_prefix="qPCR_results"):
    """
    Save calculation results to CSV files.
    
    Parameters:
    -----------
    delta_ct_df : pandas.DataFrame
        ΔCt results
    fold_change_df : pandas.DataFrame, optional
        Fold change results
    output_prefix : str
        Prefix for output filenames
    """
    delta_ct_df.to_csv(f"{output_prefix}_delta_ct.csv", index=False)
    print(f"✓ Saved: {output_prefix}_delta_ct.csv")
    
    if fold_change_df is not None:
        fold_change_df.to_csv(f"{output_prefix}_fold_change.csv", index=False)
        print(f"✓ Saved: {output_prefix}_fold_change.csv")


def plot_delta_ct(delta_ct_df, colormap='viridis', output_filename="delta_ct_plot.png"):
    """
    Create bar plot for ΔCt values with error bars.
    Uses neurodivergent-friendly colormaps.
    
    Parameters:
    -----------
    delta_ct_df : pandas.DataFrame
        DataFrame with Delta Ct values
    colormap : str
        Colormap to use ('viridis' or 'magma')
    output_filename : str
        Filename to save the plot
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    colors = plt.colormaps.get_cmap(colormap)(np.linspace(0, 1, len(delta_ct_df)))
    x_pos = np.arange(len(delta_ct_df))
    ax.bar(x_pos, delta_ct_df['Delta_Ct'], 
           yerr=delta_ct_df['Delta_Ct_STD'],
           capsize=5, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
    
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('ΔCt (Target - Reference)', fontsize=12, fontweight='bold')
    ax.set_title('qPCR Delta Ct Values', fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(delta_ct_df['Sample_Name'], rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_filename}")
    plt.show()


def plot_fold_change(fold_change_df, colormap='viridis', output_filename="fold_change_plot.png"):
    """
    Create bar plot for Fold Change values with error bars.
    
    Parameters:
    -----------
    fold_change_df : pandas.DataFrame
        DataFrame with Fold Change values
    colormap : str
        Colormap to use ('viridis' or 'magma')
    output_filename : str
        Filename to save the plot
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = plt.colormaps.get_cmap(colormap)(np.linspace(0, 1, len(fold_change_df)))
    x_pos = np.arange(len(fold_change_df))
    ax.bar(x_pos, fold_change_df['Fold_Change'],
           yerr=fold_change_df['Fold_Change_STD'],
           capsize=5, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
    
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Fold Change (2^-ΔΔCt)', fontsize=12, fontweight='bold')
    ax.set_title('qPCR Fold Change Values', fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(fold_change_df['Sample_Name'], rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=1, color='red', linestyle='--', linewidth=2, alpha=0.5, label='No change')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_filename}")
    plt.show()


def print_results(delta_ct_df, fold_change_df=None):
    """
    Print summary statistics.
    
    Parameters:
    -----------
    delta_ct_df : pandas.DataFrame
        ΔCt results
    fold_change_df : pandas.DataFrame, optional
        Fold change results
    """
    print("\n" + "="*70)
    print("DELTA CT CALCULATION RESULTS")
    print("="*70)
    print("\nΔCt = Target Ct - Reference Ct")
    print(delta_ct_df[['Sample_Number', 'Day', 'Sample_Name', 'Target_Ct', 'Reference_Ct', 'Delta_Ct', 'Delta_Ct_STD']].to_string(index=False))
    
    if fold_change_df is not None:
        print("\n" + "="*70)
        print("FOLD CHANGE RESULTS (2^-ΔΔCt)")
        print("="*70)
        print("\nNormalized to DMSO control per day where available")
        print(fold_change_df[['Sample_Number', 'Day', 'Sample_Name', 'Delta_Ct', 'Delta_Delta_Ct', 'Fold_Change', 'Fold_Change_STD']].to_string(index=False))
        
        print(f"\nMean Fold Change: {fold_change_df['Fold_Change'].mean():.4f} ± {fold_change_df['Fold_Change_STD'].mean():.4f}")


def main():
    """
    Main function - run complete qPCR analysis pipeline
    """
    # Configuration
    INPUT_FILE = "230913 qPCR and SD.txt"
    RAW_SAMPLE_FILE = "230913 qPCR and SD2.txt"
    COLORMAP = 'viridis'  # Change to 'magma' for alternative colormap
    OUTPUT_PREFIX = "qPCR_results"
    
    print("\n" + "="*70)
    print("qPCR ANALYSIS PIPELINE")
    print("="*70)
    
    # Step 1: Load data
    print("\n[1/4] Loading qPCR data...")
    df = load_qpcr_data(INPUT_FILE, raw_sample_file=RAW_SAMPLE_FILE)
    target_df, reference_df = separate_genes(df)
    print(f"✓ Loaded {len(df)} samples")
    
    # Step 2: Calculate ΔCt
    print("\n[2/4] Calculating ΔCt...")
    delta_ct_df = calculate_delta_ct(target_df, reference_df)
    print(f"✓ ΔCt calculated for {len(delta_ct_df)} samples")
    
    # Step 3: Calculate Fold Change
    print("\n[3/4] Calculating Fold Change...")
    fold_change_df = calculate_fold_change(delta_ct_df)
    print(f"✓ Fold change calculated (normalized to sample {fold_change_df['Sample_Number'].iloc[0]})")
    
    # Step 4: Print and save results
    print("\n[4/4] Saving results...")
    print_results(delta_ct_df, fold_change_df)
    save_results(delta_ct_df, fold_change_df, OUTPUT_PREFIX)
    
    # Create visualizations
    print("\n[5/4] Creating visualizations...")
    plot_delta_ct(delta_ct_df, colormap=COLORMAP, 
                  output_filename=f"{OUTPUT_PREFIX}_delta_ct_{COLORMAP}.png")
    plot_fold_change(fold_change_df, colormap=COLORMAP,
                     output_filename=f"{OUTPUT_PREFIX}_fold_change_{COLORMAP}.png")
    
    print("\n" + "="*70)
    print("✓ ANALYSIS COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
