import os
import datetime
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from fpdf import FPDF
from qPCR_pdf_layout import QPCRPDFReport, create_qpcr_summary_pdf

def get_output_directory(output_dir=None):
    """Gibt das Ausgabeverzeichnis für qPCR-Ergebnisse und Plots zurück."""
    if output_dir is None:
        output_dir = Path.home() / "Desktop" / "qPCR_Ergebnisse"
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def save_figure_to_output(fig, output_filename, output_dir=None):
    """Speichert eine Matplotlib-Figur im qPCR-Ausgabeverzeichnis."""
    target_path = Path(output_filename)
    if not target_path.is_absolute():
        target_path = get_output_directory(output_dir) / target_path

    target_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {target_path}")
    plt.close(fig)
    return target_path


def load_results_csv(csv_file):
    """Lädt zuvor berechnete Ergebnisse aus einer CSV-Datei."""
    return pd.read_csv(csv_file)

def create_melting_curve_plot(raw_file, output_path):
    """
    Liest die Rohdaten ein, filtert Programm 3 (Schmelzkurve) und 
    erstellt einen sauberen Plot der Fluoreszenz über die Temperatur.
    """
    try:
        # 1. Rohdaten einlesen (automatische Trennzeichen-Erkennung oder Tabulator)
        if str(raw_file).endswith('.csv'):
            df = pd.read_csv(raw_file, sep=None, engine='python')
        else:
            df = pd.read_csv(raw_file, sep='\t')

        # 2. Ausschließlich die Schmelzkurve (Programm 3) extrahieren
        if 'Prog#' in df.columns and 'SampleName' in df.columns:
            df_melt = df[df['Prog#'] == 3].copy()
        else:
            # Fallback falls Spalten anders heißen oder das Format abweicht
            df_melt = df.copy()

        if df_melt.empty:
            print("Schmelzkurve übersprungen: Keine Daten für Prog# == 3 gefunden.")
            return False

        # 3. Schmelzkurven visualisieren
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Iteriere durch jede vorhandene Probe und plotte sie
        sample_col = 'SampleName' if 'SampleName' in df_melt.columns else df_melt.columns[0]
        temp_col = 'Temp' if 'Temp' in df_melt.columns else next((c for c in df_melt.columns if 'temp' in c.lower()), None)
        signal_col = '465-510' if '465-510' in df_melt.columns else next((c for c in df_melt.columns if '465' in c or 'fluores' in c.lower()), None)

        if not temp_col or not signal_col:
            print(f"Schmelzkurve: Benötigte Spalten (Temp / Signal) nicht gefunden in {list(df_melt.columns)}")
            return False

        for sample in df_melt[sample_col].unique():
            sample_data = df_melt[df_melt[sample_col] == sample]
            ax.plot(sample_data[temp_col], sample_data[signal_col], label=str(sample))

        ax.set_title('Schmelzkurve (Fluoreszenz vs. Temperatur)', fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel('Temperatur (°C)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Fluoreszenz (465-510 nm)', fontsize=12, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
        ax.grid(True, linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        
        # Speichern für das PDF-Layout
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        # Optional: Direkt als JSON exportieren (z. B. für externe Pipelines / Blender-Animationen)
        json_output_path = Path(output_path).with_name('schmelzkurve_gefiltert.json')
        df_melt.to_json(json_output_path, orient='records', indent=4)
        
        return True

    except Exception as e:
        print(f"Fehler beim Erstellen der Schmelzkurve: {e}")
        return False

def create_grouped_boxplot(df, output_path, value_col='Delta_Ct', group_col='Day', title='Grouped qPCR Boxplot'):
    """Create a grouped boxplot for qPCR data (e.g. by Day)."""
    if df is None or df.empty:
        return False

    if group_col not in df.columns:
        if 'Day' in df.columns:
            group_col = 'Day'
        else:
            return False

    if value_col not in df.columns:
        return False

    plot_df = df[[group_col, value_col]].dropna().copy()
    if plot_df.empty:
        return False

    fig, ax = plt.subplots(figsize=(10, 6))
    grouped = [values.tolist() for _, values in plot_df.groupby(group_col)[value_col]]
    labels = [str(label) for label in plot_df[group_col].drop_duplicates().tolist()]

    if not grouped:
        return False

    ax.boxplot(grouped, labels=labels, patch_artist=True)
    ax.set_xlabel(group_col, fontsize=12, fontweight='bold')
    ax.set_ylabel(value_col, fontsize=12, fontweight='bold')
    ax.set_title(f'{title} by {group_col}', fontsize=14, fontweight='bold')
    ax.grid(axis='y', linestyle='--', alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return True


def create_volcano_plot(fold_change_df, output_path, title='Volcano Plot'):
    """Create a volcano plot using fold change and p-value-like significance proxies."""
    if fold_change_df is None or fold_change_df.empty:
        return False

    plot_df = fold_change_df.dropna(subset=['Fold_Change']).copy()
    if plot_df.empty:
        return False

    plot_df['log2_fc'] = np.log2(plot_df['Fold_Change'].astype(float))

    if 'P_Value' in plot_df.columns:
        p_values = pd.to_numeric(plot_df['P_Value'], errors='coerce')
    elif 'p_value' in plot_df.columns:
        p_values = pd.to_numeric(plot_df['p_value'], errors='coerce')
    else:
        uncertainty = pd.to_numeric(plot_df.get('Fold_Change_STD', pd.Series(1.0, index=plot_df.index)), errors='coerce').fillna(1.0)
        z_scores = np.abs(plot_df['log2_fc']) / uncertainty.replace(0, 1.0)
        p_values = np.exp(-z_scores)

    plot_df['neg_log10_p'] = -np.log10(p_values.clip(lower=1e-300))
    plot_df = plot_df.dropna(subset=['log2_fc', 'neg_log10_p'])

    if plot_df.empty:
        return False

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(plot_df['log2_fc'], plot_df['neg_log10_p'], s=75, color='#4C78A8', edgecolor='black', alpha=0.8)

    for _, row in plot_df.iterrows():
        label = row.get('Sample_Name') or row.get('Target_Sample_Name') or row.get('Reference_Sample_Name') or ''
        if label:
            ax.annotate(str(label), (row['log2_fc'], row['neg_log10_p']), xytext=(3, 3), textcoords='offset points', fontsize=7)

    ax.axvline(0, color='gray', linestyle='--', linewidth=1)
    ax.axvline(1, color='red', linestyle='--', linewidth=1, alpha=0.7)
    ax.axvline(-1, color='red', linestyle='--', linewidth=1, alpha=0.7)
    ax.axhline(-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.7)
    ax.set_xlabel('log2(Fold Change)', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return True


def create_violin_plot(df, output_path, value_col='Delta_Ct', group_col='Day', title='Delta Ct Distribution by Day'):
    """Create a violin plot for qPCR values grouped by a categorical variable."""
    if df is None or df.empty:
        return False

    if group_col not in df.columns or value_col not in df.columns:
        return False

    plot_df = df[[group_col, value_col]].dropna().copy()
    if plot_df.empty:
        return False

    labels = [str(label) for label in plot_df[group_col].drop_duplicates().tolist()]
    grouped_values = [values.tolist() for _, values in plot_df.groupby(group_col)[value_col]]

    if not grouped_values:
        return False

    fig, ax = plt.subplots(figsize=(8, 6))
    parts = ax.violinplot(grouped_values, showmeans=False, showmedians=True)
    for body in parts['bodies']:
        body.set_facecolor('#4C78A8')
        body.set_edgecolor('#2F4B7C')
        body.set_alpha(0.8)

    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels if labels else ['Value'])
    ax.set_ylabel(value_col, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(axis='y', linestyle='--', alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return True


def plot_delta_ct_with_colormap(delta_ct_df, colormap='viridis', output_filename=None, sample_name=None, output_dir=None):
    """Erstellt ein Balkendiagramm für Delta Ct-Werte."""
    if output_filename is None:
        output_filename = f"qPCR_results_delta_ct_{colormap}.png"
    
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = plt.colormaps.get_cmap(colormap)(np.linspace(0, 1, len(delta_ct_df)))
    x_pos = np.arange(len(delta_ct_df))
    
    yerr = delta_ct_df['Delta_Ct_STD'] if 'Delta_Ct_STD' in delta_ct_df.columns else None
    
    ax.bar(x_pos, delta_ct_df['Delta_Ct'], 
           yerr=yerr,
           capsize=5, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
    
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
    return save_figure_to_output(fig, output_filename, output_dir=output_dir)


def plot_fold_change_with_colormap(fold_change_df, colormap='viridis', output_filename=None, sample_name=None, output_dir=None):
    """Erstellt ein Balkendiagramm für Fold Change-Werte."""
    if output_filename is None:
        output_filename = f"qPCR_results_fold_change_{colormap}.png"
    
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = plt.colormaps.get_cmap(colormap)(np.linspace(0, 1, len(fold_change_df)))
    x_pos = np.arange(len(fold_change_df))
    
    yerr = fold_change_df['Fold_Change_STD'] if 'Fold_Change_STD' in fold_change_df.columns else None
    
    ax.bar(x_pos, fold_change_df['Fold_Change'],
           yerr=yerr,
           capsize=5, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
    
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
    ax.axhline(y=1, color='red', linestyle='--', linewidth=2, alpha=0.5, label='Kontrolle (1.0)')
    ax.legend()
    
    plt.tight_layout()
    return save_figure_to_output(fig, output_filename, output_dir=output_dir)


def save_combined_publication_quality_plots(delta_ct_df, fold_change_df, colormap='viridis', output_filename=None, sample_name=None, output_dir=None):
    """Erstellt ein 2-Panel Kombi-Diagramm in Publikationsqualität (300 DPI)."""
    if output_filename is None:
        output_filename = f"qPCR_publication_combined_{colormap}.png"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Panel A: Delta Ct
    colors1 = plt.colormaps.get_cmap(colormap)(np.linspace(0, 1, len(delta_ct_df)))
    x_pos1 = np.arange(len(delta_ct_df))
    yerr1 = delta_ct_df['Delta_Ct_STD'] if 'Delta_Ct_STD' in delta_ct_df.columns else None

    ax1.bar(x_pos1, delta_ct_df['Delta_Ct'], yerr=yerr1, capsize=4,
            color=colors1, edgecolor='black', linewidth=1.2, alpha=0.85)
    ax1.set_xlabel('Sample', fontsize=11, fontweight='bold')
    ax1.set_ylabel('ΔCt (Target - Reference)', fontsize=11, fontweight='bold')
    ax1.set_title('A: ΔCt Values', fontsize=12, fontweight='bold', loc='left')
    ax1.set_xticks(x_pos1)
    labels1 = delta_ct_df['Sample_Name'] if 'Sample_Name' in delta_ct_df.columns else [f"S{int(n)}" for n in delta_ct_df['Sample_Number']]
    ax1.set_xticklabels(labels1, rotation=45, ha='right', fontsize=9)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.7)

    # Panel B: Fold Change
    colors2 = plt.colormaps.get_cmap(colormap)(np.linspace(0, 1, len(fold_change_df)))
    x_pos2 = np.arange(len(fold_change_df))
    yerr2 = fold_change_df['Fold_Change_STD'] if 'Fold_Change_STD' in fold_change_df.columns else None

    ax2.bar(x_pos2, fold_change_df['Fold_Change'], yerr=yerr2, capsize=4,
            color=colors2, edgecolor='black', linewidth=1.2, alpha=0.85)
    ax2.set_xlabel('Sample', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Fold Change (2^-ΔΔCt)', fontsize=11, fontweight='bold')
    ax2.set_title('B: Relative Expression (Fold Change)', fontsize=12, fontweight='bold', loc='left')
    ax2.set_xticks(x_pos2)
    labels2 = fold_change_df['Sample_Name'] if 'Sample_Name' in fold_change_df.columns else [f"S{int(n)}" for n in fold_change_df['Sample_Number']]
    ax2.set_xticklabels(labels2, rotation=45, ha='right', fontsize=9)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    ax2.axhline(y=1, color='red', linestyle='--', linewidth=1.2, alpha=0.7, label='Kontrolle (1.0)')
    ax2.legend(loc='upper right', frameon=True)

    if sample_name:
        fig.suptitle(f'Publication Figure: {sample_name}', fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()
    return save_figure_to_output(fig, output_filename, output_dir=output_dir)

### Layout was here: