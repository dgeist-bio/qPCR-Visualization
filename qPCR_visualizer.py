import os
import datetime
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from fpdf import FPDF


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


class QPCRPDFReport(FPDF):
    """FPDF-Klasse für ein professionelles qPCR-Berichtslayout."""

    def __init__(self, sample_name=None):
        super().__init__(orientation='P', unit='mm', format='A4')
        self.set_auto_page_break(auto=True, margin=15)
        self.sample_name = sample_name or "Nicht angegeben"
        
        self.c_primary = (31, 106, 165)     # Steel Blue
        self.c_secondary = (47, 165, 114)   # Emerald Green
        self.c_text = (40, 40, 40)          # Charcoal Text
        self.c_bg_light = (248, 249, 250)   # Light Gray
        self.c_border = (220, 224, 230)     # Border

    def sanitize_text(self, text: str) -> str:
        """Ersetzt Unicode-Sonderzeichen für Standard-Helvetica."""
        if not isinstance(text, str):
            text = str(text)
        return (text.replace('Δ', 'Delta ')
                    .replace('±', '+/- ')
                    .replace('α', 'alpha')
                    .replace('–', '-')
                    .replace('—', '-'))

    def header(self):
        self.set_fill_color(*self.c_primary)
        self.rect(0, 0, 210, 4, 'F')

        self.set_font("Helvetica", "B", 8)
        self.set_text_color(130, 130, 130)
        self.set_xy(12, 8)
        self.cell(100, 5, "qPCR ANALYSIS REPORT", align="L")
        self.set_xy(100, 8)
        self.cell(98, 5, self.sanitize_text(f"Probe: {self.sample_name}"), align="R")
        
        self.set_draw_color(*self.c_border)
        self.set_line_width(0.3)
        self.line(12, 14, 198, 14)
        self.ln(8)

    def footer(self):
        self.set_y(-12)
        self.set_font("Helvetica", "I", 8)
        self.set_text_color(150, 150, 150)
        page_str = self.sanitize_text(f"Seite {self.page_no()} / {{nb}}")
        self.cell(0, 8, page_str, align="C")

    def section_title(self, title: str):
        self.ln(3)
        self.set_font("Helvetica", "B", 13)
        self.set_text_color(*self.c_primary)
        self.cell(0, 7, self.sanitize_text(title), ln=True, align="L")
        self.set_draw_color(*self.c_primary)
        self.set_line_width(0.5)
        self.line(12, self.get_y(), 198, self.get_y())
        self.set_line_width(0.2)
        self.ln(4)

    def draw_summary_cards(self, sample_count, mean_delta_ct, mean_fc, stats_summary):
        y_pos = self.get_y()
        card_w = 43.5
        card_h = 18
        spacing = 3.5
        x_start = 12

        cards = [
            ("Proben Anzahl", str(sample_count), (245, 247, 250)),
            ("Mean Delta Ct", f"{mean_delta_ct:.3f}" if not np.isnan(mean_delta_ct) else "N/A", (235, 243, 250)),
            ("Mean Fold Change", f"{mean_fc:.3f}" if not np.isnan(mean_fc) else "N/A", (235, 248, 242)),
            ("Statistik Status", stats_summary, (253, 242, 242) if "YES" in stats_summary or "Signifikant" in stats_summary else (245, 247, 250))
        ]

        for i, (label, val, bg_color) in enumerate(cards):
            x = x_start + i * (card_w + spacing)
            self.set_xy(x, y_pos)
            self.set_fill_color(*bg_color)
            self.set_draw_color(*self.c_border)
            self.rect(x, y_pos, card_w, card_h, 'DF')

            self.set_xy(x, y_pos + 2.5)
            self.set_font("Helvetica", "B", 7.5)
            self.set_text_color(110, 110, 110)
            self.cell(card_w, 4, self.sanitize_text(label), align="C")

            self.set_xy(x, y_pos + 8)
            self.set_font("Helvetica", "B", 11)
            self.set_text_color(*self.c_text)
            self.cell(card_w, 6, self.sanitize_text(val), align="C")

        self.set_xy(12, y_pos + card_h + 6)


def create_qpcr_summary_pdf(delta_ct_df, fold_change_df, stats_report, output_path, sample_name=None, plot_paths=None):
    """Erstellt das qPCR-Zusammenfassungs-PDF."""
    pdf = QPCRPDFReport(sample_name=sample_name)
    pdf.alias_nb_pages()
    pdf.add_page()

    # Title
    pdf.set_font("Helvetica", "B", 18)
    pdf.set_text_color(*pdf.c_primary)
    pdf.cell(0, 8, "qPCR Executive Summary", ln=True, align="L")
    
    pdf.set_font("Helvetica", "", 9)
    pdf.set_text_color(120, 120, 120)
    now_str = datetime.datetime.now().strftime("%d.%m.%Y um %H:%M Uhr")
    pdf.cell(0, 5, pdf.sanitize_text(f"Erstellt am: {now_str}"), ln=True, align="L")
    pdf.ln(3)

    # Metrics
    sample_count = len(delta_ct_df) if delta_ct_df is not None else 0
    mean_d_ct = float(delta_ct_df['Delta_Ct'].mean()) if (delta_ct_df is not None and 'Delta_Ct' in delta_ct_df) else np.nan
    mean_fc = float(fold_change_df['Fold_Change'].mean()) if (fold_change_df is not None and 'Fold_Change' in fold_change_df) else np.nan
    
    stats_flag = "Signifikant" if (stats_report and "p-value" in stats_report.lower() and "significant (alpha=0.05): yes" in stats_report.lower()) else "Unauffällig"

    pdf.draw_summary_cards(sample_count, mean_d_ct, mean_fc, stats_flag)

    # Table 1: Delta Ct
    if delta_ct_df is not None and not delta_ct_df.empty:
        pdf.section_title("1. Delta Ct Ergebnisse")
        
        headers = ["Sample #", "Sample Name", "Target Ct", "Ref Ct", "Delta Ct", "Std. Dev"]
        widths = [20, 50, 28, 28, 30, 30]

        pdf.set_font("Helvetica", "B", 8.5)
        pdf.set_fill_color(*pdf.c_primary)
        pdf.set_text_color(255, 255, 255)
        for h, w in zip(headers, widths):
            pdf.cell(w, 6, pdf.sanitize_text(h), border=0, align="C", fill=True)
        pdf.ln()

        pdf.set_font("Helvetica", "", 8)
        pdf.set_text_color(*pdf.c_text)
        fill = False

        for _, row in delta_ct_df.iterrows():
            pdf.set_fill_color(248, 249, 250) if fill else pdf.set_fill_color(255, 255, 255)
            
            vals = [
                str(row.get('Sample_Number', '')),
                str(row.get('Sample_Name', ''))[:26],
                f"{row.get('Target_Ct', 0):.3f}" if pd.notnull(row.get('Target_Ct')) else "N/A",
                f"{row.get('Reference_Ct', 0):.3f}" if pd.notnull(row.get('Reference_Ct')) else "N/A",
                f"{row.get('Delta_Ct', 0):.3f}" if pd.notnull(row.get('Delta_Ct')) else "N/A",
                f"{row.get('Delta_Ct_STD', 0):.3f}" if pd.notnull(row.get('Delta_Ct_STD')) else "N/A"
            ]
            
            for v, w in zip(vals, widths):
                pdf.cell(w, 5.5, pdf.sanitize_text(v), border="B", align="C", fill=fill)
            pdf.ln()
            fill = not fill

    # Table 2: Fold Change
    if fold_change_df is not None and not fold_change_df.empty:
        pdf.ln(4)
        pdf.section_title("2. Fold Change Ergebnisse (2^-Delta Delta Ct)")
        
        headers_fc = ["Sample #", "Sample Name", "Delta Ct", "dDelta Ct", "Fold Change", "FC StdDev"]
        widths_fc = [20, 50, 28, 28, 30, 30]

        pdf.set_font("Helvetica", "B", 8.5)
        pdf.set_fill_color(47, 165, 114)
        pdf.set_text_color(255, 255, 255)
        for h, w in zip(headers_fc, widths_fc):
            pdf.cell(w, 6, pdf.sanitize_text(h), border=0, align="C", fill=True)
        pdf.ln()

        pdf.set_font("Helvetica", "", 8)
        pdf.set_text_color(*pdf.c_text)
        fill = False

        for _, row in fold_change_df.iterrows():
            pdf.set_fill_color(248, 249, 250) if fill else pdf.set_fill_color(255, 255, 255)
            
            vals = [
                str(row.get('Sample_Number', '')),
                str(row.get('Sample_Name', ''))[:26],
                f"{row.get('Delta_Ct', 0):.3f}" if pd.notnull(row.get('Delta_Ct')) else "N/A",
                f"{row.get('Delta_Delta_Ct', 0):.3f}" if pd.notnull(row.get('Delta_Delta_Ct')) else "N/A",
                f"{row.get('Fold_Change', 0):.3f}" if pd.notnull(row.get('Fold_Change')) else "N/A",
                f"{row.get('Fold_Change_STD', 0):.3f}" if pd.notnull(row.get('Fold_Change_STD')) else "N/A"
            ]
            
            for v, w in zip(vals, widths_fc):
                pdf.cell(w, 5.5, pdf.sanitize_text(v), border="B", align="C", fill=fill)
            pdf.ln()
            fill = not fill

    # Page 2: Statistics Report
    if stats_report:
        pdf.add_page()
        pdf.section_title("3. Statistische Analyse & Auswertung")
        
        pdf.set_font("Courier", "", 7.5)
        pdf.set_fill_color(248, 249, 250)
        pdf.set_draw_color(*pdf.c_border)
        
        cleaned_stats = pdf.sanitize_text(stats_report)
        pdf.multi_cell(186, 3.8, cleaned_stats, border=1, align="L", fill=True)

    # Output PDF
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pdf.output(str(output_path))
    print(f"✓ PDF erfolgreich im FPDF-Layout gespeichert: {output_path}")
    return output_path