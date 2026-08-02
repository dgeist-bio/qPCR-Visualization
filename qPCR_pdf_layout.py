import os
import datetime
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from fpdf import FPDF

### layout of the PDF summary


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


def create_qpcr_summary_pdf(delta_ct_df, fold_change_df, stats_report, output_path, sample_name=None, plot_paths=None, melting_plot=None, delta_ct_plot=None, fold_change_plot=None, grouped_boxplot=None):
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
        
        headers_fc = ["Sample #", "Sample Name", "Delta Ct", "Delta Delta Ct", "Fold Change", "FC StdDev"]
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

# Page 3: Plots & Visualisierungen
    if delta_ct_plot or fold_change_plot or plot_paths or melting_plot or grouped_boxplot:
        pdf.add_page()
        pdf.section_title("4. Plots & Visualisierungen")
        
        # Falls plot_paths noch genutzt wird (z.B. für ältere Aufrufe)
        if plot_paths:
            for p_path in plot_paths:
                if p_path and os.path.exists(p_path):
                    pdf.image(str(p_path), w=170)
                    pdf.ln(5)

        # 1. Delta Ct Plot
        if delta_ct_plot and os.path.exists(delta_ct_plot):
            pdf.set_font("Helvetica", "B", 11)
            pdf.set_text_color(*pdf.c_primary)
            pdf.cell(0, 8, "Delta Ct Distribution", ln=True)
            pdf.image(str(delta_ct_plot), w=170)
            pdf.ln(8)

        # 2. Fold Change Plot
        if fold_change_plot and os.path.exists(fold_change_plot):
            pdf.set_font("Helvetica", "B", 11)
            pdf.set_text_color(*pdf.c_primary)
            pdf.cell(0, 8, "Fold Change Analysis", ln=True)
            pdf.image(str(fold_change_plot), w=170)
            pdf.ln(8)

        # 3. Grouped Boxplot
        if grouped_boxplot and os.path.exists(grouped_boxplot):
            pdf.set_font("Helvetica", "B", 11)
            pdf.set_text_color(*pdf.c_primary)
            pdf.cell(0, 8, "Grouped Delta Ct Boxplot (by Day)", ln=True)
            pdf.image(str(grouped_boxplot), w=170)
            pdf.ln(8)

    # Page 4: Schmelzkurve (auf neuer Seite für mehr Platz)
    if melting_plot and os.path.exists(melting_plot):
        pdf.add_page()
        pdf.section_title("5. Melting Curve Analysis")
        
        pdf.set_font("Helvetica", "", 9)
        pdf.set_text_color(100, 100, 100)
        pdf.cell(0, 5, "Visualisierung der Dissoziationskurve zur Ueberpruefung der Primer-Spezifitaet.", ln=True)
        pdf.ln(5)
        
        pdf.image(str(melting_plot), w=170)
        pdf.ln(5)

    # Output PDF
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # In FPDF wird das Dokument mit output() generiert (nicht mit doc.build)
    pdf.output(str(output_path))
    print(f"✓ PDF erfolgreich im FPDF-Layout gespeichert: {output_path}")
    return output_path