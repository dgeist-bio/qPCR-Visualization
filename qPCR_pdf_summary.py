import os
import json
import re
from datetime import datetime
from tkinter import messagebox
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu

# Import qPCR modules
from qPCR_data_loader import load_qpcr_data
from qPCR_delta_ct_calculation import calculate_delta_ct, calculate_fold_change
from qPCR_visualizer import (
    plot_delta_ct_with_colormap,
    plot_fold_change_with_colormap,
    create_qpcr_summary_pdf,
    get_output_directory,
    create_melting_curve_plot,
    create_grouped_boxplot,
    create_volcano_plot,
    create_violin_plot,
)
from qPCR_pdf_layout import create_qpcr_summary_pdf

class qPCRSummaryMixin:
    """Mix-in class containing calculation and export logic for qPCRAnalyzerApp."""
    
    def calculate_data(self):
        """Load files and calculate ΔCt and Fold Change."""
        if not self.target_file or not self.reference_file:
            raise ValueError("Please select both target and reference files")

        self.update_status("Loading data...")
        self.target_data = load_qpcr_data(self.target_file, raw_sample_file=self.raw_file)
        self.reference_data = load_qpcr_data(self.reference_file, raw_sample_file=self.raw_file)


        self.update_status("Calculating ΔCt...")
        self.delta_ct_data = calculate_delta_ct(self.target_data, self.reference_data)

        self.update_status("Calculating Fold Change...")
        self.fold_change_data = calculate_fold_change(self.delta_ct_data)

        self.update_status("Preparing statistical report...")
        self.stats_report = self.build_statistical_report()

        return self.stats_report

    def build_statistical_report(self):
        """Build a statistical analysis report for the current qPCR data."""
        if self.delta_ct_data is None or self.fold_change_data is None:
            return ""

        stats_report = f"""
{'='*80}
STATISTICAL ANALYSIS REPORT
{'='*80}
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

{'='*80}
1. DESCRIPTIVE STATISTICS - Delta Ct
{'='*80}
Mean:                 {self.delta_ct_data['Delta_Ct'].mean():.4f}
Std Deviation:        {self.delta_ct_data['Delta_Ct'].std():.4f}
Median:               {self.delta_ct_data['Delta_Ct'].median():.4f}
Min:                  {self.delta_ct_data['Delta_Ct'].min():.4f}
Max:                  {self.delta_ct_data['Delta_Ct'].max():.4f}
Range:                {self.delta_ct_data['Delta_Ct'].max() - self.delta_ct_data['Delta_Ct'].min():.4f}
Coefficient of Variation: {(self.delta_ct_data['Delta_Ct'].std() / self.delta_ct_data['Delta_Ct'].mean() * 100):.2f}%

{'='*80}
2. DESCRIPTIVE STATISTICS - Fold Change
{'='*80}
Mean:                 {self.fold_change_data['Fold_Change'].mean():.4f}
Std Deviation:        {self.fold_change_data['Fold_Change'].std():.4f}
Median:               {self.fold_change_data['Fold_Change'].median():.4f}
Min:                  {self.fold_change_data['Fold_Change'].min():.4f}
Max:                  {self.fold_change_data['Fold_Change'].max():.4f}
Coefficient of Variation: {(self.fold_change_data['Fold_Change'].std() / self.fold_change_data['Fold_Change'].mean() * 100):.2f}%

{'='*80}
3. TARGET vs REFERENCE CT COMPARISON
{'='*80}
"""
        target_vals = self.delta_ct_data['Target_Ct'].dropna()
        reference_vals = self.delta_ct_data['Reference_Ct'].dropna()

        if len(target_vals) > 1 and len(reference_vals) > 1:
            t_stat, p_val = ttest_ind(target_vals, reference_vals)
            stats_report += f"""
Target Ct Mean:       {target_vals.mean():.4f} ± {target_vals.std():.4f}
Reference Ct Mean:    {reference_vals.mean():.4f} ± {reference_vals.std():.4f}

Independent t-test:
  t-statistic:        {t_stat:.4f}
  p-value:            {p_val:.6f}
  Significant (α=0.05): {'YES' if p_val < 0.05 else 'NO'}

Mann-Whitney U Test (non-parametric alternative):
"""
            u_stat, p_val_mw = mannwhitneyu(target_vals, reference_vals)
            stats_report += f"""  U-statistic:        {u_stat:.4f}
  p-value:            {p_val_mw:.6f}
  Significant (α=0.05): {'YES' if p_val_mw < 0.05 else 'NO'}
"""

        return stats_report

    def build_analysis_json(self):
        """Build the raw plate metadata JSON from the original qPCR raw export."""
        if self.raw_file and os.path.exists(self.raw_file):
            from qPCR_data_loader import parse_raw_plate_data
            samples = parse_raw_plate_data(self.raw_file)
            return {
                "generated_at": datetime.now().isoformat(),
                "samples": [
                    {
                        "position": row.get('position'),
                        "well": row.get('well'),
                        "name": row.get('name'),
                        "sample_name": row.get('sample_name'),
                        "cp": row.get('cp'),
                        "cp_values": {
                            "cp": row.get('cp'),
                        },
                        "concentration": row.get('concentration'),
                        "day": row.get('day'),
                        "day_number": row.get('day_number'),
                    }
                    for row in samples
                ],
            }

        records = []
        if self.target_data is not None and not self.target_data.empty:
            for _, row in self.target_data.iterrows():
                cp = row.get('MeanCp')
                records.append({
                    "position": row.get('Well1'),
                    "well": row.get('Well1'),
                    "name": row.get('Sample_Name') or row.get('Samples'),
                    "sample_name": row.get('Sample_Name') or row.get('Samples'),
                    "cp": cp,
                    "cp_values": {"cp": cp},
                    "concentration": row.get('Mean_conc'),
                    "day": row.get('Day'),
                    "day_number": int(re.search(r'(\d+)', str(row.get('Day') or '')).group(1)) if isinstance(row.get('Day'), str) and re.search(r'(\d+)', row.get('Day') or '') else None,
                })
        if self.reference_data is not None and not self.reference_data.empty:
            for _, row in self.reference_data.iterrows():
                cp = row.get('MeanCp')
                records.append({
                    "position": row.get('Well1'),
                    "well": row.get('Well1'),
                    "name": row.get('Sample_Name') or row.get('Samples'),
                    "sample_name": row.get('Sample_Name') or row.get('Samples'),
                    "cp": cp,
                    "cp_values": {"cp": cp},
                    "concentration": row.get('Mean_conc'),
                    "day": row.get('Day'),
                    "day_number": int(re.search(r'(\d+)', str(row.get('Day') or '')).group(1)) if isinstance(row.get('Day'), str) and re.search(r'(\d+)', row.get('Day') or '') else None,
                })

        return {
            "generated_at": datetime.now().isoformat(),
            "samples": records,
        }

    def build_summary_json(self):
        """Build a compact summary JSON with well, sample name, day, and fold change."""
        if self.fold_change_data is None or self.fold_change_data.empty:
            return {"generated_at": datetime.now().isoformat(), "samples": []}

        records = []
        for _, row in self.fold_change_data.iterrows():
            sample_name = str(row.get('Sample_Name') or '').strip() or str(row.get('Target_Sample_Name') or '').strip() or str(row.get('Reference_Sample_Name') or '').strip() or "Unknown"
            day_match = re.search(r'(D\d+)', sample_name, flags=re.IGNORECASE)
            day = day_match.group(1).upper() if day_match else str(row.get('Day') or '').strip() or "Unknown"
            day_number = int(re.search(r'(\d+)', day).group(1)) if day.startswith('D') and re.search(r'(\d+)', day) else None

            target_cp = row.get('Target_Ct')
            reference_cp = row.get('Reference_Ct')

            record = {
                "well": row.get('Target_Well') or row.get('Reference_Well') or row.get('Sample_Name') or "Unknown",
                "sample_name": sample_name,
                "day": day,
                "day_number": day_number,
                "cp_values": {
                    "target_cp": float(target_cp) if pd.notna(target_cp) else None,
                    "reference_cp": float(reference_cp) if pd.notna(reference_cp) else None,
                },
                "fold_change": float(row.get('Fold_Change')) if pd.notna(row.get('Fold_Change')) else None,
            }
            records.append(record)

        return {
            "generated_at": datetime.now().isoformat(),
            "samples": records,
        }

    def export_results(self):
        """Run the full qPCR analysis and export a single PDF summary with plots and statistics."""
        if not self.target_file or not self.reference_file:
            messagebox.showerror("Fehler", "Bitte Target- und Reference-Datei auswählen.")
            return

        try:
            self.progress.start()
            self.update_status("Führe qPCR Analyse durch...")
            
            self.calculate_data()

            output_folder = get_output_directory()
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            sample_name = self.sample_name_var.get().strip() or None

# 1. Erst Plots generieren und Pfade speichern

            plot_paths_for_pdf = []

            p1_path = os.path.join(output_folder, f"qPCR_DeltaCt_viridis_{timestamp}.png")
            plot_delta_ct_with_colormap(
                self.delta_ct_data,
                'viridis',
                p1_path,
                sample_name=sample_name,
            )

            p2_path = os.path.join(output_folder, f"qPCR_FoldChange_plasma_{timestamp}.png")
            plot_fold_change_with_colormap(
                self.fold_change_data,
                'plasma',
                p2_path,
                sample_name=sample_name,
            )

            # 3. Schmelzkurve plotten
            p3_path = os.path.join(output_folder, f"qPCR_MeltingCurve_{timestamp}.png")
            melting_plot_exists = False
        
            if self.raw_file:
                self.update_status("Generiere Schmelzkurven-Plot...")
                melting_plot_exists = create_melting_curve_plot(self.raw_file, p3_path)

            # 4. Gruppierter Boxplot nach Tag/Day hinzufügen
            grouped_boxplot_path = os.path.join(output_folder, f"qPCR_GroupedDeltaCt_Boxplot_{timestamp}.png")
            grouped_boxplot_exists = False
            if 'Day' in self.delta_ct_data.columns:
                self.update_status("Generiere gruppierten Boxplot...")
                grouped_boxplot_exists = create_grouped_boxplot(
                    self.delta_ct_data,
                    grouped_boxplot_path,
                    value_col='Delta_Ct',
                    group_col='Day',
                    title='Delta Ct by Day'
                )

            volcano_plot_path = os.path.join(output_folder, f"qPCR_VolcanoPlot_{timestamp}.png")
            volcano_plot_exists = False
            if self.fold_change_data is not None and not self.fold_change_data.empty:
                self.update_status("Generiere Volcano-Plot...")
                volcano_plot_exists = create_volcano_plot(self.fold_change_data, volcano_plot_path)

            # --- Vorher stand hier group_col='Day' ---
            violin_plot_path = os.path.join(output_folder, f"qPCR_ViolinPlot_{timestamp}.png")
            violin_plot_exists = False
            if 'Sample_Name' in self.delta_ct_data.columns: # <--- Geändert auf Sample_Name
                self.update_status("Generiere Violin-Plot...")
                violin_plot_exists = create_violin_plot(
                    self.delta_ct_data,
                    violin_plot_path,
                    value_col='Delta_Ct',
                    group_col='Sample_Name',
                    title='Delta Ct Distribution by Sample'
                )

            # 5. PDF Erstellen
            summary_pdf = os.path.join(output_folder, f"qPCR_ExecutiveSummary_{timestamp}.pdf")
            create_qpcr_summary_pdf(
                self.delta_ct_data,
                self.fold_change_data,
                self.stats_report,
                summary_pdf,
                sample_name=sample_name,
                delta_ct_plot=p1_path,
                fold_change_plot=p2_path,
                melting_plot=p3_path if melting_plot_exists else None,
                grouped_boxplot=grouped_boxplot_path if grouped_boxplot_exists else None,
                volcano_plot=volcano_plot_path if volcano_plot_exists else None,
                violin_plot=violin_plot_path if violin_plot_exists else None,
                raw_target_data=self.target_data,
                raw_reference_data=self.reference_data,
            )

            analysis_json = self.build_analysis_json()
            summary_json = self.build_summary_json()

            combined_json = {
                "generated_at": datetime.now().isoformat(),
                "analysis_export": analysis_json,
                "summary_export": summary_json,
            }

            json_path = os.path.join(output_folder, f"qPCR_AnalysisExport_{timestamp}.json")
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(combined_json, f, indent=2, ensure_ascii=False)

            summary_json_path = os.path.join(output_folder, f"qPCR_SummaryExport_{timestamp}.json")
            with open(summary_json_path, 'w', encoding='utf-8') as f:
                json.dump(summary_json, f, indent=2, ensure_ascii=False)

            if os.path.exists(p1_path):
                os.remove(p1_path)
            if os.path.exists(p2_path):
                os.remove(p2_path)
            if os.path.exists(p3_path) and melting_plot_exists:
                os.remove(p3_path)
            if os.path.exists(grouped_boxplot_path) and grouped_boxplot_exists:
                os.remove(grouped_boxplot_path)
            if os.path.exists(volcano_plot_path) and volcano_plot_exists:
                os.remove(volcano_plot_path)
            if os.path.exists(violin_plot_path) and violin_plot_exists:
                os.remove(violin_plot_path)

            self.progress.stop()
            self.progress.set(1)
            self.update_status("Erfolgreich exportiert!")
            
            messagebox.showinfo("Erfolg", f"qPCR Auswertung erfolgreich!\n\nDateien gespeichert unter:\n{output_folder}")

        except Exception as e:
            self.progress.stop()
            self.progress.set(0)
            messagebox.showerror("Fehler", f"Fehler beim Exportieren:\n{str(e)}")
            self.update_status("Fehler bei der Analyse")