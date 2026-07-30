import os
from datetime import datetime
from tkinter import messagebox
from scipy.stats import ttest_ind, mannwhitneyu

from qPCR_data_loader import load_qpcr_data
from qPCR_delta_ct_calculation import calculate_delta_ct, calculate_fold_change
from qPCR_visualizer import (
    plot_delta_ct_with_colormap,
    plot_fold_change_with_colormap,
    create_qpcr_summary_pdf,
    get_output_directory,
)

class qPCRSummaryMixin:
    """Mix-in class containing calculation and export logic for qPCRAnalyzerApp."""
    
    def calculate_data(self):
        """Load files and calculate ΔCt and Fold Change."""
        if not self.target_file or not self.reference_file:
            raise ValueError("Please select both target and reference files")

        self.update_status("Loading data...")
        self.target_data = load_qpcr_data(self.target_file)
        self.reference_data = load_qpcr_data(self.reference_file)

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

            delta_ct_file = os.path.join(output_folder, f"qPCR_DeltaCt_{timestamp}.csv")
            self.delta_ct_data.to_csv(delta_ct_file, index=False)

            if self.fold_change_data is not None:
                fold_change_file = os.path.join(output_folder, f"qPCR_FoldChange_{timestamp}.csv")
                self.fold_change_data.to_csv(fold_change_file, index=False)

            stats_file = os.path.join(output_folder, f"qPCR_StatisticalAnalysis_{timestamp}.txt")
            with open(stats_file, 'w', encoding='utf-8') as handle:
                handle.write(self.stats_report)

            summary_pdf = os.path.join(output_folder, f"qPCR_ExecutiveSummary_{timestamp}.pdf")
            create_qpcr_summary_pdf(
                self.delta_ct_data,
                self.fold_change_data,
                self.stats_report,
                summary_pdf,
                sample_name=sample_name,
            )

            plot_delta_ct_with_colormap(
                self.delta_ct_data,
                'viridis',
                os.path.join(output_folder, f"qPCR_DeltaCt_viridis_{timestamp}.png"),
                sample_name=sample_name,
            )
            plot_fold_change_with_colormap(
                self.fold_change_data,
                'plasma',
                os.path.join(output_folder, f"qPCR_FoldChange_plasma_{timestamp}.png"),
                sample_name=sample_name,
            )
            self.progress.stop()
            self.progress.set(1)
            self.update_status("Erfolgreich exportiert!")
            
            messagebox.showinfo("Erfolg", f"qPCR Auswertung erfolgreich!\n\nDateien gespeichert unter:\n{output_folder}")

        except Exception as e:
            self.progress.stop()
            self.progress.set(0)
            messagebox.showerror("Fehler", f"Fehler beim Exportieren:\n{str(e)}")
            self.update_status("Fehler bei der Analyse")