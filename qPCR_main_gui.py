"""
qPCR Analysis GUI
Main interface for qPCR data analysis with file selection and cycle visualization
"""

import sys
import os
from pathlib import Path
from datetime import datetime
import customtkinter as ctk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import ttest_ind, f_oneway, mannwhitneyu

# Import qPCR modules
from qPCR_data_loader import load_qpcr_data, separate_genes, get_sample_groups
from qPCR_delta_ct_calculation import calculate_delta_ct, calculate_fold_change
from qPCR_visualizer import plot_delta_ct_with_colormap, plot_fold_change_with_colormap, save_combined_publication_quality_plots

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


class qPCRAnalyzerApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("qPCR Analysis Suite v1.1.0")
        self.geometry("1200x900")

        # Storage for selected files and data
        self.target_file = None
        self.reference_file = None
        self.target_data = None
        self.reference_data = None
        self.delta_ct_data = None
        self.fold_change_data = None

        # --- Header ---
        self.header_label = ctk.CTkLabel(
            self,
            text="qPCR Data Analysis & Cycle Visualization",
            font=("Arial", 22, "bold")
        )
        self.header_label.pack(pady=(20, 10))

        # --- File Selection Frame ---
        self.file_frame = ctk.CTkFrame(self)
        self.file_frame.pack(pady=15, padx=20, fill="x")

        ctk.CTkLabel(
            self.file_frame,
            text="Step 1: Select qPCR Result Files",
            font=("Arial", 14, "bold")
        ).pack(anchor="w", pady=(0, 10))

        # Target Gene File Selection
        target_frame = ctk.CTkFrame(self.file_frame)
        target_frame.pack(fill="x", pady=5)

        ctk.CTkLabel(target_frame, text="Target Gene File:", font=("Arial", 11)).pack(side="left", padx=10)
        self.target_file_label = ctk.CTkLabel(
            target_frame,
            text="No file selected",
            font=("Arial", 10),
            text_color="gray"
        )
        self.target_file_label.pack(side="left", padx=10, fill="x", expand=True)

        self.target_btn = ctk.CTkButton(
            target_frame,
            text="Browse",
            command=self.select_target_file,
            width=100
        )
        self.target_btn.pack(side="right", padx=10)

        # Reference Gene File Selection
        reference_frame = ctk.CTkFrame(self.file_frame)
        reference_frame.pack(fill="x", pady=5)

        ctk.CTkLabel(reference_frame, text="Reference Gene File:", font=("Arial", 11)).pack(side="left", padx=10)
        self.reference_file_label = ctk.CTkLabel(
            reference_frame,
            text="No file selected",
            font=("Arial", 10),
            text_color="gray"
        )
        self.reference_file_label.pack(side="left", padx=10, fill="x", expand=True)

        self.reference_btn = ctk.CTkButton(
            reference_frame,
            text="Browse",
            command=self.select_reference_file,
            width=100
        )
        self.reference_btn.pack(side="right", padx=10)

        # --- Analysis Controls ---
        control_frame = ctk.CTkFrame(self)
        control_frame.pack(pady=15, padx=20, fill="x")

        ctk.CTkLabel(
            control_frame,
            text="Step 2: Analysis & Visualization",
            font=("Arial", 14, "bold")
        ).pack(anchor="w", pady=(0, 10))

        button_frame = ctk.CTkFrame(control_frame)
        button_frame.pack(fill="x")

        self.calculate_btn = ctk.CTkButton(
            button_frame,
            text="Calculate ΔCt & Fold Change",
            command=self.calculate_data,
            font=("Arial", 12, "bold"),
            height=40
        )
        self.calculate_btn.pack(side="left", padx=5)

        self.visualize_btn = ctk.CTkButton(
            button_frame,
            text="Generate Cycle Plots",
            command=self.visualize_cycles,
            font=("Arial", 12, "bold"),
            height=40,
            state="disabled"
        )
        self.visualize_btn.pack(side="left", padx=5)

        self.advanced_btn = ctk.CTkButton(
            button_frame,
            text="Advanced Visualizations",
            command=self.visualize_advanced,
            font=("Arial", 12, "bold"),
            height=40,
            state="disabled"
        )
        self.advanced_btn.pack(side="left", padx=5)

        self.stats_btn = ctk.CTkButton(
            button_frame,
            text="Statistical Tests",
            command=self.run_statistical_tests,
            font=("Arial", 12, "bold"),
            height=40,
            state="disabled"
        )
        self.stats_btn.pack(side="left", padx=5)

        self.export_btn = ctk.CTkButton(
            button_frame,
            text="Export Results",
            command=self.export_results,
            font=("Arial", 12, "bold"),
            height=40,
            state="disabled"
        )
        self.export_btn.pack(side="left", padx=5)

        # --- Progress Info ---
        self.progress_label = ctk.CTkLabel(
            control_frame,
            text="Status: Ready",
            font=("Arial", 10, "italic"),
            text_color="gray"
        )
        self.progress_label.pack(anchor="w", pady=(10, 0))

        # --- Data Display Notebook (Tabs) ---
        self.tabview = ctk.CTkTabview(self)
        self.tabview.pack(pady=15, padx=20, fill="both", expand=True)

        # Tab 1: Cycle Data
        self.cycle_tab = self.tabview.add("Cycle Data Overview")
        self.cycle_text = ctk.CTkTextbox(self.cycle_tab, font=("Courier", 10))
        self.cycle_text.pack(fill="both", expand=True, padx=10, pady=10)

        # Tab 2: Delta Ct Results
        self.delta_ct_tab = self.tabview.add("ΔCt Results")
        self.delta_ct_text = ctk.CTkTextbox(self.delta_ct_tab, font=("Courier", 9))
        self.delta_ct_text.pack(fill="both", expand=True, padx=10, pady=10)

        # Tab 3: Fold Change Results
        self.fold_change_tab = self.tabview.add("Fold Change Results")
        self.fold_change_text = ctk.CTkTextbox(self.fold_change_tab, font=("Courier", 9))
        self.fold_change_text.pack(fill="both", expand=True, padx=10, pady=10)

        # Tab 4: Visualization
        self.viz_tab = self.tabview.add("Cycle Visualization")
        self.viz_frame = ctk.CTkFrame(self.viz_tab)
        self.viz_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Tab 5: Statistical Analysis
        self.stats_tab = self.tabview.add("Statistical Analysis")
        self.stats_text = ctk.CTkTextbox(self.stats_tab, font=("Courier", 9))
        self.stats_text.pack(fill="both", expand=True, padx=10, pady=10)

        # Tab 6: Advanced Visualizations
        self.advanced_viz_tab = self.tabview.add("Advanced Plots")
        self.advanced_viz_frame = ctk.CTkFrame(self.advanced_viz_tab)
        self.advanced_viz_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # --- Status Bar ---
        self.status_bar = ctk.CTkLabel(
            self,
            text="Speicherort: Desktop/qPCR_Ergebnisse",
            font=("Arial", 9, "italic"),
            text_color="gray"
        )
        self.status_bar.pack(side="bottom", pady=10)

    def select_target_file(self):
        """Select target gene data file"""
        file_path = filedialog.askopenfilename(
            title="Select Target Gene qPCR Results File",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        if file_path:
            self.target_file = file_path
            self.target_file_label.configure(text=os.path.basename(file_path), text_color="white")
            self.update_status("Target file selected")

    def select_reference_file(self):
        """Select reference gene data file"""
        file_path = filedialog.askopenfilename(
            title="Select Reference Gene qPCR Results File",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        if file_path:
            self.reference_file = file_path
            self.reference_file_label.configure(text=os.path.basename(file_path), text_color="white")
            self.update_status("Reference file selected")

    def calculate_data(self):
        """Load files and calculate ΔCt and Fold Change"""
        if not self.target_file or not self.reference_file:
            messagebox.showerror("Error", "Please select both target and reference files")
            return

        try:
            self.update_status("Loading data...")
            
            # Load data
            self.target_data = load_qpcr_data(self.target_file)
            self.reference_data = load_qpcr_data(self.reference_file)

            # Display loaded data
            self.display_cycle_data()

            # Calculate Delta Ct
            self.update_status("Calculating ΔCt...")
            self.delta_ct_data = calculate_delta_ct(self.target_data, self.reference_data)
            self.display_delta_ct()

            # Calculate Fold Change
            self.update_status("Calculating Fold Change...")
            self.fold_change_data = calculate_fold_change(self.delta_ct_data)
            self.display_fold_change()

            # Enable visualization buttons
            self.visualize_btn.configure(state="normal")
            self.export_btn.configure(state="normal")
            self.advanced_btn.configure(state="normal")
            self.stats_btn.configure(state="normal")

            self.update_status("Analysis complete!")
            messagebox.showinfo("Success", "ΔCt and Fold Change calculated successfully!")

        except Exception as e:
            messagebox.showerror("Error", f"Error during calculation:\n{str(e)}")
            self.update_status("Error during analysis")

    def display_cycle_data(self):
        """Display loaded cycle data in overview"""
        self.cycle_text.configure(state="normal")
        self.cycle_text.delete("1.0", "end")

        overview = f"""
{'='*80}
CYCLE DATA OVERVIEW
{'='*80}

TARGET GENE DATA:
{'-'*80}
Rows: {len(self.target_data)}

{self.target_data.to_string()}

{'='*80}

REFERENCE GENE DATA:
{'-'*80}
Rows: {len(self.reference_data)}

{self.reference_data.to_string()}

{'='*80}
"""
        self.cycle_text.insert("1.0", overview)
        self.cycle_text.configure(state="disabled")

    def display_delta_ct(self):
        """Display Delta Ct calculations"""
        self.delta_ct_text.configure(state="normal")
        self.delta_ct_text.delete("1.0", "end")

        summary = f"""
{'='*80}
DELTA CT (ΔCt) ANALYSIS
{'='*80}
Formula: ΔCt = Target Ct - Reference Ct

{self.delta_ct_data.to_string()}

{'='*80}

STATISTICS:
{'-'*80}
Mean ΔCt: {self.delta_ct_data['Delta_Ct'].mean():.3f}
Std Dev ΔCt: {self.delta_ct_data['Delta_Ct'].std():.3f}
Min ΔCt: {self.delta_ct_data['Delta_Ct'].min():.3f}
Max ΔCt: {self.delta_ct_data['Delta_Ct'].max():.3f}

{'='*80}
"""
        self.delta_ct_text.insert("1.0", summary)
        self.delta_ct_text.configure(state="disabled")

    def display_fold_change(self):
        """Display Fold Change calculations"""
        self.fold_change_text.configure(state="normal")
        self.fold_change_text.delete("1.0", "end")

        summary = f"""
{'='*80}
FOLD CHANGE (2^-ΔΔCt) ANALYSIS
{'='*80}
Formula: Fold Change = 2^-(ΔCt_sample - ΔCt_control)

{self.fold_change_data.to_string()}

{'='*80}

STATISTICS:
{'-'*80}
Mean Fold Change: {self.fold_change_data['Fold_Change'].mean():.3f}
Std Dev Fold Change: {self.fold_change_data['Fold_Change'].std():.3f}
Min Fold Change: {self.fold_change_data['Fold_Change'].min():.3f}
Max Fold Change: {self.fold_change_data['Fold_Change'].max():.3f}

{'='*80}
"""
        self.fold_change_text.insert("1.0", summary)
        self.fold_change_text.configure(state="disabled")

    def visualize_cycles(self):
        """Create and display cycle visualization"""
        if self.delta_ct_data is None:
            messagebox.showerror("Error", "Please calculate data first")
            return

        try:
            self.update_status("Generating cycle visualizations...")

            # Clear previous plots
            for widget in self.viz_frame.winfo_children():
                widget.destroy()

            # Create figure with subplots
            fig = Figure(figsize=(16, 10), dpi=100)

            # Plot 1: Target vs Reference Ct
            ax1 = fig.add_subplot(2, 3, 1)
            sample_labels = [f"S{int(n)}" for n in self.delta_ct_data['Sample_Number']]
            x_pos = np.arange(len(self.delta_ct_data))
            ax1.bar(x_pos - 0.2, self.delta_ct_data['Target_Ct'], 0.4, label='Target', alpha=0.8, color='#3498DB')
            ax1.bar(x_pos + 0.2, self.delta_ct_data['Reference_Ct'], 0.4, label='Reference', alpha=0.8, color='#E74C3C')
            ax1.set_xlabel('Sample')
            ax1.set_ylabel('Ct Value')
            ax1.set_title('Target vs Reference Ct Values')
            ax1.set_xticks(x_pos)
            ax1.set_xticklabels(sample_labels, rotation=45)
            ax1.legend()
            ax1.grid(axis='y', alpha=0.3)

            # Plot 2: Delta Ct with error bars
            ax2 = fig.add_subplot(2, 3, 2)
            colors = plt.cm.viridis(np.linspace(0, 1, len(self.delta_ct_data)))
            ax2.bar(x_pos, self.delta_ct_data['Delta_Ct'], yerr=self.delta_ct_data['Delta_Ct_STD'],
                   capsize=5, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
            ax2.set_xlabel('Sample')
            ax2.set_ylabel('ΔCt')
            ax2.set_title('Delta Ct with Error Bars')
            ax2.set_xticks(x_pos)
            ax2.set_xticklabels(sample_labels, rotation=45)
            ax2.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.5)
            ax2.grid(axis='y', alpha=0.3)

            # Plot 3: Fold Change with error bars
            ax3 = fig.add_subplot(2, 3, 3)
            if 'Fold_Change' in self.fold_change_data.columns:
                colors_fc = plt.cm.plasma(np.linspace(0, 1, len(self.fold_change_data)))
                ax3.bar(x_pos, self.fold_change_data['Fold_Change'], 
                       yerr=self.fold_change_data['Fold_Change_STD'],
                       capsize=5, color=colors_fc, edgecolor='black', linewidth=1.5, alpha=0.8)
                ax3.set_xlabel('Sample')
                ax3.set_ylabel('Fold Change (2^-ΔΔCt)')
                ax3.set_title('Fold Change Analysis')
                ax3.set_xticks(x_pos)
                ax3.set_xticklabels(sample_labels, rotation=45)
                ax3.axhline(y=1, color='red', linestyle='--', linewidth=2, alpha=0.5, label='No change')
                ax3.legend()
                ax3.grid(axis='y', alpha=0.3)

            # Plot 4: Scatter plot - Target vs Reference
            ax4 = fig.add_subplot(2, 3, 4)
            scatter = ax4.scatter(self.delta_ct_data['Target_Ct'], self.delta_ct_data['Reference_Ct'],
                                 s=100, alpha=0.6, c=range(len(self.delta_ct_data)), cmap='cool')
            ax4.set_xlabel('Target Ct')
            ax4.set_ylabel('Reference Ct')
            ax4.set_title('Target vs Reference Ct Correlation')
            fig.colorbar(scatter, ax=ax4, label='Sample')
            ax4.grid(True, alpha=0.3)

            # Plot 5: Distribution of Ct values
            ax5 = fig.add_subplot(2, 3, 5)
            ax5.boxplot([self.delta_ct_data['Target_Ct'].dropna(), self.delta_ct_data['Reference_Ct'].dropna()],
                       tick_labels=['Target', 'Reference'])
            ax5.set_ylabel('Ct Value')
            ax5.set_title('Ct Value Distribution')
            ax5.grid(axis='y', alpha=0.3)

            # Plot 6: Delta Ct distribution
            ax6 = fig.add_subplot(2, 3, 6)
            ax6.hist(self.delta_ct_data['Delta_Ct'].dropna(), bins=10, color='#2ECC71', alpha=0.7, edgecolor='black')
            ax6.set_xlabel('ΔCt Value')
            ax6.set_ylabel('Frequency')
            ax6.set_title('ΔCt Distribution')
            ax6.grid(axis='y', alpha=0.3)

            fig.subplots_adjust(hspace=0.75, wspace=0.45, left=0.1, right=0.9, top=0.9, bottom=0.1)

            # Embed in tkinter
            canvas = FigureCanvasTkAgg(fig, master=self.viz_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)

            self.update_status("Cycle visualizations complete!")

        except Exception as e:
            messagebox.showerror("Error", f"Error creating visualization:\n{str(e)}")
            self.update_status("Error during visualization")

    def export_results(self):
        """Export results to CSV and generate high-quality plots"""
        if self.delta_ct_data is None:
            messagebox.showerror("Error", "Please calculate data first")
            return

        try:
            # Create output folder
            output_folder = os.path.join(os.path.expanduser("~"), "Desktop", "qPCR_Ergebnisse")
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # Export Delta Ct
            delta_ct_file = os.path.join(output_folder, f"qPCR_DeltaCt_{timestamp}.csv")
            self.delta_ct_data.to_csv(delta_ct_file, index=False)
            print(f"✓ Saved: {delta_ct_file}")

            # Export Fold Change
            if self.fold_change_data is not None:
                fold_change_file = os.path.join(output_folder, f"qPCR_FoldChange_{timestamp}.csv")
                self.fold_change_data.to_csv(fold_change_file, index=False)
                print(f"✓ Saved: {fold_change_file}")

            # Export Statistical Report
            stats_report = self.stats_text.get("1.0", "end")
            if stats_report.strip():
                stats_file = os.path.join(output_folder, f"qPCR_StatisticalAnalysis_{timestamp}.txt")
                with open(stats_file, 'w') as f:
                    f.write(stats_report)
                print(f"✓ Saved: {stats_file}")

            # Generate publication-quality plots
            plot_delta_ct_with_colormap(self.delta_ct_data, 'viridis',
                                       os.path.join(output_folder, f"qPCR_DeltaCt_viridis_{timestamp}.png"))
            plot_delta_ct_with_colormap(self.delta_ct_data, 'cividis',
                                       os.path.join(output_folder, f"qPCR_DeltaCt_cividis_{timestamp}.png"))
            plot_delta_ct_with_colormap(self.delta_ct_data, 'plasma',
                                       os.path.join(output_folder, f"qPCR_DeltaCt_plasma_{timestamp}.png"))

            if self.fold_change_data is not None:
                plot_fold_change_with_colormap(self.fold_change_data, 'plasma',
                                              os.path.join(output_folder, f"qPCR_FoldChange_plasma_{timestamp}.png"))
                plot_fold_change_with_colormap(self.fold_change_data, 'viridis',
                                              os.path.join(output_folder, f"qPCR_FoldChange_viridis_{timestamp}.png"))

            combined_png = os.path.join(output_folder, f"qPCR_Publication_Quality_AllPlots_{timestamp}.png")
            save_combined_publication_quality_plots(self.delta_ct_data, self.fold_change_data, combined_png)
            combined_pdf = os.path.join(output_folder, f"qPCR_Publication_Quality_AllPlots_{timestamp}.pdf")
            save_combined_publication_quality_plots(self.delta_ct_data, self.fold_change_data, combined_pdf)

            self.update_status(f"Results exported to {output_folder}")
            messagebox.showinfo("Success", f"Results exported successfully!\n\nFolder: {output_folder}")

        except Exception as e:
            messagebox.showerror("Error", f"Error exporting results:\n{str(e)}")

    def run_statistical_tests(self):
        """Perform statistical tests on qPCR data"""
        if self.delta_ct_data is None or self.fold_change_data is None:
            messagebox.showerror("Error", "Please calculate data first")
            return

        try:
            self.update_status("Running statistical tests...")

            self.stats_text.configure(state="normal")
            self.stats_text.delete("1.0", "end")

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
            # Paired t-test for Target vs Reference
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

            # Normality tests
            stats_report += f"""
{'='*80}
4. NORMALITY TESTS (Shapiro-Wilk)
{'='*80}
"""
            if len(target_vals) >= 3:
                w_stat_t, p_val_t = stats.shapiro(target_vals)
                stats_report += f"Target Ct:  W={w_stat_t:.4f}, p-value={p_val_t:.6f}, Normal: {'YES' if p_val_t > 0.05 else 'NO'}\n"

            if len(reference_vals) >= 3:
                w_stat_r, p_val_r = stats.shapiro(reference_vals)
                stats_report += f"Reference Ct: W={w_stat_r:.4f}, p-value={p_val_r:.6f}, Normal: {'YES' if p_val_r > 0.05 else 'NO'}\n"

            delta_ct_vals = self.delta_ct_data['Delta_Ct'].dropna()
            if len(delta_ct_vals) >= 3:
                w_stat_d, p_val_d = stats.shapiro(delta_ct_vals)
                stats_report += f"Delta Ct:   W={w_stat_d:.4f}, p-value={p_val_d:.6f}, Normal: {'YES' if p_val_d > 0.05 else 'NO'}\n"

            # Correlation analysis
            stats_report += f"""
{'='*80}
5. CORRELATION ANALYSIS
{'='*80}
Target vs Reference Ct:
"""
            if len(target_vals) > 2 and len(reference_vals) > 2:
                pearson_r, pearson_p = stats.pearsonr(target_vals, reference_vals)
                spearman_r, spearman_p = stats.spearmanr(target_vals, reference_vals)
                stats_report += f"""
  Pearson r:          {pearson_r:.4f}
  p-value:            {pearson_p:.6f}
  Spearman ρ:         {spearman_r:.4f}
  p-value:            {spearman_p:.6f}
"""

            # Outlier detection using IQR
            stats_report += f"""
{'='*80}
6. OUTLIER DETECTION (IQR Method)
{'='*80}
"""
            for data_name, data_vals in [('Target Ct', target_vals), ('Reference Ct', reference_vals), ('Delta Ct', delta_ct_vals)]:
                Q1 = data_vals.quantile(0.25)
                Q3 = data_vals.quantile(0.75)
                IQR = Q3 - Q1
                lower_bound = Q1 - 1.5 * IQR
                upper_bound = Q3 + 1.5 * IQR
                outliers = data_vals[(data_vals < lower_bound) | (data_vals > upper_bound)]
                stats_report += f"\n{data_name}:\n"
                stats_report += f"  IQR: {IQR:.4f}, Bounds: [{lower_bound:.4f}, {upper_bound:.4f}]\n"
                stats_report += f"  Outliers: {len(outliers)} detected\n"
                if len(outliers) > 0:
                    stats_report += f"  Outlier values: {outliers.values}\n"

            stats_report += f"""
{'='*80}
NOTES:
- p-value < 0.05 indicates statistical significance (95% confidence)
- Normality test: if p > 0.05, data is normally distributed
- Outliers identified using Interquartile Range (IQR) method
{'='*80}
"""

            self.stats_text.insert("1.0", stats_report)
            self.stats_text.configure(state="disabled")

            self.update_status("Statistical tests complete!")

        except Exception as e:
            messagebox.showerror("Error", f"Error running statistical tests:\n{str(e)}")
            self.update_status("Error during statistical analysis")

    def visualize_advanced(self):
        """Create advanced visualizations with multiple styles"""
        if self.delta_ct_data is None:
            messagebox.showerror("Error", "Please calculate data first")
            return

        try:
            self.update_status("Generating advanced visualizations...")

            # Clear previous plots
            for widget in self.advanced_viz_frame.winfo_children():
                widget.destroy()

            # Create figure with advanced plots
            fig = Figure(figsize=(17, 11), dpi=100)

            sample_labels = [f"S{int(n)}" for n in self.delta_ct_data['Sample_Number']]
            x_pos = np.arange(len(self.delta_ct_data))

            # Plot 1: Violin plot
            ax1 = fig.add_subplot(2, 3, 1)
            target_vals = self.delta_ct_data['Target_Ct'].values
            reference_vals = self.delta_ct_data['Reference_Ct'].values
            parts = ax1.violinplot([target_vals, reference_vals], positions=[1, 2], showmeans=True, showmedians=True)
            ax1.set_xticks([1, 2])
            ax1.set_xticklabels(['Target', 'Reference'])
            ax1.set_ylabel('Ct Value')
            ax1.set_title('Violin Plot: Ct Distribution')
            ax1.grid(axis='y', alpha=0.3)

            # Plot 2: Delta Ct with different colormap
            ax2 = fig.add_subplot(2, 3, 2)
            colors_plasma = plt.cm.plasma(np.linspace(0, 1, len(self.delta_ct_data)))
            bars = ax2.bar(x_pos, self.delta_ct_data['Delta_Ct'], color=colors_plasma, edgecolor='black', linewidth=1.5, alpha=0.8)
            ax2.set_xlabel('Sample')
            ax2.set_ylabel('ΔCt')
            ax2.set_title('ΔCt (Plasma Colormap)')
            ax2.set_xticks(x_pos)
            ax2.set_xticklabels(sample_labels, rotation=45)
            ax2.axhline(y=self.delta_ct_data['Delta_Ct'].mean(), color='red', linestyle='--', linewidth=2, alpha=0.7, label='Mean')
            ax2.legend()
            ax2.grid(axis='y', alpha=0.3)

            # Plot 3: Error distribution
            ax3 = fig.add_subplot(2, 3, 3)
            colors_cool = plt.cm.cool(np.linspace(0, 1, len(self.delta_ct_data)))
            ax3.bar(x_pos, self.delta_ct_data['Delta_Ct_STD'], color=colors_cool, edgecolor='black', alpha=0.8)
            ax3.set_xlabel('Sample')
            ax3.set_ylabel('Standard Deviation')
            ax3.set_title('ΔCt Error Distribution')
            ax3.set_xticks(x_pos)
            ax3.set_xticklabels(sample_labels, rotation=45)
            ax3.grid(axis='y', alpha=0.3)

            # Plot 4: Q-Q plot for normality check
            ax4 = fig.add_subplot(2, 3, 4)
            delta_ct_clean = self.delta_ct_data['Delta_Ct'].dropna()
            if len(delta_ct_clean) >= 3:
                stats.probplot(delta_ct_clean, dist="norm", plot=ax4)
                ax4.set_title('Q-Q Plot: ΔCt Normality Check')
                ax4.grid(True, alpha=0.3)
            else:
                ax4.text(0.5, 0.5, 'Insufficient data for Q-Q plot', ha='center', va='center')

            # Plot 5: Fold Change comparison (if available)
            ax5 = fig.add_subplot(2, 3, 5)
            if 'Fold_Change' in self.fold_change_data.columns:
                colors_magma = plt.cm.magma(np.linspace(0, 1, len(self.fold_change_data)))
                ax5.bar(x_pos, self.fold_change_data['Fold_Change'], color=colors_magma, edgecolor='black', linewidth=1.5, alpha=0.8)
                ax5.set_xlabel('Sample')
                ax5.set_ylabel('Fold Change')
                ax5.set_title('Fold Change (Magma Colormap)')
                ax5.set_xticks(x_pos)
                ax5.set_xticklabels(sample_labels, rotation=45)
                ax5.axhline(y=1, color='red', linestyle='--', linewidth=2, alpha=0.7, label='No change (FC=1)')
                ax5.legend()
                ax5.grid(axis='y', alpha=0.3)

            # Plot 6: Ct values with confidence intervals
            ax6 = fig.add_subplot(2, 3, 6)
            # Calculate 95% CI for means
            target_mean = self.delta_ct_data['Target_Ct'].mean()
            reference_mean = self.delta_ct_data['Reference_Ct'].mean()
            target_se = self.delta_ct_data['Target_Ct'].std() / np.sqrt(len(self.delta_ct_data))
            reference_se = self.delta_ct_data['Reference_Ct'].std() / np.sqrt(len(self.delta_ct_data))
            ci = 1.96  # 95% CI

            means = [target_mean, reference_mean]
            errors = [ci * target_se, ci * reference_se]
            colors_viridis = ['#3498DB', '#E74C3C']

            ax6.bar([0, 1], means, yerr=errors, capsize=10, color=colors_viridis, alpha=0.7, edgecolor='black', linewidth=1.5)
            ax6.set_xticks([0, 1])
            ax6.set_xticklabels(['Target', 'Reference'])
            ax6.set_ylabel('Ct Value')
            ax6.set_title('Mean Ct ± 95% CI')
            ax6.grid(axis='y', alpha=0.3)

            fig.subplots_adjust(hspace=0.75, wspace=0.45, left=0.1, right=0.9, top=0.9, bottom=0.1)

            # Embed in tkinter
            canvas = FigureCanvasTkAgg(fig, master=self.advanced_viz_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)

            self.update_status("Advanced visualizations complete!")

        except Exception as e:
            messagebox.showerror("Error", f"Error creating advanced visualization:\n{str(e)}")
            self.update_status("Error during advanced visualization")

    def update_status(self, message):
        """Update status label"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.progress_label.configure(text=f"Status: {message} ({timestamp})")
        self.update()

class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.running = True # <--- Das ist dein Lebenszeichen-Flag
        self.protocol("WM_DELETE_WINDOW", self.on_closing) # <--- Abfangen des Schließen-Buttons
        self.update_loop()

    def on_closing(self):
        self.running = False # Stoppt die Schleife
        self.destroy() # Schließt das Fenster sauber

    def update_loop(self):
        if self.running:
            # Dein Code, der die GUI aktualisiert
            ...
            self.after(100, self.update_loop) # Erst hier wird der nächste Loop geplant

if __name__ == "__main__":
    app = qPCRAnalyzerApp()
    app.mainloop()
