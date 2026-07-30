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
from qPCR_pdf_summary import qPCRSummaryMixin
from qPCR_delta_ct_calculation import calculate_delta_ct, calculate_fold_change
from qPCR_visualizer import (
    plot_delta_ct_with_colormap,
    plot_fold_change_with_colormap,
    save_combined_publication_quality_plots,
    create_qpcr_summary_pdf,
    get_output_directory,
)

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


class qPCRAnalyzerApp(ctk.CTk, qPCRSummaryMixin):
    def __init__(self):
        super().__init__()

        self.title("qPCR Analysis Suite v1.3.0 - Professional Edition")
        self.geometry("800x600")
        self.configure(fg_color="#1A1A1A") # Tiefdunkler Hintergrund wie beim MTT Analyzer

        # Storage for selected files and data
        self.target_file = None
        self.reference_file = None
        self.target_data = None
        self.reference_data = None
        self.delta_ct_data = None
        self.fold_change_data = None
        self.stats_report = ""
        self.sample_name_var = ctk.StringVar(value="")

        # --- Header ---
        self.header_label = ctk.CTkLabel(
            self,
            text="qPCR Data Analysis & Visualization",
            font=("Helvetica", 24, "bold"),
            text_color="#FFFFFF"
        )
        self.header_label.pack(pady=(30, 15))

        # --- Main Card (Zentrierter Container) ---
        main_card = ctk.CTkFrame(self, corner_radius=15, fg_color="#242424")
        main_card.pack(pady=10, padx=60, fill="x", expand=False)

        # --- Input Card (Dateien & Parameter) ---
        input_card = ctk.CTkFrame(main_card, corner_radius=12, fg_color="#2B2B2B")
        input_card.pack(pady=20, padx=30, fill="x")

        ctk.CTkLabel(
            input_card,
            text="Schritt 1: Daten-Import & Parameter",
            font=("Helvetica", 14, "bold"),
            text_color="#E0E0E0"
        ).pack(anchor="w", padx=15, pady=(15, 10))

        # Target Gene File Selection
        target_frame = ctk.CTkFrame(input_card, fg_color="transparent")
        target_frame.pack(fill="x", pady=5, padx=15)

        ctk.CTkLabel(target_frame, text="Target Gene:", font=("Helvetica", 12), width=100, anchor="w").pack(side="left")
        self.target_file_label = ctk.CTkLabel(
            target_frame,
            text="Keine Datei ausgewählt",
            font=("Helvetica", 11, "italic"),
            text_color="#8D8D8D",
            anchor="w"
        )
        self.target_file_label.pack(side="left", padx=10, fill="x", expand=True)

        self.target_btn = ctk.CTkButton(
            target_frame,
            text="📁 Browse",
            command=self.select_target_file,
            width=110,
            fg_color="#3A3A3A",
            hover_color="#4A4A4A"
        )
        self.target_btn.pack(side="right")

        # Reference Gene File Selection
        reference_frame = ctk.CTkFrame(input_card, fg_color="transparent")
        reference_frame.pack(fill="x", pady=5, padx=15)

        ctk.CTkLabel(reference_frame, text="Reference Gene:", font=("Helvetica", 12), width=100, anchor="w").pack(side="left")
        self.reference_file_label = ctk.CTkLabel(
            reference_frame,
            text="Keine Datei ausgewählt",
            font=("Helvetica", 11, "italic"),
            text_color="#8D8D8D",
            anchor="w"
        )
        self.reference_file_label.pack(side="left", padx=10, fill="x", expand=True)

        self.reference_btn = ctk.CTkButton(
            reference_frame,
            text="📁 Browse",
            command=self.select_reference_file,
            width=110,
            fg_color="#3A3A3A",
            hover_color="#4A4A4A"
        )
        self.reference_btn.pack(side="right")

        # Sample Name
        sample_frame = ctk.CTkFrame(input_card, fg_color="transparent")
        sample_frame.pack(fill="x", pady=(15, 15), padx=15)

        ctk.CTkLabel(sample_frame, text="Sample Name:", font=("Helvetica", 12), width=100, anchor="w").pack(side="left")
        self.sample_name_entry = ctk.CTkEntry(
            sample_frame,
            textvariable=self.sample_name_var,
            placeholder_text="z. B. Wildtype_3h",
            corner_radius=8,
            height=32
        )
        self.sample_name_entry.pack(side="left", fill="x", expand=True, padx=(10, 0))

        # --- Export Button ---
        self.export_btn = ctk.CTkButton(
            main_card,
            text="Datei auswählen & analysieren",
            command=self.export_results,
            font=("Helvetica", 14, "bold"),
            height=45,
            width=280,
            corner_radius=10,
            fg_color="#1F6AA5",  # Stahlblau (MTT Style)
            hover_color="#144A75"
        )
        self.export_btn.pack(pady=(15, 20))

        # --- Progress & Status ---
        self.progress = ctk.CTkProgressBar(main_card, width=280, height=8, progress_color="#2FA572")
        self.progress.set(0)
        self.progress.pack(pady=(0, 10))

        self.progress_label = ctk.CTkLabel(
            main_card,
            text="Status: Bereit für Analyse",
            font=("Helvetica", 11, "italic"),
            text_color="#8D8D8D"
        )
        self.progress_label.pack(pady=(0, 20))

        # --- Footer ---
        self.status_bar = ctk.CTkLabel(
            self,
            text="Speicherort: Desktop/qPCR_Ergebnisse",
            font=("Helvetica", 11, "italic"),
            text_color="#7A7A7A"
        )
        self.status_bar.pack(side="bottom", pady=15)

    def select_target_file(self):
        """Select target gene data file"""
        file_path = filedialog.askopenfilename(
            title="Select Target Gene qPCR Results File",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        if file_path:
            self.target_file = file_path
            self.target_file_label.configure(text=os.path.basename(file_path), text_color="#2FA572")
            self.update_status("Target file selected")

    def select_reference_file(self):
        """Select reference gene data file"""
        file_path = filedialog.askopenfilename(
            title="Select Reference Gene qPCR Results File",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        if file_path:
            self.reference_file = file_path
            self.reference_file_label.configure(text=os.path.basename(file_path), text_color="#2FA572")
            self.update_status("Reference file selected")

### Changed into qPCR_pdf_summary.py

    def update_status(self, message):
        """Update status label"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.progress_label.configure(text=f"Status: {message} ({timestamp})")
        self.update()


if __name__ == "__main__":
    app = qPCRAnalyzerApp()
    app.mainloop()