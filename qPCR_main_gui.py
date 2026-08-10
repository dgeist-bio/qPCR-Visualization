import sys
import os
import threading
import matplotlib.pyplot as plt
import customtkinter as ctk
import pandas as pd
import numpy as np
import json
import qPCR_visualizer

from pathlib import Path
from datetime import datetime
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy import stats
from scipy.stats import ttest_ind, f_oneway, mannwhitneyu

# Import qPCR modules
from qPCR_data_loader import load_qpcr_data, separate_genes, get_sample_groups
from qPCR_pdf_summary import qPCRSummaryMixin
from qPCR_delta_ct_calculation import calculate_delta_ct, calculate_fold_change
from qPCR_pdf_layout import create_qpcr_summary_pdf, QPCRPDFReport

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


class qPCRAnalyzerApp(ctk.CTk, qPCRSummaryMixin):
    def __init__(self):
        super().__init__()

        with open('config.json', 'r', encoding='utf-8') as f:
            self.config = json.load(f)

        self.title(f"qPCR Data Analyzer v{self.config['app_meta']['version']}")
        self.geometry(f"{self.config['gui_settings']['window_size']['width']}x{self.config['gui_settings']['window_size']['height']}")

        # Storage for selected files and data
        self.target_file = None
        self.reference_file = None
        self.target_data = None
        self.reference_data = None
        self.delta_ct_data = None
        self.fold_change_data = None
        self.stats_report = ""
        self.sample_name_var = ctk.StringVar(value="")
        self._export_thread = None
        self._export_result = None

        # --- Header ---
        self.header_label = ctk.CTkLabel(
            self,
            text="qPCR Data Analysis & Visualization",
            font=("Helvetica", 24, "bold"),
            text_color="#FFFFFF"
        )
        self.header_label.pack(pady=(30, 15))

        # --- Main Card (Centralized Container) ---
        main_card = ctk.CTkFrame(self, corner_radius=15, fg_color="#242424")
        main_card.pack(pady=10, padx=60, fill="x", expand=False)

        # --- Input Card (Files & Parameter) ---
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

        self.target_remove_btn = ctk.CTkButton(
            target_frame,
            text="Remove",
            command=self.remove_target_file,
            width=90,
            fg_color="#5C2F2F",
            hover_color="#7A3B3B"
        )
        self.target_remove_btn.pack(side="right")

        self.target_btn = ctk.CTkButton(
            target_frame,
            text="📁 Browse",
            command=self.select_target_file,
            width=110,
            fg_color="#3A3A3A",
            hover_color="#4A4A4A"
        )
        self.target_btn.pack(side="right", padx=(0, 8))

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

        self.reference_remove_btn = ctk.CTkButton(
            reference_frame,
            text="Remove",
            command=self.remove_reference_file,
            width=90,
            fg_color="#5C2F2F",
            hover_color="#7A3B3B"
        )
        self.reference_remove_btn.pack(side="right")

        self.reference_btn = ctk.CTkButton(
            reference_frame,
            text="📁 Browse",
            command=self.select_reference_file,
            width=110,
            fg_color="#3A3A3A",
            hover_color="#4A4A4A"
        )
        self.reference_btn.pack(side="right", padx=(0, 8))

        # Speicherort für die Raw-Datei initialisieren
        self.raw_file = None

        # --- Raw/Melting Curve File Selection (Optional) ---
        raw_frame = ctk.CTkFrame(input_card, fg_color="transparent")
        raw_frame.pack(fill="x", pady=5, padx=15)

        ctk.CTkLabel(raw_frame, text="Raw/Melting (opt.):", font=("Helvetica", 12), width=120, anchor="w").pack(side="left")
        self.raw_file_label = ctk.CTkLabel(
            raw_frame,
            text="Keine Datei ausgewählt (Optional)",
            font=("Helvetica", 11, "italic"),
            text_color="#8D8D8D",
            anchor="w"
        )
        self.raw_file_label.pack(side="left", padx=10, fill="x", expand=True)

        self.raw_remove_btn = ctk.CTkButton(
            raw_frame,
            text="Remove",
            command=self.remove_raw_file,
            width=90,
            fg_color="#5C2F2F",
            hover_color="#7A3B3B"
        )
        self.raw_remove_btn.pack(side="right")

        self.raw_btn = ctk.CTkButton(
            raw_frame,
            text="📁 Browse",
            command=self.select_raw_file,
            width=110,
            fg_color="#3A3A3A",
            hover_color="#4A4A4A"
        )
        self.raw_btn.pack(side="right", padx=(0, 8))

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

    def remove_target_file(self):
        """Clear the selected target gene file"""
        self.target_file = None
        self.target_file_label.configure(text="Keine Datei ausgewählt", text_color="#8D8D8D")
        self.update_status("Target file removed")

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

    def remove_reference_file(self):
        """Clear the selected reference gene file"""
        self.reference_file = None
        self.reference_file_label.configure(text="Keine Datei ausgewählt", text_color="#8D8D8D")
        self.update_status("Reference file removed")

    def select_raw_file(self):
        """Select optional Raw Data / Melting Curve file"""
        file_path = filedialog.askopenfilename(
            title="Select Raw Data / Melting Curve File (Optional)",
            filetypes=[("Text Files", "*.txt"), ("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        if file_path:
            self.raw_file = file_path
            self.raw_file_label.configure(text=os.path.basename(file_path), text_color="#2FA572")
            self.update_status("Raw/Melting file selected")

    def remove_raw_file(self):
        """Clear the selected raw file"""
        self.raw_file = None
        self.raw_file_label.configure(text="Keine Datei ausgewählt (Optional)", text_color="#8D8D8D")
        self.update_status("Raw/Melting file removed")

### Changed into qPCR_pdf_summary.py

    def update_status(self, message):
        """Update status label safely from the GUI thread or a worker thread."""
        timestamp = datetime.now().strftime("%H:%M:%S")
        if not self.winfo_exists():
            return
        if threading.current_thread() is threading.main_thread():
            self.progress_label.configure(text=f"Status: {message} ({timestamp})")
            self.update()
        else:
            self.after(0, lambda: self.progress_label.configure(text=f"Status: {message} ({timestamp})"))
            self.after(0, self.update)


if __name__ == "__main__":
    app = qPCRAnalyzerApp()
    app.mainloop()