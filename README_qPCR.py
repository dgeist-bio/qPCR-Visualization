"""
qPCR Analysis Suite - Quick Start Guide
======================================

Three Python scripts for automated qPCR analysis with delta Ct calculation.

FILES CREATED:
==============

1. qPCR_data_loader.py
   - Loads and separates qPCR data
   - Parses sample information
   - Handles irregular file formatting

2. qPCR_delta_ct_calculation.py
   - Calculates ΔCt = Target Ct - Reference Ct
   - Calculates Fold Change = 2^(-ΔΔCt)
   - Generates visualizations with viridis colormap by default
   - Outputs CSV files with results

3. qPCR_visualizer.py
   - Re-generates plots with different colormaps
   - Supports: viridis, magma, plasma, inferno, cividis
   - Useful for finding the most readable colormap for your vision


WORKFLOW:
=========

Step 1: Run the calculation
   

   This will:
   ✓ Load "230913 qPCR and SD.txt"
   ✓ Separate Target (A/B wells) and Reference (C/D wells) genes
   ✓ Calculate ΔCt for each sample
   ✓ Calculate Fold Change (normalized to sample 1)
   ✓ Save results to CSV files
   ✓ Generate bar plots with viridis colormap

   OUTPUT FILES:
   - qPCR_results_delta_ct.csv          (ΔCt values)
   - qPCR_results_fold_change.csv       (Fold change values)
   - qPCR_results_delta_ct_viridis.png  (ΔCt visualization)
   - qPCR_results_fold_change_viridis.png (Fold change visualization)


Step 2: (Optional) Re-generate with different colormaps
   

   This will generate plots with ALL accessible colormaps:
   - viridis  (perceptually uniform, colorblind-friendly)
   - magma    (dark to bright, high contrast)
   - plasma   (vibrant, high contrast)
   - inferno  (warm tones, high contrast)
   - cividis (optimized for deuteranopia colorblindness)


CUSTOMIZATION:
==============

To change the colormap in qPCR_delta_ct_calculation.py:

Open the file and change line 247:
   COLORMAP = 'viridis'  # Change to 'magma', 'plasma', 'inferno', or 'cividis'

Then run the script again.


CALCULATION DETAILS:
====================

Data Structure:
  - Wells A1, B1: Sample 1, Target gene (your gene of interest)
  - Wells C1, D1: Sample 1, Reference gene (housekeeping gene)
  - Wells A2, B2: Sample 2, Target gene
  - Wells C2, D2: Sample 2, Reference gene
  ... and so on

ΔCt Calculation (what happens):
  1. For each sample, calculate: ΔCt = Target Ct - Reference Ct
  2. This normalizes for differences in reference gene expression
  3. Error bars show propagated standard deviations

Fold Change Calculation (what happens):
  1. Calculate ΔΔCt = Sample ΔCt - Control ΔCt (normalized to sample 1)
  2. Calculate Fold Change = 2^(-ΔΔCt)
  3. Values > 1 mean upregulation
  4. Values < 1 mean downregulation
  5. Error bars show propagated uncertainties


OUTPUT INTERPRETATION:
======================

CSV Files contain:
  • Sample_Number: Sample identifier
  • Sample_Name: Well positions (e.g., A1_B1)
  • Target_Ct: Ct value for your gene of interest
  • Reference_Ct: Ct value for housekeeping gene
  • Delta_Ct: Ct difference (your main normalization value)
  • Delta_Ct_STD: Standard deviation of ΔCt
  • Fold_Change: 2^(-ΔΔCt) relative to sample 1
  • Fold_Change_STD: Standard deviation of fold change


COLORMAP GUIDE FOR NEURODIVERGENT USERS:
========================================

viridis:   Purple → Blue → Green → Yellow
           Good for: General use, colorblind-friendly, widely accepted
           Best if: You like gradual, smooth transitions

magma:     Black/Purple → Red → Yellow/White
           Good for: High contrast, easier to distinguish bars
           Best if: You prefer warm-to-cool transitions or high contrast

plasma:    Purple → Red/Orange → Yellow
           Good for: Very vibrant, maximum contrast
           Best if: You need the highest visual differentiation

inferno:   Black → Dark Purple → Orange → Yellow
           Good for: Warm tones, good for print
           Best if: You prefer warm colors

cividis:   Blue → Cyan → Yellow
           Good for: Colorblind users (especially deuteranopia)
           Best if: You have color vision deficiency or want maximum accessibility


COMMON ISSUES & FIXES:
====================

Q: ModuleNotFoundError: No module named 'pandas'
A: Packages are already installed. Make sure to use the full Python path:

Q: Can't find the input file
A: Make sure "230913 qPCR and SD.txt" is in the same directory as the scripts

Q: Visualizations don't open
A: PNG files are saved automatically. They don't pop up in a window.
   Check the working directory for the .png files.

Q: I want to use a different sample as the control for fold change
A: Edit line 259 in qPCR_delta_ct_calculation.py and change:
   fold_change_df = calculate_fold_change(delta_ct_df)
   to:
   fold_change_df = calculate_fold_change(delta_ct_df, reference_sample=3)  # Use sample 3 as control


FOR MORE HELP:
=============
Check the docstrings in each Python file for detailed function documentation.
The code is well-commented and designed to be easy to modify.
"""

# This file is for documentation - not meant to be run
# Just open it in a text editor or use: cat README_qPCR.txt
