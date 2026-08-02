# Changelog

All notable changes to this project will be documented in this file.

### [1.4.0] - 2026-08-02

- **Added a "Remove" Button in the GUI for the selected files
- **Added a further line for importing raw/melting raw data (will/can be removed in the future)
- **Added a JSON export (for a Blender pipeline)
- **Added a boxplot in the generated pdf summary
- **Changed the well position to the well plates (with the help of importing raw/melting raw data)

- **Changed files:
- qPCR_main_gui.py
- qPCR_pdf_layout.py
- qPCR_pdf_summary.py
- qPCR_visualizer.py
- qPCR_delta_ct_calculations


### [1.3.0] - 2026-07-30

- **Changed the layout of the Graphical User Interface, so that it doesn't show the buttons for generating the (advanced) plots and calculations in the GUI. Now, everything is saved on the Desktop in the file qPCR_Ergebnisse. Furthermore, an executive summary of the most important values is now generated in a pdf file. Furthermore, the (advanced) plots are removed. It will be later again added. The GUI has now only one button for generating all files (pdf, png, etc.) and a status bar. After calculating and saving the file, it gives now on the output folder a csv. file, a viridis colormap of the delta CT values, the pdf summary file, a plasma colormap for the fold change values and the statistical analysis report in a .txt file.

- **Changed files:
- qPCR_main.gui.py
- qPCR_pdf_summary.py (new!)
- qPCR_visualizer.py
- qPCR_delta_ct_calculation.py

### [1.2.0] - 2026-07-18

- **Changed the output of the visual graphics into an output folder on the desktop.
- **Added a textbox for writing the sample name

- **Changed python files:
- qPCR_main_gui.py
- qPCR-visualizer.py

## [1.1.0] - 2026-06-27

### Added
- **Graphics Export**: Optimized export logic to automatically save additional analysis plots to the local filesystem.

### Fixed
- **Colorbar Handling**: Switched from `plt.colorbar` to `fig.colorbar` to resolve `UserWarning` regarding figure references.
- **Matplotlib Deprecation**: Fixed `tick_labels` warning to ensure compatibility with Matplotlib 3.9+.
- **GUI Stability**: Implemented a clean event loop handler (`running` flag and `on_closing` method) to prevent "invalid command name" crashes upon closing the application.

### Changed
- **GUI Layout**: Adjusted plot sizes and spacing within the interface to prevent graphical overlaps.