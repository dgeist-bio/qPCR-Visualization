# Changelog

All notable changes to this project will be documented in this file.

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