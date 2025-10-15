## EEG Analysis: Phase-Locking Value and Spectral Methods

This repository contains MATLAB code for analyzing EEG entrainment under different light stimulation conditions. The analyses focus on phase-locking value (PLV) between channels, frequency-domain measures (FFT/STFT), and regional network metrics across multiple days and groups (control vs active).

### Highlights
- Computes time-resolved PLV and aggregates by condition, day, and subject
- Produces group comparisons with non-parametric statistics and multiple-figure exports
- Summarizes regional effects (within- and between-region pairs) including effect sizes and p-values

## Repository Structure (selected)
- `main.m`: high-level entry point for batch analyses
- `PLV_plots.m`: end-to-end PLV computation, statistics, and figure/CSV export
- `calculate_timeresolved_plv.m`, `visualize_plv_data.m`, `compare_plv_groups.m`, `compare_network_metrics.m`: core PLV utilities
- `FFT_plots.m`, `STFT_plots.m`, `plotFFT.m`, `plot_topography.m`: spectral and visualization utilities
- `swtest.m`, `fdr_correction.m`: statistics helpers
- `sample_locs/`: electrode location templates for topography plotting

## Requirements
- MATLAB R2020a or newer
- Toolboxes: Statistics and Machine Learning Toolbox, Signal Processing Toolbox

## Data
Analyses expect preprocessed EEG arranged by subject, day, and stimulation type. Internally, scripts access data as:

`data_struct.(subject).(day).(stim).data_clean`


## Outputs
- Figures: PLV matrices, group comparisons, regional time courses, effect-size heatmaps (PDF/PNG)
- Tables: descriptive statistics and inferential results (e.g., `descriptive_analysis.csv`, `statistical_analysis_results.csv`, `plv_regional_effect_sizes.csv`)


## Citation
If you use this code, please cite the accompanying paper and dataset:


## DOIs
- Paper: EEG Entrainment under Dynamic Light Stimulation Reveals Phase-Locking Patterns Across Brain Regions
PLOS ONE, 2025, https://doi.org/10.1371/journal.pone.0332310
- Dataset: Zenodo, https://doi.org/10.5281/zenodo.17019688

