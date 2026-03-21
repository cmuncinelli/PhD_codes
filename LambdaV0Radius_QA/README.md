# LambdaV0Radius QA

## **Work in Progress**
> This study is ongoing. Results and conclusions presented here are preliminary and subject to change.

Quality assurance study investigating systematic effects in Λ reconstruction as a function of the
V0 decay radius. This is a cross-check designed to validate the stability of the analysis
selections used in the main ring polarization measurement.

## Pipeline Overview

The workflow follows these steps:

1. **Download** raw AOD or derived data from the ALICE grid (Alien) using the provided shell scripts.
2. **Extract analysis results** from the O2 output using `extract_analysis_results.cxx`.
3. **Calculate efficiency** corrections with `CalculateEfficiency.py` (or the no-signal-extraction
   variant for cross-checks).
4. **Apply efficiency** corrections with `ApplyEfficiency.py`.
5. **Extract signal** from the corrected spectra using `SignalExtraction.py`.

## Contents

| File / Script | Description |
|---|---|
| `extract_analysis_results.cxx` | C++ macro to extract and format analysis results from O2 output files. |
| `CalculateEfficiency.py` | Efficiency calculation including signal extraction. |
| `CalculateEfficiency-NoSignalExtraction.py` | Efficiency calculation without signal extraction (for cross-checks). |
| `ApplyEfficiency.py` | Applies computed efficiency corrections to the data spectra. |
| `ApplyEfficiency-NoSignalExtraction.py` | Efficiency application without signal extraction. |
| `SignalExtraction.py` | Signal extraction from invariant mass spectra. |
| `download_files.sh` | Downloads files from Alien. |
| `download_files_reworkedForOO.sh` | Download script adapted for O--O data. |
| `download_files_with_timeout.sh` | Download with timeout protection for large files. |
| `expand_alien_folders.sh` | Expands Alien folder listings into file lists. |
| `gen_input_paths.sh` / `gen_paths.sh` | Generates input path lists for the analysis. |
| `merge_download_lists.sh` | Merges multiple download lists into one. |
| `run_analysis_multithreaded.sh` | Runs the analysis with multi-threading support. |
| `run_test.sh` | Quick test run script. |
| `dpl-config-Asymmetry_test_finerYbinning.json` | DPL workflow configuration for this QA. |
| `DownloadListFromAlien.txt` | Reference list of files downloaded from Alien. |
| `Some_cuts_and_runMeFiles.txt` | Notes on selection cuts and run commands. |
| `run_extract_analysis_results.txt` | Command reference for running the extraction macro. |