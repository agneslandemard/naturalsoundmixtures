# naturalsoundmixtures

## Overview

This repository contains the Matlab code used for the analysis and figure generation presented in:
Landemard A, Bimbard C, Boubenec Y. **Hierarchical encoding of natural sound mixtures in ferret auditory cortex**. _eLife_. 2025. 14:RP106628, doi: [10.7554/eLife.106628](https://doi.org/10.7554/eLife.106628)

The code is designed to process the functional ultrasound imaging (fUSI) dataset publicly available on [Zenodo](https://doi.org/10.5281/zenodo.15800440).

This repository allows users to:
*   Load the preprocessed fUSI data.
*   Perform the main analyses described in the accompanying publication.
*   Reproduce the figures from the publication.

## Author

Agnès Landemard - Ecole Normale Supérieure, Paris

## Citation

If you use this code or data in your work, please cite:

*   Landemard A, Bimbard C, Boubenec Y. **Hierarchical encoding of natural sound mixtures in ferret auditory cortex**. _eLife_. 2025. 14:RP106628, doi: [10.7554/eLife.106628](https://doi.org/10.7554/eLife.106628)
*   Agnès Landemard. (2025). **fUS imaging of ferret auditory cortex during passive listening of natural sound mixtures** [Data set]. _Zenodo_.  [10.5281/zenodo.15800440](https://doi.org/10.5281/zenodo.15800440) 


## Requirements

*   **Matlab:** Version R2022b or later.
*   **Required Toolboxes:** Signal Processing Toolbox, Statistics and Machine Learning Toolbox.
*   **Dataset:** The corresponding dataset must be downloaded separately from [Zenodo](https://doi.org/10.5281/zenodo.15800440).

## Installation / Setup

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/agneslandemard/naturalsoundmixtures.git
    ```
    or download the ZIP file from GitHub.
2.  **Add code directories to Matlab path:**
    Open Matlab, navigate to the root directory of the cloned repository, and run the setup script:
    ```matlab
    startup; 
    ```

## Data Setup

This code relies on the dataset available on [Zenodo](https://doi.org/10.5281/zenodo.15800440). Please download the dataset and organize the files as described in the dataset documentation.

The code needs to know where the dataset files are located on your system.

1.  **Option 1: Configure a path variable:**
    Open the `startup.m` file in the repository root.
    Set the `glob_path` and optionally `metrics_path` variable to the absolute path where you downloaded the dataset files.
    ```matlab
    % startup.m
    glob_path = '/path/to/folder/containing/code/';
    metrics_path = '/path/to/folder/containing/data/';
    ```
2.  **Option 2: Place data in a specific directory:**
    Create a `../data` directory relative to the code repository and place the downloaded data files inside it.

## Code Structure

The repository is organized into the following directories:

*   `main_scripts/`: Top-level scripts to run major analyses and generate figures. Start here!
*   `functions/`: Reusable Matlab functions used across different scripts.
*   `plotting/`: Scripts or functions specifically for generating figures.
*   `data_loading/`: Functions for loading and formatting data files.
*   `figures/`: Directory where scripts save generated figures.

## Getting Started 

To reproduce the figures from the publication "Hierarchical encoding of natural sound mixtures in ferret auditory cortex**, _eLife_, 2025,
follow these steps:

1.  Ensure you have downloaded the dataset and configured the `data_dir` variable in `config.m` as described in the Data Setup section.
2.  Ensure all required Matlab toolboxes are installed and on your path.
3.  The primary scripts to run analyses and plot figures are located in the `main_scripts/` directory. Each script's name indicates which figure(s) of the paper the analyses or plots are related to.
Here is a more detailed description of these scripts:
    *   `main_scripts/Responses_Fig1.m`: Loads timecourse of voxel responses to changes in different sound snippets, and plots the average of such timecourse for each category.
    *   `main_scripts/Cross_correlations_Fig1_FigS1.m`: Computes and plots cross-correlation of voxel responses to two repeats of the same sounds, or the same sounds presented in isolation vs. in mixtures.
    *   `main_scripts/Maps_ferrets_Fig2_Fig3_Fig4.m`: Plots maps for various metrics (e.g. invariance, average response) computed on ferret voxels. This allows to visualize these metrics in the original voxel space.
    *   `main_scripts/Invariances_by_ROI_Fig2_Fig4.m`: Compares values of invariance across different ROIs and subjects, for both species. Also runs statistical tests on regional differences.
    *   `main_scripts/Compare_tuning_regions_Fig3_FigS3.m`: Explores the tuning of voxels across different ROIs obtained through the spectrotemporal model fits, for each species.
    *   `main_scripts/Humans_vs_ferrets_Fig2_Fig4_FigS4.m`: Contrasts prediction accuracy or other metrics for true and predicted data across species.

    ```matlab
    run('main_scripts/Responses_Fig1.m');
    ```
    Figures will be saved to the `figures/` directory.

The proposed code is a starting step to explore the dataset. Further analyses can be carried out by building on the different metrics or formats of the data introduced in each of these scripts.

## Support and Contact

If you encounter any issues or have questions about the code, please:

*   Open an issue on the GitHub repository's issue tracker.
*   Contact Agnès Landemard at agnes.landemard@hotmail.fr.
