# MATLAB Codes for "Mortgages, House Prices, and Business Cycle Dynamics: A Medium-Run Exploration Using the Continuous Wavelet Transform"

This repository contains the MATLAB codes required to replicate the results presented in the paper:

**Title:** Mortgages, House Prices, and Business Cycle Dynamics: A Medium-Run Exploration Using the Continuous Wavelet Transform  
**Authors:** Michele Andreani, Federico Giri  
**Journal:** International Review of Economics & Finance, 2024, Volume 94

## Requirements

- **MATLAB Version:** These codes have been written and tested using MATLAB 9.13.0.2126072 (R2022b) Update 3. Please ensure you are using this version or a compatible version to avoid compatibility issues.

- **Toolboxes:** The codes require the following toolboxes:
  1. **ASToolbox** by M. Joana Soares and L. Aguiar-Conraria  
     [ASToolbox](https://sites.google.com/site/aguiarconraria/wavelets-and-economics/the-astoolbox)
  2. **Cross Wavelet and Wavelet Coherence Toolbox** by Aslak Grinsted  
     [Wavelet Coherence Toolbox](https://grinsted.github.io/wavelet-coherence/)
  3. **Multiple Wavelet Coherence** by Hu, W., and B.C. Si (2016)  
     [Multiple Wavelet Coherence Code](https://figshare.com/articles/code/Matlab_code_for_multiple_wavelet_coherence_and_partial_wavelet_coherency/13031123)

## Usage

To ensure smooth operation, use relative paths when running the scripts.

### How to Use

1. **Clone the Repository:**
   - Clone this repository to your local machine
2. **Setup in MATLAB:**
   - Open MATLAB and navigate to the directory where the repository is located.
3. **Run the Scripts:**
   - **`data.m`**: Downloads data from FRED and creates the dataset.
   - Execute the remaining `.m` files to generate and save the main figures 1-5 of the paper.
   - Alternatively, run **`MAIN.m`** to execute all scripts in the folder sequentially.
4. **Output:**
   - The codes will generate two folders:
     1. **`data_input`**: Contains data for Wavelet coherence estimation.
     2. **`figures`**: Contains the main figures of the paper, exported in **`.eps`** format.

## Citation

Please cite the paper as follows:

Andreani, M., & Giri, F. (2024). Mortgages, house prices, and business cycle dynamics: A medium-run exploration using the continuous wavelet transform. *International Review of Economics & Finance*, 103380.

