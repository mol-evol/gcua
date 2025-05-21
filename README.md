# GCUA: General Codon Usage Analysis

![Version](https://img.shields.io/badge/version-2.2.0-blue)
![Python](https://img.shields.io/badge/python-3.8%2B-green)
![License](https://img.shields.io/badge/license-MIT-orange)

## Overview

GCUA (General Codon Usage Analysis) is a comprehensive Python tool for analyzing codon usage patterns in DNA sequences. This improved implementation builds upon James McInerney's original C program from 1997, adding modern analytical capabilities, interactive visualizations, and sequence optimization features.

The program allows researchers to:
- Analyze codon usage bias in DNA sequences
- Perform multivariate analysis (Correspondence Analysis and Principal Component Analysis)
- Calculate various codon usage metrics (ENC, CAI, Fop, SCUO)
- Visualize results with interactive plots
- Optimize gene sequences based on identified codon usage patterns

## Citation

If you use GCUA in your research, please cite:

```
McInerney JO. GCUA: general codon usage analysis.
Bioinformatics. 1998;14(4):372-3.
doi: 10.1093/bioinformatics/14.4.372. PMID: 9632833.
```

## Table of Contents

- [Installation](#installation)
  - [Requirements](#requirements)
  - [Installation Steps](#installation-steps)
- [Usage](#usage)
  - [Command-Line Interface](#command-line-interface)
  - [Input Files](#input-files)
  - [Example Workflow](#example-workflow)
- [Features](#features)
  - [Codon Usage Analysis](#codon-usage-analysis)
  - [Multivariate Analysis](#multivariate-analysis)
  - [Visualization](#visualization)
  - [Sequence Optimization](#sequence-optimization)
- [Output Files](#output-files)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Requirements

GCUA requires Python 3.8 or higher and several third-party packages:

- numpy
- pandas
- biopython
- scipy
- matplotlib
- seaborn
- scikit-learn
- plotly
- ete3

### Installation Steps

#### Option 1: Using pip (recommended)

1. Clone the repository:
   ```bash
   git clone https://github.com/username/gcua.git
   cd gcua
   ```

2. Install required packages using the provided requirements.txt:
   ```bash
   pip install -r requirements.txt
   ```

3. Verify installation:
   ```bash
   python gcua.py --version
   ```

#### Option 2: Manual installation

1. Clone the repository:
   ```bash
   git clone https://github.com/username/gcua.git
   cd gcua
   ```

2. Install the required packages manually:
   ```bash
   pip install numpy pandas biopython scipy matplotlib seaborn scikit-learn plotly
   ```

#### Option 3: Using Conda

1. Clone the repository:
   ```bash
   git clone https://github.com/username/gcua.git
   cd gcua
   ```

2. Create a conda environment and install dependencies:
   ```bash
   conda create -n gcua python=3.9
   conda activate gcua
   conda install numpy pandas biopython scipy matplotlib seaborn scikit-learn plotly
   ```

## Usage

### Command-Line Interface

GCUA features an interactive command-line interface with menu-driven navigation:

```bash
python gcua.py
```

The main menu provides access to various functions:
1. Load FASTA file
2. Analysis
3. Visualization
4. Sequence Optimization
5. Export Data
6. Preferences
7. Help
Q. Quit program

### Input Files

GCUA accepts DNA sequences in FASTA format. The sequences should be protein-coding DNA sequences (CDS) with lengths divisible by 3.

Example FASTA format:
```
>gene1
ATGGCGTACTTCGATATCGATCGATCGTAGCTAGCTGATCGATCGAT
>gene2
ATGACTGACTAGCTAGCTACGATCGATCGATCGTACGTAGCTAGCAT
```

### Example Workflow

A typical workflow might include:

1. Load a FASTA file with coding sequences
2. Calculate codon usage and base composition statistics
3. Perform multivariate analysis to identify patterns
4. Calculate codon bias metrics like ENC, CAI, and Fop
5. Visualize results through interactive plots
6. Export data or optimize sequences based on findings

## Features

### Codon Usage Analysis

- **Codon Counting**: Calculates codon frequencies for each gene
- **RSCU Calculation**: Computes Relative Synonymous Codon Usage
- **Amino Acid Usage**: Analyzes amino acid usage patterns
- **Base Composition**: Provides GC content analysis (overall, by position, GC3s)

### Multivariate Analysis

- **Correspondence Analysis (CA)**: Identifies patterns in RSCU or amino acid usage
- **Principal Component Analysis (PCA)**: Alternative dimensionality reduction technique
- **Reference Gene Identification**: Automatically identifies potential highly expressed genes

### Visualization

GCUA offers interactive visualizations using Plotly:

- **Multivariate Analysis Plots**: Interactive scatter plots of CA or PCA results
- **GC Content Plots**: GC vs GC3 content visualization
- **ENC vs GC3s Plot (Wright's Plot)**: Visualizes relationship between ENC and GC3s
- **RSCU Heatmaps**: Color-coded visualization of RSCU values
- **CAI Distribution Plots**: Histograms of CAI values
- **Custom Scatter Plots**: User-defined plots of any calculated metrics

### Sequence Optimization

- **Optimal Codon Identification**: Multiple methods for identifying optimal codons:
  - Frequency-based
  - Multivariate analysis-based
  - RSCU-based
  - Reference gene-based
- **Gene Optimization**: Replace codons with optimal alternatives
- **Comparative Analysis**: Analyze differences between original and optimized sequences

### Genetic Code Support

- Universal genetic code
- Mycoplasma/Spiroplasma genetic code (where UGA codes for Tryptophan)

## Output Files

GCUA saves all analysis results and visualizations in a `gcua_outputs` directory. Output formats include:

- **TSV Files**: Tab-separated values for data analysis
- **HTML Files**: Interactive Plotly visualizations
- **FASTA Files**: Optimized sequences
- **JSON Files**: Complex data structures like optimal codon definitions

## Advanced Usage

### Working with Reference Genes

You can identify optimal codons based on a set of reference genes, which can be:
- Manually selected
- Automatically identified through multivariate analysis
- Loaded from a text file (one gene name per line)

Example reference gene file:
```
gene1
gene3
gene7
```

### Comparing Codon Usage Between Cohorts

GCUA can compare codon usage between genes at opposite ends of the primary multivariate axis, which often separates genes by expression level:

```
Main Menu > Analysis > Compare codon usage between axis cohorts
```

This analysis helps identify codons that are statistically overrepresented in potentially highly expressed genes.

### Customizing Optimal Codon Selection

GCUA provides multiple methods for identifying optimal codons:

- **Frequency-based**: Uses the most frequent codon for each amino acid
- **Multivariate**: Uses cohorts identified by multivariate analysis
- **RSCU-based**: Uses codons with highest RSCU values
- **Raw count**: Uses the most common codons

## Troubleshooting

### Common Issues

1. **Import errors**: Make sure all dependencies are installed:
   ```bash
   pip install -r requirements.txt
   ```

2. **Memory issues with large datasets**: For very large FASTA files, ensure your system has sufficient RAM.

3. **Visualization not displaying**: If HTML visualizations don't open automatically, try:
   ```bash
   import webbrowser
   webbrowser.open('file:///path/to/visualization.html')
   ```

4. **Invalid sequence length errors**: Ensure all sequences are coding sequences with lengths divisible by 3.

## Contributing

Contributions to GCUA are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## Metrics Explained

### ENC (Effective Number of Codons)
ENC measures the degree of codon usage bias, ranging from 20 (extreme bias, only one codon per amino acid) to 61 (no bias, equal use of all synonymous codons). Lower values indicate stronger codon bias.

### CAI (Codon Adaptation Index)
CAI measures how well a gene is adapted to the codon usage of highly expressed genes. Values range from 0 to 1, with higher values indicating stronger adaptation to the reference set.

### Fop (Frequency of Optimal Codons)
Fop is the ratio of optimal codons to synonymous codons in a gene. Values range from 0 to 1, with higher values indicating more frequent use of optimal codons.

### SCUO (Synonymous Codon Usage Order)
SCUO quantifies the degree of order in synonymous codon usage based on information theory. Values range from 0 to 1, with higher values indicating more ordered (non-random) codon usage.

### RSCU (Relative Synonymous Codon Usage)
RSCU is the observed frequency of a codon divided by the expected frequency if all synonymous codons for an amino acid were used equally. Values above 1 indicate codons used more frequently than expected by chance.
