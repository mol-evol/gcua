# GCUA: General Codon Usage Analysis

![Version](https://img.shields.io/badge/version-2.5.0-blue)
![Python](https://img.shields.io/badge/python-3.8%2B-green)
![License](https://img.shields.io/badge/license-MIT-orange)

## Overview

GCUA (General Codon Usage Analysis) is a comprehensive Python tool for analyzing codon usage patterns in DNA sequences. This improved implementation builds upon James McInerney's original C program from 1997, adding modern analytical capabilities, interactive visualizations, sequence optimization features, and scalable processing for large datasets.

The program allows researchers to:
- Analyze codon usage bias in DNA sequences from small to genome-scale datasets
- Perform multivariate analysis (Correspondence Analysis and Principal Component Analysis)
- Calculate various codon usage metrics (ENC, CAI, Fop, SCUO)
- Visualize results with interactive plots
- Optimize gene sequences based on identified codon usage patterns
- Process millions of genes with intelligent memory management and parallel processing

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

The main menu now features a workflow-based structure:
1. Data Management - Load and manage sequence data
2. Quick Analysis (Guided Workflow) - Automated standard analysis
3. Custom Analysis - Detailed analysis options by category
4. Visualization & Export - Unified menu for outputs
5. Advanced Tools - Clustering, optimization, outlier detection
6. Settings & Preferences - Configure analysis parameters
7. Help & Documentation - Comprehensive help system
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

**Option 1: Quick Analysis (Recommended for new users)**
1. Data Management → Load FASTA file
2. Quick Analysis → Follow guided workflow
   - Automatically calculates all metrics
   - Performs multivariate analysis
   - Creates standard visualizations
   - Exports results

**Option 2: Custom Analysis**
1. Data Management → Load FASTA file
2. Custom Analysis → Choose specific analyses:
   - Basic metrics (codon usage, amino acids, base composition)
   - Advanced analysis (multivariate CA/PCA)
   - Codon bias metrics (ENC, CAI, Fop, SCUO)
3. Visualization & Export → Create plots and export data
4. Advanced Tools → Optimize sequences or detect clusters

## Features

### Codon Usage Analysis

- **Codon Counting**: Calculates codon frequencies for each gene
- **RSCU Calculation**: Computes Relative Synonymous Codon Usage
- **Amino Acid Usage**: Analyzes amino acid usage patterns
- **Base Composition**: Provides GC content analysis (overall, by position, GC3s)
- **Extreme Genes Export**: Extract genes from distribution extremes based on any metric

### Multivariate Analysis

- **Correspondence Analysis (CA)**: Identifies patterns in RSCU or amino acid usage
  - Automatically excludes single-codon amino acids (e.g., Met, Trp) from RSCU analysis
- **Principal Component Analysis (PCA)**: Alternative dimensionality reduction technique
  - For amino acid data: Defaults to CLR transformation for proper compositional data analysis
  - Automatically converts to frequencies and offers preprocessing options
- **Reference Gene Identification**: Automatically identifies potential highly expressed genes

### Visualization

GCUA offers interactive visualizations using Plotly:

- **Multivariate Analysis Plots**: Interactive scatter plots of CA or PCA results
  - CA plots now show biplots with both genes and codons (RSCU) or amino acids (AA)
- **GC Content Plots**: GC vs GC3 content visualization
- **ENC vs GC3s Plot (Wright's Plot)**: Visualizes relationship between ENC and GC3s
- **RSCU Heatmaps**: Color-coded visualization of RSCU values
  - For small datasets (<1000 genes): Traditional gene-by-codon heatmap
  - For large datasets (≥1000 genes): Summary statistics heatmap showing mean RSCU, standard deviation, and percentiles
- **CAI Distribution Plots**: Histograms of CAI values
- **Custom Scatter Plots**: User-defined plots of any calculated metrics

**Export Formats**: Visualizations can be exported in multiple formats:
- **HTML** (default): Interactive plots that open in web browser
- **SVG**: Vector graphics for publication-quality figures
- **PNG**: Raster images for presentations
- **PDF**: Document format for reports

Note: Static image formats (SVG, PNG, PDF) require the `kaleido` package (`pip install kaleido`)

### Sequence Optimization

- **Optimal Codon Identification**: Multiple methods for identifying optimal codons:
  - Frequency-based
  - Multivariate analysis-based
  - RSCU-based
  - Reference gene-based
- **Gene Optimization**: Replace codons with optimal alternatives
- **Comparative Analysis**: Analyze differences between original and optimized sequences

### Genetic Code Support

- Supports 33 different genetic codes (NCBI translation tables)
- Including standard code, mitochondrial codes, and various alternative nuclear codes
- Essential for analyzing sequences from diverse organisms

### Performance Features (New in v2.5.0)

- **Hybrid Memory Management**: 
  - Automatically adapts to file size
  - Small files: Keep sequences in memory for speed
  - Large files: Load sequences on demand to conserve memory
  - Configurable threshold (default: 100 MB)

- **Parallel Processing**:
  - Multi-core support for faster analysis
  - Configurable number of processes
  - Batch processing for optimal performance
  - Scales linearly with available CPU cores

- **Progress Tracking and Checkpointing**:
  - Automatic checkpoint saving for long-running analyses
  - Resume interrupted processing from last checkpoint
  - Configurable checkpoint intervals
  - No lost work due to interruptions

- **Scalable to Genome-Scale Data**:
  - Tested with millions of genes
  - Efficient memory usage even for 1M+ gene datasets
  - Fast sequence indexing for quick access

- **Smart Visualizations**:
  - Automatic adaptation for large datasets
  - RSCU heatmap switches to summary statistics for >1000 genes
  - Prevents browser crashes with massive data matrices
  - Shows meaningful patterns instead of overwhelming detail

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

### Exporting Extreme Genes

GCUA can export genes from the extremes of any calculated metric distribution, which is useful for:
- Identifying candidate highly/lowly expressed genes
- Finding genes with unusual codon usage patterns
- Selecting outliers for experimental validation
- Creating gene sets for downstream analysis

To export extreme genes:

```
Main Menu > Analysis > Export extreme genes (E)
```

Features:
- **Flexible metric selection**: Export based on GC, GC3s, ENC, CAI, Fop, SCUO, or multivariate axes
- **Customizable thresholds**: Choose any percentage (e.g., top/bottom 10%, 20%)
- **Multiple export options**: Get top only, bottom only, or both extremes
- **Comprehensive output**: Includes all calculated metrics for exported genes
- **Metadata headers**: Files include selection criteria and threshold values

Example use cases:
- Export top 10% of genes by CAI to identify highly adapted genes
- Export bottom 20% by ENC to find genes with extreme codon bias
- Export genes at multivariate axis extremes to study expression-related patterns

### Customizing Visualization Output

You can set the default visualization format in Preferences:

```
Main Menu > Preferences > Visualization Format
```

This is particularly useful when:
- Working on HPC systems without GUI/browser access (use SVG or PNG)
- Preparing publication-quality figures (use SVG for vector graphics)
- Creating presentations (use PNG for easy insertion)
- Generating reports (use PDF format)

### Customizing Optimal Codon Selection

GCUA provides multiple methods for identifying optimal codons:

- **Frequency-based**: Uses the most frequent codon for each amino acid
- **Multivariate**: Uses cohorts identified by multivariate analysis
- **RSCU-based**: Uses codons with highest RSCU values
- **Raw count**: Uses the most common codons

### Configuring Performance Settings

Access performance settings through the Preferences menu:

```
Main Menu > Preferences > Performance Settings
```

Available options:
- **Parallel Processing**: Enable/disable and set number of CPU cores to use
- **Batch Size**: Adjust the number of sequences processed together (default: 1000)
- **Progress Saving**: Enable automatic checkpointing for long analyses
- **Checkpoint Interval**: How often to save progress (default: every 5000 sequences)

For large datasets (>100,000 genes), recommended settings:
- Enable parallel processing with all available cores
- Set batch size to 5000-10000 sequences
- Enable progress saving with checkpoint interval of 10000 sequences

## Troubleshooting

### Common Issues

1. **Import errors**: Make sure all dependencies are installed:
   ```bash
   pip install -r requirements.txt
   ```

2. **Memory issues with large datasets**: 
   - GCUA automatically switches to low-memory mode for files >100MB
   - Adjust the threshold in Preferences > Memory Management Settings
   - For genome-scale data (>1M genes), ensure at least 8GB RAM available

3. **Visualization not displaying**: If HTML visualizations don't open automatically, try:
   ```bash
   import webbrowser
   webbrowser.open('file:///path/to/visualization.html')
   ```

4. **Invalid sequence length errors**: Ensure all sequences are coding sequences with lengths divisible by 3.

5. **Processing interrupted**: 
   - If analysis is interrupted, GCUA will automatically resume from the last checkpoint
   - Checkpoint files are named `[filename]_checkpoint.pkl`
   - Delete checkpoint files to restart analysis from the beginning

6. **Slow processing on large files**:
   - Enable parallel processing in Preferences > Performance Settings
   - Increase batch size for better throughput
   - Check that all CPU cores are being utilized

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
