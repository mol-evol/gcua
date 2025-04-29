# GCUA: General Codon Usage Analysis - User Manual (Version 2.3.0)
![Version](https://img.shields.io/badge/version-2.3.0-blue) ![Python](https://img.shields.io/badge/python-3.8%2B-green) ![License](https://img.shields.io/badge/license-MIT-orange)

**Author:** James McInerney
**Original C Program:** 1997
**Python Implementation:** 2025

**Citation:**
McInerney JO. GCUA: general codon usage analysis. Bioinformatics. 1998;14(4):372-3. doi: 10.1093/bioinformatics/14.4.372. PMID: 9632833.

---

## 1. Introduction & Background

**GCUA (General Codon Usage Analysis)** is a powerful and versatile bioinformatics software tool specifically designed for the in-depth exploration and comprehensive analysis of codon usage patterns within protein-coding DNA sequences. Codon usage bias, the phenomenon where synonymous codons (different codons that specify the same amino acid) are used at unequal frequencies, is a fundamental feature observed across all domains of life, from bacteria to eukaryotes, to viruses, to plasmids. Understanding these patterns is crucial as they can reflect various biological pressures and mechanisms.

Analyzing codon usage can provide valuable insights into:

* **Gene Expression Levels:** Highly expressed genes often exhibit stronger codon bias towards codons recognized by abundant tRNAs, enhancing translation efficiency.
* **Translation Efficiency & Accuracy:** Codon choice can influence the speed and fidelity of protein synthesis.
* **mRNA Stability:** Codon composition can affect the half-life of messenger RNA molecules.
* **Protein Folding:** The rate of translation, modulated by codon usage, can impact co-translational protein folding.
* **Evolutionary Processes:** Differences in codon usage can reflect mutational biases (like GC pressure) or selective pressures, providing clues about genome evolution and horizontal gene transfer.
* **Heterologous Gene Expression:** Understanding native codon preferences is vital when optimizing gene sequences for expression in different host organisms.

GCUA processes coding sequences provided in the standard FASTA format. It calculates a diverse suite of metrics, ranging from basic frequencies (codon counts, amino acid usage) to sophisticated measures of bias (Relative Synonymous Codon Usage - RSCU, Effective Number of Codons - ENC, Codon Adaptation Index - CAI, Frequency of Optimal Codons - Fop, Synonymous Codon Usage Order - SCUO). Furthermore, it provides detailed base composition statistics, including overall GC content and positional GC content (GC1, GC2, GC3, GC3s).

This program represents a significant modernization of the original C version. Rewritten entirely in Python, it leverages the power of modern scientific libraries like NumPy, Pandas, SciPy, and Biopython. This transition facilitates enhanced functionality, broader and more easily maintainable support for the diverse range of known genetic codes (essential for analyzing mitochondrial genomes or organisms with non-standard nuclear codes), and the integration of rich, interactive data visualization capabilities powered by the Plotly library, allowing for dynamic exploration of results.

---

## 2. Installation & Dependencies

GCUA is distributed as a Python script (`gcua.py`) and requires a functioning Python 3 environment. Its extensive capabilities rely on several well-established external Python libraries that must be installed.

**Required Dependencies:**

* **Core Libraries:**
    * `numpy`: Fundamental package for numerical computing in Python. Used extensively for array operations, mathematical functions, and underlying calculations in other libraries.
    * `pandas`: Provides high-performance, easy-to-use data structures (primarily the DataFrame) and data analysis tools. Essential for organizing, manipulating, and exporting the various metrics calculated by GCUA.
    * `scipy`: A library for scientific and technical computing. GCUA utilizes its statistical functions (`scipy.stats` for tests like chi-squared contingency).
* **Bioinformatics:**
    * `biopython`: The cornerstone bioinformatics library for Python. GCUA uses it heavily for:
        * Parsing FASTA files (`Bio.SeqIO`).
        * Representing sequences (`Bio.Seq`).
        * Basic sequence utilities like GC content calculation (`Bio.SeqUtils`).
        * Handling genetic code translations.
* **Machine Learning:**
    * `scikit-learn`: A comprehensive machine learning library. GCUA specifically uses its implementation of Principal Component Analysis (`sklearn.decomposition.PCA`) for multivariate analysis.
* **Visualization:**
    * `matplotlib`: The foundational plotting library in Python. While not always called directly by the user's choice, it often serves as a backend for other plotting libraries like Seaborn.
    * `seaborn`: Builds upon Matplotlib to provide a high-level interface for drawing attractive and informative statistical graphics. Used for some internal plotting logic or potentially future text-based plots.
    * `plotly`: A powerful library for creating interactive, web-based visualizations. GCUA uses Plotly to generate HTML plots (scatter plots, heatmaps, histograms) that allow users to hover over data points, zoom, pan, and switch between different views (e.g., different axes in multivariate plots).

**Installation Steps:**

1.  **Verify Python Installation:** Open your terminal or command prompt. Type `python --version` or `python3 --version`. You should see a version number, preferably 3.6 or higher. If Python is not installed, download and install it from [python.org](https://www.python.org/).
2.  **Install Dependencies via pip:** `pip` is the standard Python package installer. Run the following command in your terminal. Using a virtual environment (like `venv` or `conda`) is highly recommended to avoid conflicts with other Python projects.

    ```bash
    # If using pip directly:
    pip install numpy pandas scipy biopython scikit-learn matplotlib seaborn plotly

    # If using pip3:
    pip3 install numpy pandas scipy biopython scikit-learn matplotlib seaborn plotly
    ```
    *Troubleshooting:* If you encounter errors, ensure `pip` is up-to-date (`pip install --upgrade pip`) and that you have the necessary build tools installed on your system if any packages require compilation. Check the specific error messages for clues. Using Anaconda/Miniconda can sometimes simplify installation, especially on Windows.

3.  **Obtain the GCUA Script:** Download or save the `gcua.py` script file to a convenient directory on your computer.

---

## 3. Running GCUA

GCUA is designed to be run interactively from your terminal.

1.  **Navigate to the Script Directory:** Open your terminal or command prompt window. Use the `cd` (change directory) command to navigate to the folder where you saved the `gcua.py` script. For example: `cd /path/to/your/gcua/folder`.
2.  **Execute the Script:** Launch the program by invoking the Python interpreter followed by the script name:

    ```bash
    python gcua.py
    ```
    *(Use `python3 gcua.py` if `python` defaults to Python 2 on your system)*

Upon execution, GCUA will:
* Clear the terminal screen (on compatible systems).
* Display a welcome banner showing the program name, version, author, and citation information.
* Indicate the currently selected **Genetic Code** (defaulting to Standard/Universal Code 1 initially).
* Present the **Main Menu** with numbered options.

The interface is entirely text-based; you interact by typing the number or letter corresponding to your desired menu option and pressing Enter.

---

## 4. Program Interface & Workflow

GCUA employs a hierarchical menu system for navigation. You start at the Main Menu and select options to enter submenus for specific tasks like analysis, visualization, or configuration.

**Recommended General Workflow:**

1.  **Load Data (Option 1):** This is the essential first step. Provide the path to your FASTA file containing the coding sequences you wish to analyze. The program will parse the sequences and perform initial calculations like codon counts and base composition. Progress indicators will be shown for larger files.
2.  **Set Preferences (Option 6 - CRUCIAL):** Immediately after loading data (or before, if you know the code), navigate to the Preferences menu. The most critical setting here is the **Genetic Code**. Select the NCBI translation table ID that accurately reflects the organism(s) your sequences come from (e.g., 2 for Vertebrate Mitochondrial, 11 for Bacterial). Using the wrong code will lead to incorrect translations and invalidate many analyses (RSCU, AA Usage, CAI, Fop, Optimization). **Remember: If you change the genetic code *after* loading data, you *must* reload the FASTA file (Option 1 again) for the change to take effect on the data.** You can also configure output verbosity and visualization style here.
3.  **Perform Analysis (Option 2):** Explore the Analysis menu. You can either run specific analyses individually (e.g., calculate only ENC) or use **Option 2.1 (Calculate all metrics)** to compute everything needed for most downstream tasks and visualizations. This is often convenient.
4.  **Visualize (Option 3):** Once analyses are complete, use the Visualization menu to generate plots. Select the desired plot type. GCUA will generate an HTML file (if using Plotly) and attempt to open it in your default web browser. These plots are interactive, allowing you to explore the data visually.
5.  **Optimize (Option 4 - Optional):** If your goal is to redesign sequences (e.g., for expression in a different host), use this menu. You'll need to have identified or loaded optimal codons first (via Analysis Menu Option 9).
6.  **Export (Option 5):** Save your numerical results, analysis outputs, or optimized sequences. Choose the relevant export option. Files are saved to a dedicated `gcua_outputs` subdirectory created in the same location where you ran `gcua.py`. This keeps your results organized.
7.  **Help (Option 7):** Refer to this for a quick reminder of the program's functions and the citation.
8.  **Quit (Option Q):** Exit the program gracefully.

Navigation typically involves entering the number/letter for your choice. Use 'R' in submenus to return to the previous menu level.

---

## 5. Main Menu Options

The top-level menu provides access to all major functional areas of GCUA.

* **1. Load FASTA file:** Initiates the data loading process. Essential starting point. Handles parsing and initial counting.
* **2. Analysis:** Leads to the core analytical engine, offering a wide range of codon usage metrics and statistical analyses.
* **3. Visualization:** Accesses tools to create graphical representations of your analysis results, aiding interpretation.
* **4. Sequence Optimization:** Provides functions to redesign input sequences based on calculated or loaded optimal codon tables.
* **5. Export Data:** Allows saving calculated data, analysis results, and optimized sequences to various file formats for external use or record-keeping.
* **6. Preferences:** Enables configuration of crucial program settings, most importantly the genetic code.
* **7. Help:** Displays a concise help text outlining the program's purpose and features, along with the citation details.
* **Q. Quit program:** Terminates the GCUA application.

---

## 6. Analysis Menu (Option 2)

This is the central hub for performing calculations and statistical analyses on your loaded sequence data.

* **1. Calculate all metrics:** A convenient option that runs all major analyses sequentially: Base Composition (re-check), Codon Usage & RSCU, Amino Acid Usage, Multivariate Analysis (CA on RSCU by default), ENC, Optimal Codons (frequency-based default), Fop, CAI, and SCUO. Ideal for ensuring all data is ready for subsequent visualization or export. Displays summary statistics upon completion.
* **2. Calculate codon usage:** Focuses on codon frequencies.
    * *Display codon usage for each gene:* Shows detailed tables (paginated for large datasets) with raw counts and calculated RSCU values for each codon within each gene. RSCU normalizes for amino acid frequency, highlighting preferences among synonymous codons (RSCU > 1 indicates preference, < 1 indicates avoidance).
    * *Display cumulative codon usage:* Aggregates counts across all loaded genes and calculates overall RSCU, giving a genome-wide or dataset-wide perspective.
    * *Export codon usage to file:* Saves two TSV files: one with raw codon counts per gene, and another with RSCU values per gene.
* **3. Calculate amino acid usage:** Focuses on the resulting protein composition.
    * *Display amino acid usage for each gene:* Shows raw counts of each amino acid (and stop codons) per gene (paginated).
    * *Display cumulative amino acid usage:* Shows total counts for each amino acid across the entire dataset.
    * *Export amino acid usage to file:* Saves amino acid counts per gene to a TSV file.
* **4. Calculate base composition:** Provides detailed nucleotide statistics.
    * *Display base composition for each gene:* Shows per-gene metrics including sequence length, total codons (AA count), overall GC%, GC% at 1st, 2nd, and 3rd codon positions (GC1, GC2, GC3), GC% at synonymous third positions (GC3s), and raw counts of A, T, C, G overall and specifically at the 3rd position (A3, T3, C3, G3). Paginated display.
    * *Display cumulative base composition:* Calculates and displays the average values for these composition metrics across all loaded genes.
    * *Export base composition to file:* Saves the detailed per-gene base composition data to a TSV file.
* **5. Perform multivariate analysis:** Applies dimensionality reduction techniques to explore major trends in usage patterns across many genes and codons/amino acids simultaneously.
    * *Perform Correspondence Analysis (CA) on RSCU values:* The standard and often preferred method for codon usage. It positions genes and codons in a low-dimensional space based on their association, often revealing axes related to expression level, GC bias, or other factors.
    * *Perform Principal Component Analysis (PCA) on RSCU values:* An alternative that focuses on variance. May be less intuitive for frequency data like RSCU compared to CA. Requires Z-score standardization.
    * *Perform CA on amino acid usage:* Explores trends in amino acid composition variation across genes.
    * *Perform PCA on amino acid usage:* Alternative for exploring AA composition variance.
    * *Visualize multivariate analysis results:* Generates an interactive scatter plot showing genes on the primary axes (e.g., Axis1 vs Axis2 for CA). Requires analysis to be run first.
    * *Export multivariate analysis results:* Saves key outputs: gene coordinates on the principal axes/components, the percentage of variance explained by each axis, and the codon/amino acid loadings (indicating their contribution to each axis). Saved as TSV files.
* **6. Calculate ENC values:** Focuses on the Effective Number of Codons, a measure of how far codon usage deviates from equal usage of synonymous codons.
    * *Calculate and display ENC values:* Computes ENC (Wright's method) for each gene. Values range from 20 (extreme bias, one codon per AA) to 61 (no bias, all codons used equally).
    * *Create ENC vs GC3s plot (Wright's plot):* A standard diagnostic plot comparing observed ENC to the value expected based solely on GC content at the third synonymous position (GC3s). Genes falling significantly below the expected curve suggest selection is influencing codon usage beyond mutational bias.
    * *Export ENC values to file:* Saves per-gene ENC values to a TSV file.
* **7. Calculate Fop and CAI values:** Calculates two related measures of codon usage adaptation, often used as proxies for gene expression level or translation efficiency. Both require a set of "optimal" codons to be defined.
    * *Use all genes as reference:* Defines optimal codons based on the most frequent codons across the entire dataset. Calculates Fop (proportion of optimal codons used) and CAI (geometric mean of relative adaptiveness values, weighted by codon frequency in the gene) based on this global reference.
    * *Select reference genes manually:* Allows the user to specify a subset of genes (e.g., known highly expressed genes like ribosomal proteins) by number, name, or by loading a list from a file. Optimal codons are derived *only* from this subset, and Fop/CAI are calculated relative to this specific reference set.
    * *Use multivariate analysis to identify reference genes:* Automatically identifies a potential set of highly expressed genes by taking genes lying at one extreme of Axis 1 from a CA of RSCU values (typically the end with lower average ENC). Prompts for the percentage of genes to include in this reference set (e.g., top 10%). Calculates Fop/CAI relative to this automatically derived set.
    * *Display Fop and CAI values (using current reference):* Shows the Fop and CAI values calculated using the most recently defined set of optimal codons/reference set.
    * *Create CAI distribution plot:* Generates a histogram showing the distribution of CAI values across all genes, useful for assessing the overall adaptation landscape.
    * *Export Fop and CAI values:* Saves the currently calculated Fop and CAI values for each gene to a TSV file. The filename might reflect the reference method used if applicable.
* **8. Calculate SCUO values:** Focuses on Synonymous Codon Usage Order, an information theory-based measure of bias within synonymous codon groups.
    * *Calculate and display SCUO values:* Computes and shows the SCUO value for each gene. Higher values indicate stronger bias within synonymous groups.
    * *Export SCUO values to file:* Saves per-gene SCUO values to a TSV file.
* **9. Identify optimal codons:** Determines the single "best" or "preferred" codon for each amino acid (excluding Met, Trp, and Stop) based on various criteria. This optimal set is then used for Fop, CAI, and sequence optimization.
    * *Identify optimal codons using all genes (frequency-based):* Simplest method; picks the most frequent synonymous codon across the whole dataset.
    * *Identify optimal codons using manually selected reference genes:* Defines optimal codons based on frequency *within* a user-specified gene set.
    * *Identify optimal codons using multivariate analysis:* Defines optimal codons based on frequency *within* the gene set identified at the CA Axis 1 extreme.
    * *Identify optimal codons using highest RSCU values:* Selects the codon with the highest average RSCU value across the reference set (defaults to all genes if none specified).
    * *Identify optimal codons using most common codons:* Similar to frequency-based but uses raw counts instead of relative frequency.
    * *Load reference genes from file:* Utility to load a gene list from a text file (one gene ID per line) to be used as the reference set for methods requiring one.
    * *Load optimal codons from file:* Allows bypassing calculation and directly loading an optimal codon table previously saved by GCUA (TSV or JSON) or created externally.
    * *Compare codon usage between axis cohorts:* (See Option 0 below - this likely triggers the same function).
    * *Export optimal codons to file:* Saves the *currently defined* set of optimal codons (AA -> RNA Codon -> DNA Codon) to a user-specified format (TSV or JSON).
* **0. Compare codon usage between axis cohorts:** A specialized analysis. It performs CA on RSCU, identifies gene cohorts at opposite extremes of Axis 1 (user defines percentage), calculates aggregate codon usage for each cohort, and performs Chi-squared tests for each amino acid to find codons used significantly differently between the two cohorts. It reports statistics, identifies the likely "highly expressed" cohort (based on lower ENC), and lists preferred codons in each. Useful for identifying codons potentially under translational selection. Results can be exported.
* **R. Return to main menu:** Navigates back to the main program menu.

---

## 7. Visualization Menu (Option 3)

This menu allows graphical exploration of the results obtained from the Analysis menu. Requires the relevant analyses to have been run. Uses Plotly for interactive HTML outputs.

* **1. Create multivariate analysis plot:** Generates a scatter plot of genes based on their coordinates from CA or PCA (usually Axis 1 vs Axis 2).
    * *Interpretation:* Clusters of genes may indicate shared usage patterns. The spread along axes can reveal major sources of variation (e.g., GC bias, expression level). Hovering over points shows gene IDs. Buttons allow switching to view other axis combinations (e.g., Axis 1 vs Axis 3). Outliers (genes far from the main cluster) are highlighted in red.
    * *Requires:* Multivariate analysis (Analysis Menu Option 5) must be run first.
* **2. Create GC content vs GC3 plot:** Plots GC content at the 3rd codon position (GC3%) against the overall GC content (GC%) for each gene.
    * *Interpretation:* Helps assess the influence of mutational bias (GC pressure). Points falling near the diagonal (GC=GC3) suggest uniform GC content. Deviation can indicate non-uniform pressure or selection. A regression line shows the overall trend.
    * *Requires:* Base composition data (calculated during loading).
* **3. Create ENC vs GC3s plot (Wright's plot):** Plots the Effective Number of Codons (ENC) against GC content at synonymous third positions (GC3s%). Includes the theoretical curve representing expected ENC under mutation-drift equilibrium.
    * *Interpretation:* Genes falling significantly *below* the theoretical curve are likely experiencing selective pressure favoring specific codons, leading to higher bias (lower ENC) than expected from GC content alone. Genes on or near the curve suggest codon usage is primarily driven by mutational bias.
    * *Requires:* ENC values (Analysis Menu Option 6) and GC3s (Base Composition).
* **4. Create RSCU heatmap:** Displays RSCU values in a grid format where rows are codons, columns are genes, and cell color intensity represents the RSCU value.
    * *Interpretation:* Provides a visual overview of codon preferences across the entire dataset. Patterns in rows indicate consistent preference/avoidance of specific codons. Patterns in columns highlight genes with similar usage biases. Stop codons are typically excluded.
    * *Requires:* RSCU values (calculated during loading or via Analysis Menu Option 2).
* **5. Create CAI distribution plot:** Generates a histogram showing the frequency distribution of Codon Adaptation Index (CAI) values across all genes.
    * *Interpretation:* Shows the overall landscape of gene adaptation to the reference codon usage bias. A distribution skewed towards high CAI values might indicate strong translational selection in the organism. A rug plot at the bottom shows individual gene CAI values.
    * *Requires:* CAI values (Analysis Menu Option 7).
* **6. Show codon bias overview:** Calculates and displays correlations between key metrics (ENC, CAI, Fop, SCUO, GC, GC3s) in a table format printed to the terminal. Also shows a summary table for the first few genes. Can be exported to a TSV file containing the full table and the correlation matrix.
    * *Interpretation:* Helps understand the relationships between different measures of codon usage and composition (e.g., is CAI negatively correlated with ENC? Is GC3s correlated with CAI?).
    * *Requires:* ENC, CAI, Fop, SCUO, GC, GC3s values.
* **7. Create custom scatter plot:** A flexible plotting tool. Allows the user to select any two calculated metrics (from Base Comp, ENC, CAI, Fop, SCUO, Multivariate Axes) and plot them against each other for all genes.
    * *Features:* Includes gene ID on hover, calculates and plots a linear regression line, displays the Pearson correlation coefficient (r). Offers an optional outlier analysis based on Z-score distance from the centroid, highlighting and potentially exporting outlier genes.
    * *Requires:* The two selected metrics must have been calculated.
* **R. Return to main menu:** Goes back to the main program menu.

---

## 8. Sequence Optimization Menu (Option 4)

Provides tools to redesign input DNA sequences, typically to improve expression potential in a specific host system by replacing existing codons with 'optimal' ones. This requires a set of optimal codons to be defined first (either calculated within GCUA or loaded from a file).

* **1. Optimize a single gene using all genes as reference:** Select a gene by number/name. GCUA calculates optimal codons based on the frequency across *all* loaded genes and then rewrites the selected gene's sequence using these codons. Displays original vs. optimized sequence snippets and the percentage of change. Optionally saves the optimized sequence to a FASTA file.
* **2. Optimize a single gene using multivariate analysis:** Select a gene. GCUA identifies optimal codons based on frequencies within the gene set at the CA Axis 1 extreme (user specifies percentage) and optimizes the selected gene using this set. Displays comparison and allows saving.
* **3. Optimize a single gene using manually selected reference genes:** Select a gene. GCUA prompts the user to define a reference gene set (by number, name, or file). Optimal codons are derived from this *specific* set, and the target gene is optimized accordingly. Displays comparison and allows saving.
* **4. Optimize all genes:** Automatically optimizes *all* sequences loaded in the current session using the *currently defined* set of optimal codons (ensure the desired optimal set is active). Saves all optimized sequences together in a single output FASTA file. The filename will indicate it contains all optimized sequences.
* **5. Load optimal codons from file for optimization:** First prompts for a TSV or JSON file containing an optimal codon table. After successfully loading the table, it prompts whether to optimize a single selected gene or all genes using this externally provided table.
* **R. Return to main menu:** Goes back to the main program menu.

*Important Note:* The effectiveness of sequence optimization depends heavily on choosing an appropriate set of optimal codons, often derived from highly expressed genes in the target expression host. Using optimal codons derived from the source organism itself might not be ideal for heterologous expression.

---

## 9. Export Data Menu (Option 5)

This menu centralizes options for saving the results of your analyses to files for documentation, sharing, or use in other software (e.g., spreadsheets, statistical packages, custom scripts). Files are saved to the `gcua_outputs` directory.

* **1. Export comprehensive metrics:** Creates a single, detailed TSV file. Each row represents a gene, and columns include all major calculated metrics: Length, AA_count, GC%, GC1/2/3%, GC3s%, base counts (A,T,G,C, A3,T3,G3,C3), ENC, Fop, CAI, SCUO, and coordinates on the first few multivariate axes (e.g., Axis1, Axis2, PC1, PC2). Includes metadata headers. Highly useful for downstream analysis.
* **2. Export codon usage data:** Saves two separate TSV files:
    * `*_codon_counts.tsv`: Raw counts for each of the 64 codons for every gene.
    * `*_rscu_values.tsv`: Calculated RSCU values for each codon for every gene.
* **3. Export amino acid usage data:** Saves a TSV file with counts for each of the 20 amino acids and stop signals for every gene.
* **4. Export multivariate analysis results:** Saves the primary outputs of CA or PCA to TSV files:
    * `*_coordinates.tsv`: Coordinates of each gene on the principal axes/components.
    * `*_loadings.tsv`: Contribution (loading) of each codon/amino acid to each principal axis/component.
    * `*_explained_variance.tsv`: The proportion of total variance captured by each axis/component.
* **5. Export optimized sequences:** Functionally identical to "Optimize all genes" (Option 4 in the Optimization Menu). Saves all optimized sequences into one multi-FASTA file. Requires optimal codons to be defined.
* **6. Export optimal codons:** Saves the *currently active* set of optimal codons. Prompts the user to choose between TSV format (AA, RNA_Codon, DNA_Codon columns) or JSON format (structured data with optional metadata). Essential for saving custom codon tables or those derived from specific reference sets.
* **7. Export codon usage comparison results:** If the axis cohort comparison (Analysis Menu Option 0) has been run, this option saves the detailed findings into multiple files: a summary TSV of Chi-squared results per AA, a detailed TSV of codon counts/RSCU for significantly different AAs, TSV files for the aggregate usage in each cohort (left and right), TSV files defining the preferred codons for each cohort, an optimal codon TSV based on the likely highly expressed cohort, and a comprehensive JSON file containing all results and metadata.
* **R. Return to main menu:** Goes back to the main program menu.

---

## 10. Preferences Menu (Option 6)

Allows customization of key program settings that affect calculations and output.

* **1. Genetic Code:** This is the most critical preference. It presents a numbered list of supported genetic codes (Standard, Vertebrate Mitochondrial, Yeast Mitochondrial, Bacterial, etc.) identified by their NCBI translation table ID and name. Select the number corresponding to the code used by your input sequences.
    * **WARNING:** Setting the incorrect genetic code will lead to errors in translation-dependent analyses (RSCU, AA Usage, CAI, Fop, Optimization). Always verify this setting. If you change the code *after* loading data, the program will warn you and prompt you to reload the FASTA file for the changes to apply correctly to the internal data structures.
* **2. Output Detail Level:** Controls how much information is printed to the terminal during calculations.
    * `1: Basic`: Minimal output, mostly just final results.
    * `2: Moderate`: Shows some intermediate calculation steps or progress indicators.
    * `3: Detailed`: Prints extensive information about ongoing calculations (can be very verbose for large datasets).
* **3. Visualization Method:** Determines the format for generated plots.
    * `1: Plotly` (Default): Creates interactive HTML files that open in a web browser. Recommended for data exploration.
    * `2: Text`: Prints rudimentary text-based representations of plots to the terminal (limited utility).
* **4. Toggle Beep:** Turns an audible beep sound (played on completion of some tasks) on or off.
* **R. Return to main menu:** Goes back to the main program menu.

---

## 11. Output Files

GCUA organizes its output files neatly into a dedicated subdirectory.

* **Directory:** A folder named `gcua_outputs` is automatically created in the same directory where you execute `gcua.py`. All exported files are saved here.
* **Naming Convention:** Filenames are generally constructed using:
    * The base name of the input FASTA file (e.g., `my_genes`).
    * A descriptor of the content (e.g., `codon_counts`, `rscu_values`, `enc`, `optimal_codons_multivariate`, `axis_comparison_summary`).
    * A timestamp (e.g., `20250429_124500`) to ensure uniqueness and prevent overwriting previous results.
    * The appropriate file extension (`.tsv`, `.fasta`, `.html`, `.json`, `.txt`).
    * *Example:* `gcua_outputs/my_genes_rscu_values_20250429_124500.tsv`
* **File Formats:**
    * **`.tsv`:** Standard tab-delimited text files. Easily opened in spreadsheet software (Excel, Google Sheets, LibreOffice Calc) or parsed by scripts. Metadata (like genetic code used, date, version) is often included at the top in lines beginning with `#`.
    * **`.fasta`:** Standard sequence format used for exporting optimized sequences. Contains sequence headers starting with `>` followed by the sequence data.
    * **`.html`:** Self-contained HTML files for Plotly visualizations. Can be opened offline in any modern web browser. Allows interactive zooming, panning, hovering for details, and saving static images (PNG).
    * **`.json`:** Structured data format, human-readable but also easily parsed programmatically. Used for exporting optimal codons or complex comparison results, preserving metadata alongside the data.
    * **`.txt`:** Plain text files, potentially used for simple lists like reference gene IDs.

* **Metadata:** Pay attention to the commented lines (`#`) at the beginning of TSV files or within JSON structures. They provide crucial context about how the data was generated (e.g., which genetic code was used, which reference set was employed for CAI/Fop, the date of analysis).

---

## 12. Help & Citation (Option 7)

This option provides a quick reference within the program.

* **Content:** Displays a concise summary of GCUA's purpose, the typical workflow, a list of major analysis types and visualization options.
* **Citation:** Reminds users of the correct way to cite the original GCUA publication if the software is used in research work. Proper citation is essential for acknowledging the tool's development.
* **Current Settings:** Also displays the currently active genetic code for user awareness.



