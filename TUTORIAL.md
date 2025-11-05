# GCUA Tutorial: Investigating Synonymous Codon Usage

## Introduction

Welcome to this hands-on tutorial on analyzing synonymous codon usage! By the end of this tutorial, you will understand how different organisms show preferences for certain codons over others, even when those codons encode the same amino acid. This phenomenon, called **codon usage bias**, has important implications for gene expression, protein production, and molecular evolution.

## Learning Objectives

By completing this tutorial, you will be able to:
1. Understand what synonymous codons are and why organisms show codon bias
2. Use GCUA to analyze codon usage patterns in real genomic data
3. Calculate and interpret key codon usage metrics (RSCU, ENC, CAI)
4. Visualize codon usage patterns using multivariate analysis
5. Compare codon usage between different genes and predict expression levels
6. Design an independent investigation of codon usage in a biological system

---

## Part 1: Background and Theory (30 minutes)

### What are Synonymous Codons?

The genetic code is **degenerate** - most amino acids are encoded by more than one codon (DNA triplet). For example:
- **Leucine** can be encoded by: UUA, UUG, CUU, CUC, CUA, CUG (6 codons!)
- **Glycine** can be encoded by: GGU, GGC, GGA, GGG (4 codons)
- **Methionine** is encoded by: AUG (only 1 codon - no synonyms)

These alternative codons that encode the same amino acid are called **synonymous codons**.

### Why Does Codon Usage Bias Exist?

If synonymous codons all produce the same protein, why would an organism prefer one over another? Here are the major hypotheses:

1. **Translation Efficiency**: Some codons are decoded faster because their corresponding tRNAs are more abundant
2. **Translation Accuracy**: Certain codons may be decoded more accurately, reducing errors
3. **mRNA Stability**: Codon choice can affect mRNA secondary structure and half-life
4. **Co-translational Folding**: Translation speed affects how proteins fold during synthesis
5. **Mutational Bias**: Simple compositional biases (e.g., GC-rich vs AT-rich genomes)

**Key Insight**: Highly expressed genes (like ribosomal proteins) often show stronger codon bias toward "optimal" codons than lowly expressed genes.

### Key Concepts and Metrics

**1. RSCU (Relative Synonymous Codon Usage)**
- Measures how frequently a codon is used relative to other synonymous codons
- RSCU = 1.0 means no bias (used at expected frequency)
- RSCU > 1.0 means over-represented (preferred)
- RSCU < 1.0 means under-represented (avoided)

**2. ENC (Effective Number of Codons)**
- Summarizes overall codon bias in a gene
- Range: 20 (maximum bias - one codon per amino acid) to 61 (no bias - all codons used equally)
- Lower values = stronger bias

**3. CAI (Codon Adaptation Index)**
- Measures how well a gene's codons match a reference set (typically highly expressed genes)
- Range: 0 (poor adaptation) to 1 (perfect adaptation)
- Useful for predicting expression levels and designing synthetic genes

**4. GC3s (GC content at synonymous third positions)**
- Measures GC content at wobble positions where changes are often synonymous
- Can indicate mutational bias vs. selection

---

## Part 2: Getting Started with GCUA (15 minutes)

### Installation Check

1. Open your terminal and navigate to the GCUA directory:
```bash
cd gcua
```

2. Check that GCUA is installed correctly:
```bash
python gcua.py --version
```

You should see: `GCUA version 2.5.1` (or similar)

3. Launch GCUA:
```bash
python gcua.py
```

You should see the main menu with 7 options.

### Understanding the Interface

GCUA uses an interactive menu system:
- Type a **number** to select an option
- Type **H** for help at any menu
- Type a **number followed by ?** (e.g., `1?`) for specific option help
- Type **B** to go back to the previous menu
- Type **Q** to quit

**PRO TIP**: Always read the menu carefully and use the help system (`H`) if you're unsure!

---

## Part 3: Hands-On Exercises with *Plasmodium falciparum* (90 minutes)

### About the Dataset

The included file `pfal.fas` contains protein-coding gene sequences from *Plasmodium falciparum*, the parasite that causes the most deadly form of malaria. This organism has an extremely AT-rich genome (~80% AT content), making it an interesting case study for codon usage analysis.

### Exercise 1: Loading Data and Basic Exploration

**Task**: Load the sample dataset and explore its basic properties.

1. Launch GCUA and select option **1** (Data Management)
2. Select option **1** (Load FASTA file)
3. When prompted, enter: `pfal.fas`
4. The program will load the sequences - note how many sequences were loaded

**Questions to Answer**:
- Q1.1: How many gene sequences are in the dataset?
- Q1.2: What is the mean sequence length?
- Q1.3: From the main menu, use option 1 → 2 (View data statistics) - what is the range of sequence lengths (min to max)?

**Reflection**: Why might sequence length vary so much between genes? What biological factors determine protein length?

---

### Exercise 2: Analyzing Codon Usage Patterns

**Task**: Calculate codon usage and RSCU values to identify codon preferences.

1. From the main menu, select **3** (Custom Analysis)
2. Select **1** (Basic Metrics)
3. Select **1** (Codon Usage)
4. Choose to calculate:
   - **1** (Codon Counts)
   - **2** (Codon Frequencies)
   - **3** (RSCU values)

The program will calculate these metrics and save them to the `gcua_outputs/` directory.

5. Use option **4** (Visualization & Export) → **1** (Create Visualizations) → **3** (RSCU Heatmap)
6. Open the generated HTML file in your browser

**Questions to Answer**:
- Q2.1: Look at leucine (Leu) codons. Which leucine codon has the highest RSCU? Which has the lowest?
- Q2.2: Examine the overall pattern. Do you see more AT-rich codons (ending in A or T/U) or GC-rich codons (ending in G or C) being preferred?
- Q2.3: Find alanine (Ala) which has 4 synonymous codons: GCU, GCC, GCA, GCG. List their approximate RSCU values. Is there a clear preference?
- Q2.4: Look at the two-codon amino acids (like Asp, Cys, Phe). For these, is there generally a strong preference for one codon over the other?

**Reflection**: Based on your observations, do you think *P. falciparum*'s codon usage is driven more by AT-richness (mutational bias) or by selection for specific "optimal" codons? Explain your reasoning.

---

### Exercise 3: Measuring Codon Bias with ENC

**Task**: Calculate ENC (Effective Number of Codons) for all genes to measure codon bias strength.

1. From the main menu, select **3** (Custom Analysis)
2. Select **2** (Advanced Analysis)
3. Select **1** (Calculate ENC)

The program will calculate ENC for each gene.

4. Create a visualization: **4** (Visualization & Export) → **1** (Create Visualizations) → **4** (ENC Plot)
5. Open the generated HTML file

**Questions to Answer**:
- Q3.1: What is the range of ENC values in this dataset? (Look at the x-axis of your plot)
- Q3.2: Are most genes closer to 20 (maximum bias) or 61 (no bias)?
- Q3.3: The ENC plot shows ENC vs. GC3s (GC content at synonymous sites). Describe the overall pattern. Do genes with different GC3s values show different levels of codon bias?
- Q3.4: Do you see genes that deviate from the expected curve? These are interesting outliers!

**Advanced Challenge**:
- Export the ENC data (option **4** → **2** → appropriate export option)
- Identify the 10 genes with the lowest ENC values (strongest codon bias)
- Identify the 10 genes with the highest ENC values (weakest codon bias)
- Research what these genes encode. Is there a pattern? (Hint: Use NCBI BLAST or gene annotations)

---

### Exercise 4: Multivariate Analysis - Finding Patterns

**Task**: Use Correspondence Analysis to explore patterns in codon usage across genes.

Correspondence Analysis (CA) is a statistical technique that reduces the complexity of codon usage data. It finds the major "axes of variation" and can reveal:
- Clusters of genes with similar codon usage
- The relationship between genes and specific codons
- Potential different "codon usage strategies"

**Steps**:
1. From the main menu, select **3** (Custom Analysis)
2. Select **2** (Advanced Analysis)
3. Select **4** (Correspondence Analysis)
4. When prompted, select:
   - Analysis type: **1** (Codon usage based)
   - Number of axes: **4** (or more)

5. Create visualizations: **4** (Visualization & Export) → **1** (Create Visualizations) → **1** (Multivariate Analysis)
6. View the interactive plot - you can zoom and hover over points

**Questions to Answer**:
- Q4.1: What percentage of the total variation is explained by Axis 1? By Axis 2?
- Q4.2: Do you see any distinct clusters of genes, or is there a continuous distribution?
- Q4.3: Use the plot controls to color points by GC3s content. Does GC3s correlate with any of the major axes?
- Q4.4: Try coloring by ENC. Do genes with strong codon bias (low ENC) cluster together?

**Reflection**: What do you think the major axes represent biologically? Are they capturing compositional bias (GC vs AT), expression level differences, or something else?

---

### Exercise 5: Identifying Optimal Codons and Calculating CAI

**Task**: Identify which codons are "optimal" (used preferentially in highly expressed genes) and calculate CAI.

**Background**: Highly expressed genes (like ribosomal proteins, translation factors) are under strong selection for efficient translation. Their codon usage can serve as a reference for identifying "optimal" codons.

**Steps**:
1. From the main menu, select **5** (Advanced Tools)
2. Select **2** (Identify reference genes)
3. Choose method **1** (Automatic - based on ENC)
4. Use the default top percentage (e.g., 5% or 10%)

The program will identify high-bias genes and determine optimal codons from them.

5. Now calculate CAI: **3** (Custom Analysis) → **2** (Advanced Analysis) → **2** (Calculate CAI)
6. Visualize: **4** (Visualization & Export) → **1** (Create Visualizations) → **5** (CAI Distribution)

**Questions to Answer**:
- Q5.1: How many genes were identified as the reference set (high codon bias genes)?
- Q5.2: Examine the optimal codons (saved in JSON format). For leucine, which codon(s) were identified as optimal?
- Q5.3: Compare the optimal codons to the RSCU patterns you observed in Exercise 2. Are they the same as the most frequently used codons overall?
- Q5.4: Look at the CAI distribution. What is the range of CAI values? Is the distribution skewed or roughly normal?

**Advanced Challenge**:
- Research what genes are in your reference set (the high-bias genes)
- Are they enriched for any particular functional categories (e.g., ribosomal proteins, metabolic enzymes)?
- Why might these genes show strong codon bias?

---

### Exercise 6: Comparative Analysis - Expression Level Predictions

**Task**: Use codon usage metrics to predict which genes might be highly vs. lowly expressed.

**Hypothesis**: If codon bias is driven by selection for translation efficiency, then genes with:
- Low ENC (strong bias)
- High CAI (good adaptation to optimal codons)
- High Fop (high frequency of optimal codons)

...should be more highly expressed.

**Steps**:
1. Calculate Fop: **3** (Custom Analysis) → **2** (Advanced Analysis) → **3** (Calculate Fop)
2. Export all metrics: **4** (Visualization & Export) → **2** (Export Data)
   - Export summary statistics (includes ENC, CAI, Fop, GC3s, etc.)

3. Open the exported TSV file in Excel, Google Sheets, or Python/R

**Analysis Tasks**:
- Sort genes by CAI (highest to lowest)
- Identify the top 20 genes with highest CAI
- Identify the bottom 20 genes with lowest CAI
- For each group, calculate the mean ENC and mean Fop

**Questions to Answer**:
- Q6.1: Is there a correlation between CAI and ENC? (Do high-CAI genes have low ENC?)
- Q6.2: Is there a correlation between CAI and Fop?
- Q6.3: Look up what the top 20 high-CAI genes encode. Are they enriched for any functional categories?
- Q6.4: Do the same for the low-CAI genes. What types of proteins do they encode?

**Reflection**: Do your findings support the hypothesis that codon bias is related to expression level? What evidence would you need to be more confident?

---

### Exercise 7: The ENC-GC3s Plot - Selection vs. Mutation

**Task**: Use the relationship between ENC and GC3s to distinguish selection from mutational bias.

**Background**:
- If codon usage is determined solely by mutational bias (e.g., AT-richness), genes should fall on a predictable curve on an ENC vs. GC3s plot
- Genes that fall **below** this curve show stronger codon bias than expected by mutation alone, suggesting **selection** for specific codons
- This analysis is based on a classic 1990 paper by Wright (1990)

**Steps**:
1. If you haven't already, create an ENC plot: **4** (Visualization & Export) → **1** (Create Visualizations) → **4** (ENC Plot)
2. The plot shows observed data points and the expected curve (representing pure mutational bias)

**Questions to Answer**:
- Q7.1: Do most genes fall on, above, or below the expected curve?
- Q7.2: Identify genes that are far below the expected curve. What does this suggest about the forces shaping their codon usage?
- Q7.3: Are there genes that fall above the expected curve? (This is unusual but can happen)
- Q7.4: Is there a range of GC3s values where selection appears strongest (biggest deviation from the curve)?

**Advanced Challenge**:
- Divide your genes into three groups:
  - Group A: On the expected curve (weak/no selection)
  - Group B: Below the curve (strong selection)
  - Group C: Above the curve (unusual)
- Export your data and categorize genes
- Compare the functions of genes in each group
- Write a paragraph explaining what biological processes might explain the pattern

---

## Part 4: Independent Investigation (Extended Project)

### Your Research Question

Now it's time to design and conduct your own investigation! Choose **one** of the following research questions, or propose your own:

#### Option 1: Codon Usage in Metabolic Pathways
**Question**: Do genes in the same metabolic pathway show similar codon usage patterns?

**Approach**:
- Identify genes involved in glycolysis, TCA cycle, or another pathway in *P. falciparum*
- Filter your dataset to include only these genes (or load a custom FASTA file)
- Compare their codon usage, ENC, and CAI values
- Test whether they cluster together in correspondence analysis

#### Option 2: Codon Usage and Gene Length
**Question**: Do longer genes show different codon usage patterns than shorter genes?

**Approach**:
- Divide your dataset into quantiles by sequence length (e.g., shortest 25%, middle 50%, longest 25%)
- Compare mean ENC, CAI, and Fop between groups
- Examine whether longer genes prefer different codons
- Consider why length might affect codon usage (hint: think about folding and translation time)

#### Option 3: Strand-Specific Codon Bias
**Question**: Do genes on the leading vs. lagging DNA replication strand show different codon usage?

**Approach**:
- You'll need gene annotations with strand information
- Compare codon usage between plus and minus strand genes
- This can reveal replication-associated mutational biases

#### Option 4: Design Your Own Question!

**Some ideas**:
- Compare codon usage in different organisms (you'll need to download FASTA files)
- Investigate whether genes acquired by horizontal gene transfer show unusual codon usage
- Analyze how codon usage changes in different life stages (if data available)
- Compare codon usage in essential vs. non-essential genes

### Investigation Requirements

Your investigation should include:

1. **Clear Research Question**: State a specific, testable hypothesis
2. **Methods**: Describe what analyses you performed using GCUA
3. **Results**: Present your findings with appropriate visualizations
   - Include at least 2 figures
   - Use appropriate statistical summaries
4. **Interpretation**: Explain what your results mean biologically
5. **Limitations**: Identify potential confounding factors or limitations
6. **Conclusion**: Summarize your findings and suggest future directions

### Deliverable Format

Prepare a short research report (3-5 pages) with:
- **Title** and your name
- **Introduction** (background and research question)
- **Methods** (what you did with GCUA)
- **Results** (figures and statistical findings)
- **Discussion** (biological interpretation)
- **References** (cite relevant papers)

---

## Part 5: Advanced Topics and Extensions

### Optional Advanced Exercises

#### A. Sequence Optimization

**Task**: Optimize a gene sequence for expression in *E. coli*

1. Download a *P. falciparum* gene sequence
2. Use GCUA's sequence optimization tools (**5. Advanced Tools** → **1. Optimize Sequences**)
3. Compare the original and optimized sequences:
   - What percentage of codons were changed?
   - Did the GC content change?
   - How did CAI change?

This is highly relevant for biotechnology - expressing parasite proteins in bacteria for drug development!

#### B. Comparing Multiple Organisms

**Task**: Compare codon usage between *P. falciparum* and another organism

1. Download coding sequences from NCBI for:
   - *E. coli* (bacteria)
   - *S. cerevisiae* (yeast)
   - *H. sapiens* (human)
   - Or your favorite organism!

2. Analyze each dataset separately in GCUA
3. Compare:
   - Which has the strongest codon bias (lowest mean ENC)?
   - Are optimal codons the same or different?
   - How do GC3s values compare?

**Biological Insight**: This reveals how different organisms have "solved" the translation problem in different ways!

#### C. Cluster Analysis

**Task**: Use unsupervised clustering to group genes by codon usage

1. Use GCUA's clustering tools (**5. Advanced Tools** → **3. Cluster genes**)
2. Try different numbers of clusters (k = 2, 3, 4, 5)
3. Export cluster assignments
4. Characterize each cluster:
   - Mean ENC, CAI, Fop
   - GC content
   - Functional enrichment (what types of genes?)

---

## Part 6: Discussion Questions and Thinking Points

### Conceptual Questions

1. **Why is codon bias important for biotechnology?**
   - If you wanted to express a human protein in bacteria, why would you care about codon usage?
   - How might you use GCUA to design a better expression construct?

2. **Evolution and Selection**
   - Is codon bias an example of natural selection? Why or why not?
   - What is the selective advantage of using optimal codons?
   - Why don't all genes in a genome show maximum codon bias?

3. **Mutational Bias vs. Selection**
   - How can you distinguish between codon usage patterns caused by mutation vs. selection?
   - In *P. falciparum*, which force seems more important? Defend your answer.

4. **Clinical Relevance**
   - *P. falciparum* causes malaria. Could codon usage analysis help identify drug targets?
   - Highly expressed proteins might be good targets - why?

5. **Horizontal Gene Transfer**
   - If a gene was recently transferred from another organism, how might its codon usage differ from native genes?
   - How could you use GCUA to identify horizontally transferred genes?

### Critical Thinking Challenges

1. **Limitations of RSCU**:
   - RSCU is calculated across all genes. Why might this be misleading?
   - How would you identify optimal codons without simply taking the most common ones?

2. **CAI Assumptions**:
   - CAI assumes you know which codons are optimal. How can you be sure you've identified the right reference genes?
   - Could CAI be misleading in organisms with very low codon bias?

3. **Causation vs. Correlation**:
   - You find that high-CAI genes are enriched for ribosomal proteins. Does this prove that codon bias increases expression?
   - What alternative explanations could there be?
   - What experiment would you design to test causation?

---

## Part 7: Resources and Further Reading

### Key Papers (Recommended Reading)

1. **Sharp, P.M., & Li, W.H. (1987)**. The codon adaptation index - a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281-1295.
   - The original CAI paper

2. **Wright, F. (1990)**. The 'effective number of codons' used in a gene. *Gene*, 87(1), 23-29.
   - Introduces ENC and the ENC-GC3s plot

3. **Plotkin, J.B., & Kudla, G. (2011)**. Synonymous but not the same: the causes and consequences of codon bias. *Nature Reviews Genetics*, 12(1), 32-42.
   - Excellent modern review of codon bias

4. **Hershberg, R., & Petrov, D.A. (2008)**. Selection on codon bias. *Annual Review of Genetics*, 42, 287-299.
   - Discusses selection vs. drift in codon evolution

5. **Quax, T.E., et al. (2015)**. Codon bias as a means to fine-tune gene expression. *Molecular Cell*, 59(2), 149-161.
   - Modern perspectives on codon usage and regulation

### Online Resources

- **NCBI Genetic Codes**: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
- **Codon Usage Database**: http://www.kazusa.or.jp/codon/
- **PlasmoDB** (for *P. falciparum* annotations): https://plasmodb.org/

### Data Sources for Your Own Analyses

- **NCBI RefSeq**: https://www.ncbi.nlm.nih.gov/refseq/
- **Ensembl**: https://www.ensembl.org/
- **UniProt**: https://www.uniprot.org/

---

## Appendix: Quick Reference Guide

### GCUA Menu Structure

```
Main Menu
├── 1. Data Management
│   ├── 1. Load FASTA file
│   ├── 2. View data statistics
│   ├── 3. Filter sequences
│   └── 4. Return to main menu
│
├── 2. Quick Analysis (Guided Workflow)
│
├── 3. Custom Analysis
│   ├── 1. Basic Metrics
│   │   ├── 1. Codon Usage
│   │   ├── 2. Amino Acid Usage
│   │   └── 3. Base Composition
│   └── 2. Advanced Analysis
│       ├── 1. Calculate ENC
│       ├── 2. Calculate CAI
│       ├── 3. Calculate Fop
│       └── 4. Correspondence Analysis
│
├── 4. Visualization & Export
│   ├── 1. Create Visualizations
│   └── 2. Export Data
│
├── 5. Advanced Tools
│   ├── 1. Optimize Sequences
│   ├── 2. Identify reference genes
│   └── 3. Cluster genes
│
├── 6. Settings & Preferences
│   └── 1. Set genetic code (IMPORTANT!)
│
└── 7. Help & Documentation
```

### Common Commands

- **H**: Help at any menu
- **B**: Back to previous menu
- **Q**: Quit program
- **number?**: Help for specific option (e.g., "1?")

### Important File Locations

- Input: `*.fas` or `*.fasta` files
- Output: `gcua_outputs/` directory
- All outputs timestamped for easy tracking

### Troubleshooting Tips

1. **Wrong genetic code**: Always set the correct genetic code first! (Option 6 → 1)
2. **Memory issues**: For very large datasets, GCUA automatically switches to low-memory mode
3. **Missing dependencies**: Check that all Python packages are installed: `pip install -r requirements.txt`
4. **Can't find output files**: Look in the `gcua_outputs/` directory in the same folder as `gcua.py`

---

## Assessment Rubric (For Instructors)

### Exercise Completion (40 points)
- Exercises 1-3: 15 points (basic understanding)
- Exercises 4-5: 15 points (intermediate analysis)
- Exercises 6-7: 10 points (advanced synthesis)

### Independent Investigation (40 points)
- Research question clarity: 5 points
- Methods appropriateness: 10 points
- Results presentation: 10 points
- Biological interpretation: 10 points
- Critical thinking: 5 points

### Discussion Questions (20 points)
- Conceptual understanding: 10 points
- Critical analysis: 10 points

---

## Acknowledgments

This tutorial was developed for use with GCUA (General Codon Usage Analysis), based on the original work by James McInerney (1998). The Python implementation includes modern analytical and visualization capabilities while maintaining the core functionality of the original tool.

**Citation**:
McInerney JO. GCUA: general codon usage analysis. *Bioinformatics*. 1998;14(4):372-3.

---

## Need Help?

- Type **H** in any GCUA menu for help
- Check the `manual.md` file for detailed documentation
- See the `README.md` for installation and troubleshooting

**Good luck with your investigation into synonymous codon usage!**
