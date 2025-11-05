# GCUA Tutorial - Dataset Preparation Guide

## Overview

This guide helps instructors prepare datasets for use in the GCUA tutorial. It covers how to obtain, format, and validate sequence data for different organisms and teaching scenarios.

---

## Table of Contents

1. [Quick Start: Using Existing Datasets](#quick-start-using-existing-datasets)
2. [Downloading New Datasets](#downloading-new-datasets)
3. [Data Formatting Requirements](#data-formatting-requirements)
4. [Creating Custom Subsets](#creating-custom-subsets)
5. [Validation and Quality Control](#validation-and-quality-control)
6. [Recommended Organisms for Teaching](#recommended-organisms-for-teaching)
7. [Troubleshooting Data Issues](#troubleshooting-data-issues)

---

## Quick Start: Using Existing Datasets

### Included Dataset

**`pfal.fas`** - *Plasmodium falciparum* genes
- **Size**: ~12 MB, ~5,000-6,000 genes
- **Characteristics**: Extremely AT-rich (~80% AT), moderate codon bias
- **Advantages**:
  - Medically relevant (malaria parasite)
  - Clear AT-bias patterns
  - Good for teaching compositional bias vs. selection
- **Disadvantages**:
  - Many hypothetical proteins (limited functional annotation)
  - Relatively weak overall codon bias

### When to Use Different Datasets

Consider alternatives if:
- Students want to analyze specific organisms
- You want stronger codon bias examples (*E. coli*, *S. cerevisiae*)
- You want weaker bias for contrast (*H. sapiens*)
- You want to demonstrate comparative analysis (multiple organisms)
- You need smaller datasets (faster analysis for demonstrations)

---

## Downloading New Datasets

### Option 1: NCBI RefSeq (Recommended for Complete Genomes)

**Best for**: Complete, well-annotated genomes

**Steps**:

1. Go to NCBI RefSeq: https://www.ncbi.nlm.nih.gov/refseq/

2. Search for your organism:
   - Use full scientific name (e.g., "Escherichia coli K-12")
   - Or use taxonomic ID

3. Navigate to the genome page

4. Download CDS (coding sequences):
   - Look for "CDS FASTA" or "protein-coding genes"
   - **Important**: Download **nucleotide** sequences, NOT protein sequences
   - File format: `.fna` or `.fasta`

5. Alternative: Use command-line (faster for large datasets)
```bash
# Example: Download E. coli K-12 MG1655 CDS
# First, find the RefSeq accession (e.g., NC_000913.3)
# Then use datasets or entrez-direct

# Using NCBI datasets tool
datasets download genome accession NC_000913.3 --include cds

# Or use FTP directly
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz
gunzip *.fna.gz
```

**Tip**: For bacteria, look for "reference" or "representative" genomes

### Option 2: Ensembl (Best for Eukaryotes)

**Best for**: Model eukaryotes with high-quality annotations

**Steps**:

1. Go to Ensembl: https://www.ensembl.org/

2. Select organism (e.g., Human, Mouse, Yeast)

3. Navigate to Downloads → FTP site

4. Download CDS FASTA:
```bash
# Example: Human CDS
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
gunzip Homo_sapiens.GRCh38.cds.all.fa.gz

# Example: Yeast CDS
wget ftp://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz
```

**Note**: Ensembl provides transcript-level CDS. One gene may have multiple transcripts (isoforms).

### Option 3: Organism-Specific Databases

**For specialized organisms**:

- **PlasmoDB** (*Plasmodium*): https://plasmodb.org/
- **WormBase** (*C. elegans*): https://wormbase.org/
- **FlyBase** (*Drosophila*): https://flybase.org/
- **SGD** (*S. cerevisiae*): https://www.yeastgenome.org/
- **TAIR** (*Arabidopsis*): https://www.arabidopsis.org/

These often have better annotations and additional metadata.

### Option 4: Codon Usage Database

**Best for**: Quick examples and comparisons

Website: http://www.kazusa.or.jp/codon/

- Provides pre-processed codon usage tables
- Good for validating your GCUA results
- Limited to codon frequencies (not full analysis)

---

## Data Formatting Requirements

### FASTA Format Requirements

GCUA requires standard FASTA format:

```
>sequence_id1 optional description
ATGAAACCCGGGTTT...
>sequence_id2 optional description
ATGCCCGGGAAA...
```

**Critical Requirements**:

1. **Nucleotide sequences** (A, T/U, G, C)
   - NOT protein sequences (GCUA cannot analyze protein FASTA)

2. **Divisible by 3** (complete codons)
   - Sequences must represent complete CDS (no partial codons)
   - Start codon to stop codon (or excluding stop codon)

3. **Standard FASTA format**
   - Header line starts with `>`
   - Sequence on subsequent lines
   - No line length restrictions (GCUA handles any line width)

4. **DNA or RNA alphabet**
   - A, T, G, C (DNA) or A, U, G, C (RNA)
   - Both work (T and U are equivalent)
   - Lowercase or uppercase both acceptable

### Common Format Issues and Fixes

#### Issue 1: Protein FASTA Instead of Nucleotide

**Problem**:
```
>protein1
MKKLAVACLAVF...
```

**Detection**: Presence of amino acid letters (M, K, L, etc.)

**Solution**: Download nucleotide CDS, not protein FASTA

---

#### Issue 2: Sequences Not Divisible by 3

**Problem**: Partial codons at sequence ends

**Detection**: GCUA will warn about sequences not divisible by 3

**Solution**: Clean sequences with this Python script:

```python
from Bio import SeqIO

def clean_fasta(input_file, output_file):
    """Remove sequences not divisible by 3 and trim to complete codons."""
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(input_file, 'fasta'):
            seq_len = len(record.seq)
            trimmed_len = seq_len - (seq_len % 3)
            if trimmed_len >= 300:  # Keep only sequences ≥100 codons
                record.seq = record.seq[:trimmed_len]
                SeqIO.write(record, out, 'fasta')

clean_fasta('raw_sequences.fasta', 'cleaned_sequences.fasta')
```

**Or use command-line (if you have `seqkit`)**:
```bash
seqkit seq --min-len 300 raw.fasta | \
seqkit subseq -r 1:-1 | \
awk '/^>/ {if (seq) {if (length(seq) % 3 == 0) print header ORS seq} header=$0; seq=""; next} {seq=seq $0} END {if (length(seq) % 3 == 0) print header ORS seq}' > cleaned.fasta
```

---

#### Issue 3: Sequences Contain Stop Codons

**Problem**: Some files include stop codons, some don't

**Detection**: GCUA may give warnings or unexpected results

**Solution**: Decide on policy and apply consistently

**Option A: Remove stop codons** (recommended)
```python
from Bio import SeqIO
from Bio.Seq import Seq

def remove_stop_codons(input_file, output_file, genetic_code=1):
    """Remove terminal stop codons from sequences."""
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(input_file, 'fasta'):
            # Check if last codon is stop
            last_codon = str(record.seq[-3:])
            if last_codon.upper() in ['TAA', 'TAG', 'TGA']:  # Standard code
                record.seq = record.seq[:-3]
            SeqIO.write(record, out, 'fasta')

remove_stop_codons('with_stops.fasta', 'no_stops.fasta')
```

**Option B: Document and keep consistent**

---

#### Issue 4: Non-Standard Characters

**Problem**: Sequences contain N, R, Y, or other IUPAC ambiguous bases

**Detection**: GCUA may skip or error on these sequences

**Solution**: Filter out ambiguous sequences

```python
from Bio import SeqIO

def filter_clean_sequences(input_file, output_file, max_n_percent=1.0):
    """Keep only sequences with standard bases (A,T,G,C)."""
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(input_file, 'fasta'):
            seq_str = str(record.seq).upper()
            # Count non-ATGC characters
            non_standard = sum(1 for base in seq_str if base not in 'ATGCU')
            if non_standard / len(seq_str) * 100 <= max_n_percent:
                SeqIO.write(record, out, 'fasta')

filter_clean_sequences('raw.fasta', 'clean.fasta', max_n_percent=1.0)
```

---

#### Issue 5: Very Long or Very Short Sequences

**Problem**:
- Very short sequences (<300 bp): Unreliable codon usage statistics
- Very long sequences (>15,000 bp): May be concatenated genes or chromosomes

**Solution**: Filter by length

```python
from Bio import SeqIO

def filter_by_length(input_file, output_file, min_len=300, max_len=15000):
    """Keep sequences within specified length range."""
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(input_file, 'fasta'):
            if min_len <= len(record.seq) <= max_len:
                SeqIO.write(record, out, 'fasta')

filter_by_length('all_genes.fasta', 'filtered_genes.fasta')
```

**Recommended thresholds**:
- Minimum: 300 bp (100 codons)
- Maximum: 15,000 bp (5,000 codons) or no max

---

## Creating Custom Subsets

### Subset by Gene Function

**Use case**: Create dataset of ribosomal proteins, metabolic genes, etc.

**Method 1: Using gene names (simple)**
```python
from Bio import SeqIO

def subset_by_keywords(input_file, output_file, keywords):
    """Extract sequences with keywords in header."""
    keywords = [k.lower() for k in keywords]
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(input_file, 'fasta'):
            header = record.description.lower()
            if any(keyword in header for keyword in keywords):
                SeqIO.write(record, out, 'fasta')

# Example: Extract ribosomal protein genes
subset_by_keywords('all_genes.fasta', 'ribosomal_genes.fasta',
                   keywords=['ribosomal', 'ribosome', '40S', '60S', 'RPS', 'RPL'])
```

**Method 2: Using annotation files (advanced)**

Requires GFF/GTF annotation file with functional annotations or GO terms.

```python
from Bio import SeqIO
import pandas as pd

def subset_by_annotation(fasta_file, annotation_file, gene_list, output_file):
    """Extract sequences based on gene ID list."""
    # Read gene list (from annotation processing)
    genes_to_keep = set(gene_list)

    with open(output_file, 'w') as out:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            gene_id = record.id.split('|')[0]  # Adjust parsing as needed
            if gene_id in genes_to_keep:
                SeqIO.write(record, out, 'fasta')
```

### Subset by Expression Level

**Use case**: Create high-expression vs. low-expression datasets for comparison

**Requirements**: Expression data (RNA-seq, microarray, etc.)

```python
from Bio import SeqIO
import pandas as pd

def subset_by_expression(fasta_file, expression_file, percentile, output_file):
    """Extract high or low expression genes."""
    # Read expression data (CSV with gene_id, expression columns)
    expr = pd.read_csv(expression_file)

    # Get genes above/below percentile
    threshold = expr['expression'].quantile(percentile / 100)
    if percentile > 50:
        selected_genes = set(expr[expr['expression'] >= threshold]['gene_id'])
    else:
        selected_genes = set(expr[expr['expression'] <= threshold]['gene_id'])

    # Extract sequences
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id in selected_genes:
                SeqIO.write(record, out, 'fasta')

# Example: Top 10% expressed genes
subset_by_expression('all_genes.fasta', 'expression_data.csv',
                    percentile=90, output_file='high_expression.fasta')
```

### Random Sampling (for Quick Tests)

**Use case**: Create small test datasets for demonstrations

```python
from Bio import SeqIO
import random

def random_sample(input_file, output_file, n_sequences):
    """Randomly sample n sequences."""
    records = list(SeqIO.parse(input_file, 'fasta'))
    sampled = random.sample(records, min(n_sequences, len(records)))
    SeqIO.write(sampled, output_file, 'fasta')

# Example: Sample 100 genes for quick demo
random_sample('pfal.fas', 'pfal_small.fas', n_sequences=100)
```

---

## Validation and Quality Control

### Pre-Analysis Checks

**Before giving dataset to students**, verify:

1. **File format is correct**:
```bash
head -20 dataset.fasta
# Should see FASTA format with > headers and ATGC sequences
```

2. **Sequences are nucleotides, not proteins**:
```python
from Bio import SeqIO
record = next(SeqIO.parse('dataset.fasta', 'fasta'))
seq = str(record.seq).upper()
if any(aa in seq for aa in 'EFIJLOPQXZ'):
    print("ERROR: This is a protein sequence!")
else:
    print("OK: Nucleotide sequence")
```

3. **Sequences are divisible by 3**:
```python
from Bio import SeqIO
for record in SeqIO.parse('dataset.fasta', 'fasta'):
    if len(record.seq) % 3 != 0:
        print(f"WARNING: {record.id} length not divisible by 3")
```

4. **No empty sequences**:
```bash
grep -c "^>" dataset.fasta  # Count sequences
# Should match number of records
```

5. **Reasonable sequence lengths**:
```python
from Bio import SeqIO
lengths = [len(record.seq) for record in SeqIO.parse('dataset.fasta', 'fasta')]
print(f"Min: {min(lengths)}, Max: {max(lengths)}, Mean: {sum(lengths)/len(lengths):.0f}")
```

### Test Run with GCUA

**Always test dataset before distributing**:

1. Load into GCUA
2. Check how many sequences load successfully
3. Run basic codon usage analysis
4. Verify output files are generated
5. Check visualizations render correctly

**Red flags**:
- Many sequences rejected
- Crashes or errors
- Empty output files
- Visualizations show extreme outliers (may indicate data quality issues)

---

## Recommended Organisms for Teaching

### For Strong Codon Bias Examples

**1. *Escherichia coli* K-12**
- **Codon bias**: Moderate to strong
- **GC content**: ~50%
- **Advantages**: Well-studied, clear optimal codons, good reference set
- **Dataset size**: ~4,300 genes
- **Download**: NCBI RefSeq NC_000913.3

**2. *Saccharomyces cerevisiae* (Yeast)**
- **Codon bias**: Moderate
- **GC content**: ~40%
- **Advantages**: Model eukaryote, excellent annotations, lots of expression data available
- **Dataset size**: ~6,000 genes
- **Download**: Ensembl or SGD

**3. *Bacillus subtilis***
- **Codon bias**: Moderate to strong
- **GC content**: ~43%
- **Advantages**: Gram-positive model, industrial applications
- **Dataset size**: ~4,100 genes

### For Weak Codon Bias Examples

**1. *Homo sapiens* (Human)**
- **Codon bias**: Weak
- **GC content**: ~40-42%
- **Advantages**: Medical relevance, isochores (GC variation), large genome
- **Dataset size**: ~20,000 genes (may want to subset)
- **Download**: Ensembl

**2. *Arabidopsis thaliana***
- **Codon bias**: Weak to moderate
- **GC content**: ~36%
- **Advantages**: Plant model, weak but detectable bias
- **Dataset size**: ~27,000 genes (subset recommended)

### For Extreme Compositional Bias

**1. *Plasmodium falciparum*** (included)
- **AT-rich**: ~80% AT
- **Codon bias**: Weak overall, moderate in subset
- **Advantages**: Extreme composition, medical relevance

**2. *Borrelia burgdorferi* (Lyme disease)**
- **GC content**: ~28% (AT-rich)
- **Advantages**: Extreme AT bias, small genome
- **Dataset size**: ~1,500 genes

**3. *Streptomyces coelicolor***
- **GC content**: ~72% (GC-rich)
- **Advantages**: Extreme GC bias, antibiotic producer
- **Dataset size**: ~8,000 genes

### For Comparative Analysis

**Suggested combinations**:

1. **Bias comparison**:
   - *E. coli* (strong bias) vs. *H. sapiens* (weak bias)

2. **Composition comparison**:
   - *P. falciparum* (AT-rich) vs. *Streptomyces* (GC-rich)

3. **Phylogenetic comparison**:
   - Three bacteria: *E. coli*, *B. subtilis*, *M. tuberculosis*
   - Three eukaryotes: *S. cerevisiae*, *C. elegans*, *H. sapiens*

---

## Troubleshooting Data Issues

### Issue: GCUA rejects many sequences

**Possible causes**:
1. Sequences not divisible by 3 → Use cleaning script
2. Protein sequences loaded by mistake → Download nucleotide CDS
3. Non-standard characters → Filter with script
4. Very short sequences → Apply length filter

### Issue: Results seem wrong (all RSCU ≈ 1.0)

**Possible causes**:
1. Dataset too small → Need at least 50-100 genes
2. Very weak codon bias organism → Expected (e.g., some eukaryotes)
3. Heterogeneous dataset → Check if combining different species/genetic codes

### Issue: Correspondence Analysis fails

**Possible causes**:
1. Too few sequences → Need at least 30-50 for meaningful CA
2. All sequences very similar → No variation to analyze
3. Dataset includes stop codons inconsistently → Clean and standardize

### Issue: CAI calculation gives strange results

**Possible causes**:
1. Wrong genetic code selected → Check Settings menu
2. Reference set too small → Need at least 50 genes for reference
3. Reference set unrepresentative → Manually curate reference genes

### Issue: Students get different results from each other

**Check**:
1. Are they using same input file?
2. Are they using same genetic code setting?
3. Did they filter sequences differently?
4. Are they analyzing same dataset but at different times (if data updated)?

---

## Creating Organism-Specific Documentation

For courses focusing on specific organisms, create supplementary documentation:

### Template: Organism-Specific Guide

```markdown
# GCUA Tutorial: [Organism Name]

## Background

[Brief intro to organism, why it's interesting, medical/industrial relevance]

## Expected Patterns

Based on literature:
- **GC content**: [XX]%
- **Codon bias strength**: [Weak/Moderate/Strong]
- **Optimal codons**: [List if known]
- **Key references**: [Citations]

## Dataset Information

- **Source**: [Database and accession]
- **Number of genes**: [N]
- **Mean gene length**: [XXX] bp
- **Genetic code**: [Standard/Mitochondrial/Other]

## Genetic Code Settings

**IMPORTANT**: Set genetic code to [X] in GCUA Settings menu!

## Expected Results

### Exercise 1
[What students should find...]

### Exercise 2
[What students should find...]

[etc.]

## Biological Context

[Interpret findings in context of organism biology]

## Additional Resources

- [Organism database]
- [Relevant papers]
```

---

## Advanced: Creating Synthetic Datasets

### For Teaching Specific Concepts

**Use case**: Create artificial datasets with known properties to teach concepts

**Example: Extreme codon bias dataset**

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random

def create_synthetic_dataset(n_genes, output_file, bias_strength='strong'):
    """Create synthetic genes with specified codon bias."""

    # Define optimal codons (for "strong bias" example)
    optimal_codons = {
        'Ala': 'GCG', 'Arg': 'CGC', 'Asn': 'AAC', 'Asp': 'GAC',
        'Cys': 'TGC', 'Gln': 'CAG', 'Glu': 'GAG', 'Gly': 'GGC',
        'His': 'CAC', 'Ile': 'ATC', 'Leu': 'CTG', 'Lys': 'AAG',
        'Met': 'ATG', 'Phe': 'TTC', 'Pro': 'CCG', 'Ser': 'TCG',
        'Thr': 'ACC', 'Trp': 'TGG', 'Tyr': 'TAC', 'Val': 'GTG'
    }

    amino_acids = list(optimal_codons.keys())
    records = []

    for i in range(n_genes):
        # Generate random amino acid sequence
        aa_seq = [random.choice(amino_acids) for _ in range(300)]  # 300 aa = 900 bp

        # Convert to codons
        if bias_strength == 'strong':
            # Always use optimal codon
            codon_seq = ''.join([optimal_codons[aa] for aa in aa_seq])
        elif bias_strength == 'weak':
            # Use random synonymous codon (not implemented here - would need codon table)
            codon_seq = ''.join([optimal_codons[aa] for aa in aa_seq])

        # Create record
        record = SeqRecord(
            Seq(codon_seq),
            id=f"synthetic_gene_{i+1}",
            description=f"Synthetic gene with {bias_strength} codon bias"
        )
        records.append(record)

    SeqIO.write(records, output_file, 'fasta')

# Create dataset
create_synthetic_dataset(100, 'synthetic_strong_bias.fasta', bias_strength='strong')
```

**Educational use**: Students can compare synthetic (perfect bias) vs. real organisms

---

## Quick Reference: Data Preparation Checklist

- [ ] Download CDS nucleotide sequences (not proteins)
- [ ] Check file is valid FASTA format
- [ ] Remove or trim sequences not divisible by 3
- [ ] Filter very short sequences (<300 bp)
- [ ] Remove sequences with ambiguous bases (if needed)
- [ ] Standardize stop codon inclusion/exclusion
- [ ] Test load in GCUA
- [ ] Run basic analysis to verify
- [ ] Document source and any processing steps
- [ ] Note correct genetic code setting for students

---

## Conclusion

Proper dataset preparation is crucial for successful GCUA tutorials. Well-formatted, validated datasets ensure:
- Students can focus on biology, not troubleshooting
- Results are interpretable and educational
- Analysis runs smoothly without errors

**Key principle**: Test everything yourself before distributing to students!

---

## Resources

### Tools

- **BioPython**: Sequence manipulation (https://biopython.org/)
- **seqkit**: Fast sequence filtering (https://bioinf.shenwei.me/seqkit/)
- **NCBI Datasets**: Command-line genome download (https://www.ncbi.nlm.nih.gov/datasets/)

### Documentation

- **FASTA format**: https://en.wikipedia.org/wiki/FASTA_format
- **NCBI genetic codes**: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
- **GCUA manual**: See `manual.md` in repository

---

## Contact

For questions about dataset preparation for GCUA tutorial, please consult the GCUA repository or contact your bioinformatics support team.

**Happy teaching!**
