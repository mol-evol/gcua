"""
GCUA: General Codon Usage Analysis
Improved Python implementation of GCUA for analyzing codon usage in DNA sequences.

Original C program by James McInerney (1997)
Python implementation by James McInerney (2025)
Enhanced version with additional features including full genetic code support
"""

import argparse
import os
import sys
import math
import numpy as np
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
import pickle
import json
from datetime import datetime
import threading
try:
    from Bio.SeqUtils import GC
except ImportError:
    from Bio.SeqUtils import gc_fraction as GC
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
try:
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans, DBSCAN
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import StandardScaler
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    PCA = None
    KMeans = None
    DBSCAN = None
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import webbrowser
from scipy.spatial.distance import cdist
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Optional, Union, Any, Callable
from scipy.stats import chi2_contingency, chi2
import json

# Program constants
VERSION = "2.5.0"  # Added parallel processing, checkpointing, and batch processing
PROGRAM_NAME = "GCUA: General Codon Usage Analysis"
AUTHOR = "James McInerney"

# Standard Genetic code definitions (NCBI translation table IDs)
GENETIC_CODES = {
    1: "Standard (Universal)",
    2: "Vertebrate Mitochondrial",
    3: "Yeast Mitochondrial",
    4: "Mold, Protozoan, Coelenterate Mitochondrial & Mycoplasma/Spiroplasma",
    5: "Invertebrate Mitochondrial",
    6: "Ciliate, Dasycladacean and Hexamita Nuclear",
    9: "Echinoderm and Flatworm Mitochondrial",
    10: "Euplotid Nuclear",
    11: "Bacterial, Archaeal, and Plant Plastid",
    12: "Alternative Yeast Nuclear",
    13: "Ascidian Mitochondrial",
    14: "Alternative Flatworm Mitochondrial",
    16: "Chlorophycean Mitochondrial",
    21: "Trematode Mitochondrial",
    22: "Scenedesmus obliquus Mitochondrial",
    23: "Thraustochytrium Mitochondrial",
    24: "Pterobranchia Mitochondrial",
    25: "Candidate Division SR1 and Gracilibacteria",
    26: "Pachysolen tannophilus Nuclear",
    27: "Karyorelict Nuclear",
    28: "Condylostoma Nuclear",
    29: "Mesodinium Nuclear",
    30: "Peritrich Nuclear",
    31: "Blastocrithidia Nuclear",
    33: "Cephalodiscidae Mitochondrial UAA-Tyr"
}

# Codon tables
CODONS = [
    'UUU', 'UUC', 'UUA', 'UUG', 'UCU', 'UCC', 'UCA', 'UCG',
    'UAU', 'UAC', 'UAA', 'UAG', 'UGU', 'UGC', 'UGA', 'UGG',
    'CUU', 'CUC', 'CUA', 'CUG', 'CCU', 'CCC', 'CCA', 'CCG',
    'CAU', 'CAC', 'CAA', 'CAG', 'CGU', 'CGC', 'CGA', 'CGG',
    'AUU', 'AUC', 'AUA', 'AUG', 'ACU', 'ACC', 'ACA', 'ACG',
    'AAU', 'AAC', 'AAA', 'AAG', 'AGU', 'AGC', 'AGA', 'AGG',
    'GUU', 'GUC', 'GUA', 'GUG', 'GCU', 'GCC', 'GCA', 'GCG',
    'GAU', 'GAC', 'GAA', 'GAG', 'GGU', 'GGC', 'GGA', 'GGG'
]

# DNA equivalents of the RNA codons
DNA_CODONS = [codon.replace('U', 'T') for codon in CODONS]

# =====================================================
# GENETIC CODE TRANSLATION TABLES
# =====================================================

# Standard (Universal) code (NCBI translation table 1)
AA_MAP_1 = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

# Vertebrate Mitochondrial code (NCBI translation table 2)
AA_MAP_2 = AA_MAP_1.copy()
AA_MAP_2.update({
    'AUA': 'Met',   # Ile -> Met
    'UGA': 'Trp',   # STOP -> Trp
    'AGA': 'STOP',  # Arg -> STOP
    'AGG': 'STOP'   # Arg -> STOP
})

# Yeast Mitochondrial code (NCBI translation table 3)
AA_MAP_3 = AA_MAP_1.copy()
AA_MAP_3.update({
    'AUA': 'Met',   # Ile -> Met
    'UGA': 'Trp',   # STOP -> Trp
    'CUU': 'Thr',   # Leu -> Thr
    'CUC': 'Thr',   # Leu -> Thr
    'CUA': 'Thr',   # Leu -> Thr
    'CUG': 'Thr'    # Leu -> Thr
})

# Mold, Protozoan, Coelenterate Mitochondrial & Mycoplasma/Spiroplasma
# (NCBI translation table 4)
AA_MAP_4 = AA_MAP_1.copy()
AA_MAP_4.update({'UGA': 'Trp'})

# Invertebrate Mitochondrial code (NCBI translation table 5)
AA_MAP_5 = AA_MAP_1.copy()
AA_MAP_5.update({
    'AUA': 'Met',   # Ile -> Met
    'UGA': 'Trp',   # STOP -> Trp
    'AGA': 'Ser',   # Arg -> Ser
    'AGG': 'Ser'    # Arg -> Ser
})

# Ciliate, Dasycladacean and Hexamita Nuclear (NCBI translation table 6)
AA_MAP_6 = AA_MAP_1.copy()
AA_MAP_6.update({'UAA': 'Gln', 'UAG': 'Gln'})

# Echinoderm and Flatworm Mitochondrial (NCBI translation table 9)
AA_MAP_9 = AA_MAP_1.copy()
AA_MAP_9.update({
    'AAA': 'Asn', 'AGA': 'Ser', 'AGG': 'Ser',
    'UGA': 'Trp', 'AUA': 'Ile'
})

# Euplotid Nuclear (NCBI translation table 10)
AA_MAP_10 = AA_MAP_1.copy()
AA_MAP_10.update({'UGA': 'Cys'})

# Bacterial, Archaeal, and Plant Plastid (NCBI translation table 11)
AA_MAP_11 = AA_MAP_1.copy()
# Same as standard

# Alternative Yeast Nuclear (NCBI translation table 12)
AA_MAP_12 = AA_MAP_1.copy()
AA_MAP_12.update({'CUG': 'Ser'})  # Leu -> Ser

# Ascidian Mitochondrial (NCBI translation table 13)
AA_MAP_13 = AA_MAP_1.copy()
AA_MAP_13.update({
    'AUA': 'Met',   # Ile -> Met
    'UGA': 'Trp',   # STOP -> Trp
    'AGA': 'Gly',   # Arg -> Gly
    'AGG': 'Gly'    # Arg -> Gly
})

# Alternative Flatworm Mitochondrial (NCBI translation table 14)
AA_MAP_14 = AA_MAP_1.copy()
AA_MAP_14.update({
    'AAA': 'Asn', 'AGA': 'Ser', 'AGG': 'Ser',
    'UAA': 'Tyr', 'UGA': 'Trp'
})

# Chlorophycean Mitochondrial (NCBI translation table 16)
AA_MAP_16 = AA_MAP_1.copy()
AA_MAP_16.update({'UAG': 'Leu'})

# Trematode Mitochondrial (NCBI translation table 21)
AA_MAP_21 = AA_MAP_1.copy()
AA_MAP_21.update({
    'UGA': 'Trp', 'AUA': 'Met',
    'AAA': 'Asn', 'AGA': 'Ser', 'AGG': 'Ser'
})

# Scenedesmus obliquus Mitochondrial (NCBI translation table 22)
AA_MAP_22 = AA_MAP_1.copy()
AA_MAP_22.update({
    'UCA': 'STOP', 'UAG': 'Leu'
})

# Thraustochytrium Mitochondrial (NCBI translation table 23)
AA_MAP_23 = AA_MAP_1.copy()
AA_MAP_23.update({'UUA': 'STOP'})

# Pterobranchia Mitochondrial (NCBI translation table 24)
AA_MAP_24 = AA_MAP_1.copy()
AA_MAP_24.update({
    'AGA': 'Ser', 'AGG': 'Lys',
    'UGA': 'Trp'
})

# Candidate Division SR1 and Gracilibacteria (NCBI translation table 25)
AA_MAP_25 = AA_MAP_1.copy()
AA_MAP_25.update({'UGA': 'Gly'})

# Pachysolen tannophilus Nuclear (NCBI translation table 26)
AA_MAP_26 = AA_MAP_1.copy()
AA_MAP_26.update({'CUG': 'Ala'})

# Karyorelict Nuclear (NCBI translation table 27)
AA_MAP_27 = AA_MAP_1.copy()
AA_MAP_27.update({
    'UAA': 'Gln',   # STOP -> Gln
    'UAG': 'Gln',   # STOP -> Gln
    'UGA': 'Trp'    # STOP -> Trp
})

# Condylostoma Nuclear (NCBI translation table 28)
AA_MAP_28 = AA_MAP_1.copy()
AA_MAP_28.update({
    'UAA': 'Gln',   # STOP -> Gln
    'UAG': 'Gln',   # STOP -> Gln
    'UGA': 'Trp'    # STOP -> Trp
})

# Mesodinium Nuclear (NCBI translation table 29)
AA_MAP_29 = AA_MAP_1.copy()
AA_MAP_29.update({
    'UAA': 'Tyr', 'UAG': 'Tyr'
})

# Peritrich Nuclear (NCBI translation table 30)
AA_MAP_30 = AA_MAP_1.copy()
AA_MAP_30.update({
    'UAA': 'Glu', 'UAG': 'Glu'
})

# Blastocrithidia Nuclear (NCBI translation table 31)
AA_MAP_31 = AA_MAP_1.copy()
AA_MAP_31.update({'UGA': 'Trp', 'UAG': 'Glu', 'UAA': 'Glu'})

# Cephalodiscidae Mitochondrial UAA-Tyr (NCBI translation table 33)
AA_MAP_33 = AA_MAP_1.copy()
AA_MAP_33.update({
    'UAA': 'Tyr',   # STOP -> Tyr
    'UGA': 'Trp',   # STOP -> Trp
    'AGA': 'Ser',   # Arg -> Ser
    'AGG': 'Lys'    # Arg -> Lys
})

# Create a lookup table for all genetic code mappings
GENETIC_CODE_MAP = {
    1: AA_MAP_1,
    2: AA_MAP_2,
    3: AA_MAP_3,
    4: AA_MAP_4,
    5: AA_MAP_5,
    6: AA_MAP_6,
    9: AA_MAP_9,
    10: AA_MAP_10,
    11: AA_MAP_11,
    12: AA_MAP_12,
    13: AA_MAP_13,
    14: AA_MAP_14,
    16: AA_MAP_16,
    21: AA_MAP_21,
    22: AA_MAP_22,
    23: AA_MAP_23,
    24: AA_MAP_24,
    25: AA_MAP_25,
    26: AA_MAP_26,
    27: AA_MAP_27,
    28: AA_MAP_28,
    29: AA_MAP_29,
    30: AA_MAP_30,
    31: AA_MAP_31,
    33: AA_MAP_33
}

# Create DNA versions of the amino acid mappings
DNA_GENETIC_CODE_MAP = {}
for code_id, aa_map in GENETIC_CODE_MAP.items():
    DNA_GENETIC_CODE_MAP[code_id] = {
        codon.replace(
            'U',
            'T'): aa for codon,
        aa in aa_map.items()}

# Function to create a mapping of amino acids to their codons for any
# genetic code


def create_aa_to_codons_map(aa_map):
    aa_to_codons = defaultdict(list)
    for codon, aa in aa_map.items():
        aa_to_codons[aa].append(codon)
    return aa_to_codons


# Create a dictionary to store AA_TO_CODONS mappings for each genetic code
AA_TO_CODONS_MAP = {}
for code_id, aa_map in GENETIC_CODE_MAP.items():
    AA_TO_CODONS_MAP[code_id] = create_aa_to_codons_map(aa_map)

# Also create DNA_AA_TO_CODONS mappings for each genetic code
DNA_AA_TO_CODONS_MAP = {}
for code_id, dna_aa_map in DNA_GENETIC_CODE_MAP.items():
    DNA_AA_TO_CODONS_MAP[code_id] = create_aa_to_codons_map(dna_aa_map)


def process_sequence_batch(args):
    """
    Process a batch of sequences in parallel.
    This function is designed to be called by multiprocessing.
    """
    # Get process ID for debugging
    import os
    pid = os.getpid()
    
    # Unpack arguments - now includes file path and sequence info
    batch_info, genetic_code = args
    
    if isinstance(batch_info[0], str):  # If we received sequences directly
        sequences = batch_info
        gene_names = None
    else:  # If we received (gene_names, file_path, sequence_index)
        gene_names, file_path, sequence_indices = batch_info
        sequences = []
        
        # Open file once for the entire batch
        try:
            with open(file_path, 'r') as f:
                # Load sequences from file using indices
                for gene_name in gene_names:
                    if gene_name in sequence_indices:
                        try:
                            f.seek(sequence_indices[gene_name]['start'])
                            seq_data = ""
                            end_pos = sequence_indices[gene_name]['end']
                            
                            while f.tell() < end_pos:
                                line = f.readline().strip()
                                if line and not line.startswith('>'):
                                    seq_data += line
                            
                            sequences.append(seq_data)
                        except Exception as e:
                            sequences.append("")  # Empty sequence on error
                    else:
                        sequences.append("")  # Empty sequence if not in index
        except Exception as e:
            print(f"[Process {pid}] Error opening file: {e}")
            sequences = [""] * len(gene_names)
    
    # If gene_names not provided, use indices
    if gene_names is None:
        gene_names = [f"seq_{i}" for i in range(len(sequences))]
    
    # Get the appropriate amino acid mapping
    aa_map = GENETIC_CODE_MAP.get(genetic_code, GENETIC_CODE_MAP[1])
    
    # Initialize results dictionaries
    codon_counts = {}
    base_composition = {}
    
    for seq, gene_name in zip(sequences, gene_names):
        # Skip empty sequences
        if not seq:
            continue
            
        # Process each sequence
        dna_seq = seq.upper()
        length = len(dna_seq) - (len(dna_seq) % 3)
        
        # Initialize counts for this gene
        gene_codon_counts = defaultdict(int)
        base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        base_counts_1st = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        base_counts_2nd = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        base_counts_3rd = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        base_counts_3rd_4fold = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        
        codon_count = 0
        
        # Define 4-fold degenerate codons (3rd position doesn't change AA)
        fourfold_codons = set(['CUN', 'CCN', 'CGN', 'ACN', 'GUN', 'GCN', 'GGN'])
        
        # Count codons and bases
        for j in range(0, length, 3):
            if j + 3 <= length:
                codon_dna = dna_seq[j:j + 3]
                
                if 'N' in codon_dna or len(codon_dna) < 3:
                    continue
                
                # Update base counts
                for pos, base in enumerate(codon_dna):
                    if base in 'ATGC':
                        base_counts[base] += 1
                        if pos == 0:
                            base_counts_1st[base] += 1
                        elif pos == 1:
                            base_counts_2nd[base] += 1
                        elif pos == 2:
                            base_counts_3rd[base] += 1
                            # Check if this is a 4-fold degenerate codon
                            codon_pattern = codon_dna[:2] + 'N'
                            if codon_pattern.replace('T', 'U') in fourfold_codons:
                                base_counts_3rd_4fold[base] += 1
                
                # Convert to RNA codon
                codon_rna = codon_dna.replace('T', 'U')
                if codon_rna in CODONS:
                    gene_codon_counts[codon_rna] += 1
                    codon_count += 1
        
        # Calculate base composition stats
        base_comp = {
            'Length': length,
            'AA_count': codon_count,
            'A': base_counts['A'],
            'T': base_counts['T'],
            'G': base_counts['G'],
            'C': base_counts['C'],
            'A3': base_counts_3rd['A'],
            'T3': base_counts_3rd['T'],
            'G3': base_counts_3rd['G'],
            'C3': base_counts_3rd['C']
        }
        
        # Calculate GC content
        total_bases = sum(base_counts.values())
        if total_bases > 0:
            base_comp['GC'] = (base_counts['G'] + base_counts['C']) / total_bases * 100
        else:
            base_comp['GC'] = 0
            
        # Calculate positional GC content
        for pos, counts, label in [(1, base_counts_1st, 'GC1'), 
                                   (2, base_counts_2nd, 'GC2'), 
                                   (3, base_counts_3rd, 'GC3')]:
            total = sum(counts.values())
            if total > 0:
                base_comp[label] = (counts['G'] + counts['C']) / total * 100
            else:
                base_comp[label] = 0
        
        # Calculate GC3s (GC at 4-fold degenerate sites)
        total_3rd_4fold = sum(base_counts_3rd_4fold.values())
        if total_3rd_4fold > 0:
            base_comp['GC3s'] = (base_counts_3rd_4fold['G'] + base_counts_3rd_4fold['C']) / total_3rd_4fold * 100
        else:
            base_comp['GC3s'] = 0
        
        # Store results
        codon_counts[gene_name] = dict(gene_codon_counts)
        base_composition[gene_name] = base_comp
    
    return codon_counts, base_composition


class Config:
    def __init__(self):
        self.genetic_code = 1  # Default: Standard (Universal)
        self.detailed_output = 1  # 1=basic, 2=some details, 3=all details
        self.beep_on = False
        self.visualization_method = "plotly"  # "plotly" or "text"
        self.visualization_format = "html"  # "html", "svg", "png", "pdf", etc.
        
        # Memory management settings
        self.memory_mode = "auto"  # "auto", "low_memory", or "full_memory"
        self.memory_threshold_mb = 100  # File size threshold in MB for auto mode
        self.enable_sequence_indexing = True  # Enable fast disk-based sequence access
        
        # Performance settings
        self.enable_parallel_processing = True  # Use multiprocessing for large datasets
        self.num_processes = mp.cpu_count() - 1  # Number of processes (leave 1 for system)
        self.batch_size = 1000  # Number of sequences per batch
        self.enable_progress_saving = True  # Save progress for resuming
        self.checkpoint_interval = 5000  # Save checkpoint every N sequences

    def get_aa_map(self):
        """Return the appropriate amino acid mapping based on the genetic code setting."""
        if self.genetic_code in GENETIC_CODE_MAP:
            return GENETIC_CODE_MAP[self.genetic_code]
        else:
            print(
                f"Warning: Genetic code {self.genetic_code} not found. Using Standard code (1).")
            return GENETIC_CODE_MAP[1]

    def get_dna_aa_map(self):
        """Return the appropriate DNA amino acid mapping based on the genetic code setting."""
        if self.genetic_code in DNA_GENETIC_CODE_MAP:
            return DNA_GENETIC_CODE_MAP[self.genetic_code]
        else:
            print(
                f"Warning: Genetic code {self.genetic_code} not found. Using Standard code (1).")
            return DNA_GENETIC_CODE_MAP[1]

    def get_genetic_code_name(self):
        """Return the name of the current genetic code."""
        if self.genetic_code in GENETIC_CODES:
            return GENETIC_CODES[self.genetic_code]
        else:
            return "Unknown"


class Visualization(ABC):
    @abstractmethod
    def create(self, data, output_path=None):
        """Create visualization and optionally save to file"""
        pass

    def open_in_browser(self, path):
        """Open a saved visualization in the browser"""
        try:
            webbrowser.open('file://' + os.path.abspath(path), new=2)
            print(f"\nVisualization opened in browser.")
        except Exception as e:
            print(f"\nFailed to open visualization in browser: {e}")


class PlotlyVisualization(Visualization):
    """Base class for Plotly visualizations"""

    def create(self, data, output_path=None):
        fig = self._create_figure(data)

        if output_path:
            ext = os.path.splitext(output_path)[1].lower()
            
            # Check if this is a static image format
            if ext in ['.svg', '.png', '.pdf', '.eps', '.jpeg', '.jpg', '.webp']:
                # Check if kaleido is available for static image export
                try:
                    import kaleido
                    has_kaleido = True
                except ImportError:
                    has_kaleido = False
                    
                if has_kaleido:
                    # Export as static image
                    try:
                        fig.write_image(output_path)
                        print(f"Visualization saved to {output_path}")
                        return output_path
                    except Exception as e:
                        print(f"Error exporting to {ext}: {e}")
                        print("Falling back to HTML export...")
                        ext = '.html'  # Fall back to HTML
                else:
                    print(f"Warning: kaleido package not installed. Cannot export to {ext} format.")
                    print("Install with: pip install kaleido")
                    print("Falling back to HTML export...")
                    ext = '.html'
            
            # Default to HTML export
            if ext != '.html':
                # If no extension or unrecognized extension, default to .html
                html_path = os.path.splitext(output_path)[0] + '.html'
            else:
                html_path = output_path
                
            fig.write_html(
                html_path,
                include_plotlyjs='cdn',
                full_html=True,
                config={'responsive': True}
            )
            print(f"Visualization saved to {html_path}")
            return html_path

        return None

    @abstractmethod
    def _create_figure(self, data):
        """Create the Plotly figure object"""
        pass


class MultivariateVisualization(PlotlyVisualization):
    """Create multivariate analysis visualizations"""

    def _create_figure(self, data):
        # Extract data
        result = data
        coords = result['coordinates']
        explained_var = result['explained_variance']
        analysis_type = result.get('analysis_type', 'CA')
        data_type = result.get('data_type', 'RSCU')
        
        # Get selected axes or use defaults
        selected_axes = result.get('selected_axes', (0, 1))
        x_axis, y_axis = selected_axes

        # Calculate outliers based on distance from center using selected axes
        selected_coords = coords.iloc[:, [x_axis, y_axis]]
        centroid = selected_coords.mean().values.reshape(1, -1)
        distances = cdist(selected_coords.values,
                          centroid, 'euclidean').flatten()
        threshold = np.mean(distances) + 2 * np.std(distances)
        is_outlier = distances > threshold

        # Create the figure
        fig = go.Figure()

        # Add main cluster points
        fig.add_trace(go.Scatter(
            x=coords.iloc[~is_outlier, x_axis],
            y=coords.iloc[~is_outlier, y_axis],
            mode='markers',
            marker=dict(
                # Steel blue with transparency
                color='rgba(70, 130, 180, 0.5)',
                size=5
            ),
            text=coords.index[~is_outlier],
            hovertemplate='<b>%{text}</b><extra></extra>',
            name='Main cluster'
        ))

        # Add outlier points
        fig.add_trace(go.Scatter(
            x=coords.iloc[is_outlier, x_axis],
            y=coords.iloc[is_outlier, y_axis],
            mode='markers',
            marker=dict(
                color='rgba(255, 0, 0, 0.7)',  # Red with transparency
                size=8,
                line=dict(width=1, color='DarkSlateGrey')
            ),
            text=coords.index[is_outlier],
            hovertemplate='<b>%{text}</b><extra></extra>',
            name='Outliers'
        ))
        
        # Add codons/amino acids to the plot (only for CA)
        if analysis_type == 'CA' and 'loadings' in result:
            loadings = result['loadings']
            
            # Determine if we're showing codons or amino acids
            if data_type == 'AA':
                item_name = 'Amino acid'
                items_name = 'Amino acids'
            else:  # RSCU
                item_name = 'Codon'
                items_name = 'Codons'
            
            # Get axis names for hover template
            x_axis_name = coords.columns[x_axis]
            y_axis_name = coords.columns[y_axis]
            
            # Add codon/amino acid points
            fig.add_trace(go.Scatter(
                x=loadings.iloc[:, x_axis],
                y=loadings.iloc[:, y_axis],
                mode='markers+text',
                marker=dict(
                    color='rgba(0, 128, 0, 0.8)',  # Green
                    size=10,
                    symbol='diamond',
                    line=dict(width=1, color='darkgreen')
                ),
                text=loadings.index,
                textposition="top center",
                textfont=dict(size=8),
                hovertemplate=f'<b>{item_name}: %{{text}}</b><br>{x_axis_name}: %{{x:.3f}}<br>{y_axis_name}: %{{y:.3f}}<extra></extra>',
                name=items_name,
                visible=True  # Can be toggled in the legend
            ))

        # Create buttons for axis switching
        buttons = []
        axes_combinations = []

        # Get all valid combinations of axes (up to 5)
        max_dim = min(5, coords.shape[1])
        for i in range(max_dim):
            for j in range(i + 1, max_dim):
                axes_combinations.append((i, j))

        # Create buttons for all axes combinations
        for i, (dim1, dim2) in enumerate(axes_combinations):
            # If these dimensions exist in our data
            if dim1 < coords.shape[1] and dim2 < coords.shape[1]:
                dim1_name = coords.columns[dim1]
                dim2_name = coords.columns[dim2]

                # Create a new set of traces for this combination
                main_cluster_trace = go.Scatter(
                    x=coords.iloc[~is_outlier, dim1],
                    y=coords.iloc[~is_outlier, dim2],
                    mode='markers',
                    marker=dict(
                        color='rgba(70, 130, 180, 0.5)',
                        size=5
                    ),
                    text=coords.index[~is_outlier],
                    hovertemplate='<b>%{text}</b><extra></extra>',
                    name='Main cluster'
                )

                outlier_trace = go.Scatter(
                    x=coords.iloc[is_outlier, dim1],
                    y=coords.iloc[is_outlier, dim2],
                    mode='markers',
                    marker=dict(
                        color='rgba(255, 0, 0, 0.7)',
                        size=8,
                        line=dict(width=1, color='DarkSlateGrey')
                    ),
                    text=coords.index[is_outlier],
                    hovertemplate='<b>%{text}</b><extra></extra>',
                    name='Outliers'
                )
                
                # Prepare button args - include codons/amino acids if CA
                if analysis_type == 'CA' and 'loadings' in result:
                    loadings = result['loadings']
                    
                    # Determine if we're showing codons or amino acids
                    if data_type == 'AA':
                        item_name = 'Amino acid'
                        items_name = 'Amino acids'
                    else:  # RSCU
                        item_name = 'Codon'
                        items_name = 'Codons'
                    
                    codon_trace = go.Scatter(
                        x=loadings.iloc[:, dim1],
                        y=loadings.iloc[:, dim2],
                        mode='markers+text',
                        marker=dict(
                            color='rgba(0, 128, 0, 0.8)',
                            size=10,
                            symbol='diamond',
                            line=dict(width=1, color='darkgreen')
                        ),
                        text=loadings.index,
                        textposition="top center",
                        textfont=dict(size=8),
                        hovertemplate=f'<b>{item_name}: %{{text}}</b><br>' + dim1_name + f': %{{x:.3f}}<br>' + dim2_name + f': %{{y:.3f}}<extra></extra>',
                        name=items_name
                    )
                    x_data = [main_cluster_trace.x, outlier_trace.x, codon_trace.x]
                    y_data = [main_cluster_trace.y, outlier_trace.y, codon_trace.y]
                else:
                    x_data = [main_cluster_trace.x, outlier_trace.x]
                    y_data = [main_cluster_trace.y, outlier_trace.y]

                # Add to the button list
                button = dict(
                    method='update',
                    label=f"{dim1_name} vs {dim2_name}",
                    args=[
                        {'x': x_data, 'y': y_data},
                        {'xaxis.title': f"{dim1_name} ({explained_var[dim1]:.2%})",
                         'yaxis.title': f"{dim2_name} ({explained_var[dim2]:.2%})"}
                    ]
                )
                buttons.append(button)

        # Improved layout
        biplot_text = " (Biplot)" if analysis_type == 'CA' and 'loadings' in result else ""
        fig.update_layout(
            title=f"{analysis_type} of {data_type} data{biplot_text}",
            xaxis=dict(
                title=f"{coords.columns[x_axis]} ({explained_var[x_axis]:.2%})",
                showgrid=True,
                gridwidth=1,
                gridcolor='lightgray',
                zeroline=True,
                zerolinewidth=1,
                zerolinecolor='gray',
                showline=True,
                linewidth=1,
                linecolor='black'
            ),
            yaxis=dict(
                title=f"{coords.columns[y_axis]} ({explained_var[y_axis]:.2%})",
                showgrid=True,
                gridwidth=1,
                gridcolor='lightgray',
                zeroline=True,
                zerolinewidth=1,
                zerolinecolor='gray',
                showline=True,
                linewidth=1,
                linecolor='black'
            ),
            plot_bgcolor='white',
            width=1000,
            height=800,
            margin=dict(l=80, r=80, t=100, b=80),
            updatemenus=[dict(
                buttons=buttons,
                direction="right",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.5,
                xanchor="center",
                y=1.15,
                yanchor="top"
            )],
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.2,
                xanchor="center",
                x=0.5
            ),
            hoverlabel=dict(
                bgcolor="white",
                font_size=12,
                font_family="Arial"
            )
        )

        # Add info about number of genes and outliers
        info_text = f"Total genes: {len(coords)}<br>Outliers: {sum(is_outlier)}"
        if analysis_type == 'CA' and 'loadings' in result:
            if data_type == 'AA':
                info_text += f"<br>Amino acids: {len(result['loadings'])}"
            else:
                info_text += f"<br>Codons: {len(result['loadings'])}"
        
        fig.add_annotation(
            x=0.98,
            y=0.02,
            xref="paper",
            yref="paper",
            text=info_text,
            showarrow=False,
            font=dict(size=12),
            align="right",
            bordercolor="black",
            borderwidth=1,
            bgcolor="white",
            opacity=0.8
        )
        
        # Add variance explained annotation
        total_var = sum(explained_var[:2]) * 100
        var_text = f"<b>Variance Explained</b><br>Total (2D): {total_var:.1f}%"
        if len(explained_var) >= 3:
            total_var_3d = sum(explained_var[:3]) * 100
            var_text += f"<br>Total (3D): {total_var_3d:.1f}%"
        
        fig.add_annotation(
            x=0.02,
            y=0.98,
            xref="paper",
            yref="paper",
            text=var_text,
            showarrow=False,
            font=dict(size=14),
            align="left",
            bordercolor="darkblue",
            borderwidth=2,
            bgcolor="lightblue",
            opacity=0.9
        )

        return fig


class CustomScatterPlotVisualization(PlotlyVisualization):
    """Create custom scatter plots with Plotly."""

    def _create_figure(self, data):
        # Extract the x and y descriptions from column names if available,
        # otherwise use defaults
        x_desc = data.get('x_desc', 'X Value')
        y_desc = data.get('y_desc', 'Y Value')
        title = data.get('title', f"{y_desc} vs {x_desc}")

        # Create scatter plot
        fig = go.Figure()

        # Add scatter plot points
        fig.add_trace(go.Scatter(
            x=data['x'],
            y=data['y'],
            mode='markers',
            marker=dict(
                color='rgba(0, 100, 200, 0.7)',
                size=8
            ),
            text=data['gene'],
            hovertemplate='<b>%{text}</b><br>X: %{x:.4f}<br>Y: %{y:.4f}<extra></extra>'
        ))

        # Add regression line if we have enough data points
        if len(data) > 2:  # Need at least 3 points for regression
            slope, intercept = np.polyfit(data['x'], data['y'], 1)
            x_range = np.linspace(data['x'].min(), data['x'].max(), 100)
            y_fit = slope * x_range + intercept

            fig.add_trace(go.Scatter(
                x=x_range,
                y=y_fit,
                mode='lines',
                line=dict(color='red', dash='dash'),
                name=f'Regression (y = {slope:.4f}x + {intercept:.4f})'
            ))

            # Calculate correlation
            correlation = data['x'].corr(data['y'])

            # Add correlation annotation
            fig.add_annotation(
                x=0.05,
                y=0.95,
                xref="paper",
                yref="paper",
                text=f"r = {correlation:.4f}",
                showarrow=False,
                font=dict(size=14),
                bgcolor="white",
                bordercolor="black",
                borderwidth=1
            )

        # Update layout
        fig.update_layout(
            title=title,
            xaxis=dict(
                title=x_desc,
                showgrid=True,
                zeroline=True,
                zerolinecolor='gray',
                showline=True,
                linewidth=1,
                linecolor='black'
            ),
            yaxis=dict(
                title=y_desc,
                showgrid=True,
                zeroline=True,
                zerolinecolor='gray',
                showline=True,
                linewidth=1,
                linecolor='black'
            ),
            plot_bgcolor='white',
            width=900,
            height=700,
            margin=dict(l=80, r=80, t=100, b=80),
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99
            ),
            hovermode='closest'
        )

        # Add color-coded outliers if available
        if 'outlier' in data.columns:
            # Create a new figure with color-coded points
            fig_with_outliers = go.Figure()

            # Add non-outlier points
            non_outliers = data[~data['outlier']]
            fig_with_outliers.add_trace(
                go.Scatter(
                    x=non_outliers['x'],
                    y=non_outliers['y'],
                    mode='markers',
                    marker=dict(
                        color='rgba(0, 100, 200, 0.7)',
                        size=8),
                    text=non_outliers['gene'],
                    hovertemplate='<b>%{text}</b><br>X: %{x:.4f}<br>Y: %{y:.4f}<extra></extra>',
                    name='Regular points'))

            # Add outlier points
            outliers = data[data['outlier']]
            if len(outliers) > 0:
                fig_with_outliers.add_trace(
                    go.Scatter(
                        x=outliers['x'],
                        y=outliers['y'],
                        mode='markers',
                        marker=dict(
                            color='rgba(255, 0, 0, 0.7)',
                            size=10,
                            line=dict(
                                width=1,
                                color='black')),
                        text=outliers['gene'],
                        hovertemplate='<b>%{text}</b><br>X: %{x:.4f}<br>Y: %{y:.4f}<br>Distance: %{customdata:.4f}<extra></extra>',
                        customdata=outliers['distance'] if 'distance' in outliers.columns else None,
                        name='Outliers'))

                # Use the layout from the original figure
                fig_with_outliers.update_layout(fig.layout)

                # Add annotation about outliers
                fig_with_outliers.add_annotation(
                    x=0.05,
                    y=0.89,
                    xref="paper",
                    yref="paper",
                    text=f"Outliers: {len(outliers)} points",
                    showarrow=False,
                    font=dict(size=14, color='red'),
                    bgcolor="white",
                    bordercolor="red",
                    borderwidth=1
                )

                return fig_with_outliers

        return fig


class GCContentVisualization(PlotlyVisualization):
    """Create GC content visualization"""

    def _create_figure(self, data):
        # Extract GC and GC3 content data
        gc = data.base_composition['GC'].values
        gc3 = data.base_composition['GC3'].values

        # Create plot
        fig = go.Figure()

        # Add scatter plot
        fig.add_trace(go.Scatter(
            x=gc,
            y=gc3,
            mode='markers',
            marker=dict(
                color='blue',
                size=8,
                opacity=0.7
            ),
            text=data.gene_names,
            hovertemplate='<b>%{text}</b><br>GC: %{x:.2f}%<br>GC3: %{y:.2f}%<extra></extra>'
        ))

        # Add theoretical line (GC = GC3)
        x_range = np.linspace(min(gc), max(gc), 100)
        fig.add_trace(go.Scatter(
            x=x_range,
            y=x_range,
            mode='lines',
            line=dict(color='red', dash='dash'),
            name='GC = GC3'
        ))

        # Add regression line
        slope, intercept = np.polyfit(gc, gc3, 1)
        y_fit = slope * x_range + intercept
        fig.add_trace(go.Scatter(
            x=x_range,
            y=y_fit,
            mode='lines',
            line=dict(color='green'),
            name=f'Regression (y = {slope:.2f}x + {intercept:.2f})'
        ))

        # Update layout
        fig.update_layout(
            title='GC3 vs GC Content',
            xaxis=dict(
                title='GC Content (%)',
                showgrid=True,
                zeroline=True
            ),
            yaxis=dict(
                title='GC3 Content (%)',
                showgrid=True,
                zeroline=True
            ),
            plot_bgcolor='white',
            width=800,
            height=600,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=0.01
            )
        )

        return fig


class ENCPlotVisualization(PlotlyVisualization):
    """Create ENC vs GC3s plot visualization"""

    def _create_figure(self, data):
        # Calculate ENC values if not already present
        enc_values = data.calculate_enc()

        # Extract GC3s values
        gc3s_values = data.base_composition['GC3s'].values

        # Create theoretical ENC curve
        # Wright's formula: ENC = 2 + s + 29/(s^2 + (1-s)^2)
        # where s is the fraction of GC3s (0-1)
        x_range = np.linspace(0, 100, 100)  # GC3s from 0% to 100%
        y_range = []
        for gc3s in x_range:
            s = gc3s / 100  # Convert percentage to fraction
            if 0 <= s <= 1:
                enc = 2 + s + 29 / ((s**2) + ((1 - s)**2))
                y_range.append(min(enc, 61))
            else:
                y_range.append(np.nan)

        # Create plot
        fig = go.Figure()

        # Add scatter plot of actual values
        fig.add_trace(go.Scatter(
            x=gc3s_values,
            y=[enc_values[gene] for gene in data.gene_names],
            mode='markers',
            marker=dict(
                color='blue',
                size=8,
                opacity=0.7
            ),
            text=data.gene_names,
            hovertemplate='<b>%{text}</b><br>GC3s: %{x:.2f}%<br>ENC: %{y:.2f}<extra></extra>',
            name='Observed ENC'
        ))

        # Add theoretical curve
        fig.add_trace(go.Scatter(
            x=x_range,
            y=y_range,
            mode='lines',
            line=dict(color='red'),
            name='Expected ENC'
        ))

        # Update layout
        fig.update_layout(
            title="ENC vs GC3s (Wright's Plot)",
            xaxis=dict(
                title='GC3s (%)',
                showgrid=True,
                zeroline=True,
                range=[0, 100]
            ),
            yaxis=dict(
                title='ENC',
                showgrid=True,
                zeroline=True,
                range=[0, 62]
            ),
            plot_bgcolor='white',
            width=800,
            height=600,
            legend=dict(
                yanchor="bottom",
                y=0.01,
                xanchor="right",
                x=0.99
            )
        )

        return fig


class RSCUHeatmapVisualization(PlotlyVisualization):
    """Create RSCU heatmap visualization"""

    def _create_figure(self, data):
        # Get the RSCU values
        rscu_values = data.get_rscu_values()
        
        # Check dataset size and choose appropriate visualization
        num_genes = len(rscu_values)
        
        if num_genes > 1000:
            # For large datasets, show summary statistics
            return self._create_summary_heatmap(rscu_values, data.config)
        else:
            # For small datasets, show traditional gene-level heatmap
            return self._create_gene_heatmap(rscu_values, data.config)
    
    def _create_summary_heatmap(self, rscu_values, config):
        """Create summary statistics heatmap for large datasets"""
        # Filter out stop codons
        aa_map = config.get_aa_map()
        non_stop_codons = [c for c in CODONS if aa_map[c] != 'STOP']
        rscu_filtered = rscu_values[non_stop_codons]
        
        # Calculate summary statistics
        mean_rscu = rscu_filtered.mean()
        std_rscu = rscu_filtered.std()
        
        # Group codons by amino acid
        aa_to_codons = {}
        for codon in non_stop_codons:
            aa = aa_map[codon]
            if aa not in aa_to_codons:
                aa_to_codons[aa] = []
            aa_to_codons[aa].append(codon)
        
        # Prepare data for grouped heatmap
        codon_order = []
        aa_labels = []
        for aa in sorted(aa_to_codons.keys()):
            codons = sorted(aa_to_codons[aa])
            codon_order.extend(codons)
            aa_labels.extend([aa] * len(codons))
        
        # Create summary data matrix
        summary_data = []
        summary_labels = []
        
        # Mean RSCU
        summary_data.append([mean_rscu[codon] for codon in codon_order])
        summary_labels.append("Mean RSCU")
        
        # Standard deviation
        summary_data.append([std_rscu[codon] for codon in codon_order])
        summary_labels.append("Std Dev")
        
        # Percentiles
        for pct in [10, 25, 50, 75, 90]:
            pct_values = rscu_filtered[codon_order].quantile(pct/100)
            summary_data.append(list(pct_values))
            summary_labels.append(f"{pct}th percentile")
        
        # Create custom hover text
        hover_text = []
        for i, label in enumerate(summary_labels):
            row_hover = []
            for j, codon in enumerate(codon_order):
                aa = aa_labels[j]
                value = summary_data[i][j]
                if label == "Mean RSCU":
                    # Add interpretation for mean RSCU
                    if value > 1.5:
                        interpretation = "Strongly preferred"
                    elif value > 1.2:
                        interpretation = "Preferred"
                    elif value < 0.5:
                        interpretation = "Strongly avoided"
                    elif value < 0.8:
                        interpretation = "Avoided"
                    else:
                        interpretation = "Neutral"
                    row_hover.append(f"{aa} - {codon}<br>{label}: {value:.3f}<br>{interpretation}")
                else:
                    row_hover.append(f"{aa} - {codon}<br>{label}: {value:.3f}")
            hover_text.append(row_hover)
        
        # Create the heatmap
        fig = go.Figure(data=go.Heatmap(
            z=summary_data,
            x=[f"{aa}-{codon}" for aa, codon in zip(aa_labels, codon_order)],
            y=summary_labels,
            colorscale='RdBu_r',  # Red (low) to Blue (high)
            colorbar=dict(title='RSCU Value'),
            hovertext=hover_text,
            hovertemplate='%{hovertext}<extra></extra>',
            zmid=1.0  # Center color scale at RSCU=1
        ))
        
        # Update layout
        fig.update_layout(
            title=f'RSCU Summary Statistics<br><sub>Based on {len(rscu_values):,} genes</sub>',
            xaxis=dict(
                title='Codon (grouped by amino acid)',
                tickangle=-45,
                showgrid=False
            ),
            yaxis=dict(
                title='Statistic',
                showgrid=False
            ),
            width=1400,
            height=600
        )
        
        # Add vertical lines to separate amino acids
        aa_boundaries = []
        current_pos = -0.5
        for aa in sorted(aa_to_codons.keys()):
            current_pos += len(aa_to_codons[aa])
            aa_boundaries.append(current_pos)
        
        for boundary in aa_boundaries[:-1]:  # Don't add line after last group
            fig.add_vline(x=boundary, line_dash="dash", line_color="gray", opacity=0.5)
        
        # Add annotation about dataset size
        fig.add_annotation(
            xref="paper", yref="paper",
            x=0.02, y=0.98,
            text=f"Large dataset mode: Showing summary statistics",
            showarrow=False,
            font=dict(size=12, color="gray"),
            align="left"
        )
        
        return fig
    
    def _create_gene_heatmap(self, rscu_values, config):
        """Create traditional gene-level heatmap for small datasets"""
        # Filter out stop codons
        aa_map = config.get_aa_map()
        non_stop_codons = [c for c in CODONS if aa_map[c] != 'STOP']
        rscu_filtered = rscu_values[non_stop_codons]

        # Create the heatmap
        fig = go.Figure(data=go.Heatmap(
            z=rscu_filtered.values.T,
            x=rscu_filtered.index,
            y=non_stop_codons,
            colorscale='RdBu_r',
            colorbar=dict(title='RSCU Value'),
            hovertemplate='Gene: %{x}<br>Codon: %{y}<br>RSCU: %{z:.2f}<extra></extra>',
            zmid=1.0
        ))

        # Update layout
        fig.update_layout(
            title=f'RSCU Values Heatmap<br><sub>{len(rscu_values)} genes × {len(non_stop_codons)} codons</sub>',
            xaxis=dict(
                title='Gene',
                tickangle=-45
            ),
            yaxis=dict(
                title='Codon',
                categoryorder='trace'
            ),
            width=1200,
            height=800
        )

        return fig


class CAIDistributionVisualization(PlotlyVisualization):
    """Create CAI distribution visualization"""

    def _create_figure(self, data):
        # Get CAI values
        if isinstance(data, dict):
            cai_values = data
        else:
            # Calculate CAI if not already done
            if not hasattr(data, 'cai_values') or not data.cai_values:
                print("Calculating Codon Adaptation Index (CAI)...")
                cai_values = data.calculate_cai()
            else:
                cai_values = data.cai_values

        # Convert dictionary to DataFrame for easier plotting
        cai_df = pd.DataFrame(
            list(
                cai_values.items()),
            columns=[
                'Gene',
                'CAI'])

        # Create the figure
        fig = go.Figure()

        # Add histogram
        fig.add_trace(go.Histogram(
            x=cai_df['CAI'],
            nbinsx=20,
            marker_color='rgba(0, 100, 200, 0.7)',
            name='CAI Distribution'
        ))

        # Add individual points as a rug plot
        fig.add_trace(go.Scatter(
            x=cai_df['CAI'],
            y=[0] * len(cai_df),
            mode='markers',
            marker=dict(
                symbol='line-ns',
                color='black',
                size=8,
                line=dict(width=1)
            ),
            name='Individual Genes',
            hovertext=cai_df['Gene'],
            hovertemplate='Gene: %{hovertext}<br>CAI: %{x:.4f}<extra></extra>'
        ))

        # Update layout
        fig.update_layout(
            title='Distribution of CAI Values',
            xaxis=dict(
                title='CAI',
                range=[0, 1],
                showgrid=True
            ),
            yaxis=dict(
                title='Frequency',
                showgrid=True
            ),
            plot_bgcolor='white',
            width=800,
            height=500,
            barmode='overlay',
            bargap=0.1
        )

        return fig


class TextVisualization(Visualization):
    """Base class for text-based visualizations"""

    def create(self, data, output_path=None):
        """Create a text visualization and optionally save to file"""
        text_content = self._create_text_visualization(data)
        print(text_content)

        if output_path:
            text_path = os.path.splitext(output_path)[0] + '.txt'
            with open(text_path, 'w') as f:
                f.write(text_content)
            print(f"Text visualization saved to {text_path}")
            return text_path

        return None

    @abstractmethod
    def _create_text_visualization(self, data):
        """Create the text visualization"""
        pass


class ScreePlotVisualization(PlotlyVisualization):
    """Create scree plot for multivariate analysis showing variance explained by each component."""
    
    def _create_figure(self, data):
        result = data.multivariate_results
        if not result:
            raise ValueError("No multivariate analysis results available")
        
        explained_var = result['explained_variance']
        analysis_type = result['analysis_type']
        data_type = result['data_type']
        
        # Calculate cumulative variance
        cumulative_var = np.cumsum(explained_var)
        
        # Create figure with secondary y-axis
        fig = make_subplots(
            rows=1, cols=1,
            specs=[[{"secondary_y": True}]]
        )
        
        # Add bar plot for individual variance explained
        fig.add_trace(
            go.Bar(
                x=[f'Component {i+1}' for i in range(len(explained_var))],
                y=explained_var * 100,  # Convert to percentage
                name='Individual Variance',
                marker_color='steelblue',
                hovertemplate='Component %{x}<br>Variance: %{y:.2f}%<extra></extra>'
            ),
            secondary_y=False
        )
        
        # Add line plot for cumulative variance
        fig.add_trace(
            go.Scatter(
                x=[f'Component {i+1}' for i in range(len(explained_var))],
                y=cumulative_var * 100,  # Convert to percentage
                name='Cumulative Variance',
                mode='lines+markers',
                line=dict(color='red', width=2),
                marker=dict(size=8),
                hovertemplate='Component %{x}<br>Cumulative: %{y:.2f}%<extra></extra>'
            ),
            secondary_y=True
        )
        
        # Update layout
        fig.update_xaxes(title_text="Principal Components")
        fig.update_yaxes(title_text="Variance Explained (%)", secondary_y=False)
        fig.update_yaxes(title_text="Cumulative Variance (%)", secondary_y=True)
        
        # Add 80% threshold line
        fig.add_hline(
            y=80,
            line_dash="dash",
            line_color="gray",
            annotation_text="80% threshold",
            annotation_position="right",
            secondary_y=True
        )
        
        fig.update_layout(
            title=f"Scree Plot - {analysis_type} of {data_type} data",
            plot_bgcolor='white',
            width=900,
            height=600,
            margin=dict(l=80, r=80, t=100, b=80),
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.2,
                xanchor="center",
                x=0.5
            ),
            hovermode='x'
        )
        
        # Add annotation with key statistics
        total_components = len(explained_var)
        if any(cumulative_var >= 0.8):
            components_80 = np.argmax(cumulative_var >= 0.8) + 1
            components_80_text = f"Components for 80% variance: {components_80}"
        else:
            max_var = cumulative_var[-1] * 100
            components_80_text = f"Components for 80% variance: >{total_components}<br>Max available: {max_var:.1f}%"
        
        fig.add_annotation(
            x=0.98,
            y=0.02,
            xref="paper",
            yref="paper",
            text=f"Total components: {total_components}<br>" +
                 f"{components_80_text}<br>" +
                 f"First component: {explained_var[0]*100:.1f}%",
            showarrow=False,
            font=dict(size=12),
            align="right",
            bordercolor="black",
            borderwidth=1,
            bgcolor="white",
            opacity=0.8
        )
        
        return fig


class Multivariate3DVisualization(PlotlyVisualization):
    """Create 3D scatter plot for multivariate analysis results."""
    
    def _create_figure(self, data):
        result = data.multivariate_results
        if not result:
            raise ValueError("No multivariate analysis results available")
        
        coords = result['coordinates']
        explained_var = result['explained_variance']
        analysis_type = result['analysis_type']
        data_type = result['data_type']
        
        # Check if we have at least 3 components
        if coords.shape[1] < 3:
            raise ValueError("Need at least 3 components for 3D visualization")
        
        # Create 3D scatter plot
        fig = go.Figure()
        
        # Add gene points
        fig.add_trace(go.Scatter3d(
            x=coords.iloc[:, 0],
            y=coords.iloc[:, 1],
            z=coords.iloc[:, 2],
            mode='markers',
            marker=dict(
                size=5,
                color=coords.iloc[:, 2],  # Color by third component
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(
                    title=f"{coords.columns[2]}<br>({explained_var[2]:.2%})"
                )
            ),
            text=coords.index,
            hovertemplate='<b>%{text}</b><br>' +
                         f'{coords.columns[0]}: %{{x:.3f}}<br>' +
                         f'{coords.columns[1]}: %{{y:.3f}}<br>' +
                         f'{coords.columns[2]}: %{{z:.3f}}<extra></extra>'
        ))
        
        # Update layout
        fig.update_layout(
            title=f"3D {analysis_type} of {data_type} data",
            scene=dict(
                xaxis=dict(
                    title=f"{coords.columns[0]} ({explained_var[0]:.2%})",
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray'
                ),
                yaxis=dict(
                    title=f"{coords.columns[1]} ({explained_var[1]:.2%})",
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray'
                ),
                zaxis=dict(
                    title=f"{coords.columns[2]} ({explained_var[2]:.2%})",
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray'
                ),
                bgcolor='white'
            ),
            width=1000,
            height=800,
            margin=dict(l=0, r=0, t=50, b=0)
        )
        
        # Add annotation with variance info
        total_var_3d = sum(explained_var[:3]) * 100
        fig.add_annotation(
            x=0.02,
            y=0.98,
            xref="paper",
            yref="paper",
            text=f"Total variance (3D): {total_var_3d:.1f}%",
            showarrow=False,
            font=dict(size=14),
            bgcolor="white",
            bordercolor="black",
            borderwidth=1
        )
        
        return fig


class VisualizationFactory:
    @staticmethod
    def create_visualization(vis_type, vis_method="plotly"):
        """Create appropriate visualization based on type and method"""
        if vis_method == "plotly":
            if vis_type == "multivariate":
                return MultivariateVisualization()
            elif vis_type == "gc_content":
                return GCContentVisualization()
            elif vis_type == "enc":
                return ENCPlotVisualization()
            elif vis_type == "rscu_heatmap":
                return RSCUHeatmapVisualization()
            elif vis_type == "cai_distribution":
                return CAIDistributionVisualization()
            elif vis_type == "scree_plot":
                return ScreePlotVisualization()
            elif vis_type == "multivariate_3d":
                return Multivariate3DVisualization()
            elif vis_type == "custom_scatter":
                return CustomScatterPlotVisualization()
        else:  # text visualizations
            # Implement text visualization classes in the next version
            pass

        raise ValueError(
            f"Unsupported visualization type: {vis_type} with method: {vis_method}")


class FileManager:
    @staticmethod
    def get_output_path(data, base_filename, output_type, extension=None):
        """Generate output file paths for various output types."""
        # Create a directory for outputs if it doesn't exist
        output_dir = "gcua_outputs"
        os.makedirs(output_dir, exist_ok=True)

        # Get base name from input file if none provided
        if not base_filename:
            base_filename = os.path.splitext(
                os.path.basename(data.file_path))[0]

        # Generate paths based on output type
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")

        base_path = os.path.join(
            output_dir,
            f"{base_filename}_{output_type}_{timestamp}")
        
        # Add extension if provided
        if extension:
            return base_path + extension
        return base_path

    @staticmethod
    def save_data_to_file(
            data,
            output_path,
            header=None,
            delimiter="\t",
            metadata=None):
        """Save data to a file with specified format and optional metadata"""
        with open(output_path, 'w') as f:
            # Add metadata if provided
            if metadata:
                f.write("# GCUA Metadata\n")
                for key, value in metadata.items():
                    f.write(f"# {key}: {value}\n")
                f.write("#" + "-" * 50 + "\n\n")

            if header:
                f.write(header + "\n")

            if isinstance(data, pd.DataFrame):
                data.to_csv(f, sep=delimiter, header=True, index=True)
            elif isinstance(data, dict):
                for key, value in data.items():
                    if isinstance(value, (list, tuple)):
                        f.write(
                            f"{key}{delimiter}" +
                            delimiter.join(
                                str(v) for v in value) +
                            "\n")
                    else:
                        f.write(f"{key}{delimiter}{value}\n")
            elif isinstance(data, str):
                f.write(data)
            else:
                f.write(str(data))

        print(f"Data saved to {output_path}")
        return output_path


class GCUAData:
    def __init__(self, config):
        self.config = config
        self.sequences = []  # List of Bio.SeqRecord objects
        self.gene_names = []  # List of gene names
        self.codon_counts = pd.DataFrame()  # Codon counts for each gene
        self.rscu_values = pd.DataFrame()  # RSCU values for each gene
        self.amino_acid_usage = pd.DataFrame()  # Amino acid counts for each gene
        self.base_composition = pd.DataFrame()  # Base composition for each gene
        self.file_path = None  # Path to the loaded FASTA file
        self.multivariate_results = None  # Store multivariate analysis results
        self.enc_values = None  # Store ENC values
        self.fop_values = None  # Store Fop values
        self.cai_values = None  # Store CAI values
        self.optimal_codons = None  # Store optimal codons
        self.visualization_factory = VisualizationFactory()
        
        # Memory management attributes
        self.memory_mode = None  # Will be set based on file size or config
        self.sequence_index = {}  # Maps gene_id to file position for fast access
        self.file_size_mb = 0  # Store file size for decision making
        
        # Analysis status tracking
        self.analysis_status = {
            'basic_metrics': False,
            'codon_usage': False,
            'amino_acid_usage': False, 
            'base_composition': False,
            'rscu_calculated': False,
            'enc_calculated': False,
            'cai_calculated': False,
            'fop_calculated': False,
            'scuo_calculated': False,
            'optimal_codons': False,
            'multivariate_performed': False,
            'clusters_detected': False,
            'axis_cohorts_compared': False
        }
        
        # Progress tracking
        self.checkpoint_file = None  # Path to checkpoint file
        self.processed_genes = set()  # Track which genes have been processed
        
        # Check parallel computing capabilities
        self._check_parallel_backends()
    
    def _check_parallel_backends(self):
        """Check available parallel computing backends."""
        self.parallel_info = {
            'numpy_threads': 1,
            'blas_library': 'unknown',
            'has_mkl': False,
            'has_openblas': False
        }
        
        try:
            # Check NumPy configuration
            import numpy as np
            numpy_config = np.show_config(mode='dicts')
            
            # Try to detect BLAS library
            if hasattr(np, '__config__'):
                config_info = str(np.__config__.show())
                if 'mkl' in config_info.lower():
                    self.parallel_info['blas_library'] = 'MKL'
                    self.parallel_info['has_mkl'] = True
                elif 'openblas' in config_info.lower():
                    self.parallel_info['blas_library'] = 'OpenBLAS'
                    self.parallel_info['has_openblas'] = True
            
            # Check current thread settings
            import os
            self.parallel_info['numpy_threads'] = int(os.environ.get('OMP_NUM_THREADS', mp.cpu_count()))
            
        except:
            pass  # Silently fail if we can't detect configuration

    def load_fasta(self, file_path):
        """
        Load sequences from a FASTA file with hybrid memory management.
        
        Automatically decides whether to keep sequences in memory based on:
        - File size and configured threshold
        - User-specified memory mode
        """
        try:
            # Store file path for later access
            self.file_path = file_path
            
            # Get file size
            file_size_bytes = os.path.getsize(file_path)
            self.file_size_mb = file_size_bytes / (1024 * 1024)
            
            # Determine memory mode
            if self.config.memory_mode == "auto":
                if self.file_size_mb <= self.config.memory_threshold_mb:
                    self.memory_mode = "full_memory"
                    print(f"\nFile size ({self.file_size_mb:.1f} MB) is below threshold ({self.config.memory_threshold_mb} MB).")
                    print("Using full memory mode for optimal performance.")
                else:
                    self.memory_mode = "low_memory"
                    print(f"\nFile size ({self.file_size_mb:.1f} MB) exceeds threshold ({self.config.memory_threshold_mb} MB).")
                    print("Using low memory mode. Sequences will be loaded on demand.")
            else:
                self.memory_mode = self.config.memory_mode
                print(f"\nUsing {self.memory_mode} mode as configured.")
            
            # First pass: validation and counting
            print(f"\nLoading sequences from {file_path}...")
            sequence_count = sum(1 for _ in SeqIO.parse(file_path, "fasta"))
            
            if sequence_count == 0:
                print(f"\nNo sequences found in {file_path}")
                return False
            
            print(f"Total sequences detected: {sequence_count}")
            
            # Clear any existing data
            self.sequences = []
            self.gene_names = []
            self.sequence_index = {}
            
            # Load sequences based on memory mode
            if self.memory_mode == "full_memory":
                # Original behavior - keep all sequences in memory
                print("Loading sequences into memory...", end=" ", flush=True)
                progress_interval = max(1, sequence_count // 10)
                
                # Track duplicate IDs
                id_count = {}
                
                for i, record in enumerate(SeqIO.parse(file_path, "fasta"), 1):
                    base_id = record.id
                    
                    # Handle duplicate IDs
                    if base_id in id_count:
                        id_count[base_id] += 1
                        record.id = f"{base_id}_{id_count[base_id]}"
                    else:
                        id_count[base_id] = 0
                    
                    self.sequences.append(record)
                    self.gene_names.append(record.id)
                    
                    if i % progress_interval == 0:
                        print(f"{i/sequence_count*100:.0f}%", end=" ", flush=True)
                
                print(" 100% Done!")
            else:
                # Low memory mode - only extract gene names
                print("Extracting gene names...", end=" ", flush=True)
                progress_interval = max(1, sequence_count // 10)
                
                # Track duplicate IDs to match index building
                id_count = {}
                
                for i, record in enumerate(SeqIO.parse(file_path, "fasta"), 1):
                    base_id = record.id
                    
                    # Handle duplicate IDs the same way as index building
                    if base_id in id_count:
                        id_count[base_id] += 1
                        gene_id = f"{base_id}_{id_count[base_id]}"
                    else:
                        id_count[base_id] = 0
                        gene_id = base_id
                    
                    self.gene_names.append(gene_id)
                    
                    if i % progress_interval == 0:
                        print(f"{i/sequence_count*100:.0f}%", end=" ", flush=True)
                
                print(" 100% Done!")
                
                # Build sequence index AFTER extracting gene names
                if self.config.enable_sequence_indexing:
                    print("Building sequence index for fast access...")
                    self._build_sequence_index()
            
            # Process sequences to generate codon counts
            self._process_sequences()
            
            # Display comprehensive loading summary
            self._display_loading_summary(file_path)
            
            return True
            
        except Exception as e:
            print(f"\nError loading FASTA file: {e}")
            return False
    
    def _display_loading_summary(self, file_path):
        """Display comprehensive summary statistics after loading sequences."""
        print("\n" + "="*60)
        print("LOADING SUMMARY")
        print("="*60)
        
        # Basic file info
        print(f"File: {os.path.basename(file_path)}")
        print(f"File size: {self.file_size_mb:.1f} MB")
        print(f"Memory mode: {self.memory_mode}")
        if self.memory_mode == "low_memory":
            print("  Note: Sequences stored on disk, loaded on demand")
        else:
            print("  Note: All sequences loaded into memory")
        
        # Sequence statistics
        print(f"\nSequence Statistics:")
        print(f"  Total sequences: {len(self.gene_names):,}")
        
        # Codon statistics
        total_codons = int(self.codon_counts.sum().sum())
        print(f"\nCodon Statistics:")
        print(f"  Total codons counted: {total_codons:,}")
        if len(self.gene_names) > 0:
            avg_codons = total_codons / len(self.gene_names)
            print(f"  Average codons per sequence: {avg_codons:.1f}")
            print(f"  Average sequence length: {avg_codons * 3:.1f} bp")
        
        # Base composition summary
        if not self.base_composition.empty:
            gc_values = self.base_composition['GC']
            print(f"\nBase Composition Summary:")
            print(f"  GC content range: {gc_values.min():.1f}% - {gc_values.max():.1f}%")
            print(f"  Average GC content: {gc_values.mean():.1f}%")
            # Note: GC3s will be shown if already calculated, but not triggered
            if 'GC3s' in self.base_composition.columns:
                gc3s_values = self.base_composition['GC3s']
                print(f"  GC3s range: {gc3s_values.min():.1f}% - {gc3s_values.max():.1f}%")
                print(f"  Average GC3s: {gc3s_values.mean():.1f}%")
        
        # Check for potential issues
        issues = []
        empty_seqs = (self.codon_counts.sum(axis=1) == 0).sum()
        if empty_seqs > 0:
            issues.append(f"{empty_seqs} sequences with no valid codons")
        
        # Count sequences with very short lengths
        if len(self.gene_names) > 0:
            codon_sums = self.codon_counts.sum(axis=1)
            short_seqs = (codon_sums < 10).sum()
            if short_seqs > 0:
                issues.append(f"{short_seqs} sequences with < 10 codons")
        
        if issues:
            print(f"\nPotential Issues:")
            for issue in issues:
                print(f"  ⚠️  {issue}")
        else:
            print(f"\n✓ All sequences processed successfully")
        
        # Duplicate ID handling info
        if self.memory_mode == "low_memory" and hasattr(self, 'sequence_index'):
            unique_ids = len(set(gene.split('_')[0] for gene in self.gene_names))
            if unique_ids < len(self.gene_names):
                print(f"\nDuplicate ID Handling:")
                print(f"  Original unique IDs: {unique_ids:,}")
                print(f"  Total sequences (with suffixes): {len(self.gene_names):,}")
        
        print("="*60)

    def _build_sequence_index(self):
        """Build an index mapping gene IDs to file positions for fast random access."""
        try:
            # Count sequences while building index
            sequences_found = 0
            id_count = {}  # Track how many times each ID appears
            
            with open(self.file_path, 'r') as f:
                current_id = None
                header_pos = 0
                
                while True:
                    current_pos = f.tell()
                    line = f.readline()
                    
                    if not line:
                        # End of file - save last sequence
                        if current_id is not None:
                            self.sequence_index[current_id]['end'] = current_pos
                        break
                    
                    if line.startswith('>'):
                        # Save previous sequence end position
                        if current_id is not None:
                            self.sequence_index[current_id]['end'] = current_pos
                        
                        # Extract gene ID (first word after >)
                        # Use same method as BioPython for consistency
                        header_text = line[1:].strip()
                        if header_text:
                            base_id = header_text.split()[0]
                            
                            # Handle duplicate IDs
                            if base_id in id_count:
                                id_count[base_id] += 1
                                current_id = f"{base_id}_{id_count[base_id]}"
                            else:
                                id_count[base_id] = 0
                                current_id = base_id
                        else:
                            continue
                        sequences_found += 1
                        
                        # Store position after the header
                        header_pos = f.tell()
                        self.sequence_index[current_id] = {
                            'start': header_pos,
                            'header': line.strip()
                        }
            
            print(f"Built index for {len(self.sequence_index)} sequences (found {sequences_found} headers).")
            
            # Verify index completeness
            if len(self.sequence_index) != len(self.gene_names):
                print(f"Warning: Index has {len(self.sequence_index)} entries but {len(self.gene_names)} gene names.")
                # Print sample of missing genes for debugging
                missing = [g for g in self.gene_names if g not in self.sequence_index]
                if missing:
                    print(f"Sample of missing genes: {missing[:5]}")
            
        except Exception as e:
            print(f"Error building sequence index: {e}")
            import traceback
            traceback.print_exc()
            self.sequence_index = {}
    
    def _load_sequence_from_file(self, gene_id):
        """Load a specific sequence from file using the index."""
        if not self.file_path:
            return None
            
        # If we have an index, use it for fast access
        if gene_id in self.sequence_index:
            try:
                with open(self.file_path, 'r') as f:
                    # Seek to the start position
                    f.seek(self.sequence_index[gene_id]['start'])
                    
                    # Read until the end position
                    seq_data = ""
                    end_pos = self.sequence_index[gene_id]['end']
                    
                    while f.tell() < end_pos:
                        line = f.readline().strip()
                        if line and not line.startswith('>'):
                            seq_data += line
                    
                    # Create a SeqRecord object
                    from Bio.SeqRecord import SeqRecord
                    seq_record = SeqRecord(
                        Seq(seq_data),
                        id=gene_id,
                        description=self.sequence_index[gene_id]['header']
                    )
                    return seq_record
                    
            except Exception as e:
                print(f"Error loading sequence {gene_id} from file: {e}")
                return None
        else:
            # Fall back to parsing the entire file (slower)
            try:
                for record in SeqIO.parse(self.file_path, "fasta"):
                    if record.id == gene_id:
                        return record
                return None
            except Exception as e:
                print(f"Error loading sequence {gene_id} from file: {e}")
                return None
    
    def get_sequence(self, gene_id):
        """Get a sequence by ID, from memory or disk depending on mode."""
        if self.memory_mode == "full_memory":
            # Find sequence in memory
            for seq in self.sequences:
                if seq.id == gene_id:
                    return seq
            return None
        else:
            # Load from file
            return self._load_sequence_from_file(gene_id)
    
    def _save_checkpoint(self):
        """Save current processing state to checkpoint file."""
        if not self.config.enable_progress_saving or not self.checkpoint_file:
            return
            
        checkpoint_data = {
            'processed_genes': list(self.processed_genes),
            'codon_counts': self.codon_counts.to_dict(),
            'base_composition': self.base_composition.to_dict(),
            'timestamp': datetime.now().isoformat(),
            'total_genes': len(self.gene_names),
            'genes_processed': len(self.processed_genes)
        }
        
        try:
            with open(self.checkpoint_file, 'wb') as f:
                pickle.dump(checkpoint_data, f)
            if self.config.detailed_output >= 2:
                print(f"\nCheckpoint saved ({len(self.processed_genes)} genes processed)")
        except Exception as e:
            print(f"\nWarning: Failed to save checkpoint: {e}")
    
    def _load_checkpoint(self):
        """Load processing state from checkpoint file."""
        if not self.config.enable_progress_saving or not self.checkpoint_file:
            return False
            
        if not os.path.exists(self.checkpoint_file):
            return False
            
        try:
            with open(self.checkpoint_file, 'rb') as f:
                checkpoint_data = pickle.load(f)
            
            self.processed_genes = set(checkpoint_data['processed_genes'])
            self.codon_counts = pd.DataFrame(checkpoint_data['codon_counts'])
            self.base_composition = pd.DataFrame(checkpoint_data['base_composition'])
            
            print(f"\nResuming from checkpoint: {len(self.processed_genes)} genes already processed")
            print(f"Checkpoint timestamp: {checkpoint_data['timestamp']}")
            return True
            
        except Exception as e:
            print(f"\nWarning: Failed to load checkpoint: {e}")
            return False
    
    def _clear_checkpoint(self):
        """Remove checkpoint file after successful completion."""
        if self.checkpoint_file and os.path.exists(self.checkpoint_file):
            try:
                os.remove(self.checkpoint_file)
                if self.config.detailed_output >= 2:
                    print("\nCheckpoint file removed (processing complete)")
            except:
                pass

    def _process_sequences(self):
        """Process all loaded sequences to generate codon counts and other statistics."""
        # Set up checkpoint file if enabled
        if self.config.enable_progress_saving:
            base_name = os.path.splitext(os.path.basename(self.file_path))[0]
            self.checkpoint_file = f"{base_name}_checkpoint.pkl"
        
        # Try to load from checkpoint
        resumed = self._load_checkpoint()
        
        if not resumed:
            # Initialize DataFrames with explicit data types
            self.codon_counts = pd.DataFrame(
                0, index=self.gene_names, columns=CODONS, dtype=float)
            self.base_composition = pd.DataFrame(
                0,
                index=self.gene_names,
                columns=[
                    'Length',
                    'AA_count',
                    'GC',
                    'GC1',
                    'GC2',
                    'GC3',
                    'GC3s',
                    'C',
                    'T',
                    'A',
                    'G',
                    'C3',
                    'T3',
                    'A3',
                    'G3'],
                dtype=float)
            self.processed_genes = set()

        # Get unprocessed genes
        unprocessed_genes = [g for g in self.gene_names if g not in self.processed_genes]
        
        if not unprocessed_genes:
            print("All sequences already processed!")
            return
            
        print(f"Processing {len(unprocessed_genes)} sequences...")
        if resumed:
            print(f"(Resuming from {len(self.processed_genes)} already processed)")
            
        # Decide whether to use parallel processing
        use_parallel = (
            self.config.enable_parallel_processing and 
            len(unprocessed_genes) > self.config.batch_size * 2 and
            self.config.num_processes > 1
        )
        
        if use_parallel:
            self._process_sequences_parallel(unprocessed_genes)
        else:
            self._process_sequences_sequential(unprocessed_genes)
        
        # Clear checkpoint after successful completion
        self._clear_checkpoint()
        
        # Check total codon counts
        total_codons = self.codon_counts.sum().sum()
        print(f"\nTotal codons counted across all sequences: {int(total_codons):,}")
        
        # Update analysis status flags
        self.analysis_status['basic_metrics'] = True
        self.analysis_status['codon_usage'] = True
        self.analysis_status['amino_acid_usage'] = True
        self.analysis_status['base_composition'] = True
        
        # Post-processing calculations are now done on-demand to save time
        # RSCU, amino acid usage, and GC3s will be calculated when first accessed

    def _process_sequences_parallel(self, unprocessed_genes):
        """Process sequences in parallel using multiprocessing."""
        print(f"Using parallel processing with {self.config.num_processes} processes...")
        print(f"Total genes to process: {len(unprocessed_genes)}")
        print(f"Batch size: {self.config.batch_size}")
        print("Progress: ", end="", flush=True)
        
        total_genes = len(unprocessed_genes)
        genes_processed = 0
        last_percent = -1
        
        # Create batches
        batches = []
        for i in range(0, total_genes, self.config.batch_size):
            batch_genes = unprocessed_genes[i:i + self.config.batch_size]
            
            if self.memory_mode == "full_memory":
                # Get sequences from memory
                batch_seqs = []
                for gene_id in batch_genes:
                    for seq in self.sequences:
                        if seq.id == gene_id:
                            batch_seqs.append(str(seq.seq))
                            break
                    else:
                        batch_seqs.append("")  # Empty if not found
                
                batches.append((batch_seqs, self.config.genetic_code))
            else:
                # Pass file info for workers to load sequences themselves
                # Create a subset of the index for this batch
                batch_index = {gene: self.sequence_index[gene] 
                              for gene in batch_genes if gene in self.sequence_index}
                
                batches.append(((batch_genes, self.file_path, batch_index), self.config.genetic_code))
        
        print(f"\nCreated {len(batches)} batches for parallel processing")
        print(f"Launching {self.config.num_processes} worker processes...")
        
        # Process batches in parallel
        with ProcessPoolExecutor(max_workers=self.config.num_processes) as executor:
            # Submit all batches at once
            futures = [executor.submit(process_sequence_batch, batch) for batch in batches]
            print(f"Submitted {len(futures)} batch jobs to process pool")
            
            for future in as_completed(futures):
                try:
                    codon_counts, base_composition = future.result()
                    
                    # Update DataFrames with results
                    for gene_name, counts in codon_counts.items():
                        for codon, count in counts.items():
                            self.codon_counts.at[gene_name, codon] = count
                        self.processed_genes.add(gene_name)
                    
                    for gene_name, base_comp in base_composition.items():
                        for col, value in base_comp.items():
                            if col in self.base_composition.columns:
                                self.base_composition.at[gene_name, col] = value
                    
                    genes_processed += len(codon_counts)
                    
                    # Update progress
                    current_percent = int((genes_processed / total_genes) * 100)
                    if current_percent > last_percent:
                        if current_percent % 10 == 0:
                            print(f"{current_percent}%", end="", flush=True)
                        else:
                            print(".", end="", flush=True)
                        last_percent = current_percent
                    
                    # Save checkpoint if needed
                    if genes_processed % self.config.checkpoint_interval == 0:
                        self._save_checkpoint()
                        
                except Exception as e:
                    print(f"\nError processing batch: {e}")
        
        print(" 100%")
        
        # Report on parallel processing performance
        if self.config.detailed_output >= 2:
            print(f"\nParallel processing summary:")
            print(f"  Batches processed: {len(batches)}")
            print(f"  Sequences per batch: {self.config.batch_size}")
            print(f"  Worker processes used: {self.config.num_processes}")
    
    def _process_sequences_sequential(self, unprocessed_genes):
        """Process sequences sequentially (original method)."""
        print("Progress: ", end="", flush=True)
        
        self.last_percent = -1
        total_seqs = len(unprocessed_genes)
        
        for i, seq_id in enumerate(unprocessed_genes):
            if self.memory_mode == "full_memory":
                # Find sequence in memory
                seq_record = None
                for seq in self.sequences:
                    if seq.id == seq_id:
                        seq_record = seq
                        break
                if seq_record:
                    dna_seq = str(seq_record.seq).upper()
                    self._process_single_sequence(i, seq_id, dna_seq, total_seqs)
            else:
                # Load sequence from file
                seq_record = self._load_sequence_from_file(seq_id)
                if seq_record:
                    dna_seq = str(seq_record.seq).upper()
                    self._process_single_sequence(i, seq_id, dna_seq, total_seqs)
                else:
                    print(f"\nWarning: Could not load sequence for {seq_id}")
            
            # Mark as processed
            self.processed_genes.add(seq_id)
            
            # Save checkpoint periodically
            if (i + 1) % self.config.checkpoint_interval == 0:
                self._save_checkpoint()
        
        print(" 100%")
    
    def _process_single_sequence(self, i, seq_id, dna_seq, total_seqs):
        """Process a single sequence to extract codon counts and base composition."""
        # Update progress indicator
        current_percent = int((i / total_seqs) * 100)
        if current_percent > self.last_percent:
            if current_percent % 10 == 0 and current_percent != 0:
                print(f"{current_percent}%", end="", flush=True)
            else:
                print(".", end="", flush=True)
            self.last_percent = current_percent
        
        # Print details for the first few sequences
        if i < 3 and self.config.detailed_output >= 2:
            print(f"\nSequence {i+1}: ID={seq_id}, Length={len(dna_seq)}")
            print(f"First 30 bases: {dna_seq[:30]}")

        # Ensure length is divisible by 3
        length = len(dna_seq) - (len(dna_seq) % 3)
        self.base_composition.at[seq_id, 'Length'] = length

        # Initialize base counts
        base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        base_counts_1st = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        base_counts_2nd = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        base_counts_3rd = {'A': 0, 'T': 0, 'G': 0, 'C': 0}

        # Count codons and bases
        codon_count = 0
        for j in range(0, length, 3):
            if j + 3 <= length:  # Make sure we don't go beyond sequence length
                codon_dna = dna_seq[j:j + 3]

                # Skip if codon contains invalid characters
                if 'N' in codon_dna or len(codon_dna) < 3:
                    continue

                # Update base counts
                for pos, base in enumerate(codon_dna):
                    if base in 'ATGC':
                        base_counts[base] += 1
                        if pos == 0:
                            base_counts_1st[base] += 1
                        elif pos == 1:
                            base_counts_2nd[base] += 1
                        elif pos == 2:
                            base_counts_3rd[base] += 1

                # Convert DNA to RNA codon
                codon_rna = codon_dna.replace('T', 'U')

                # Check if this is a valid codon
                if codon_rna in CODONS:
                    # Increase codon count
                    self.codon_counts.at[seq_id, codon_rna] += 1
                    codon_count += 1

        # Store base composition statistics
        self.base_composition.at[seq_id, 'AA_count'] = codon_count

        # Store base counts
        for base in 'ATGC':
            self.base_composition.at[seq_id, base] = base_counts[base]

        # Store 3rd position base counts
        for base in 'ATGC':
            self.base_composition.at[seq_id,
                                     f"{base}3"] = base_counts_3rd[base]

        # Calculate GC content for all positions
        total_bases = sum(base_counts.values())
        if total_bases > 0:
            gc_content = (
                base_counts['G'] + base_counts['C']) / total_bases * 100
            self.base_composition.at[seq_id, 'GC'] = gc_content

        # Calculate GC content for 1st positions
        total_1st = sum(base_counts_1st.values())
        if total_1st > 0:
            gc1_content = (
                base_counts_1st['G'] + base_counts_1st['C']) / total_1st * 100
            self.base_composition.at[seq_id, 'GC1'] = gc1_content

        # Calculate GC content for 2nd positions
        total_2nd = sum(base_counts_2nd.values())
        if total_2nd > 0:
            gc2_content = (
                base_counts_2nd['G'] + base_counts_2nd['C']) / total_2nd * 100
            self.base_composition.at[seq_id, 'GC2'] = gc2_content

        # Calculate GC content for 3rd positions
        total_3rd = sum(base_counts_3rd.values())
        if total_3rd > 0:
            gc3_content = (
                base_counts_3rd['G'] + base_counts_3rd['C']) / total_3rd * 100
            self.base_composition.at[seq_id, 'GC3'] = gc3_content

    def _calculate_rscu(self):
        """Calculate RSCU (Relative Synonymous Codon Usage) values."""
        print("  Calculating RSCU values...", end=" ", flush=True)
        
        aa_map = self.config.get_aa_map()
        # Get the correct AA_TO_CODONS mapping for current genetic code
        aa_to_codons = AA_TO_CODONS_MAP.get(
            self.config.genetic_code, AA_TO_CODONS_MAP[1])

        self.rscu_values = pd.DataFrame(
            0.0, index=self.gene_names, columns=CODONS)
        
        # Use vectorized operations for better performance
        # First, calculate amino acid totals for each gene
        aa_totals = {}
        for aa, codons in aa_to_codons.items():
            if aa != 'STOP':
                aa_totals[aa] = self.codon_counts[codons].sum(axis=1)
        
        # Calculate RSCU for each codon
        for codon in CODONS:
            aa = aa_map[codon]
            if aa != 'STOP':
                synonymous_codons = len(aa_to_codons[aa])
                # Vectorized calculation for all genes at once
                mask = aa_totals[aa] > 0
                self.rscu_values.loc[mask, codon] = (
                    self.codon_counts.loc[mask, codon] * synonymous_codons / aa_totals[aa][mask]
                )
        
        print("Done!")

    def _calculate_aa_usage(self):
        """Calculate amino acid usage."""
        print("  Calculating amino acid usage...", end=" ", flush=True)
        
        aa_map = self.config.get_aa_map()
        # Get unique amino acids excluding stop codons
        unique_aas = sorted(
            set([aa for aa in aa_map.values() if aa != 'STOP']) | {'STOP'})

        self.amino_acid_usage = pd.DataFrame(
            0, index=self.gene_names, columns=unique_aas)
        
        # Use vectorized operations - sum codons for each amino acid
        for aa in unique_aas:
            # Get all codons that code for this amino acid
            codons_for_aa = [codon for codon in CODONS if aa_map[codon] == aa]
            if codons_for_aa:
                # Sum all codons for this amino acid using vectorized operation
                self.amino_acid_usage[aa] = self.codon_counts[codons_for_aa].sum(axis=1)
        
        print("Done!")

    def get_rscu_values(self):
        """Get RSCU values, calculating them if needed."""
        if self.rscu_values.empty and not self.codon_counts.empty:
            print("\nCalculating RSCU values (first time)...")
            self._calculate_rscu()
        return self.rscu_values
    
    def get_amino_acid_usage(self):
        """Get amino acid usage, calculating it if needed."""
        if self.amino_acid_usage.empty and not self.codon_counts.empty:
            print("\nCalculating amino acid usage (first time)...")
            self._calculate_aa_usage()
        return self.amino_acid_usage

    def calculate_enc(self):
        """Calculate ENC (Effective Number of Codons)."""
        # Implementation of Wright's method
        if self.enc_values is not None:
            return self.enc_values

        enc_values = {}
        print(f"\nCalculating ENC values for {len(self.gene_names)} genes...")
        print("Progress: ", end="", flush=True)

        last_percent = -1
        total_genes = len(self.gene_names)

        for i, gene in enumerate(self.gene_names):
            # Update progress indicator
            current_percent = int((i / total_genes) * 100)
            if current_percent > last_percent:
                if current_percent % 10 == 0 and current_percent != 0:
                    print(f"{current_percent}%", end="", flush=True)
                else:
                    print(".", end="", flush=True)
                last_percent = current_percent

        # Get the current genetic code mappings
        aa_map = self.config.get_aa_map()
        aa_to_codons = AA_TO_CODONS_MAP.get(
            self.config.genetic_code, AA_TO_CODONS_MAP[1])

        for gene in self.gene_names:
            # Group by degeneracy class (2-fold, 3-fold, 4-fold, 6-fold)
            by_degeneracy = {2: [], 3: [], 4: [], 6: []}

            # Determine degeneracy class for each AA
            for aa, codons in aa_to_codons.items():
                if aa == 'STOP':
                    continue  # Skip stop codons

                degeneracy = len(codons)
                if degeneracy > 1:  # Only consider degenerate AAs
                    # Calculate frequencies within each AA
                    total_aa = sum(
                        self.codon_counts.at[gene, codon] for codon in codons)
                    if total_aa > 0:
                        freqs = [self.codon_counts.at[gene, codon] /
                                 total_aa for codon in codons]
                        # Calculate homozygosity
                        homozygosity = sum(f * f for f in freqs)

                        # Store by degeneracy class
                        if degeneracy == 2:
                            by_degeneracy[2].append(homozygosity)
                        elif degeneracy == 3:
                            by_degeneracy[3].append(homozygosity)
                        elif degeneracy == 4:
                            by_degeneracy[4].append(homozygosity)
                        elif degeneracy == 6:
                            by_degeneracy[6].append(homozygosity)

            # Calculate average homozygosity for each class
            F_values = {}
            for degeneracy, values in by_degeneracy.items():
                if values:
                    F_values[degeneracy] = sum(values) / len(values)
                else:
                    F_values[degeneracy] = 0

            # Calculate ENC
            enc = 2  # Start with non-degenerate AAs (Met, Trp)

            # Add contribution from 2-fold degenerate AAs
            if F_values[2] > 0:
                enc += 9 / F_values[2]

            # Add contribution from 3-fold degenerate AAs (Ile)
            if F_values[3] > 0:
                enc += 1 / F_values[3]

            # Add contribution from 4-fold degenerate AAs
            if F_values[4] > 0:
                enc += 5 / F_values[4]

            # Add contribution from 6-fold degenerate AAs (Arg, Leu)
            if F_values[6] > 0:
                enc += 3 / F_values[6]

            enc_values[gene] = min(enc, 61)  # Cap at 61

        print(" Done!")
        self.enc_values = enc_values
        return enc_values

    def calculate_scuo(self):
        """Calculate SCUO (Synonymous Codon Usage Order)."""
        # Implementation of SCUO based on information theory
        scuo_values = {}
        aa_map = self.config.get_aa_map()
        aa_to_codons = AA_TO_CODONS_MAP.get(
            self.config.genetic_code, AA_TO_CODONS_MAP[1])

        for gene in self.gene_names:
            # Calculate SCUO for each amino acid
            aa_scuo = {}
            for aa, codons in aa_to_codons.items():
                if aa == 'STOP' or len(
                        codons) <= 1:  # Skip non-degenerate amino acids and stop codons
                    continue

                # Calculate total usage of this amino acid
                total_aa = sum(
                    self.codon_counts.at[gene, codon] for codon in codons)
                if total_aa == 0:
                    continue

                # Calculate normalized frequencies
                freqs = [self.codon_counts.at[gene, codon] /
                         total_aa for codon in codons]

                # Calculate entropy
                entropy = -sum(f * math.log(f) if f > 0 else 0 for f in freqs)

                # Calculate maximum entropy for this amino acid
                max_entropy = math.log(len(codons))

                # Calculate normalized complexity (SCUO)
                if max_entropy > 0:
                    aa_scuo[aa] = 1.0 - (entropy / max_entropy)
                else:
                    aa_scuo[aa] = 0.0

            # Calculate overall SCUO as weighted average
            if aa_scuo:
                total_weight = sum(sum(
                    self.codon_counts.at[gene, codon] for codon in aa_to_codons[aa]) for aa in aa_scuo)
                if total_weight > 0:
                    weighted_scuo = sum(
                        aa_scuo[aa] * sum(self.codon_counts.at[gene, codon]
                                          for codon in aa_to_codons[aa])
                        for aa in aa_scuo
                    ) / total_weight
                    scuo_values[gene] = weighted_scuo
                else:
                    scuo_values[gene] = 0.0
            else:
                scuo_values[gene] = 0.0

        return scuo_values

    def calculate_optimal_codons(self,reference_genes=None,method="frequency",percentage=10):
        """
        Identify optimal codons using specified method.

        Parameters:
        reference_genes (list): List of gene names to use as reference (for frequency and multivariate)
        method (str): Method to use for selecting optimal codons:
           "frequency" - most frequent codon (default)
           "multivariate" - multivariate analysis
           "rscu" - highest RSCU value
           "raw_count" - most common codon
        percentage (float): Percentage of genes to use for multivariate analysis

        Returns:
        dict: Dictionary mapping amino acids to their optimal codons
        """
        aa_map = self.config.get_aa_map()
        optimal_codons = {}
        aa_to_codons = AA_TO_CODONS_MAP.get(
            self.config.genetic_code, AA_TO_CODONS_MAP[1])

        print(f"Identifying optimal codons using {method} method...")

        # Get reference genes based on method
        if method == "multivariate":
            reference_genes = self.identify_reference_genes_multivariate(percentage)
            if not reference_genes:
                print("Multivariate analysis failed. Using all genes.")
                reference_genes = self.gene_names
        elif reference_genes is None:
            reference_genes = self.gene_names

        # Calculate reference data
        ref_counts = self.codon_counts.loc[reference_genes].sum()

        if method == "rscu":
            # Calculate mean RSCU across reference genes
            rscu_values = self.get_rscu_values()
            mean_rscu = rscu_values.loc[reference_genes].mean()

        # Determine optimal codon for each amino acid
        for aa in set(aa_map.values()):
            if aa == 'STOP':
                continue  # Skip stop codons

            # Get all codons for this amino acid based on current genetic code
            aa_codons = aa_to_codons.get(aa, [])
            if not aa_codons:
                continue

            # Select optimal codon based on method
            if method == "rscu":
                aa_rscu = {codon: mean_rscu[codon] for codon in aa_codons}
                optimal_codon = max(aa_rscu, key=aa_rscu.get)
            else:
                optimal_codon = max(aa_codons, key=lambda c: ref_counts[c])

            optimal_codons[aa] = optimal_codon

        # Store the optimal codons
        self.optimal_codons = optimal_codons
        return optimal_codons

    def calculate_fop(self, reference_genes=None):
        """
        Calculate Frequency of Optimal Codons (Fop) for each gene.
        Fop is the ratio of optimal codons to synonymous codons in a gene.

        Parameters:
        reference_genes (list): List of gene names to use as reference. If None, all genes are used.

        Returns:
        dict: Dictionary mapping gene names to Fop values
        """
        # Identify optimal codons if not already done
        if self.optimal_codons is None:
            self.calculate_optimal_codons(reference_genes)

        fop_values = {}
        aa_map = self.config.get_aa_map()
        aa_to_codons = AA_TO_CODONS_MAP.get(
            self.config.genetic_code, AA_TO_CODONS_MAP[1])

        for gene in self.gene_names:
            optimal_codon_count = 0
            synonymous_codon_count = 0

            # Count optimal and synonymous codons
            for codon in CODONS:
                aa = aa_map[codon]
                if aa == 'STOP':
                    continue  # Skip stop codons

                # Get count of this codon in the gene
                codon_count = self.codon_counts.at[gene, codon]

                # Only count if this amino acid has synonymous codons
                if len(aa_to_codons[aa]) > 1:
                    synonymous_codon_count += codon_count

                    # Check if this is an optimal codon
                    if aa in self.optimal_codons and self.optimal_codons[aa] == codon:
                        optimal_codon_count += codon_count

            # Calculate Fop
            if synonymous_codon_count > 0:
                fop_values[gene] = optimal_codon_count / synonymous_codon_count
            else:
                fop_values[gene] = 0.0

        # Store the values
        self.fop_values = fop_values
        return fop_values

    def calculate_cai(self, reference_genes=None):
        """
        Calculate Codon Adaptation Index (CAI) for each gene.
        CAI measures how well a gene is adapted to the codon usage bias of highly expressed genes.

        Parameters:
        reference_genes (list): List of gene names to use as reference. If None, all genes are used.

        Returns:
        dict: Dictionary mapping gene names to CAI values
        """
        # Generate codon weights based on reference genes
        aa_map = self.config.get_aa_map()
        aa_to_codons = AA_TO_CODONS_MAP.get(
            self.config.genetic_code, AA_TO_CODONS_MAP[1])

        # If no reference genes provided, use all genes
        if reference_genes is None:
            reference_genes = self.gene_names

        # Calculate reference codon counts
        ref_counts = self.codon_counts.loc[reference_genes].sum()

        # Calculate weights for each codon
        codon_weights = {}
        for aa in set(aa_map.values()):
            if aa == 'STOP':
                continue  # Skip stop codons

            # Get all codons for this amino acid
            aa_codons = aa_to_codons.get(aa, [])
            if not aa_codons:
                continue

            # Find the most used codon for this amino acid in reference genes
            max_codon = max(aa_codons, key=lambda c: ref_counts[c])
            max_count = ref_counts[max_codon]

            # Calculate relative adaptiveness for each codon
            for codon in aa_codons:
                count = ref_counts[codon]
                if count > 0 and max_count > 0:
                    codon_weights[codon] = count / max_count
                else:
                    codon_weights[codon] = 0.001  # Small non-zero value

        # Calculate CAI for each gene
        cai_values = {}
        for gene in self.gene_names:
            weighted_sum = 0
            codon_count = 0

            for codon in CODONS:
                aa = aa_map[codon]
                if aa == 'STOP':
                    continue  # Skip stop codons

                # Get count of this codon in the gene
                count = self.codon_counts.at[gene, codon]
                if count > 0:
                    if codon in codon_weights and codon_weights[codon] > 0:
                        weighted_sum += count * math.log(codon_weights[codon])
                        codon_count += count

            # Calculate CAI
            if codon_count > 0:
                cai = math.exp(weighted_sum / codon_count)
                cai_values[gene] = cai
            else:
                cai_values[gene] = 0.0

        # Store the values
        self.cai_values = cai_values
        return cai_values

    def perform_multivariate_analysis(self, analysis_type='CA', data_type='RSCU'):
        """Perform multivariate analysis on the data with progress indicator."""
        if not self.gene_names:
            print("No data loaded. Please load a FASTA file first.")
            return None

        print(f"\nPerforming {analysis_type} on {data_type} data...")
        
        # Initialize preprocessing parameters tracking
        self.last_preprocessing_params = {
            'analysis_type': analysis_type,
            'data_type': data_type
        }
        
        # Check if we can use multiple cores
        use_parallel = self.config.enable_parallel_processing and self.config.num_processes > 1
        if use_parallel:
            print(f"Using parallel processing with up to {self.config.num_processes} cores")

        # Prepare data matrix
        print("Preparing data matrix...", end="", flush=True)
        if data_type == 'RSCU':
            rscu_values = self.get_rscu_values()
            data_matrix = rscu_values.copy()
            # Filter out stop codons
            aa_map = self.config.get_aa_map()
            stop_codons = [codon for codon in CODONS if aa_map[codon] == 'STOP']
            for codon in stop_codons:
                if codon in data_matrix.columns:
                    data_matrix = data_matrix.drop(codon, axis=1)
            
            # Also filter out single-codon amino acids (e.g., Met, Trp)
            # These have no synonymous alternatives and don't inform about codon usage bias
            aa_to_codons = AA_TO_CODONS_MAP.get(self.config.genetic_code, AA_TO_CODONS_MAP[1])
            single_codon_aas = []
            single_codons_removed = []
            
            for aa, codons in aa_to_codons.items():
                if len(codons) == 1 and aa != 'STOP':
                    single_codon_aas.append(aa)
                    codon = codons[0]
                    if codon in data_matrix.columns:
                        data_matrix = data_matrix.drop(codon, axis=1)
                        single_codons_removed.append(f"{codon} ({aa})")
            
            if single_codons_removed:
                print(f" (excluded {len(single_codons_removed)} single-codon amino acids: {', '.join(single_codons_removed)})", end="")
        elif data_type == 'AA':
            amino_acid_usage = self.get_amino_acid_usage()
            data_matrix = amino_acid_usage.copy()
            if 'STOP' in data_matrix.columns:
                data_matrix = data_matrix.drop('STOP', axis=1)
                
            # Convert to frequencies for PCA to avoid length effects
            if analysis_type == 'PCA':
                row_sums = data_matrix.sum(axis=1)
                # Avoid division by zero for genes with no amino acids
                row_sums[row_sums == 0] = 1
                data_matrix = data_matrix.div(row_sums, axis=0)
                print(" (converted to frequencies)", end="")
                
                # Ask user about CLR transformation for compositional data analysis
                print("\n\nCompositional data preprocessing options:")
                print("1. Standard (frequencies only)")
                print("2. CLR (Centered Log-Ratio) transformation - recommended for compositional data")
                print("3. Log transformation with pseudocount")
                
                # Default to CLR for amino acid compositional data
                if data_type == 'AA':
                    default_choice = '2'
                    prompt_text = "\nSelect preprocessing method (1-3, default=2 CLR for amino acids): "
                    print("\nNote: CLR is the mathematically appropriate default for compositional data.")
                else:
                    default_choice = '1'
                    prompt_text = "\nSelect preprocessing method (1-3, default=1): "
                
                preprocess_choice = input(prompt_text).strip() or default_choice
                
                if preprocess_choice == '2':
                    print("\nApplying CLR transformation...")
                    # Add small pseudocount to avoid log(0)
                    pseudocount = 1e-6
                    data_matrix_pseudo = data_matrix + pseudocount
                    # Renormalize after adding pseudocount
                    data_matrix_pseudo = data_matrix_pseudo.div(data_matrix_pseudo.sum(axis=1), axis=0)
                    
                    # Calculate geometric mean for each row
                    log_data = np.log(data_matrix_pseudo)
                    geometric_means = log_data.mean(axis=1)
                    
                    # Apply CLR transformation
                    data_matrix = log_data.sub(geometric_means, axis=0)
                    print("CLR transformation applied successfully.")
                    self.last_preprocessing_params['transformation'] = 'CLR'
                    self.last_preprocessing_params['pseudocount'] = pseudocount
                    
                elif preprocess_choice == '3':
                    print("\nApplying log transformation...")
                    pseudocount = float(input("Enter pseudocount value (default=0.001): ").strip() or "0.001")
                    data_matrix = np.log(data_matrix + pseudocount)
                    print(f"Log transformation with pseudocount {pseudocount} applied.")
                    self.last_preprocessing_params['transformation'] = 'Log'
                    self.last_preprocessing_params['pseudocount'] = pseudocount
                else:
                    print("\nUsing standard frequency-based analysis.")
                    self.last_preprocessing_params['transformation'] = 'Frequency only'
        else:
            print(f"Unsupported data type: {data_type}")
            return None

        # Filter out columns with zero variance
        data_matrix = data_matrix.loc[:, data_matrix.var() > 0]
        print(" Done!")

        if data_matrix.empty or data_matrix.shape[0] < 2 or data_matrix.shape[1] < 2:
            print("Not enough data for multivariate analysis. Need at least 2 sequences and 2 variables with variance.")
            return None

        print(f"Analysis matrix dimensions: {data_matrix.shape[0]} genes × {data_matrix.shape[1]} variables")
        
        # Outlier detection option
        if data_matrix.shape[0] > 10000:
            print(f"\nNote: Large dataset ({data_matrix.shape[0]} genes). Outlier detection may take time.")
            print("For compositional data (amino acids), outliers are often less meaningful.")
        
        outlier_option = input("\nPerform outlier detection? (y/n, default=n): ").strip().lower()
        if outlier_option == 'y':
            print("\nDetecting outliers...")
            
            # Convert to numeric array for distance calculation
            data_array = data_matrix.values
            
            # For large datasets or compositional data, use simpler outlier detection
            if data_matrix.shape[0] > 50000 or data_type == 'AA':
                print("Using simplified outlier detection for large/compositional dataset...")
                # Use z-score based outlier detection on first few PCs
                from sklearn.decomposition import PCA as PCA_outlier
                pca_temp = PCA_outlier(n_components=min(5, data_matrix.shape[1]))
                pca_scores = pca_temp.fit_transform(stats.zscore(data_array, axis=0))
                
                # Calculate distance from origin in PC space
                distances = np.sqrt(np.sum(pca_scores**2, axis=1))
                threshold = np.mean(distances) + 3 * np.std(distances)
                outliers = distances > threshold
                mahal_dist = distances  # For consistency in reporting
                print(f"Using PCA-based distance (threshold: {threshold:.2f})")
            else:
                # Original Mahalanobis distance calculation
                print("Using Mahalanobis distance...")
                
                # Calculate mean and covariance
                mean = np.mean(data_array, axis=0)
                
                # Use robust covariance estimation for better outlier detection
                # Check if CLR transformation was applied
                use_robust = True
                if hasattr(self, 'last_preprocessing_params') and self.last_preprocessing_params.get('transformation') == 'CLR':
                    use_robust = False  # CLR creates singular covariance matrix
                    print("(Using standard covariance for CLR-transformed data)")
                
                try:
                    if use_robust and data_matrix.shape[0] < 10000:  # Only use robust for smaller datasets
                        # Try to use MinCovDet for robust covariance estimation
                        from sklearn.covariance import MinCovDet
                        # Suppress warnings
                        import warnings
                        with warnings.catch_warnings():
                            warnings.filterwarnings('ignore', category=RuntimeWarning)
                            robust_cov = MinCovDet(random_state=42, support_fraction=0.9).fit(data_array)
                            mahal_dist = robust_cov.mahalanobis(data_array)
                        print("Using robust covariance estimation")
                    else:
                        raise Exception("Use standard covariance")
                except:
                    # Fall back to standard covariance
                    cov = np.cov(data_array.T)
                    # Add regularization to avoid singular matrix
                    regularization = np.trace(cov) / data_matrix.shape[1] * 1e-4
                    cov += np.eye(cov.shape[0]) * regularization
                    
                    try:
                        inv_cov = np.linalg.inv(cov)
                        # Calculate Mahalanobis distance for each point
                        diff = data_array - mean
                        mahal_dist = np.sqrt(np.sum(diff @ inv_cov * diff, axis=1))
                        print("Using standard covariance estimation")
                    except:
                        # If inversion still fails, use SVD-based approach
                        print("Using SVD-based pseudo-inverse for singular covariance")
                        U, s, Vt = np.linalg.svd(cov)
                        # Keep only non-zero singular values
                        s_inv = np.zeros_like(s)
                        s_inv[s > 1e-10] = 1.0 / s[s > 1e-10]
                        inv_cov = Vt.T @ np.diag(s_inv) @ U.T
                        diff = data_array - mean
                        mahal_dist = np.sqrt(np.sum(diff @ inv_cov * diff, axis=1))
                
                # Identify outliers using chi-square distribution
                threshold = np.sqrt(chi2.ppf(0.975, data_matrix.shape[1]))  # 97.5% confidence
                outliers = mahal_dist > threshold
            
            n_outliers = np.sum(outliers)
            
            print(f"\nOutliers detected: {n_outliers} ({n_outliers/len(data_matrix)*100:.1f}%)")
            
            if n_outliers > 0:
                # Show some outlier information
                outlier_genes = data_matrix.index[outliers].tolist()
                print(f"Outlier genes: {', '.join(outlier_genes[:10])}", end="")
                if n_outliers > 10:
                    print(f" ... and {n_outliers-10} more")
                else:
                    print()
                
                remove_option = input("\nRemove outliers before analysis? (y/n, default=n): ").strip().lower()
                if remove_option == 'y':
                    data_matrix = data_matrix.loc[~outliers]
                    print(f"Removed {n_outliers} outliers. New matrix: {data_matrix.shape[0]} genes × {data_matrix.shape[1]} variables")
                    self.last_preprocessing_params['outliers_removed'] = n_outliers
                    self.last_preprocessing_params['outlier_method'] = 'Mahalanobis distance'
                else:
                    print("Keeping all genes including outliers.")
                    self.last_preprocessing_params['outliers_removed'] = 0
            else:
                print("No significant outliers detected.")
                self.last_preprocessing_params['outliers_removed'] = 0
        
        print("\nRunning analysis calculations...")
        print("Progress: ", end="", flush=True)

        if analysis_type == 'PCA':
            if not HAS_SKLEARN:
                print("\nPCA requires scikit-learn, which is not installed.")
                print("Please install it with: pip install scikit-learn")
                return None
            print("10%", end="", flush=True)
            
            # Enable parallel computation for z-score normalization
            scaled_data = stats.zscore(data_matrix, axis=0)
            print("...30%", end="", flush=True)

            # Use randomized SVD for large datasets, which can be parallelized
            n_samples, n_features = data_matrix.shape
            if n_samples > 5000 and n_features > 50:
                # Use randomized SVD for better performance on large datasets
                pca = PCA(n_components=min(50, n_features, n_samples - 1), 
                         svd_solver='randomized', 
                         random_state=42)
                print("...40% (using randomized SVD)", end="", flush=True)
            else:
                pca = PCA()
            print("...50%", end="", flush=True)

            pca_result = pca.fit_transform(scaled_data)
            print("...80%", end="", flush=True)

            n_components = min(5, data_matrix.shape[1], data_matrix.shape[0] - 1)

            result = {
                'coordinates': pd.DataFrame(
                    pca_result[:, :n_components],
                    index=data_matrix.index,
                    columns=[f'PC{i+1}' for i in range(n_components)]
                ),
                'explained_variance': pca.explained_variance_ratio_[:n_components],
                'loadings': pd.DataFrame(
                    pca.components_.T[:, :n_components],
                    index=data_matrix.columns,
                    columns=[f'PC{i+1}' for i in range(n_components)]
                ),
                'analysis_type': analysis_type,
                'data_type': data_type
            }
            print("...100% Done!")

        elif analysis_type == 'CA':
            data_array = data_matrix.values
            print("...20%", end="", flush=True)

            if data_array.size == 0 or data_array.shape[0] < 2 or data_array.shape[1] < 2:
                print("\nNot enough data for correspondence analysis.")
                return None

            # Use parallel computation for large matrices
            if use_parallel and data_array.shape[0] > 1000:
                # Set number of threads for NumPy operations
                import os
                original_threads = os.environ.get('OMP_NUM_THREADS', '1')
                os.environ['OMP_NUM_THREADS'] = str(self.config.num_processes)
                os.environ['MKL_NUM_THREADS'] = str(self.config.num_processes)
                os.environ['OPENBLAS_NUM_THREADS'] = str(self.config.num_processes)

            row_sums = data_array.sum(axis=1, keepdims=True)
            col_sums = data_array.sum(axis=0, keepdims=True)
            total = data_array.sum()
            print("...30%", end="", flush=True)

            if total == 0:
                print("\nSum of data matrix is zero, cannot perform correspondence analysis.")
                return None

            P = data_array / total
            print("...40%", end="", flush=True)

            r = row_sums / total
            c = col_sums / total
            print("...50%", end="", flush=True)

            # For very large matrices, use sparse diagonal matrices
            n_rows, n_cols = data_array.shape
            if n_rows > 10000 or n_cols > 1000:
                # Use vectorized operations instead of creating large diagonal matrices
                r_inv_sqrt = 1.0 / np.sqrt(r.flatten())
                c_inv_sqrt = 1.0 / np.sqrt(c.flatten())
                
                # Compute S using broadcasting
                rc_outer = np.outer(r, c)
                S = (P - rc_outer) * r_inv_sqrt[:, np.newaxis] * c_inv_sqrt[np.newaxis, :]
                print("...60% (using optimized computation)", end="", flush=True)
            else:
                Dr_inv_sqrt = np.diag(1.0 / np.sqrt(r.flatten()))
                Dc_inv_sqrt = np.diag(1.0 / np.sqrt(c.flatten()))
                print("...60%", end="", flush=True)
                S = Dr_inv_sqrt @ (P - np.outer(r, c)) @ Dc_inv_sqrt
            
            print("...70%", end="", flush=True)

            # Use appropriate SVD solver based on matrix size
            if n_rows > 5000 and n_cols > 50:
                # For large matrices, compute only the first few components
                from scipy.sparse.linalg import svds
                try:
                    # Compute top k singular values/vectors
                    k = min(50, min(n_rows, n_cols) - 1)
                    u, s, vh = svds(S, k=k)
                    # svds returns singular values in ascending order, so reverse
                    idx = np.argsort(s)[::-1]
                    u = u[:, idx]
                    s = s[idx]
                    vh = vh[idx, :]
                    print("...80% (using sparse SVD)", end="", flush=True)
                except:
                    # Fall back to standard SVD if sparse fails
                    u, s, vh = np.linalg.svd(S, full_matrices=False)
                    print("...80%", end="", flush=True)
            else:
                u, s, vh = np.linalg.svd(S, full_matrices=False)
                print("...80%", end="", flush=True)
            
            # Restore original thread settings
            if use_parallel and data_array.shape[0] > 1000:
                os.environ['OMP_NUM_THREADS'] = original_threads
                if 'MKL_NUM_THREADS' in locals():
                    del os.environ['MKL_NUM_THREADS']
                if 'OPENBLAS_NUM_THREADS' in locals():
                    del os.environ['OPENBLAS_NUM_THREADS']

            n_components = min(5, len(s))

            # Calculate coordinates based on which method was used
            if n_rows > 10000 or n_cols > 1000:
                # Use broadcasting for large matrices
                row_coords = u[:, :n_components] * s[:n_components] * r_inv_sqrt[:, np.newaxis]
                col_coords = vh.T[:, :n_components] * s[:n_components] * c_inv_sqrt[:, np.newaxis]
            else:
                row_coords = Dr_inv_sqrt @ u[:, :n_components] * s[:n_components]
                col_coords = Dc_inv_sqrt @ vh.T[:, :n_components] * s[:n_components]

            result = {
                'coordinates': pd.DataFrame(
                    row_coords,
                    index=data_matrix.index,
                    columns=[f'Axis{i+1}' for i in range(n_components)]
                ),
                'explained_variance': (s**2 / np.sum(s**2))[:n_components],
                'loadings': pd.DataFrame(
                    col_coords,
                    index=data_matrix.columns,
                    columns=[f'Axis{i+1}' for i in range(n_components)]
                ),
                'analysis_type': analysis_type,
                'data_type': data_type
            }
            print("...100% Done!")
        else:
            print(f"Unsupported analysis type: {analysis_type}")
            return None

        print("\nExplained variance by dimension:")
        cum_var = 0
        for i, var in enumerate(result['explained_variance']):
            cum_var += var
            print(f"  {result['coordinates'].columns[i]}: {var:.4f} ({var*100:.2f}%) - Cumulative: {cum_var*100:.2f}%")

        # Store the result with metadata about what type of analysis was performed
        self.multivariate_results = result
        self.multivariate_history = getattr(self, 'multivariate_history', [])
        self.multivariate_history.append({
            'analysis_type': analysis_type,
            'data_type': data_type,
            'result': result,
            'timestamp': datetime.now().isoformat()
        })
        return result

    def detect_clusters(self, method='kmeans', n_clusters=None, min_samples=5, eps=None, max_k=10):
        """
        Detect clusters in multivariate analysis results.
        
        Parameters:
        - method: 'kmeans', 'dbscan', or 'auto'
        - n_clusters: Number of clusters for k-means (None for auto-detection)
        - min_samples: Minimum samples for DBSCAN
        - eps: Epsilon parameter for DBSCAN (None for auto)
        - max_k: Maximum k to test for auto k-means
        
        Returns:
        - Dictionary with cluster assignments and statistics
        """
        if not hasattr(self, 'multivariate_results') or not self.multivariate_results:
            print("No multivariate analysis results found. Please perform multivariate analysis first.")
            return None
        
        # Get coordinates (use first 2-3 dimensions)
        coords = self.multivariate_results['coordinates']
        n_dims = min(3, coords.shape[1])
        X = coords.iloc[:, :n_dims].values
        
        # Standardize the data
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        cluster_labels = None
        cluster_method = method
        
        if method == 'kmeans' or (method == 'auto' and len(X) < 50000):
            if n_clusters is None:
                # Find optimal k using silhouette score
                print("Finding optimal number of clusters...")
                silhouette_scores = []
                K_range = range(2, min(max_k + 1, len(X) // 5))  # Don't test more than n/5 clusters
                
                for k in K_range:
                    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
                    labels = kmeans.fit_predict(X_scaled)
                    score = silhouette_score(X_scaled, labels)
                    silhouette_scores.append(score)
                    print(f"  k={k}: silhouette score = {score:.3f}")
                
                # Find k with highest silhouette score
                optimal_k = K_range[np.argmax(silhouette_scores)]
                n_clusters = optimal_k
                print(f"\nOptimal number of clusters: {n_clusters}")
            
            # Perform k-means clustering
            print(f"\nPerforming k-means clustering with k={n_clusters}...")
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            cluster_labels = kmeans.fit_predict(X_scaled)
            cluster_centers = scaler.inverse_transform(kmeans.cluster_centers_)
            cluster_method = f'k-means (k={n_clusters})'
            
        elif method == 'dbscan' or (method == 'auto' and len(X) >= 50000):
            if eps is None:
                # Estimate eps using k-distance graph (k = min_samples)
                from sklearn.neighbors import NearestNeighbors
                neighbors = NearestNeighbors(n_neighbors=min_samples)
                neighbors_fit = neighbors.fit(X_scaled)
                distances, indices = neighbors_fit.kneighbors(X_scaled)
                distances = np.sort(distances[:, min_samples-1], axis=0)
                
                # Use the knee/elbow point as eps (simplified: 90th percentile)
                eps = np.percentile(distances, 90)
                print(f"Estimated eps parameter: {eps:.3f}")
            
            # Perform DBSCAN clustering
            print(f"\nPerforming DBSCAN clustering (eps={eps:.3f}, min_samples={min_samples})...")
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            cluster_labels = dbscan.fit_predict(X_scaled)
            n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
            n_noise = np.sum(cluster_labels == -1)
            print(f"Found {n_clusters} clusters and {n_noise} noise points")
            cluster_method = f'DBSCAN (eps={eps:.3f})'
        
        if cluster_labels is None:
            return None
        
        # Calculate cluster statistics
        cluster_stats = self._calculate_cluster_statistics(cluster_labels, coords, X)
        
        # Store results
        self.cluster_results = {
            'method': cluster_method,
            'labels': cluster_labels,
            'n_clusters': len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0),
            'coordinates': coords.columns[:n_dims].tolist(),
            'statistics': cluster_stats,
            'analysis_type': self.multivariate_results.get('analysis_type'),
            'data_type': self.multivariate_results.get('data_type'),
            'gene_indices': coords.index.tolist()  # Store which genes were clustered
        }
        
        return self.cluster_results
    
    def _calculate_cluster_statistics(self, labels, coords, X):
        """Calculate statistics for each cluster."""
        unique_labels = sorted(set(labels))
        stats = {}
        
        for label in unique_labels:
            if label == -1:  # Skip noise points for DBSCAN
                continue
                
            mask = labels == label
            cluster_coords = X[mask]
            cluster_genes = coords.index[mask].tolist()
            
            # Calculate centroid
            centroid = np.mean(cluster_coords, axis=0)
            
            # Calculate within-cluster variance
            distances = np.sqrt(np.sum((cluster_coords - centroid)**2, axis=1))
            
            stats[f'Cluster_{label}'] = {
                'size': np.sum(mask),
                'genes': cluster_genes,
                'centroid': centroid.tolist(),
                'mean_distance': np.mean(distances),
                'std_distance': np.std(distances),
                'compactness': np.mean(distances) / (np.std(distances) + 1e-6)  # Lower is more compact
            }
        
        return stats
    
    def characterize_clusters(self):
        """Characterize clusters by their codon usage patterns."""
        if not hasattr(self, 'cluster_results') or not self.cluster_results:
            print("No cluster results found. Please run cluster detection first.")
            return None
        
        labels = self.cluster_results['labels']
        unique_labels = sorted(set(labels))
        data_type = self.cluster_results.get('data_type', 'RSCU')
        
        print(f"\nCharacterizing clusters based on {data_type} patterns...")
        
        # Get the appropriate data matrix
        if data_type == 'RSCU':
            data_matrix = self.get_rscu_values()
        else:  # AA
            data_matrix = self.get_amino_acid_usage()
            # Convert to frequencies
            row_sums = data_matrix.sum(axis=1)
            data_matrix = data_matrix.div(row_sums[row_sums > 0], axis=0)
        
        # Get the genes that were used in multivariate analysis (and clustering)
        coords = self.multivariate_results['coordinates']
        clustered_genes = coords.index
        
        # Filter data matrix to only include clustered genes
        data_matrix = data_matrix.loc[clustered_genes]
        
        characterization = {}
        
        for label in unique_labels:
            if label == -1:  # Skip noise
                continue
                
            cluster_name = f'Cluster_{label}'
            mask = labels == label
            cluster_data = data_matrix[mask]
            
            # Calculate mean usage for cluster
            cluster_mean = cluster_data.mean()
            
            # Calculate mean usage for all other genes
            other_mean = data_matrix[~mask].mean()
            
            # Find most distinctive features (largest absolute differences)
            differences = cluster_mean - other_mean
            abs_diff = np.abs(differences)
            
            # Get top distinctive features
            n_top = min(10, len(abs_diff))
            top_indices = abs_diff.nlargest(n_top).index
            
            distinctive_features = []
            for feature in top_indices:
                diff = differences[feature]
                cluster_val = cluster_mean[feature]
                other_val = other_mean[feature]
                
                if data_type == 'RSCU':
                    # For RSCU, show actual values
                    distinctive_features.append({
                        'codon': feature,
                        'cluster_rscu': cluster_val,
                        'other_rscu': other_val,
                        'difference': diff,
                        'enrichment': cluster_val / (other_val + 1e-6)
                    })
                else:
                    # For AA, show as percentages
                    distinctive_features.append({
                        'amino_acid': feature,
                        'cluster_freq': cluster_val * 100,
                        'other_freq': other_val * 100,
                        'difference': diff * 100,
                        'enrichment': cluster_val / (other_val + 1e-6)
                    })
            
            # Calculate additional metrics
            if hasattr(self, 'enc_values') and self.enc_values:
                cluster_enc = [self.enc_values[gene] for gene in self.cluster_results['statistics'][cluster_name]['genes'] if gene in self.enc_values]
                mean_enc = np.mean(cluster_enc) if cluster_enc else None
            else:
                mean_enc = None
            
            # Calculate GC content for clustered genes only
            cluster_genes = self.cluster_results['statistics'][cluster_name]['genes']
            # Find intersection with base_composition index
            available_genes = [g for g in cluster_genes if g in self.base_composition.index]
            if available_genes:
                cluster_gc = self.base_composition.loc[available_genes, 'GC'].mean()
                cluster_gc3s = self.base_composition.loc[available_genes, 'GC3s'].mean()
            else:
                cluster_gc = None
                cluster_gc3s = None
            
            characterization[cluster_name] = {
                'size': self.cluster_results['statistics'][cluster_name]['size'],
                'distinctive_features': distinctive_features,
                'mean_GC': cluster_gc,
                'mean_GC3s': cluster_gc3s,
                'mean_ENC': mean_enc,
                'compactness': self.cluster_results['statistics'][cluster_name]['compactness']
            }
        
        self.cluster_characterization = characterization
        return characterization

    def identify_reference_genes_multivariate(self, percentage=10):
        """
        Use multivariate analysis to identify potential reference genes for optimal codon determination.

        Parameters:
        percentage (float): Percentage of genes at the extreme of Axis1 to consider as reference

        Returns:
        list: List of gene names identified as reference genes
        """
        print("Performing multivariate analysis to identify reference genes...")

        # Perform CA on RSCU values if not already done
        if not self.multivariate_results or self.multivariate_results.get('analysis_type') != 'CA':
            result = self.perform_multivariate_analysis('CA', 'RSCU')
        else:
            result = self.multivariate_results

        if not result:
            print("Multivariate analysis failed.")
            return None

        # Sort genes by Axis1 position
        coords = result['coordinates']
        axis1_sorted = coords.sort_values(by=coords.columns[0])
        num_genes = len(axis1_sorted)
        num_extreme_genes = max(2, int(num_genes * percentage / 100))

        # Get genes at both extremes
        left_extreme_genes = axis1_sorted.index[:num_extreme_genes].tolist()
        right_extreme_genes = axis1_sorted.index[-num_extreme_genes:].tolist()

        # Calculate ENC values for both extremes to determine which group contains highly expressed genes
        print("Calculating ENC values to identify highly expressed genes...")
        enc_values = self.calculate_enc()

        left_avg_enc = sum(enc_values[gene] for gene in left_extreme_genes) / len(left_extreme_genes)
        right_avg_enc = sum(enc_values[gene] for gene in right_extreme_genes) / len(right_extreme_genes)

        # Calculate GC3s averages for both extremes
        # GC3s is already calculated during sequence processing
        left_avg_gc3s = sum(self.base_composition.at[gene, 'GC3s'] for gene in left_extreme_genes) / len(left_extreme_genes)
        right_avg_gc3s = sum(self.base_composition.at[gene, 'GC3s'] for gene in right_extreme_genes) / len(right_extreme_genes)

        print(f"Left extreme - Avg ENC: {left_avg_enc:.2f}, Avg GC3s: {left_avg_gc3s:.2f}")
        print(f"Right extreme - Avg ENC: {right_avg_enc:.2f}, Avg GC3s: {right_avg_gc3s:.2f}")

        # Lower ENC and higher GC3s → likely highly expressed
        left_score = (61 - left_avg_enc) + (left_avg_gc3s / 2)
        right_score = (61 - right_avg_enc) + (right_avg_gc3s / 2)

        if left_score > right_score:
            print("Left extreme genes likely contain highly expressed genes.")
            reference_genes = left_extreme_genes
        else:
            print("Right extreme genes likely contain highly expressed genes.")
            reference_genes = right_extreme_genes

        return reference_genes

    def calculate_optimal_codons_multivariate(self, percentage=10):
        """
        Identify optimal codons based on multivariate analysis.

        Parameters:
        percentage (float): Percentage of genes at the extreme of Axis1 to consider as reference

        Returns:
        dict: Dictionary mapping amino acids to their optimal codons
        """
        # Identify reference genes through multivariate analysis
        reference_genes = self.identify_reference_genes_multivariate(
            percentage)

        if not reference_genes:
            print("Failed to identify reference genes. Using all genes instead.")
            reference_genes = self.gene_names

        # Use the reference genes to calculate optimal codons
        return self.calculate_optimal_codons(reference_genes)

    def calculate_fop_multivariate(self, percentage=10):
        """
        Calculate Frequency of Optimal Codons (Fop) using multivariate-selected reference genes.

        Parameters:
        percentage (float): Percentage of genes to use from multivariate analysis

        Returns:
        dict: Dictionary mapping gene names to Fop values
        """
        # Identify optimal codons using multivariate analysis
        if self.optimal_codons is None:
            self.calculate_optimal_codons_multivariate(percentage)

        # Now calculate Fop using these optimal codons
        return self.calculate_fop()

    def calculate_cai_multivariate(self, percentage=10):
        """
        Calculate Codon Adaptation Index (CAI) using multivariate-selected reference genes.

        Parameters:
        percentage (float): Percentage of genes to use from multivariate analysis

        Returns:
        dict: Dictionary mapping gene names to CAI values
        """
        # Identify reference genes through multivariate analysis
        reference_genes = self.identify_reference_genes_multivariate(
            percentage)

        if not reference_genes:
            print("Failed to identify reference genes. Using all genes instead.")
            reference_genes = self.gene_names

        # Calculate CAI using these reference genes
        return self.calculate_cai(reference_genes)

    def get_comprehensive_metrics(self):
        """Generate a comprehensive DataFrame with all metrics for each gene."""
        # Make sure all metrics are calculated
        self.perform_multivariate_analysis() if self.multivariate_results is None else None
        self.calculate_enc() if self.enc_values is None else None
        self.calculate_fop() if self.fop_values is None else None
        self.calculate_cai() if self.cai_values is None else None
        scuo_values = self.calculate_scuo()

        # Create a complete metrics DataFrame
        metrics = pd.DataFrame(index=self.gene_names)

        # Add base composition data
        for col in self.base_composition.columns:
            metrics[col] = self.base_composition[col]

        # Add multivariate coordinates (first 4 axes)
        if self.multivariate_results:
            coords = self.multivariate_results['coordinates']
            for col in coords.columns[:min(4, len(coords.columns))]:
                metrics[col] = coords[col]

        # Add ENC, Fop, CAI, and SCUO values
        metrics['ENC'] = pd.Series(self.enc_values)
        metrics['Fop'] = pd.Series(self.fop_values)
        metrics['CAI'] = pd.Series(self.cai_values)
        metrics['SCUO'] = pd.Series(scuo_values)

        return metrics

    def optimize_gene_sequence(self,gene_id,method="frequency",reference_genes=None,percentage=10,output_format='fasta'):
        """
        Rewrite gene sequence using optimal codons determined by the specified method.

        Parameters:
        gene_id (str): ID of the gene to optimize
        method (str): Method to use for optimal codon selection
           "frequency" - most frequent codon (default)
           "multivariate" - multivariate analysis
           "rscu" - highest RSCU value
           "raw_count" - most common codon
        reference_genes (list): Reference genes for frequency method
        percentage (float): Percentage of genes for multivariate analysis
        output_format (str): Output format ('fasta', 'sequence', or 'both')

        Returns:
        str or tuple: Optimized sequence in requested format
        """
        # Determine optimal codons if not already done
        if self.optimal_codons is None:
            self.calculate_optimal_codons(reference_genes, method, percentage)

        # Get the gene sequence using unified access method
        gene_seq = self.get_sequence(gene_id)
        
        if gene_seq is None:
            print(f"Gene {gene_id} not found in loaded sequences.")
            return None

        # Get the appropriate DNA-amino acid mapping for current genetic code
        dna_aa_map = self.config.get_dna_aa_map()

        # Create a mapping of amino acids to their optimal DNA codons
        optimal_dna_codons = {}
        for aa, rna_codon in self.optimal_codons.items():
            dna_codon = rna_codon.replace('U', 'T')
            optimal_dna_codons[aa] = dna_codon

        # Translate the original sequence and rebuild with optimal codons
        dna_seq = str(gene_seq.seq).upper()
        optimized_seq = ""

        i = 0
        while i + 3 <= len(dna_seq):
            codon = dna_seq[i:i + 3]

            # Skip if codon contains invalid characters
            if 'N' in codon or len(codon) < 3:
                optimized_seq += codon
                i += 3
                continue

            # Get amino acid for this codon
            aa = dna_aa_map.get(codon, 'X')  # Use 'X' for unknown amino acids

            if aa == 'STOP':
                # Keep original stop codon
                optimized_seq += codon
            elif aa in optimal_dna_codons:
                # Replace with optimal codon
                optimized_seq += optimal_dna_codons[aa]
            else:
                # Keep original codon if we don't have an optimal one
                optimized_seq += codon

            i += 3

        # Add any remaining partial codon at the end (if any)
        if i < len(dna_seq):
            optimized_seq += dna_seq[i:]

        # Return the optimized sequence in the requested format
        if output_format == 'fasta':
            return f">{gene_id}_optimized\n{optimized_seq}"
        elif output_format == 'sequence':
            return optimized_seq
        elif output_format == 'both':
            return gene_seq, optimized_seq
        else:
            return optimized_seq

    def optimize_all_genes(self, output_path=None):
        """Rewrite all gene sequences using optimal codons."""
        if self.optimal_codons is None:
            self.calculate_optimal_codons()

        optimized_sequences = []

        for gene_id in self.gene_names:
            optimized_fasta = self.optimize_gene_sequence(gene_id, 'fasta')
            if optimized_fasta:
                optimized_sequences.append(optimized_fasta)

        # Join all optimized sequences
        combined_fasta = "\n".join(optimized_sequences)

        # Save to file if output path provided
        if output_path:
            with open(output_path, 'w') as f:
                f.write(combined_fasta)
            print(f"Optimized sequences saved to {output_path}")

        return combined_fasta

    def load_reference_genes_from_file(self, file_path):
        """
        Load reference gene names from a file.
        File should contain one gene name per line.

        Parameters:
        file_path (str): Path to the file containing reference gene names

        Returns:
        list: List of valid gene names found in the file
        """
        print(f"Loading reference genes from {file_path}...")
        try:
            with open(file_path, 'r') as f:
                genes = [line.strip() for line in f if line.strip()]

            # Validate gene names
            valid_genes = [g for g in genes if g in self.gene_names]

            if not valid_genes:
                print("No valid gene names found in the file.")
                return None

            if len(valid_genes) != len(genes):
                print(f"Warning: {len(genes) - len(valid_genes)} gene names from the file were not found.")
            print(f"Proceeding with {len(valid_genes)} valid genes.")

            return valid_genes
        except Exception as e:
            print(f"Error reading file: {e}")
            return None

    def visualize(self, vis_type, data=None, output_path=None):
        """
        Create and display a visualization of the specified type.

        Parameters:
        vis_type (str): Type of visualization to generate
        data: Optional data to visualize; defaults to the current object
        output_path (str): Optional path to save the visualization

        Returns:
        str or None: Path to the saved visualization, if any
        """
        try:
            # Instantiate the appropriate visualization object
            vis = self.visualization_factory.create_visualization(
                vis_type, self.config.visualization_method
            )

            # Determine data to use
            vis_data = data if data is not None else self

            # Generate and optionally save the visualization
            result_path = vis.create(vis_data, output_path)

            # Open the visualization if saved to a path (only for HTML files)
            if result_path and result_path.lower().endswith('.html'):
                vis.open_in_browser(result_path)

            return result_path

        except Exception as e:
            print(f"Error creating visualization '{vis_type}': {e}")
            return None

    def save_optimal_codons_to_file(self,output_path,optimal_codons=None,reference_genes=None,format="tsv"):
        """
        Save optimal codons to a file in either TSV or JSON format.

        Parameters:
        output_path (str): Path where the file will be saved
        optimal_codons (dict): Dictionary of optimal codons (if None, will use self.optimal_codons)
        reference_genes (list): List of reference genes used (optional)
        format (str): Format to save in - "tsv" or "json"

        Returns:
        str: Path to saved file
        """
        # Calculate if not provided
        if not optimal_codons:
            if not self.optimal_codons:
                print("\nIdentifying optimal codons...")
                optimal_codons = self.calculate_optimal_codons(reference_genes)
            else:
                optimal_codons = self.optimal_codons

        # Create metadata about the genetic code
        metadata = {
            "date": pd.Timestamp.now().isoformat(),
            "genetic_code": self.config.genetic_code,
            "genetic_code_name": self.config.get_genetic_code_name(),
            "reference_genes": reference_genes if reference_genes else [],
            "total_genes": len(self.gene_names),
            "gcua_version": VERSION
        }

        if format.lower() == "json":
            # Create a JSON structure with metadata
            data = {
                "optimal_codons": {},
                "metadata": metadata
            }

            # Add the optimal codons
            for aa in sorted(optimal_codons.keys()):
                rna_codon = optimal_codons[aa]
                dna_codon = rna_codon.replace('U', 'T')
                data["optimal_codons"][aa] = {
                    "rna_codon": rna_codon,
                    "dna_codon": dna_codon
                }

            # Write to JSON file
            with open(output_path, 'w') as f:
                json.dump(data, f, indent=4)

        else:  # TSV format
            with open(output_path, 'w') as f:
                # Write metadata as comments
                f.write("# GCUA Optimal Codons\n")
                for key, value in metadata.items():
                    f.write(f"# {key}: {value}\n")
                f.write("#" + "-" * 50 + "\n\n")

                # Write data
                f.write("AA\tRNA_Codon\tDNA_Codon\n")
                for aa in sorted(optimal_codons.keys()):
                    rna_codon = optimal_codons[aa]
                    dna_codon = rna_codon.replace('U', 'T')
                    f.write(f"{aa}\t{rna_codon}\t{dna_codon}\n")

        print(f"\nOptimal codons saved to {output_path}")
        return output_path

    def load_optimal_codons_from_file(self, file_path):
        """
        Load optimal codons from a file (TSV or JSON).
        Parameters:
        file_path (str): Path to the file containing optimal codons
        Returns:
        dict: Dictionary mapping amino acids to their optimal codons
        """
        print(f"Loading optimal codons from {file_path}...")
        try:
            # Determine file type from extension
            ext = os.path.splitext(file_path)[1].lower()

            if ext == '.json':
                # JSON format
                with open(file_path, 'r') as f:
                    data = json.load(f)

                # Extract optimal codons
                if "optimal_codons" in data:
                    optimal_codons = {}
                    for aa, codon_info in data["optimal_codons"].items():
                        optimal_codons[aa] = codon_info["rna_codon"]

                    # Print metadata if available
                    if "metadata" in data:
                        print(f"Metadata from file:")
                        for key, value in data["metadata"].items():
                            print(f"  {key}: {value}")
                else:
                    # Try a simpler format
                    optimal_codons = {
                        aa: codon_info["rna_codon"] for aa,
                        codon_info in data.items()}
            else:  # TSV format
                optimal_codons = {}
                with open(file_path, 'r') as f:
                    header = f.readline().strip().split('\t')

                    # Find column indices
                    aa_idx = header.index('AA') if 'AA' in header else 0
                    rna_idx = header.index(
                        'RNA_Codon') if 'RNA_Codon' in header else 1

                    # Read data
                    for line in f:
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) > max(aa_idx, rna_idx):
                                aa = parts[aa_idx]
                                rna_codon = parts[rna_idx]
                                optimal_codons[aa] = rna_codon

            # Store and return
            self.optimal_codons = optimal_codons
            print(f"Loaded {len(optimal_codons)} optimal codons from {file_path}")
            return optimal_codons

        except Exception as e:
            print(f"Error loading optimal codons: {e}")
            return None

    def optimize_with_external_optimal_codons(self, file_path, gene_id=None, output_format='fasta'):
        """
        Optimize sequences using optimal codons loaded from an external file.

        Parameters:
        file_path (str): Path to the file containing optimal codons
        gene_id (str): ID of a specific gene to optimize (if None, all genes will be optimized)
        output_format (str): Output format - 'fasta', 'sequence', or 'both'

        Returns:
        str or dict: Optimized sequence(s) in the requested format
        """
        # Load optimal codons
        optimal_codons = self.load_optimal_codons_from_file(file_path)
        if not optimal_codons:
            return None

        # Store the optimal codons
        self.optimal_codons = optimal_codons

        # Optimize a single gene or all genes
        if gene_id:
            return self.optimize_gene_sequence(gene_id, output_format=output_format)
        else:
            return self.optimize_all_genes()

    def compare_codon_usage_axis_cohorts(self,percentage=10,chi_squared_threshold=0.05):
        """
        Compare codon usage between cohorts at opposite ends of Axis 1 from multivariate analysis.
        Uses chi-squared test to identify codons with significant usage differences.

        Parameters:
            percentage (float): Percentage of genes at each end of Axis 1 to use
            chi_squared_threshold (float): P-value threshold for significance

        Returns:
            dict: Results of the comparison including significant codons and their statistics
        """
        print(f"Comparing codon usage between gene cohorts at opposite ends of Axis 1...")

        # Perform multivariate analysis if not already done
        if not self.multivariate_results or 'CA' not in str(self.multivariate_results):
            result = self.perform_multivariate_analysis('CA', 'RSCU')
        else:
            result = self.multivariate_results

        if not result:
            print("Multivariate analysis failed.")
            return None

        # Sort genes by Axis1 position
        coords = result['coordinates']
        axis1_sorted = coords.sort_values(by=coords.columns[0])
        num_genes = len(axis1_sorted)
        num_cohort_genes = max(2, int(num_genes * percentage / 100))

        # Get genes at both extremes
        left_extreme_genes = axis1_sorted.index[:num_cohort_genes].tolist()
        right_extreme_genes = axis1_sorted.index[-num_cohort_genes:].tolist()

        print(f"Left extreme cohort: {len(left_extreme_genes)} genes")
        print(f"Right extreme cohort: {len(right_extreme_genes)} genes")

        # Calculate ENC values if not already done
        if not hasattr(self, 'enc_values') or not self.enc_values:
            print("Calculating ENC values...")
            self.calculate_enc()

        # Calculate average ENC for both cohorts
        left_enc_values = [self.enc_values[gene] for gene in left_extreme_genes if gene in self.enc_values]
        right_enc_values = [self.enc_values[gene] for gene in right_extreme_genes if gene in self.enc_values]

        left_avg_enc = sum(left_enc_values) / len(left_enc_values) if left_enc_values else float('nan')
        right_avg_enc = sum(right_enc_values) / len(right_enc_values) if right_enc_values else float('nan')

        # Calculate GC3s averages for both cohorts
        left_avg_gc3s = sum(self.base_composition.at[gene, 'GC3s'] for gene in left_extreme_genes) / len(left_extreme_genes)
        right_avg_gc3s = sum(self.base_composition.at[gene, 'GC3s'] for gene in right_extreme_genes) / len(right_extreme_genes)

        # Identify which cohort likely contains highly expressed genes (lower ENC)
        if left_avg_enc < right_avg_enc:
            highly_expressed_cohort = "left"
            print(f"\nLeft cohort likely contains highly expressed genes (Lower ENC values)")
        else:
            highly_expressed_cohort = "right"
            print(f"\nRight cohort likely contains highly expressed genes (Lower ENC values)")

        print(f"Left cohort - Avg ENC: {left_avg_enc:.2f}, Avg GC3s: {left_avg_gc3s:.2f}%")
        print(f"Right cohort - Avg ENC: {right_avg_enc:.2f}, Avg GC3s: {right_avg_gc3s:.2f}%")

        # Calculate cumulative codon usage for both cohorts
        left_codon_usage = self.codon_counts.loc[left_extreme_genes].sum()
        right_codon_usage = self.codon_counts.loc[right_extreme_genes].sum()

        # Calculate cumulative RSCU for both cohorts
        left_data = GCUAData(self.config)
        left_data.gene_names = ["Left_Cohort"]
        left_data.codon_counts = pd.DataFrame([left_codon_usage], index=["Left_Cohort"])

        right_data = GCUAData(self.config)
        right_data.gene_names = ["Right_Cohort"]
        right_data.codon_counts = pd.DataFrame([right_codon_usage], index=["Right_Cohort"])

        aa_map = self.config.get_aa_map()
        aa_to_codons = AA_TO_CODONS_MAP.get(self.config.genetic_code, AA_TO_CODONS_MAP[1])

        results = {
            "left_cohort": left_extreme_genes,
            "right_cohort": right_extreme_genes,
            "left_cohort_usage": {codon: left_codon_usage[codon] for codon in CODONS},
            "right_cohort_usage": {codon: right_codon_usage[codon] for codon in CODONS},
            "left_cohort_rscu": {codon: left_data.get_rscu_values().at["Left_Cohort", codon] for codon in CODONS},
            "right_cohort_rscu": {codon: right_data.get_rscu_values().at["Right_Cohort", codon] for codon in CODONS},
            "left_avg_enc": left_avg_enc,
            "right_avg_enc": right_avg_enc,
            "left_avg_gc3s": left_avg_gc3s,
            "right_avg_gc3s": right_avg_gc3s,
            "highly_expressed_cohort": highly_expressed_cohort,
            "significant_differences": {},
            "chi_squared_results": {},
            "preferred_codons": {"left": {}, "right": {}}
        }

        print("\nStatistical comparison of codon usage between cohorts:")
        print("\nAmino Acid\tChi-squared\tP-value\tSignificant")

        for aa, codons in aa_to_codons.items():
            if aa == 'STOP' or len(codons) < 2:
                continue

            contingency = []
            for codon in codons:
                left_count = left_codon_usage[codon]
                right_count = right_codon_usage[codon]
                contingency.append([left_count, right_count])

            if sum(sum(row) for row in contingency) == 0:
                continue

            row_sums = [sum(row) for row in contingency]
            col_sums = [sum(row[i] for row in contingency) for i in range(2)]

            if 0 in row_sums or 0 in col_sums:
                print(f"{aa}\tSkipped\tSkipped\tNo (zero counts)")
                continue

            min_count = min(count for row in contingency for count in row if count > 0)
            if min_count < 5:
                contingency = [[count + 0.5 for count in row] for row in contingency]

            try:
                chi2, p, dof, expected = chi2_contingency(contingency)
                results["chi_squared_results"][aa] = {
                    "chi2": chi2,
                    "p_value": p,
                    "significant": p < chi_squared_threshold,
                    "codons": codons,
                    "contingency": contingency
                }

                significant = "Yes" if p < chi_squared_threshold else "No"
                print(f"{aa}\t{chi2:.4f}\t{p:.4f}\t{significant}")

                if p < chi_squared_threshold:
                    results["significant_differences"][aa] = {
                        "codons": codons,
                        "contingency": contingency,
                        "chi2": chi2,
                        "p_value": p
                    }

                    left_preferred = max(codons, key=lambda c: left_codon_usage[c])
                    right_preferred = max(codons, key=lambda c: right_codon_usage[c])

                    results["preferred_codons"]["left"][aa] = left_preferred
                    results["preferred_codons"]["right"][aa] = right_preferred

                    codon_details = []
                    for i, codon in enumerate(codons):
                        left_count = contingency[i][0]
                        right_count = contingency[i][1]
                        total_aa_left = sum(row[0] for row in contingency)
                        total_aa_right = sum(row[1] for row in contingency)

                        left_rscu = left_count * len(codons) / total_aa_left if total_aa_left > 0 else 0
                        right_rscu = right_count * len(codons) / total_aa_right if total_aa_right > 0 else 0

                        codon_details.append({
                            "codon": codon,
                            "left_count": left_count,
                            "right_count": right_count,
                            "left_rscu": left_rscu,
                            "right_rscu": right_rscu,
                            "preferred_in": "left" if left_count > right_count else "right"
                        })

                    results["significant_differences"][aa]["codon_details"] = codon_details

            except Exception as e:
                print(f"{aa}\tError\tError\tNo (test failed: {str(e)})")
                continue

        sig_count = sum(1 for res in results["chi_squared_results"].values() if res["significant"])
        print(f"\nFound significant differences in codon usage for {sig_count} out of {len(results['chi_squared_results'])} amino acids")

        print("\nComprehensive Codon Usage Comparison:\n")
        print(f"Left cohort (Avg ENC: {left_avg_enc:.2f}) vs Right cohort (Avg ENC: {right_avg_enc:.2f})")
        print(f"Likely highly expressed cohort: {highly_expressed_cohort.upper()} (Lower ENC values)")
        print("\n   AA  Codon   Left Count   Left RSCU   Right Count   Right RSCU   Significant")
        print("  -------------------------------------------------------------------------")

        significant_codons = set()
        for aa, result in results["chi_squared_results"].items():
            if result["significant"]:
                significant_codons.update(result["codons"])

        first_letters = ['U', 'C', 'A', 'G']
        second_letters = ['U', 'C', 'A', 'G']
        third_letters = ['U', 'C', 'A', 'G']

        for first in first_letters:
            for second in second_letters:
                for third in third_letters:
                    codon = f"{first}{second}{third}"
                    if codon in CODONS:
                        aa = aa_map[codon]
                        left_count = results["left_cohort_usage"][codon]
                        right_count = results["right_cohort_usage"][codon]
                        left_rscu = results["left_cohort_rscu"][codon]
                        right_rscu = results["right_cohort_rscu"][codon]

                        is_significant = codon in significant_codons
                        sig_marker = "*" if is_significant else " "

                        left_preferred = results["preferred_codons"]["left"].get(aa) == codon
                        right_preferred = results["preferred_codons"]["right"].get(aa) == codon

                        left_marker = "+" if left_preferred and is_significant else " "
                        right_marker = "+" if right_preferred and is_significant else " "

                        highly_exp_marker = ""
                        if is_significant:
                            if highly_expressed_cohort == "left" and left_preferred:
                                highly_exp_marker = " ◄"
                            elif highly_expressed_cohort == "right" and right_preferred:
                                highly_exp_marker = " ►"

                        print(f"  {aa:3s} {codon:5s} {left_count:10.0f} {left_marker}{left_rscu:9.2f} {right_count:12.0f} {right_marker}{right_rscu:9.2f}   {sig_marker}{highly_exp_marker}")

        print("\nLegend:")
        print("  * - Indicates amino acid with significant differences in codon usage")
        print("  + - Indicates preferred codon in that cohort")
        print("  ◄ - Indicates preferred codon in left cohort (if left is highly expressed)")
        print("  ► - Indicates preferred codon in right cohort (if right is highly expressed)")
        print(f"\nOptimal codons from the {highly_expressed_cohort.upper()} cohort (likely highly expressed genes) can be used for gene optimization.")

        return results

    def save_codon_comparison_results(self, results, output_path):
        """
        Save codon usage comparison results to files.

        Parameters:
        results (dict): Results from compare_codon_usage_axis_cohorts
        output_path (str): Base path for output files

        Returns:
        list: Paths to saved files
        """
        saved_files = []

        # Save summary results
        summary_path = f"{output_path}_summary.tsv"
        with open(summary_path, 'w') as f:
            # Add metadata
            f.write(f"# GCUA Codon Usage Comparison\n")
            f.write(f"# Version: {VERSION}\n")
            f.write(
                f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
            f.write(f"# Left Cohort: {len(results['left_cohort'])} genes\n")
            f.write(f"# Right Cohort: {len(results['right_cohort'])} genes\n")
            f.write(f"# Left Avg ENC: {results['left_avg_enc']:.2f}\n")
            f.write(f"# Right Avg ENC: {results['right_avg_enc']:.2f}\n")
            f.write(f"# Highly Expressed Cohort: {results['highly_expressed_cohort']}\n")
            f.write(f"#{'=' * 50}\n\n")

            # Write data
            f.write("Amino_Acid\tChi_squared\tP_value\tSignificant\tLeft_preferred\tRight_preferred\n")
            for aa, res in results["chi_squared_results"].items():
                significant = "Yes" if res["significant"] else "No"
                left_preferred = results["preferred_codons"]["left"].get(aa, "N/A")
                right_preferred = results["preferred_codons"]["right"].get(aa, "N/A")
                f.write(
                    f"{aa}\t{res['chi2']:.4f}\t{res['p_value']:.4f}\t{significant}\t{left_preferred}\t{right_preferred}\n"
                )

        saved_files.append(summary_path)
        print(f"Summary results saved to {summary_path}")

        # Save detailed results for significant differences
        if results["significant_differences"]:
            detailed_path = f"{output_path}_detailed.tsv"
            with open(detailed_path, 'w') as f:
                # Add metadata headers
                f.write(f"# GCUA Detailed Codon Usage Comparison\n")
                f.write(f"# Version: {VERSION}\n")
                f.write(
                    f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
                f.write(f"#{'=' * 50}\n\n")

                f.write("Amino_Acid\tCodon\tLeft_count\tRight_count\tLeft_RSCU\tRight_RSCU\tPreferred_in\n")
                for aa, res in results["significant_differences"].items():
                    for codon_detail in res["codon_details"]:
                        f.write(
                            f"{aa}\t{codon_detail['codon']}\t{codon_detail['left_count']}\t{codon_detail['right_count']}\t"
                            f"{codon_detail['left_rscu']:.4f}\t{codon_detail['right_rscu']:.4f}\t{codon_detail['preferred_in']}\n"
                        )

            saved_files.append(detailed_path)
            print(f"Detailed results saved to {detailed_path}")

        # Save cohort-specific codon usage
        left_usage_path = f"{output_path}_left_cohort_usage.tsv"
        with open(left_usage_path, 'w') as f:
            f.write("Codon\tAmino_Acid\tCount\tRSCU\n")
            for codon in CODONS:
                aa = self.config.get_aa_map()[codon]
                count = results["left_cohort_usage"][codon]
                rscu = results["left_cohort_rscu"][codon]
                f.write(f"{codon}\t{aa}\t{count}\t{rscu:.4f}\n")

        saved_files.append(left_usage_path)
        print(f"Left cohort codon usage saved to {left_usage_path}")

        right_usage_path = f"{output_path}_right_cohort_usage.tsv"
        with open(right_usage_path, 'w') as f:
            f.write("Codon\tAmino_Acid\tCount\tRSCU\n")
            for codon in CODONS:
                aa = self.config.get_aa_map()[codon]
                count = results["right_cohort_usage"][codon]
                rscu = results["right_cohort_rscu"][codon]
                f.write(f"{codon}\t{aa}\t{count}\t{rscu:.4f}\n")

        saved_files.append(right_usage_path)
        print(f"Right cohort codon usage saved to {right_usage_path}")

        # Save preferred codons for both cohorts in optimal codons format
        if results["preferred_codons"]["left"]:
            left_preferred_path = f"{output_path}_left_preferred_codons.tsv"
            self.save_optimal_codons_to_file(
                left_preferred_path,
                results["preferred_codons"]["left"]
            )
            saved_files.append(left_preferred_path)

        if results["preferred_codons"]["right"]:
            right_preferred_path = f"{output_path}_right_preferred_codons.tsv"
            self.save_optimal_codons_to_file(
                right_preferred_path,
                results["preferred_codons"]["right"]
            )
            saved_files.append(right_preferred_path)

        # Save optimized versions of the highly expressed preferred codons
        highly_exp_preferred = results["preferred_codons"][results["highly_expressed_cohort"]]
        if highly_exp_preferred:
            optimal_path = f"{output_path}_optimal_codons.tsv"
            self.save_optimal_codons_to_file(optimal_path, highly_exp_preferred)
            saved_files.append(optimal_path)
            print(f"Optimal codons from highly expressed cohort saved to {optimal_path}")

        # Save JSON with all results - with improved serialization to handle numpy types
        json_path = f"{output_path}_complete.json"
        # (Assuming there's missing JSON write code here)

        return saved_files


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)


class GCUAInterface:
    """Command-line interface for the GCUA tool."""

    def __init__(self, config=None):
        """Initialize the interface with optional config."""
        self.config = config or Config()
        self.data = GCUAData(self.config)
        self.file_manager = FileManager()
        self.show_help_hint = True  # Show help hint on first menu entry
        self._init_help_content()
    
    def _init_help_content(self):
        """Initialize help content for all menu items."""
        self.help_content = {
            # Main menu help
            'main_menu': {
                '1': {
                    'brief': 'Load and manage sequence data files',
                    'detailed': '''DATA MANAGEMENT
                    Load FASTA files, view data summaries, filter sequences, and manage your dataset.
                    Start here to load your sequences before running any analyses.
                    
                    Supports: Standard FASTA format with protein-coding DNA sequences (length divisible by 3)
                    Memory modes: Automatic switching between full and low-memory modes for large files'''
                },
                '2': {
                    'brief': 'Automated workflow for standard codon usage analysis',
                    'detailed': '''QUICK ANALYSIS (GUIDED WORKFLOW)
                    Perfect for new users or standard analyses. Automatically:
                    1. Calculates all basic metrics (codon usage, amino acids, base composition, RSCU)
                    2. Performs multivariate analysis (CA on RSCU)
                    3. Identifies optimal codons
                    4. Creates standard visualizations
                    5. Exports all results
                    
                    Complete analysis in one step - no need to navigate multiple menus!'''
                },
                '3': {
                    'brief': 'Choose specific analyses and customize parameters',
                    'detailed': '''CUSTOM ANALYSIS
                    For experienced users who want control over specific analyses.
                    Access individual metrics, advanced statistical analyses, and comparative tools.
                    Allows you to run only what you need and customize parameters.'''
                },
                '4': {
                    'brief': 'Create plots and export results',
                    'detailed': '''VISUALIZATION & EXPORT
                    Generate publication-quality plots and export data in various formats.
                    Includes: Multivariate plots, heatmaps, scatter plots, and more.
                    Export formats: TSV (data), HTML (interactive plots), SVG/PNG/PDF (static images)'''
                },
                '5': {
                    'brief': 'Advanced features for specialized analyses',
                    'detailed': '''ADVANCED TOOLS
                    Specialized analyses including:
                    - Sequence optimization for heterologous expression
                    - Cluster detection in multivariate space
                    - Outlier analysis to find unusual genes
                    - Export extreme genes for further study
                    - Statistical comparison of gene cohorts'''
                },
                '6': {
                    'brief': 'Configure genetic codes and program settings',
                    'detailed': '''SETTINGS & PREFERENCES
                    Critical settings:
                    - GENETIC CODE: Must be set correctly for your organism!
                    - Visualization formats and methods
                    - Performance settings for large datasets
                    - Output detail levels
                    
                    ⚠️  Always set genetic code BEFORE loading sequences!'''
                },
                '7': {
                    'brief': 'View documentation and program information',
                    'detailed': '''HELP & DOCUMENTATION
                    Access comprehensive help, including:
                    - Program overview and citation
                    - Available analyses and their purposes
                    - Metric explanations (ENC, CAI, Fop, SCUO)
                    - Current genetic code information'''
                }
            },
            # Custom Analysis menu help
            'analysis_menu': {
                '1': {
                    'brief': 'View or export raw codon counts',
                    'detailed': '''CODON USAGE DETAILS
                    Access the raw codon count data for your sequences.
                    - View counts for individual genes
                    - Export full codon count matrix
                    - Check data quality and completeness
                    
                    Useful for: Manual inspection, external analysis, quality control'''
                },
                '2': {
                    'brief': 'Analyze amino acid composition patterns',
                    'detailed': '''AMINO ACID USAGE DETAILS
                    Examine amino acid frequencies in your sequences.
                    - Identify compositional biases
                    - Compare protein functional classes
                    - Detect species-specific preferences
                    
                    Note: Calculated automatically during sequence loading'''
                },
                '3': {
                    'brief': 'Examine nucleotide composition and GC content',
                    'detailed': '''BASE COMPOSITION DETAILS
                    Analyze nucleotide frequencies and GC content:
                    - Overall GC content
                    - Positional GC (GC1, GC2, GC3)
                    - GC3s (synonymous sites only)
                    
                    GC3s is particularly important for gene expression studies'''
                },
                '4': {
                    'brief': 'Calculate Relative Synonymous Codon Usage',
                    'detailed': '''RSCU CALCULATION
                    RSCU removes amino acid composition effects to reveal true codon preferences.
                    Values >1 indicate preferred codons, <1 indicate avoided codons.
                    
                    Essential for: Multivariate analysis, identifying selection pressures'''
                },
                '5': {
                    'brief': 'Reduce data complexity to identify major patterns',
                    'detailed': '''MULTIVARIATE ANALYSIS (CA/PCA)
                    Reduces high-dimensional codon usage data to 2-3 major axes:
                    
                    CA (Correspondence Analysis):
                    - Best for RSCU data (count-based)
                    - Handles sparse data well
                    - Standard for codon usage studies
                    
                    PCA (Principal Component Analysis):
                    - Best for amino acid frequencies
                    - Use CLR transformation for compositional data
                    
                    Axis 1 often correlates with gene expression level!'''
                },
                '6': {
                    'brief': 'Calculate codon bias metrics (ENC, CAI, Fop, SCUO)',
                    'detailed': '''CODON BIAS METRICS
                    Suite of metrics to quantify codon usage bias:
                    
                    ENC (Effective Number of Codons): 20-61
                    - Lower values = stronger bias
                    - Species-independent measure
                    
                    CAI (Codon Adaptation Index): 0-1
                    - Higher values = more adapted to highly expressed genes
                    - Requires reference gene set
                    
                    Fop (Frequency of Optimal Codons): 0-1
                    - Proportion of optimal codons used
                    
                    SCUO (Synonymous Codon Usage Order): 0-1
                    - Information theory-based measure'''
                },
                '7': {
                    'brief': 'Determine preferred codons for each amino acid',
                    'detailed': '''IDENTIFY OPTIMAL CODONS
                    Find the preferred codon for each amino acid using various methods:
                    - Frequency-based (most common)
                    - Multivariate axis extremes (expression-correlated)
                    - RSCU-based (highest RSCU values)
                    - Reference gene set (user-defined)
                    
                    Required for: CAI/Fop calculation, sequence optimization'''
                },
                '8': {
                    'brief': 'Compare codon usage between gene groups',
                    'detailed': '''COMPARE AXIS COHORTS
                    Statistical comparison of genes at opposite ends of the primary multivariate axis.
                    Identifies codons significantly associated with high vs low expression.
                    
                    Uses chi-squared tests to find differential codon usage.
                    Results can define optimal codons for sequence optimization.'''
                },
                '9': {
                    'brief': 'Run all analyses comprehensively',
                    'detailed': '''CALCULATE ALL METRICS
                    Comprehensive analysis that calculates:
                    - All basic metrics (if not already done)
                    - RSCU values
                    - Multivariate analysis (CA on RSCU)
                    - All codon bias metrics (ENC, CAI, Fop, SCUO)
                    - Optimal codon identification
                    
                    One-stop option for complete analysis!'''
                }
            },
            # Advanced tools help
            'advanced_tools': {
                '1': {
                    'brief': 'Redesign sequences using optimal codons',
                    'detailed': '''SEQUENCE OPTIMIZATION
                    Replace codons with optimal alternatives for expression in a host.
                    - Uses identified optimal codons
                    - Maintains amino acid sequence
                    - Can optimize single genes or entire dataset
                    - Compare original vs optimized sequences
                    
                    Applications: Heterologous expression, synthetic biology'''
                },
                '2': {
                    'brief': 'Find groups of similar genes',
                    'detailed': '''CLUSTER DETECTION
                    Identify groups of genes with similar codon usage patterns:
                    
                    K-means: For spherical, well-separated clusters
                    - Auto-detects optimal K using silhouette score
                    
                    DBSCAN: For arbitrary shapes and outlier detection
                    - Density-based, finds natural groupings
                    
                    Useful for: Functional classification, expression classes'''
                },
                '3': {
                    'brief': 'Identify genes with unusual patterns',
                    'detailed': '''OUTLIER ANALYSIS
                    Detect abnormal genes using multiple methods:
                    
                    Statistical: Z-score based on single metrics
                    Multivariate: Mahalanobis distance in PC space
                    Domain-specific: Extreme codon bias patterns
                    Combined: Consensus outliers from all methods
                    
                    Applications: Quality control, horizontal transfer detection'''
                },
                '4': {
                    'brief': 'Export genes from distribution extremes',
                    'detailed': '''EXPORT EXTREME GENES
                    Extract genes from the tails of any metric distribution:
                    - Top/bottom X% by GC, ENC, CAI, etc.
                    - Multivariate axis extremes
                    - Customizable thresholds
                    
                    Useful for: Creating reference sets, validation studies'''
                },
                '5': {
                    'brief': 'Statistical comparison of gene groups',
                    'detailed': '''COMPARE AXIS COHORTS
                    Compares codon usage between genes at opposite multivariate extremes.
                    Statistical tests identify significantly different codon usage.
                    
                    Often reveals expression-related codon preferences.'''
                }
            },
            # Data management menu help
            'data_management': {
                '1': {
                    'brief': 'Load sequences from a FASTA file',
                    'detailed': '''LOAD FASTA FILE
                    Load protein-coding DNA sequences for analysis.
                    Requirements:
                    - Standard FASTA format
                    - DNA sequences (A, T, G, C)
                    - Length divisible by 3 (complete codons)
                    
                    Large files automatically use memory-efficient loading.'''
                },
                '2': {
                    'brief': 'View summary of loaded data',
                    'detailed': '''VIEW DATA SUMMARY
                    Display statistics about your loaded sequences:
                    - Number of sequences
                    - Length distribution
                    - Basic composition
                    - Memory usage
                    
                    Useful for verifying data loaded correctly.'''
                },
                '3': {
                    'brief': 'Filter sequences by various criteria',
                    'detailed': '''FILTER SEQUENCES
                    Remove sequences based on:
                    - Length (min/max)
                    - GC content
                    - Ambiguous nucleotides
                    - Custom criteria
                    
                    Creates a subset for focused analysis.'''
                },
                '4': {
                    'brief': 'Save current dataset',
                    'detailed': '''SAVE CURRENT DATA
                    Export the current working dataset:
                    - After filtering
                    - With calculated metrics
                    - In FASTA format
                    
                    Preserves your preprocessed data.'''
                }
            },
            # Visualization & Export menu help
            'visualization_export': {
                '1': {
                    'brief': 'Interactive scatter plots of multivariate results',
                    'detailed': '''MULTIVARIATE PLOTS
                    Create interactive visualizations of CA/PCA results:
                    - Gene positions in reduced dimensional space
                    - Biplot showing codon/AA contributions (CA only)
                    - Outlier highlighting
                    - Hover for gene details
                    
                    Requires multivariate analysis to be run first.'''
                },
                '2': {
                    'brief': 'Visualize codon preferences across genes',
                    'detailed': '''RSCU HEATMAP
                    Color-coded visualization of RSCU values:
                    - Small datasets: Full gene × codon matrix
                    - Large datasets: Summary statistics
                    - Red = underused, Blue = overused
                    - Grouped by amino acid
                    
                    Quickly identify codon preferences.'''
                },
                '3': {
                    'brief': 'Various bias and composition plots',
                    'detailed': '''OTHER PLOTS
                    Additional visualization options:
                    - ENC vs GC3s (Wright plot)
                    - GC vs GC3 scatter
                    - CAI distribution histogram
                    - Custom scatter plots
                    
                    Publication-quality outputs in multiple formats.'''
                },
                '4': {
                    'brief': 'Save analysis results as data files',
                    'detailed': '''EXPORT DATA
                    Export calculated metrics in TSV format:
                    - Comprehensive metrics (all in one file)
                    - Individual analyses
                    - Optimal codon tables
                    - Multivariate results
                    
                    Compatible with Excel, R, Python, etc.'''
                }
            },
            # Sequence optimization menu help
            'sequence_optimization': {
                '1': {
                    'brief': 'Optimize individual sequences',
                    'detailed': '''SINGLE GENE OPTIMIZATION
                    Replace codons in one gene with optimal alternatives:
                    - Select gene by name or number
                    - Uses current optimal codon set
                    - Shows before/after comparison
                    - Maintains amino acid sequence
                    
                    Perfect for optimizing specific genes of interest.'''
                },
                '2': {
                    'brief': 'Optimize all sequences at once',
                    'detailed': '''BATCH OPTIMIZATION
                    Optimize entire dataset in one operation:
                    - Uses current optimal codon set
                    - Processes all genes
                    - Saves to single FASTA file
                    - Progress tracking for large datasets
                    
                    Efficient for whole-genome optimization.'''
                },
                '3': {
                    'brief': 'Choose optimal codon identification method',
                    'detailed': '''SELECT OPTIMIZATION METHOD
                    Define which codons are "optimal":
                    - Frequency-based (most common)
                    - Multivariate extremes (expression-correlated)
                    - Reference gene set
                    - Custom codon table
                    
                    Critical for optimization success.'''
                },
                '4': {
                    'brief': 'Import optimal codons from file',
                    'detailed': '''LOAD CODON TABLE
                    Use a pre-defined optimal codon set:
                    - TSV or JSON format
                    - From previous analysis
                    - From literature
                    - Species-specific tables
                    
                    Ensures reproducible optimization.'''
                }
            }
        }

    def display_banner(self):
        """Display program banner with a clean, minimalist layout."""
        banner = f"""
    *******************************************************
    GCUA: General Codon Usage Analysis {VERSION}
    by {AUTHOR}

    Cite: McInerney JO. GCUA: general codon usage analysis.
    Bioinformatics. 1998;14(4):372-3.
    doi: 10.1093/bioinformatics/14.4.372
    *******************************************************

    Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}
    """
        print(banner)

    def main_menu(self):
        """Display main menu and handle user inputs."""
        # Add breadcrumb tracking
        self.menu_breadcrumb = ["Main Menu"]
        
        menu_options = [
            ("Data Management", "1", self._data_management_menu),
            ("Quick Analysis (Guided Workflow)", "2", self._quick_analysis_workflow),
            ("Custom Analysis", "3", self._analysis_menu),
            ("Visualization & Export", "4", self._visualization_export_menu),
            ("Advanced Tools", "5", self._advanced_tools_menu),
            ("Settings & Preferences", "6", self._preferences_menu),
            ("Help & Documentation", "7", self._help),
            ("Quit program", "Q", lambda: False)
        ]

        while True:
            self._clear_screen()
            self.display_banner()

            # Show data status
            if self.data.gene_names:
                self._display_data_status()
            else:
                print("\nNo data loaded.")
                print(
                    f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

            # Show analysis status if data is loaded
            if self.data.gene_names:
                self._display_analysis_status()

            # Display menu with improved formatting
            print("\n+----------------------+")
            print("|      MAIN MENU       |")
            print("+----------------------+")

            # Calculate padding for alignment
            max_option_length = max(len(option)
                                    for option, _, _ in menu_options)

            for option, key, _ in menu_options:
                # Check if option is available
                available = self._is_option_available(option)
                if not available:
                    print(f"  {key}. {option} (requires data)")
                else:
                    padding = " " * (max_option_length - len(option))
                    print(f"  {key}. {option}{padding}")
            
            print("\n  H. Show help for all menu items")
            
            # Show help hint on first menu display
            if self.show_help_hint:
                print("\n💡 Tip: Type a number followed by '?' (e.g., '1?') for help on that option")
                self.show_help_hint = False

            # Show recent actions
            if hasattr(self, 'recent_actions') and self.recent_actions:
                print("\nRecent actions:")
                for i, action in enumerate(self.recent_actions[-3:], 1):
                    print(f"  • {action}")

            choice = input("\nEnter your choice: ").strip().upper()

            # Check for help requests
            if choice == 'H':
                self._show_menu_help('main_menu', menu_options)
                self._pause()
                continue
            elif choice.endswith('?'):
                option_key = choice[:-1]
                self._show_quick_help('main_menu', option_key)
                self._pause()
                continue

            for option, key, handler in menu_options:
                if choice == key:
                    if key == "Q":  # only exit program on Q
                        return
                    if self._is_option_available(option):
                        # Track action
                        if not hasattr(self, 'recent_actions'):
                            self.recent_actions = []
                        self.recent_actions.append(option)
                        handler()
                    else:
                        print("\nThis option requires data to be loaded first.")
                        self._pause()
                    break
            else:
                if choice != 'Q':  # 'Q' is already handled in the loop
                    print("\nUnrecognized command.")
                    self._pause()

    def _display_data_status(self):
        """Display current data status."""
        print(f"\nCurrent data: {len(self.data.gene_names)} sequences")
        print(f"Source: {self.data.file_path}")
        print(f"Genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        if self.data.memory_mode == "low_memory":
            print(f"Memory mode: Low memory")

    def _display_analysis_status(self):
        """Display analysis completion status."""
        print("\nAnalysis Status:")
        
        # Group statuses
        basic_complete = all([
            self.data.analysis_status.get('codon_usage', False),
            self.data.analysis_status.get('amino_acid_usage', False),
            self.data.analysis_status.get('base_composition', False)
        ])
        
        bias_complete = all([
            self.data.analysis_status.get('enc_calculated', False),
            self.data.analysis_status.get('cai_calculated', False),
            self.data.analysis_status.get('fop_calculated', False),
            self.data.analysis_status.get('scuo_calculated', False)
        ])
        
        # Display grouped status
        print(f"  {'✓' if basic_complete else '✗'} Basic metrics " + 
              f"{'(cached)' if basic_complete and hasattr(self.data, 'codon_counts') and not self.data.codon_counts.empty else ''}")
        print(f"  {'✓' if self.data.analysis_status.get('rscu_calculated', False) else '✗'} RSCU values")
        print(f"  {'✓' if bias_complete else '✗'} Codon bias metrics (ENC, CAI, Fop, SCUO)")
        print(f"  {'✓' if self.data.analysis_status.get('multivariate_performed', False) else '✗'} Multivariate analysis" +
              f" ({len(getattr(self.data, 'multivariate_history', []))} performed)" if self.data.analysis_status.get('multivariate_performed', False) else "")
        print(f"  {'✓' if self.data.analysis_status.get('clusters_detected', False) else '✗'} Cluster detection")
        print(f"  {'✓' if self.data.analysis_status.get('optimal_codons', False) else '✗'} Optimal codons identified")

    def _is_option_available(self, option_name):
        """Check if a menu option is available based on current state."""
        # Options that always work
        always_available = ["Settings & Preferences", "Help & Documentation", "Quit program", "Data Management"]
        
        if option_name in always_available:
            return True
            
        # All other options require data
        return bool(self.data.gene_names)

    def _display_menu(self, title, options, menu_key=None):
        """
        Enhanced menu display with a cleaner, more minimalist layout and help system.
        """
        while True:
            self._clear_screen()

            # Title with underline
            print(title)
            print("=" * len(title))

            # Show data status if applicable
            if self.data.gene_names:
                print(f"\nCurrent Data: {len(self.data.gene_names)} sequences")
                print(f"Source: {self.data.file_path or 'N/A'}")
                if self.data.memory_mode == "low_memory":
                    print(f"Memory mode: Low memory")

            print("\nOptions:")
            # Calculate padding for alignment
            max_option_length = max(len(f"{option_key}. {option_text}") for option_text, option_key, _ in options)

            for option_text, option_key, _ in options:
                # Right-align the key and add padding to align all options
                full_option = f"{option_key}. {option_text}"
                print(full_option)

            print("\nR. Return to previous menu")
            print("H. Show help for all menu items")
            
            # Show help hint on first menu display
            if self.show_help_hint:
                print("\n💡 Tip: Type a number followed by '?' (e.g., '1?') for help on that option")
                self.show_help_hint = False

            choice = input("\nSelect an option: ").strip().upper()

            # Check for help requests
            if choice == 'H':
                self._show_menu_help(menu_key, options)
                self._pause()
                continue
            elif choice.endswith('?'):
                option_key = choice[:-1]
                self._show_quick_help(menu_key, option_key)
                self._pause()
                continue

            for option_text, option_key, handler in options:
                if choice == option_key:
                    result = handler()
                    if result is not None and not result:
                        return False
                    break
            else:
                if choice == 'R':
                    return False
                else:
                    print("\nInvalid option. Please try again.")
                    self._pause()

    def _export_multivariate(self, result=None, analysis_type=None, data_type=None):
        """Export multivariate analysis results to a file."""
        if result is None:
            # Check for multivariate history
            if not hasattr(self.data, 'multivariate_history') or not self.data.multivariate_history:
                print("\nNo multivariate analysis results found.")
                print("Please perform a multivariate analysis first.")
                self._pause()
                return True
            
            # If only one analysis, use it
            if len(self.data.multivariate_history) == 1:
                analysis = self.data.multivariate_history[0]
                result = analysis['result']
                analysis_type = analysis['analysis_type']
                data_type = analysis['data_type']
            else:
                # Multiple analyses, let user choose
                print("\nAvailable analyses to export:")
                for i, analysis in enumerate(self.data.multivariate_history, 1):
                    print(f"{i}. {analysis['analysis_type']} of {analysis['data_type']} data")
                
                choice = input("\nWhich analysis would you like to export? (number): ").strip()
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(self.data.multivariate_history):
                        analysis = self.data.multivariate_history[idx]
                        result = analysis['result']
                        analysis_type = analysis['analysis_type']
                        data_type = analysis['data_type']
                    else:
                        print("Invalid choice.")
                        self._pause()
                        return True
                except ValueError:
                    print("Invalid input.")
                    self._pause()
                    return True

        if not result:
            print("No results available to export.")
            self._pause()
            return True

        # Determine analysis type and data type if not provided
        if analysis_type is None:
            analysis_type = result.get('analysis_type', 'CA')
        if data_type is None:
            data_type = result.get('data_type', 'RSCU')

        # Get output path
        output_path = self.file_manager.get_output_path(
            self.data, None, f"{analysis_type}_{data_type}_results") + '.tsv'

        # Extract coordinates and explained variance
        coords = result['coordinates']
        explained_var = result['explained_variance']

        with open(output_path, 'w') as f:
            # Add metadata
            f.write(f"# GCUA {analysis_type} of {data_type} Results\n")
            f.write(f"# Version: {VERSION}\n")
            f.write(f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
            f.write(f"# Total Genes: {len(self.data.gene_names)}\n")
            f.write(f"# Genes in Analysis: {len(coords)}\n")
            f.write(f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            
            # Add preprocessing information if available
            if hasattr(self.data, 'last_preprocessing_params'):
                f.write("\n# Preprocessing Parameters:\n")
                params = self.data.last_preprocessing_params
                if 'outliers_removed' in params:
                    f.write(f"# Outliers removed: {params['outliers_removed']}\n")
                if 'transformation' in params:
                    f.write(f"# Transformation: {params['transformation']}\n")
                if 'pseudocount' in params:
                    f.write(f"# Pseudocount: {params['pseudocount']}\n")

            # Write explained variance with cumulative values
            f.write("\n# Explained variance by dimension:\n")
            cumulative_var = 0
            for i, var in enumerate(explained_var):
                cumulative_var += var
                f.write(f"# {coords.columns[i]}: {var:.4f} ({var*100:.2f}%) - Cumulative: {cumulative_var*100:.2f}%\n")
            
            # Add summary statistics
            f.write(f"\n# Summary Statistics:\n")
            f.write(f"# Total variance explained (first 2 axes): {sum(explained_var[:2])*100:.2f}%\n")
            if len(explained_var) >= 3:
                f.write(f"# Total variance explained (first 3 axes): {sum(explained_var[:3])*100:.2f}%\n")
            
            # Find number of components for 80% variance
            cumulative = np.cumsum(explained_var)
            if any(cumulative >= 0.8):
                components_80 = np.argmax(cumulative >= 0.8) + 1
                f.write(f"# Components needed for 80% variance: {components_80}\n")
            else:
                # If we can't reach 80% with available components, report what we can achieve
                max_var = cumulative[-1] * 100
                f.write(f"# Components needed for 80% variance: >{len(explained_var)} (max available: {max_var:.1f}%)\n")
            
            f.write("#" + "-" * 50 + "\n\n")

            # Write coordinates
            coords.to_csv(f, sep='\t', float_format='%.4f')

        print(f"\nMultivariate analysis results saved to {output_path}")

        # Ask if user wants to export loadings as well
        if 'loadings' in result:
            export_loadings = input("\nDo you want to export loadings as well? (y/n): ").strip().lower()
            if export_loadings == 'y':
                loadings_path = output_path.replace('.tsv', '_loadings.tsv')
                with open(loadings_path, 'w') as f:
                    f.write(f"# GCUA {analysis_type} of {data_type} Loadings\n")
                    f.write(f"# Version: {VERSION}\n")
                    f.write(f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
                    f.write("#" + "-" * 50 + "\n\n")
                    result['loadings'].to_csv(f, sep='\t', float_format='%.4f')
                print(f"Loadings saved to {loadings_path}")

        # Ask if user wants to export scree plot data
        export_scree = input("\nDo you want to export scree plot data? (y/n): ").strip().lower()
        if export_scree == 'y':
            scree_path = output_path.replace('.tsv', '_scree.tsv')
            with open(scree_path, 'w') as f:
                f.write(f"# GCUA {analysis_type} of {data_type} Scree Plot Data\n")
                f.write(f"# Version: {VERSION}\n")
                f.write(f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("#" + "-" * 50 + "\n\n")
                
                # Create scree plot dataframe
                scree_df = pd.DataFrame({
                    'Component': [f'{coords.columns[i]}' for i in range(len(explained_var))],
                    'Component_Number': list(range(1, len(explained_var) + 1)),
                    'Variance_Explained': explained_var,
                    'Variance_Percent': explained_var * 100,
                    'Cumulative_Variance': np.cumsum(explained_var),
                    'Cumulative_Percent': np.cumsum(explained_var) * 100
                })
                scree_df.to_csv(f, sep='\t', index=False, float_format='%.4f')
            print(f"Scree plot data saved to {scree_path}")

        self._pause()
        return True

    def _input_fasta(self):
        """Handle FASTA file input."""
        file_path = input("\nName of the FASTA-formatted file: ").strip()
        if file_path:
            self.data.load_fasta(file_path)
            self._pause()
        return True

    def _data_management_menu(self):
        """Handle data loading and management."""
        options = [
            ("Load FASTA file", "1", self._input_fasta),
            ("View data summary", "2", self._view_data_summary),
            ("Filter sequences", "3", self._filter_sequences),
            ("Clear current data", "4", self._clear_data),
            ("Return to main menu", "R", lambda: False)
        ]
        
        return self._display_menu("Data Management", options, menu_key='data_management')
    
    def _quick_analysis_workflow(self):
        """Guided workflow for common analysis tasks."""
        if not self._check_data():
            return True
            
        print("\nQuick Analysis Workflow")
        print("=" * 30)
        print("\nThis will guide you through a standard analysis:")
        print("1. Calculate basic metrics")
        print("2. Perform multivariate analysis")
        print("3. Identify optimal codons")
        print("4. Create visualizations")
        print("5. Export results")
        
        proceed = input("\nProceed with workflow? (y/n): ").strip().lower()
        if proceed != 'y':
            return True
            
        # Step 1: Basic metrics
        print("\n[Step 1/5] Calculating basic metrics...")
        print("  - Base composition: Already calculated")
        print("  - Codon usage: Already calculated")
        print("  - Amino acid usage: Already calculated")
        print("  - Calculating RSCU values...")
        self.data.get_rscu_values()  # This will calculate if needed
        
        # Update status
        self.data.analysis_status['basic_metrics'] = True
        self.data.analysis_status['codon_usage'] = True
        self.data.analysis_status['amino_acid_usage'] = True
        self.data.analysis_status['base_composition'] = True
        self.data.analysis_status['rscu_calculated'] = True
        
        # Step 2: Multivariate analysis
        print("\n[Step 2/5] Performing multivariate analysis...")
        self.data.perform_multivariate_analysis()
        self.data.analysis_status['multivariate_performed'] = True
        
        # Step 3: Optimal codons
        print("\n[Step 3/5] Identifying optimal codons...")
        self.data.calculate_optimal_codons()
        self.data.analysis_status['optimal_codons'] = True
        
        # Also calculate ENC, CAI, Fop for comprehensive analysis
        print("  - Calculating ENC values...")
        self.data.calculate_enc()
        self.data.analysis_status['enc_calculated'] = True
        
        print("  - Calculating CAI values...")
        self.data.calculate_cai()
        self.data.analysis_status['cai_calculated'] = True
        
        print("  - Calculating Fop values...")
        self.data.calculate_fop()
        self.data.analysis_status['fop_calculated'] = True
        
        # Step 4: Visualizations
        print("\n[Step 4/5] Creating visualizations...")
        self._create_standard_plots()
        
        # Step 5: Export
        print("\n[Step 5/5] Exporting results...")
        self._export_workflow_results()
        
        print("\nWorkflow complete! Check gcua_outputs folder for results.")
        self._pause()
        return True
    
    def _visualization_export_menu(self):
        """Unified menu for visualization and export."""
        if not self._check_data():
            return True
            
        options = [
            ("Quick Export All Results", "1", self._quick_export_all),
            ("Visualizations", "2", self._visualization_submenu),
            ("Export Data Tables", "3", self._export_submenu),
            ("Generate Report", "4", self._generate_report),
            ("Return to main menu", "R", lambda: False)
        ]
        
        return self._display_menu("Visualization & Export", options, menu_key='visualization_export')
    
    def _advanced_tools_menu(self):
        """Menu for advanced analysis tools."""
        if not self._check_data():
            return True
            
        options = [
            ("Sequence Optimization", "1", self._sequence_optimization_menu),
            ("Cluster Detection", "2", self._detect_clusters_menu),
            ("Outlier Analysis", "3", self._outlier_analysis_menu),
            ("Export Extreme Genes", "4", self._export_extreme_genes_menu),
            ("Compare Axis Cohorts", "5", self._compare_codon_usage_axis_cohorts),
            ("Return to main menu", "R", lambda: False)
        ]
        
        return self._display_menu("Advanced Tools", options, menu_key='advanced_tools')
    
    def _view_data_summary(self):
        """Display a summary of loaded data."""
        if not self._check_data():
            return True
            
        print("\nData Summary")
        print("=" * 40)
        print(f"File: {self.data.file_path}")
        print(f"Total sequences: {len(self.data.gene_names)}")
        print(f"Genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        
        if hasattr(self.data, 'base_composition') and not self.data.base_composition.empty:
            print(f"\nSequence length statistics:")
            lengths = self.data.base_composition['Length']
            print(f"  Min: {lengths.min()} bp")
            print(f"  Max: {lengths.max()} bp")
            print(f"  Mean: {lengths.mean():.1f} bp")
            print(f"  Median: {lengths.median():.1f} bp")
            
        print("\nFirst 10 sequences:")
        for i, gene in enumerate(self.data.gene_names[:10]):
            print(f"  {i+1}. {gene}")
            
        if len(self.data.gene_names) > 10:
            print(f"  ... and {len(self.data.gene_names) - 10} more")
            
        self._pause()
        return True
    
    def _filter_sequences(self):
        """Filter sequences based on criteria."""
        if not self._check_data():
            return True
            
        print("\nFilter Sequences")
        print("=" * 40)
        print("1. Filter by length")
        print("2. Filter by GC content")
        print("3. Filter by gene name pattern")
        print("R. Return")
        
        choice = input("\nSelect filter type: ").strip().upper()
        
        if choice == "1":
            min_len = input("Minimum length (press Enter to skip): ").strip()
            max_len = input("Maximum length (press Enter to skip): ").strip()
            # Implementation would go here
            print("\nLength filtering not yet implemented.")
        elif choice == "2":
            min_gc = input("Minimum GC% (press Enter to skip): ").strip()
            max_gc = input("Maximum GC% (press Enter to skip): ").strip()
            # Implementation would go here
            print("\nGC filtering not yet implemented.")
        elif choice == "3":
            pattern = input("Gene name pattern (regex): ").strip()
            # Implementation would go here
            print("\nPattern filtering not yet implemented.")
            
        self._pause()
        return True
    
    def _clear_data(self):
        """Clear all loaded data."""
        if not self.data.gene_names:
            print("\nNo data loaded.")
            self._pause()
            return True
            
        confirm = input("\nAre you sure you want to clear all data? (y/n): ").strip().lower()
        if confirm == 'y':
            # Reset data
            self.data = GCUAData(self.config)
            print("\nAll data cleared.")
        else:
            print("\nOperation cancelled.")
            
        self._pause()
        return True
    
    def _create_standard_plots(self):
        """Create standard visualization plots non-interactively."""
        print("\nCreating standard plots...")
        
        # Multivariate plot - use the most recent analysis
        if hasattr(self.data, 'multivariate_history') and self.data.multivariate_history:
            print("  - Creating multivariate analysis plot...")
            analysis = self.data.multivariate_history[-1]  # Use most recent
            result = analysis['result']
            output_path = self.file_manager.get_output_path(
                self.data, None, f"multivariate_{analysis['analysis_type']}_{analysis['data_type']}", 
                f".{self.config.visualization_format}")
            self.data.visualize('multivariate', result, output_path)
            
        # GC content plot
        print("  - Creating GC vs GC3 plot...")
        gc_output_path = self.file_manager.get_output_path(
            self.data, None, "gc_content_plot", f".{self.config.visualization_format}")
        self.data.visualize('gc_content', None, gc_output_path)
        
        # ENC plot if available
        if hasattr(self.data, 'enc_values') and self.data.enc_values:
            print("  - Creating ENC vs GC3s plot...")
            enc_output_path = self.file_manager.get_output_path(
                self.data, None, "enc_plot", f".{self.config.visualization_format}")
            self.data.visualize('enc', None, enc_output_path)
            
        print("\nPlots created in gcua_outputs folder.")
        
    def _export_workflow_results(self):
        """Export all results from workflow."""
        print("\nExporting all results...")
        
        # Export multivariate results
        if hasattr(self.data, 'multivariate_results') and self.data.multivariate_results:
            self._export_multivariate()
            
        # Export metrics
        self._export_comprehensive_metrics()
        
        # Export optimal codons
        if hasattr(self.data, 'optimal_codons') and self.data.optimal_codons:
            output_path = self.file_manager.get_output_path(
                self.data, None, "optimal_codons") + ".tsv"
            self.data.save_optimal_codons_to_file(output_path, self.data.optimal_codons)
            
        print("\nAll results exported to gcua_outputs folder.")
        
    def _quick_export_all(self):
        """Quick export of all available results."""
        print("\nQuick Export - All Results")
        print("=" * 40)
        
        exported = []
        
        # Export basic metrics
        if hasattr(self.data, 'base_composition') and not self.data.base_composition.empty:
            self._export_base_composition()
            exported.append("Base composition")
            
        # Export RSCU
        if hasattr(self.data, 'rscu_values') and not self.data.rscu_values.empty:
            self._export_rscu()
            exported.append("RSCU values")
            
        # Export multivariate
        if hasattr(self.data, 'multivariate_history') and self.data.multivariate_history:
            for analysis in self.data.multivariate_history:
                self._export_multivariate(analysis['result'], 
                                        analysis['analysis_type'], 
                                        analysis['data_type'])
            exported.append(f"Multivariate results ({len(self.data.multivariate_history)} analyses)")
            
        # Export comprehensive metrics
        if any(hasattr(self.data, attr) for attr in ['enc_values', 'cai_values', 'fop_values']):
            self._export_comprehensive_metrics()
            exported.append("Comprehensive metrics")
            
        if exported:
            print("\nExported:")
            for item in exported:
                print(f"  ✓ {item}")
        else:
            print("\nNo results available to export.")
            print("Please run analyses first.")
            
        self._pause()
        return True
    
    def _visualization_submenu(self):
        """Submenu for visualizations."""
        return self._visualization_menu()
    
    def _export_submenu(self):
        """Submenu for data exports."""
        return self._export_menu()
    
    def _generate_report(self):
        """Generate a comprehensive analysis report."""
        print("\nGenerating comprehensive report...")
        
        # This would generate a markdown or HTML report
        print("\nReport generation not yet implemented.")
        print("This feature will create a complete analysis report with:")
        print("  - Data summary")
        print("  - All calculated metrics")
        print("  - Key findings")
        print("  - Visualizations")
        
        self._pause()
        return True
    
    def _sequence_optimization_menu(self):
        """Handle sequence optimization options."""
        if not self._check_data():
            return True
            
        if not hasattr(self.data, 'optimal_codons') or not self.data.optimal_codons:
            print("\nOptimal codons not yet identified.")
            print("Please identify optimal codons first (Analysis -> Identify optimal codons).")
            self._pause()
            return True
            
        options = [
            ("Optimize a single gene", "1", self._optimize_single_gene),
            ("Optimize all genes", "2", self._optimize_all_genes),
            ("Compare optimized sequences", "3", self._compare_optimized_sequences),
            ("Export optimized sequences", "4", self._export_optimized_sequences),
            ("Return to advanced tools", "R", lambda: False)
        ]
        
        return self._display_menu("Sequence Optimization", options, menu_key='sequence_optimization')
    
    def _outlier_analysis_menu(self):
        """Handle outlier analysis options."""
        if not self._check_data():
            return True
            
        print("\nOutlier Analysis")
        print("=" * 40)
        
        # Check if multivariate analysis has been performed
        has_multivariate = hasattr(self.data, 'multivariate_results') and self.data.multivariate_results
        
        print("\nSelect outlier detection method:")
        print("1. Statistical outliers (Z-score based)")
        print("2. Multivariate outliers (Mahalanobis distance)" + (" [Requires multivariate analysis]" if not has_multivariate else ""))
        print("3. Domain-specific outliers (extreme codon bias)")
        print("4. Combined approach (all methods)")
        print("R. Return")
        
        choice = input("\nSelect method (1-4, R): ").strip().upper()
        
        if choice == 'R':
            return True
        elif choice == '1':
            self._detect_statistical_outliers()
        elif choice == '2':
            if not has_multivariate:
                print("\nMultivariate analysis required. Please perform CA or PCA first.")
                self._pause()
            else:
                self._detect_multivariate_outliers()
        elif choice == '3':
            self._detect_domain_outliers()
        elif choice == '4':
            self._detect_combined_outliers()
        else:
            print("\nInvalid choice.")
            self._pause()
            
        return True
    
    def _detect_statistical_outliers(self):
        """Detect outliers using Z-score method on various metrics."""
        print("\nStatistical Outlier Detection")
        print("=" * 40)
        
        # Select metric
        print("\nSelect metric for outlier detection:")
        print("1. GC content")
        print("2. GC3s content")
        print("3. Gene length")
        print("4. ENC values" + (" [Not calculated]" if not hasattr(self.data, 'enc_values') or not self.data.enc_values else ""))
        print("5. CAI values" + (" [Not calculated]" if not hasattr(self.data, 'cai_values') or not self.data.cai_values else ""))
        
        metric_choice = input("\nSelect metric (1-5): ").strip()
        
        # Z-score threshold
        threshold_input = input("\nZ-score threshold (default=3): ").strip()
        threshold = float(threshold_input) if threshold_input else 3.0
        
        # Get the data based on choice
        if metric_choice == '1':
            data = self.data.base_composition['GC']
            metric_name = "GC content"
        elif metric_choice == '2':
            data = self.data.base_composition['GC3s']
            metric_name = "GC3s content"
        elif metric_choice == '3':
            data = self.data.base_composition['Length']
            metric_name = "gene length"
        elif metric_choice == '4' and hasattr(self.data, 'enc_values') and self.data.enc_values:
            data = pd.Series(self.data.enc_values)
            metric_name = "ENC"
        elif metric_choice == '5' and hasattr(self.data, 'cai_values') and self.data.cai_values:
            data = pd.Series(self.data.cai_values)
            metric_name = "CAI"
        else:
            print("\nInvalid choice or metric not available.")
            self._pause()
            return
        
        # Calculate Z-scores
        z_scores = np.abs(stats.zscore(data))
        outliers = data[z_scores > threshold]
        
        print(f"\nFound {len(outliers)} outliers (|Z-score| > {threshold}) for {metric_name}:")
        print(f"\nTop 10 outliers:")
        print(f"{'Gene':30} {metric_name:15} Z-score")
        print("-" * 60)
        
        # Sort by absolute Z-score
        outlier_z = z_scores[z_scores > threshold].sort_values(ascending=False)
        for i, (gene, z) in enumerate(outlier_z.head(10).items()):
            value = data[gene]
            print(f"{gene:30} {value:15.3f} {z:7.3f}")
        
        if len(outliers) > 10:
            print(f"\n... and {len(outliers) - 10} more outliers")
        
        # Ask if user wants to export
        export = input("\nExport outlier list? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, f"outliers_{metric_name.replace(' ', '_')}") + '.tsv'
            
            with open(output_path, 'w') as f:
                f.write(f"# Statistical outliers for {metric_name}\n")
                f.write(f"# Z-score threshold: {threshold}\n")
                f.write(f"# Total outliers: {len(outliers)}\n")
                f.write(f"#{'=' * 50}\n\n")
                f.write(f"Gene\t{metric_name}\tZ_score\n")
                
                for gene in outlier_z.index:
                    f.write(f"{gene}\t{data[gene]:.4f}\t{outlier_z[gene]:.4f}\n")
            
            print(f"\nOutliers saved to {output_path}")
        
        self._pause()
    
    def _detect_multivariate_outliers(self):
        """Detect outliers in multivariate space using Mahalanobis distance."""
        print("\nMultivariate Outlier Detection")
        print("=" * 40)
        
        # Get multivariate coordinates
        coords = self.data.multivariate_results['coordinates']
        
        # Select number of dimensions to use
        max_dims = min(10, coords.shape[1])
        dims_input = input(f"\nNumber of dimensions to use (1-{max_dims}, default=2): ").strip()
        n_dims = int(dims_input) if dims_input else 2
        n_dims = min(max(1, n_dims), max_dims)
        
        # Get data for selected dimensions
        X = coords.iloc[:, :n_dims].values
        
        # Calculate Mahalanobis distances
        print("\nCalculating Mahalanobis distances...")
        
        # Calculate mean and covariance
        mean = np.mean(X, axis=0)
        cov = np.cov(X.T)
        
        # Calculate Mahalanobis distance for each point
        inv_cov = np.linalg.inv(cov)
        distances = []
        
        for i in range(len(X)):
            diff = X[i] - mean
            dist = np.sqrt(diff.dot(inv_cov).dot(diff))
            distances.append(dist)
        
        distances = pd.Series(distances, index=coords.index)
        
        # Determine threshold using chi-square distribution
        # For multivariate normal data, Mahalanobis distance squared follows chi-square distribution
        threshold_pct = 99.5  # Default to 99.5th percentile
        threshold = np.sqrt(chi2.ppf(threshold_pct/100, n_dims))
        
        outliers = distances[distances > threshold]
        
        print(f"\nFound {len(outliers)} multivariate outliers (distance > {threshold:.2f}):")
        print(f"\nTop 10 outliers:")
        print(f"{'Gene':30} {'Mahalanobis Distance':20}")
        print("-" * 55)
        
        for gene, dist in outliers.nlargest(10).items():
            print(f"{gene:30} {dist:20.3f}")
        
        if len(outliers) > 10:
            print(f"\n... and {len(outliers) - 10} more outliers")
        
        # Visualize outliers
        vis = input("\nVisualize outliers on multivariate plot? (y/n): ").strip().lower()
        if vis == 'y':
            self._visualize_multivariate_outliers(outliers.index.tolist(), coords, n_dims)
        
        # Export outliers
        export = input("\nExport outlier list? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, "multivariate_outliers") + '.tsv'
            
            with open(output_path, 'w') as f:
                f.write(f"# Multivariate outliers\n")
                f.write(f"# Dimensions used: {n_dims}\n")
                f.write(f"# Threshold: {threshold:.3f} (chi-square {threshold_pct}th percentile)\n")
                f.write(f"# Total outliers: {len(outliers)}\n")
                f.write(f"#{'=' * 50}\n\n")
                f.write("Gene\tMahalanobis_Distance\n")
                
                for gene, dist in outliers.sort_values(ascending=False).items():
                    f.write(f"{gene}\t{dist:.4f}\n")
            
            print(f"\nOutliers saved to {output_path}")
        
        self._pause()
    
    def _detect_domain_outliers(self):
        """Detect outliers based on extreme codon usage patterns."""
        print("\nDomain-specific Outlier Detection")
        print("=" * 40)
        
        # Calculate RSCU if needed
        if not hasattr(self.data, 'rscu_values') or self.data.rscu_values.empty:
            print("\nCalculating RSCU values...")
            self.data.get_rscu_values()  # This will calculate if needed
        
        print("\nDetecting genes with extreme codon usage patterns...")
        
        # Method 1: Genes with many extremely biased codons
        rscu_values = self.data.rscu_values
        
        # Count codons with extreme RSCU (>2 or <0.5) per gene
        extreme_high = (rscu_values > 2).sum(axis=1)
        extreme_low = (rscu_values < 0.5).sum(axis=1)
        extreme_total = extreme_high + extreme_low
        
        # Find outliers (genes with many extreme codons)
        threshold_pct = 95
        threshold = np.percentile(extreme_total, threshold_pct)
        outliers = extreme_total[extreme_total > threshold]
        
        print(f"\nFound {len(outliers)} genes with extreme codon bias:")
        print(f"(>{int(threshold)} codons with RSCU >2 or <0.5)")
        print(f"\nTop 10 genes:")
        print(f"{'Gene':30} {'Extreme codons':15} {'High RSCU':12} {'Low RSCU':12}")
        print("-" * 75)
        
        for gene in outliers.nlargest(10).index:
            print(f"{gene:30} {extreme_total[gene]:15} {extreme_high[gene]:12} {extreme_low[gene]:12}")
        
        # Method 2: Check for stop codon readthrough
        print("\n\nChecking for potential stop codon readthrough...")
        stop_codons = ['UAA', 'UAG', 'UGA']
        aa_map = self.config.get_aa_map()
        
        # Find which stop codons are actually stops in this genetic code
        true_stops = [codon for codon in stop_codons if aa_map.get(codon) == 'STOP']
        
        if true_stops:
            # Check if any genes use stop codons (shouldn't happen in normal CDS)
            stop_usage = self.data.codon_counts[true_stops].sum(axis=1)
            readthrough_genes = stop_usage[stop_usage > 0]
            
            if len(readthrough_genes) > 0:
                print(f"\nWarning: Found {len(readthrough_genes)} genes with stop codons!")
                print("These may be pseudogenes or have readthrough:")
                for gene, count in readthrough_genes.head(10).items():
                    print(f"  {gene}: {int(count)} stop codons")
        
        # Export results
        export = input("\nExport outlier list? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, "domain_outliers") + '.tsv'
            
            with open(output_path, 'w') as f:
                f.write(f"# Domain-specific outliers\n")
                f.write(f"# Threshold: {threshold_pct}th percentile\n")
                f.write(f"# Total outliers: {len(outliers)}\n")
                f.write(f"#{'=' * 50}\n\n")
                f.write("Gene\tExtreme_codons\tHigh_RSCU\tLow_RSCU\n")
                
                for gene in outliers.sort_values(ascending=False).index:
                    f.write(f"{gene}\t{extreme_total[gene]}\t{extreme_high[gene]}\t{extreme_low[gene]}\n")
            
            print(f"\nOutliers saved to {output_path}")
        
        self._pause()
    
    def _detect_combined_outliers(self):
        """Detect outliers using multiple methods and find consensus."""
        print("\nCombined Outlier Detection")
        print("=" * 40)
        print("\nThis will run all outlier detection methods and find consensus outliers.")
        print("This may take a moment...")
        
        outlier_sets = {}
        
        # Statistical outliers on GC content
        print("\n1. Detecting statistical outliers (GC content)...")
        gc_data = self.data.base_composition['GC']
        z_scores = np.abs(stats.zscore(gc_data))
        outlier_sets['GC_outliers'] = set(gc_data[z_scores > 3].index)
        
        # Statistical outliers on GC3s
        print("2. Detecting statistical outliers (GC3s)...")
        gc3s_data = self.data.base_composition['GC3s']
        z_scores = np.abs(stats.zscore(gc3s_data))
        outlier_sets['GC3s_outliers'] = set(gc3s_data[z_scores > 3].index)
        
        # Multivariate outliers if available
        if hasattr(self.data, 'multivariate_results') and self.data.multivariate_results:
            print("3. Detecting multivariate outliers...")
            coords = self.data.multivariate_results['coordinates']
            X = coords.iloc[:, :2].values
            mean = np.mean(X, axis=0)
            cov = np.cov(X.T)
            inv_cov = np.linalg.inv(cov)
            
            distances = []
            for i in range(len(X)):
                diff = X[i] - mean
                dist = np.sqrt(diff.dot(inv_cov).dot(diff))
                distances.append(dist)
            
            distances = pd.Series(distances, index=coords.index)
            threshold = np.sqrt(chi2.ppf(0.995, 2))
            outlier_sets['Multivariate_outliers'] = set(distances[distances > threshold].index)
        
        # Domain-specific outliers
        print("4. Detecting domain-specific outliers...")
        if not hasattr(self.data, 'rscu_values') or self.data.rscu_values.empty:
            self.data.get_rscu_values()  # This will calculate if needed
        
        rscu_values = self.data.rscu_values
        extreme_total = (rscu_values > 2).sum(axis=1) + (rscu_values < 0.5).sum(axis=1)
        threshold = np.percentile(extreme_total, 95)
        outlier_sets['Codon_bias_outliers'] = set(extreme_total[extreme_total > threshold].index)
        
        # Find consensus outliers
        all_outliers = set()
        for outliers in outlier_sets.values():
            all_outliers.update(outliers)
        
        # Count how many methods flagged each gene
        outlier_counts = {}
        for gene in all_outliers:
            count = sum(1 for outliers in outlier_sets.values() if gene in outliers)
            outlier_counts[gene] = count
        
        # Sort by number of methods
        sorted_outliers = sorted(outlier_counts.items(), key=lambda x: x[1], reverse=True)
        
        print(f"\n\nSummary:")
        print(f"Total unique outliers found: {len(all_outliers)}")
        for method, outliers in outlier_sets.items():
            print(f"  {method}: {len(outliers)} genes")
        
        # Show consensus outliers
        consensus = [gene for gene, count in sorted_outliers if count >= 2]
        print(f"\nConsensus outliers (flagged by 2+ methods): {len(consensus)} genes")
        
        if consensus:
            print("\nTop consensus outliers:")
            print(f"{'Gene':30} {'Methods':10} {'Details'}")
            print("-" * 80)
            
            for gene, count in sorted_outliers[:10]:
                if count >= 2:
                    methods = [m.replace('_outliers', '') for m, s in outlier_sets.items() if gene in s]
                    print(f"{gene:30} {count:10} {', '.join(methods)}")
        
        # Export results
        export = input("\nExport combined outlier analysis? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, "combined_outliers") + '.tsv'
            
            with open(output_path, 'w') as f:
                f.write(f"# Combined outlier analysis\n")
                f.write(f"# Total genes analyzed: {len(self.data.gene_names)}\n")
                f.write(f"# Total outliers: {len(all_outliers)}\n")
                f.write(f"# Consensus outliers (2+ methods): {len(consensus)}\n")
                f.write(f"#{'=' * 50}\n\n")
                f.write("Gene\tNum_Methods\tMethods\n")
                
                for gene, count in sorted_outliers:
                    methods = [m.replace('_outliers', '') for m, s in outlier_sets.items() if gene in s]
                    f.write(f"{gene}\t{count}\t{', '.join(methods)}\n")
            
            print(f"\nResults saved to {output_path}")
        
        self._pause()
    
    def _visualize_multivariate_outliers(self, outlier_genes, coords, n_dims):
        """Visualize outliers on multivariate plot."""
        fig = go.Figure()
        
        # Plot non-outliers
        non_outliers = [g for g in coords.index if g not in outlier_genes]
        
        if n_dims >= 3:
            # 3D plot
            fig.add_trace(go.Scatter3d(
                x=coords.loc[non_outliers].iloc[:, 0],
                y=coords.loc[non_outliers].iloc[:, 1],
                z=coords.loc[non_outliers].iloc[:, 2],
                mode='markers',
                marker=dict(size=4, color='lightblue', opacity=0.5),
                text=non_outliers,
                name='Normal genes',
                hovertemplate='<b>%{text}</b><br>Axis1: %{x:.3f}<br>Axis2: %{y:.3f}<br>Axis3: %{z:.3f}<extra></extra>'
            ))
            
            fig.add_trace(go.Scatter3d(
                x=coords.loc[outlier_genes].iloc[:, 0],
                y=coords.loc[outlier_genes].iloc[:, 1],
                z=coords.loc[outlier_genes].iloc[:, 2],
                mode='markers',
                marker=dict(size=8, color='red', symbol='diamond'),
                text=outlier_genes,
                name='Outliers',
                hovertemplate='<b>%{text}</b><br>Axis1: %{x:.3f}<br>Axis2: %{y:.3f}<br>Axis3: %{z:.3f}<extra></extra>'
            ))
        else:
            # 2D plot
            fig.add_trace(go.Scatter(
                x=coords.loc[non_outliers].iloc[:, 0],
                y=coords.loc[non_outliers].iloc[:, 1],
                mode='markers',
                marker=dict(size=6, color='lightblue', opacity=0.5),
                text=non_outliers,
                name='Normal genes',
                hovertemplate='<b>%{text}</b><br>%{xaxis.title.text}: %{x:.3f}<br>%{yaxis.title.text}: %{y:.3f}<extra></extra>'
            ))
            
            fig.add_trace(go.Scatter(
                x=coords.loc[outlier_genes].iloc[:, 0],
                y=coords.loc[outlier_genes].iloc[:, 1],
                mode='markers',
                marker=dict(size=10, color='red', symbol='diamond'),
                text=outlier_genes,
                name='Outliers',
                hovertemplate='<b>%{text}</b><br>%{xaxis.title.text}: %{x:.3f}<br>%{yaxis.title.text}: %{y:.3f}<extra></extra>'
            ))
        
        # Update layout
        explained_var = self.data.multivariate_results['explained_variance']
        if n_dims >= 3:
            fig.update_layout(
                title="Multivariate Outliers (3D)",
                scene=dict(
                    xaxis_title=f"{coords.columns[0]} ({explained_var[0]:.2%})",
                    yaxis_title=f"{coords.columns[1]} ({explained_var[1]:.2%})",
                    zaxis_title=f"{coords.columns[2]} ({explained_var[2]:.2%})"
                ),
                width=1000,
                height=800
            )
        else:
            fig.update_layout(
                title="Multivariate Outliers",
                xaxis_title=f"{coords.columns[0]} ({explained_var[0]:.2%})",
                yaxis_title=f"{coords.columns[1]} ({explained_var[1]:.2%})",
                plot_bgcolor='white',
                width=1000,
                height=800
            )
        
        # Save and display
        output_path = self.file_manager.get_output_path(
            self.data, None, "multivariate_outliers_plot", f".{self.config.visualization_format}")
        
        if self.config.visualization_format == 'html':
            fig.write_html(output_path)
            try:
                webbrowser.open(f"file:///{os.path.abspath(output_path)}")
                print(f"\nOutlier visualization opened in browser.")
            except:
                print(f"\nOutlier visualization saved to {output_path}")
        else:
            fig.write_image(output_path, format=self.config.visualization_format)
            print(f"\nOutlier visualization saved to {output_path}")
    
    def _codon_bias_submenu(self):
        """Submenu for codon bias metrics."""
        if not self._check_data():
            return True
            
        options = [
            ("Calculate ENC values", "1", self._enc_menu),
            ("Calculate CAI & Fop values", "2", self._fop_cai_menu),
            ("Calculate SCUO values", "3", self._scuo_menu),
            ("Calculate all bias metrics", "4", self._calculate_all_bias_metrics),
            ("Return to analysis menu", "R", lambda: False)
        ]
        
        return self._display_menu("Codon Bias Metrics", options)
    
    def _quick_calculate_basic(self):
        """Quick calculation of basic metrics."""
        print("\nCalculating all basic metrics...")
        print("  - Codon usage: Already calculated")
        print("  - Amino acid usage: Already calculated")
        print("  - Base composition: Already calculated")
        print("  - Calculating RSCU values...")
        
        # Calculate RSCU which is the main additional metric needed
        self.data.get_rscu_values()  # This will calculate if needed
        
        # Update analysis status
        self.data.analysis_status['basic_metrics'] = True
        self.data.analysis_status['codon_usage'] = True
        self.data.analysis_status['amino_acid_usage'] = True
        self.data.analysis_status['base_composition'] = True
        self.data.analysis_status['rscu_calculated'] = True
        
        print("\nBasic metrics calculated successfully!")
        self._pause()
        return True
    
    def _quick_multivariate_workflow(self):
        """Quick multivariate analysis workflow."""
        print("\nQuick Multivariate Workflow")
        print("=" * 40)
        
        # Calculate RSCU if needed
        if not hasattr(self.data, 'rscu_values') or self.data.rscu_values.empty:
            print("\nCalculating RSCU values...")
            self.data.get_rscu_values()  # This will calculate if needed
            
        # Perform CA on RSCU
        print("\nPerforming Correspondence Analysis on RSCU...")
        result = self.data.perform_multivariate_analysis('CA', 'RSCU')
        
        if result:
            self._show_multivariate_summary(result, 'CA', 'RSCU')
            
            # Create visualization
            create_viz = input("\nCreate visualization? (y/n): ").strip().lower()
            if create_viz == 'y':
                self._visualize_multivariate(result, 'CA', 'RSCU')
                
            # Export results
            export = input("\nExport results? (y/n): ").strip().lower()
            if export == 'y':
                self._export_multivariate(result, 'CA', 'RSCU')
                
        return True
    
    def _calculate_all_bias_metrics(self):
        """Calculate all codon bias metrics."""
        print("\nCalculating all codon bias metrics...")
        
        # ENC
        print("  - Calculating ENC values...")
        self.data.calculate_enc()
        
        # Optimal codons (needed for CAI/Fop)
        if not hasattr(self.data, 'optimal_codons') or not self.data.optimal_codons:
            print("  - Identifying optimal codons...")
            self.data.calculate_optimal_codons()
            
        # CAI
        print("  - Calculating CAI values...")
        self.data.calculate_cai()
        
        # Fop
        print("  - Calculating Fop values...")
        self.data.calculate_fop()
        
        # SCUO
        print("  - Calculating SCUO values...")
        scuo_values = self.data.calculate_scuo()
        
        print("\nAll codon bias metrics calculated successfully!")
        
        # Show summary
        show_summary = input("\nShow summary? (y/n): ").strip().lower()
        if show_summary == 'y':
            self._codon_bias_overview()
        else:
            self._pause()
            
        return True

    def _visualization_menu(self):
        """Handle visualization options with enhanced plotting capabilities."""
        if not self._check_data():
            return True

        options = [
            ("View multivariate analysis results (CA/PCA plots)", "1", self._multivariate_plot),
            ("Create scree plot for multivariate analysis", "2", self._scree_plot),
            ("Create 3D multivariate plot", "3", self._multivariate_3d_plot),
            ("Create GC content vs GC3 plot", "4", self._gc_content_plot),
            ("Create ENC vs GC3s plot (Wright's plot)", "5", self._enc_plot),
            ("Create RSCU heatmap", "6", self._rscu_heatmap),
            ("Create CAI distribution plot", "7", self._cai_distribution_plot),
            ("Show codon bias overview", "8", self._codon_bias_overview),
            ("Create custom scatter plot", "9", self._custom_scatter_plot),
            ("Return to main menu", "R", lambda: False)
        ]

        return self._display_menu("Visualization Menu", options)

    def _preferences_menu(self):
        """Handle program preferences menu."""
        options = [
            ("Genetic Code", "1", self._genetic_code_menu),
            ("Output Detail Level", "2", self._detail_level_menu),
            ("Visualization Method", "3", self._visualization_method_menu),
            ("Visualization Format", "4", self._visualization_format_menu),
            ("Toggle Beep", "5", self._toggle_beep),
            ("Memory Management Settings", "6", self._memory_settings_menu),
            ("Performance Settings", "7", self._performance_settings_menu),
            ("Return to main menu", "R", lambda: False)
        ]

        return self._display_menu("Preferences Menu", options)

    def _toggle_beep(self):
        """Toggle beep setting."""
        self.config.beep_on = not self.config.beep_on
        print(f"\nBeep {'enabled' if self.config.beep_on else 'disabled'}.")
        self._pause()
        return True
    
    def _memory_settings_menu(self):
        """Handle memory management settings menu."""
        while True:
            self._clear_screen()
            print("\nMemory Management Settings")
            print("=" * 40)
            print(f"Current mode: {self.config.memory_mode}")
            print(f"File size threshold: {self.config.memory_threshold_mb} MB")
            print(f"Sequence indexing: {'Enabled' if self.config.enable_sequence_indexing else 'Disabled'}")
            if hasattr(self.data, 'memory_mode'):
                print(f"\nCurrent file memory mode: {self.data.memory_mode}")
                print(f"Current file size: {self.data.file_size_mb:.1f} MB")
            print("\nOptions:")
            print("1. Change memory mode (auto/low_memory/full_memory)")
            print("2. Set file size threshold for auto mode")
            print("3. Toggle sequence indexing")
            print("R. Return to preferences menu")
            
            choice = input("\nSelect option: ").strip().upper()
            
            if choice == '1':
                print("\nSelect memory mode:")
                print("1. auto - Automatically decide based on file size")
                print("2. low_memory - Always use low memory mode")
                print("3. full_memory - Always keep sequences in memory")
                mode_choice = input("\nSelect mode (1-3): ").strip()
                
                mode_map = {'1': 'auto', '2': 'low_memory', '3': 'full_memory'}
                if mode_choice in mode_map:
                    self.config.memory_mode = mode_map[mode_choice]
                    print(f"\nMemory mode set to: {self.config.memory_mode}")
                    if hasattr(self.data, 'sequences'):
                        print("Note: This will take effect when loading the next file.")
                else:
                    print("\nInvalid choice.")
                self._pause()
                
            elif choice == '2':
                try:
                    threshold = float(input("\nEnter file size threshold in MB (current: {}): ".format(
                        self.config.memory_threshold_mb)))
                    if threshold > 0:
                        self.config.memory_threshold_mb = threshold
                        print(f"\nThreshold set to: {threshold} MB")
                    else:
                        print("\nThreshold must be positive.")
                except ValueError:
                    print("\nInvalid input. Please enter a number.")
                self._pause()
                
            elif choice == '3':
                self.config.enable_sequence_indexing = not self.config.enable_sequence_indexing
                print(f"\nSequence indexing {'enabled' if self.config.enable_sequence_indexing else 'disabled'}.")
                self._pause()
                
            elif choice == 'R':
                return True
            else:
                print("\nInvalid choice. Please try again.")
                self._pause()
    
    def _performance_settings_menu(self):
        """Handle performance settings menu."""
        while True:
            self._clear_screen()
            print("\nPerformance Settings")
            print("=" * 40)
            print(f"Parallel processing: {'Enabled' if self.config.enable_parallel_processing else 'Disabled'}")
            print(f"Number of processes: {self.config.num_processes}")
            print(f"Batch size: {self.config.batch_size} sequences")
            print(f"Progress saving: {'Enabled' if self.config.enable_progress_saving else 'Disabled'}")
            print(f"Checkpoint interval: {self.config.checkpoint_interval} sequences")
            print(f"\nCPU cores available: {mp.cpu_count()}")
            print("\nOptions:")
            print("1. Toggle parallel processing")
            print("2. Set number of processes")
            print("3. Set batch size")
            print("4. Toggle progress saving")
            print("5. Set checkpoint interval")
            print("R. Return to preferences menu")
            
            choice = input("\nSelect option: ").strip().upper()
            
            if choice == '1':
                self.config.enable_parallel_processing = not self.config.enable_parallel_processing
                print(f"\nParallel processing {'enabled' if self.config.enable_parallel_processing else 'disabled'}.")
                self._pause()
                
            elif choice == '2':
                try:
                    num_proc = int(input(f"\nEnter number of processes (1-{mp.cpu_count()}, current: {self.config.num_processes}): "))
                    if 1 <= num_proc <= mp.cpu_count():
                        self.config.num_processes = num_proc
                        print(f"\nNumber of processes set to: {num_proc}")
                    else:
                        print(f"\nNumber must be between 1 and {mp.cpu_count()}.")
                except ValueError:
                    print("\nInvalid input. Please enter a number.")
                self._pause()
                
            elif choice == '3':
                try:
                    batch_size = int(input(f"\nEnter batch size (current: {self.config.batch_size}): "))
                    if batch_size > 0:
                        self.config.batch_size = batch_size
                        print(f"\nBatch size set to: {batch_size}")
                    else:
                        print("\nBatch size must be positive.")
                except ValueError:
                    print("\nInvalid input. Please enter a number.")
                self._pause()
                
            elif choice == '4':
                self.config.enable_progress_saving = not self.config.enable_progress_saving
                print(f"\nProgress saving {'enabled' if self.config.enable_progress_saving else 'disabled'}.")
                self._pause()
                
            elif choice == '5':
                try:
                    interval = int(input(f"\nEnter checkpoint interval (current: {self.config.checkpoint_interval}): "))
                    if interval > 0:
                        self.config.checkpoint_interval = interval
                        print(f"\nCheckpoint interval set to: {interval}")
                    else:
                        print("\nInterval must be positive.")
                except ValueError:
                    print("\nInvalid input. Please enter a number.")
                self._pause()
                
            elif choice == 'R':
                return True
            else:
                print("\nInvalid choice. Please try again.")
                self._pause()

    def _genetic_code_menu(self):
        """Handle genetic code selection menu."""
        # Create dynamic options from all available genetic codes
        options = []

        # Sort genetic code IDs numerically
        code_ids = sorted(GENETIC_CODES.keys())

        # Build the options list
        for code_id in code_ids:
            code_name = GENETIC_CODES[code_id]
            selected_str = "[SELECTED]" if self.config.genetic_code == code_id else ""

            # Create a lambda function that captures the current code_id
            def set_code_func(
                code=code_id): return self._set_genetic_code(code)

            options.append(
                (f"{code_id}: {code_name} {selected_str}",
                 str(code_id),
                    set_code_func))

        # Add return option
        options.append(("Return to preferences menu", "R", lambda: False))

        self._clear_screen()
        self.display_banner()
        print("\nWARNING: If you change the genetic code, you must re-load the datafile.")
        print("\nChoose the genetic code that matches your sequence data:")

        return self._display_menu("Genetic Code Selection", options)

    def _set_genetic_code(self, code):
        """Set the genetic code and inform the user."""
        # Verify the code exists
        if code in GENETIC_CODES:
            self.config.genetic_code = code
            code_name = GENETIC_CODES[code]
            print(f"\nGenetic code [{code}] {code_name} selected.")

            # If sequences are already loaded, warn about reloading
            if self.data.gene_names:
                print(
                    "\nWARNING: You have changed the genetic code while sequences are loaded.")
                print(
                    "You should reload your sequences to ensure proper codon translation.")
                reload = input(
                    "\nDo you want to reload your sequences now? (y/n): ").strip().lower()
                if reload == 'y' and self.data.file_path:
                    self.data.load_fasta(self.data.file_path)
        else:
            print(f"\nInvalid genetic code: {code}. Keeping current setting.")

        self._pause()
        return True

    def _detail_level_menu(self):
        """Handle detail level selection menu."""
        options = [
            (f"Basic (results only) [{self.config.detailed_output == 1 and 'SELECTED' or ''}]",
             "1",
             lambda: self._set_detail_level(1)),
            (f"Moderate (some calculations) [{self.config.detailed_output == 2 and 'SELECTED' or ''}]",
             "2",
             lambda: self._set_detail_level(2)),
            (f"Detailed (all calculations) [{self.config.detailed_output == 3 and 'SELECTED' or ''}]",
             "3",
             lambda: self._set_detail_level(3)),
            ("Return to preferences menu",
             "R",
             lambda: False)]

        return self._display_menu("Output Detail Level", options)

    def _set_detail_level(self, level):
        """Set the detail level and inform the user."""
        self.config.detailed_output = level
        level_names = {1: "Basic", 2: "Moderate", 3: "Detailed"}
        print(f"\n{level_names[level]} output level selected.")
        self._pause()
        return True

    def _visualization_method_menu(self):
        """Handle visualization method selection menu."""
        options = [
            (f"Plotly (interactive) [{self.config.visualization_method == 'plotly' and 'SELECTED' or ''}]",
             "1",
             lambda: self._set_visualization_method('plotly')),
            (f"Text - not currently implemented [{self.config.visualization_method == 'text' and 'SELECTED' or ''}]",
             "2",
             lambda: self._set_visualization_method('plotly')), #this is a HACK placeholder until text-based vizusualization is implemented.
            ("Return to preferences menu",
             "R",
             lambda: False)]

        return self._display_menu("Visualization Method", options)

    def _set_visualization_method(self, method):
        """Set the visualization method and inform the user."""
        self.config.visualization_method = method
        print(f"\n{method.capitalize()} visualization method selected.")
        self._pause()
        return True
    
    def _visualization_format_menu(self):
        """Handle visualization format selection menu."""
        current = self.config.visualization_format
        
        # Check if kaleido is available
        try:
            import kaleido
            has_kaleido = True
        except ImportError:
            has_kaleido = False
        
        options = [
            (f"HTML (interactive, opens in browser) [{current == 'html' and 'SELECTED' or ''}]",
             "1",
             lambda: self._set_visualization_format('html')),
             (f"SVG (vector graphics) [{current == 'svg' and 'SELECTED' or ''}]" + (" [Requires kaleido]" if not has_kaleido else ""),
             "2",
             lambda: self._set_visualization_format('svg')),
            (f"PNG (raster image) [{current == 'png' and 'SELECTED' or ''}]" + (" [Requires kaleido]" if not has_kaleido else ""),
             "3",
             lambda: self._set_visualization_format('png')),
            (f"PDF (document format) [{current == 'pdf' and 'SELECTED' or ''}]" + (" [Requires kaleido]" if not has_kaleido else ""),
             "4",
             lambda: self._set_visualization_format('pdf')),
            ("Return to preferences menu",
             "R",
             lambda: False)]
        
        if not has_kaleido and current != 'html':
            print("\nWARNING: kaleido package not installed. Static image formats require kaleido.")
            print("Install with: pip install kaleido")
            print("")
        
        return self._display_menu("Visualization Format", options)
    
    def _set_visualization_format(self, format):
        """Set the visualization format and inform the user."""
        # Check if kaleido is available for non-HTML formats
        if format != 'html':
            try:
                import kaleido
            except ImportError:
                print(f"\nError: {format.upper()} format requires the kaleido package.")
                print("Install with: pip install kaleido")
                print(f"Keeping current format: {self.config.visualization_format}")
                self._pause()
                return True
        
        self.config.visualization_format = format
        print(f"\n{format.upper()} visualization format selected.")
        if format != 'html':
            print(f"Visualizations will be saved as .{format} files.")
            print("Note: Static image files will not open automatically in browser.")
        self._pause()
        return True

    def _multivariate_plot(self):
        """Create and display a multivariate analysis plot."""
        # Check if any multivariate analyses have been performed
        if not hasattr(self.data, 'multivariate_history') or not self.data.multivariate_history:
            print("\nNo multivariate analysis results found.")
            print("Please perform a multivariate analysis first:")
            print("  Main Menu > Analysis > Perform multivariate analysis")
            self._pause()
            return True
        
        # If only one analysis exists, use it
        if len(self.data.multivariate_history) == 1:
            analysis = self.data.multivariate_history[0]
            result = analysis['result']
            print(f"\nVisualizing {analysis['analysis_type']} of {analysis['data_type']} data...")
        else:
            # Multiple analyses exist, let user choose
            print("\nAvailable multivariate analyses:")
            for i, analysis in enumerate(self.data.multivariate_history, 1):
                print(f"{i}. {analysis['analysis_type']} of {analysis['data_type']} data")
                print(f"   Performed at: {analysis['timestamp'][:19]}")
            
            choice = input("\nWhich analysis would you like to visualize? (number): ").strip()
            try:
                idx = int(choice) - 1
                if 0 <= idx < len(self.data.multivariate_history):
                    analysis = self.data.multivariate_history[idx]
                    result = analysis['result']
                    print(f"\nVisualizing {analysis['analysis_type']} of {analysis['data_type']} data...")
                else:
                    print("Invalid choice.")
                    self._pause()
                    return True
            except ValueError:
                print("Invalid input.")
                self._pause()
                return True

        # Ask user which axes to plot
        coords = result['coordinates']
        max_axes = min(5, coords.shape[1])  # Show up to 5 axes
        
        print("\nAvailable axes:")
        for i in range(max_axes):
            axis_name = coords.columns[i]
            variance = result['explained_variance'][i]
            print(f"{i+1}. {axis_name} ({variance:.1f}% variance)")
        
        # Get axis selection
        print("\nWhich axes would you like to plot?")
        x_axis = input("X-axis (default=1): ").strip()
        x_axis = int(x_axis) - 1 if x_axis else 0
        
        y_axis = input("Y-axis (default=2): ").strip()
        y_axis = int(y_axis) - 1 if y_axis else 1
        
        # Validate axes
        if not (0 <= x_axis < max_axes and 0 <= y_axis < max_axes):
            print("Invalid axis selection. Using defaults (1 vs 2).")
            x_axis, y_axis = 0, 1
        elif x_axis == y_axis:
            print("Cannot plot the same axis against itself. Using defaults (1 vs 2).")
            x_axis, y_axis = 0, 1
        
        # Store selected axes in result for visualization
        result['selected_axes'] = (x_axis, y_axis)
        
        # Create filename that includes axis information
        axis_suffix = f"_axis{x_axis+1}_vs_axis{y_axis+1}"
        output_path = self.file_manager.get_output_path(
            self.data, None, f"multivariate_plot{axis_suffix}", f".{self.config.visualization_format}")
        
        print(f"\nCreating plot: {coords.columns[x_axis]} vs {coords.columns[y_axis]}...")
        
        # Create and display the plot
        self.data.visualize('multivariate', result, output_path)
        self._pause()
        return True

    def _gc_content_plot(self):
        """Create and display GC content plot."""
        print("\nCreating GC content vs GC3 plot...")
        output_path = self.file_manager.get_output_path(
            self.data, None, "gc_content_plot", f".{self.config.visualization_format}")
        self.data.visualize('gc_content', output_path=output_path)
        self._pause()
        return True

    def _enc_plot(self):
        """Create and display ENC vs GC3s plot."""
        print("\nCreating ENC vs GC3s plot (Wright's plot)...")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

        # Make sure ENC values are calculated
        if not hasattr(self.data, 'enc_values') or not self.data.enc_values:
            print("Calculating ENC values...")
            self.data.calculate_enc()

        output_path = self.file_manager.get_output_path(
            self.data, None, "enc_plot", f".{self.config.visualization_format}")
        self.data.visualize('enc', output_path=output_path)
        self._pause()
        return True

    def _rscu_heatmap(self):
        """Create and display RSCU heatmap."""
        print("\nCreating RSCU heatmap...")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        
        # Get RSCU values to check dataset size
        rscu_values = self.data.get_rscu_values()
        num_genes = len(rscu_values)
        
        if num_genes > 1000:
            print(f"\nLarge dataset detected ({num_genes:,} genes)")
            print("Creating summary statistics visualization instead of gene-level heatmap...")
            
            # Calculate and display summary statistics
            aa_map = self.config.get_aa_map()
            non_stop_codons = [c for c in CODONS if aa_map[c] != 'STOP']
            rscu_filtered = rscu_values[non_stop_codons]
            
            mean_rscu = rscu_filtered.mean()
            
            # Find most over- and under-represented codons
            print("\nRSCU Summary Statistics:")
            print(f"Total genes analyzed: {num_genes:,}")
            
            # Sort codons by mean RSCU
            sorted_codons = mean_rscu.sort_values(ascending=False)
            
            print("\nMost overrepresented codons (mean RSCU > 1.5):")
            overrep = sorted_codons[sorted_codons > 1.5]
            if len(overrep) > 0:
                for codon, value in overrep.items():
                    aa = aa_map[codon]
                    std = rscu_filtered[codon].std()
                    print(f"  {codon} ({aa}): {value:.3f} ± {std:.3f}")
            else:
                print("  None")
            
            print("\nMost underrepresented codons (mean RSCU < 0.5):")
            underrep = sorted_codons[sorted_codons < 0.5]
            if len(underrep) > 0:
                for codon, value in underrep.items():
                    aa = aa_map[codon]
                    std = rscu_filtered[codon].std()
                    print(f"  {codon} ({aa}): {value:.3f} ± {std:.3f}")
            else:
                print("  None")
        else:
            print(f"\nDataset size: {num_genes} genes")
            print("Creating traditional gene-level heatmap...")
        
        output_path = self.file_manager.get_output_path(
            self.data, None, "rscu_heatmap", f".{self.config.visualization_format}")
        self.data.visualize('rscu_heatmap', output_path=output_path)
        self._pause()
        return True

    def _cai_distribution_plot(self):
        """Create and display a CAI distribution plot."""
        print("\nCreating CAI distribution plot...")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

        # Calculate CAI if not already done
        if not hasattr(self.data, 'cai_values') or not self.data.cai_values:
            print("Calculating Codon Adaptation Index (CAI)...")
            cai_values = self.data.calculate_cai()
        else:
            cai_values = self.data.cai_values

        output_path = self.file_manager.get_output_path(
            self.data, None, "cai_distribution", f".{self.config.visualization_format}")
        self.data.visualize('cai_distribution', cai_values, output_path)
        self._pause()
        return True

    def _scree_plot(self):
        """Create and display a scree plot for multivariate analysis."""
        # Check if any multivariate analyses have been performed
        if not hasattr(self.data, 'multivariate_history') or not self.data.multivariate_history:
            print("\nNo multivariate analysis results found.")
            print("Please perform a multivariate analysis first:")
            print("  Main Menu > Analysis > Perform multivariate analysis")
            self._pause()
            return True
        
        # If only one analysis exists, use it
        if len(self.data.multivariate_history) == 1:
            analysis = self.data.multivariate_history[0]
            self.data.multivariate_results = analysis['result']
            print(f"\nCreating scree plot for {analysis['analysis_type']} of {analysis['data_type']} data...")
        else:
            # Multiple analyses exist, let user choose
            print("\nAvailable multivariate analyses:")
            for i, analysis in enumerate(self.data.multivariate_history, 1):
                print(f"{i}. {analysis['analysis_type']} of {analysis['data_type']} data")
                print(f"   Performed at: {analysis['timestamp'][:19]}")
            
            choice = input("\nWhich analysis would you like to create a scree plot for? (number): ").strip()
            try:
                idx = int(choice) - 1
                if 0 <= idx < len(self.data.multivariate_history):
                    analysis = self.data.multivariate_history[idx]
                    self.data.multivariate_results = analysis['result']
                    print(f"\nCreating scree plot for {analysis['analysis_type']} of {analysis['data_type']} data...")
                else:
                    print("Invalid choice.")
                    self._pause()
                    return True
            except ValueError:
                print("Invalid input.")
                self._pause()
                return True

        # Create and display the scree plot
        output_path = self.file_manager.get_output_path(
            self.data, None, "scree_plot", f".{self.config.visualization_format}")
        self.data.visualize('scree_plot', output_path=output_path)
        self._pause()
        return True

    def _multivariate_3d_plot(self):
        """Create and display a 3D multivariate analysis plot."""
        # Check if any multivariate analyses have been performed
        if not hasattr(self.data, 'multivariate_history') or not self.data.multivariate_history:
            print("\nNo multivariate analysis results found.")
            print("Please perform a multivariate analysis first:")
            print("  Main Menu > Analysis > Perform multivariate analysis")
            self._pause()
            return True
        
        # If only one analysis exists, use it
        if len(self.data.multivariate_history) == 1:
            analysis = self.data.multivariate_history[0]
            result = analysis['result']
            # Check if we have at least 3 components
            if result['coordinates'].shape[1] < 3:
                print(f"\nError: {analysis['analysis_type']} of {analysis['data_type']} data has only {result['coordinates'].shape[1]} components.")
                print("3D visualization requires at least 3 components.")
                self._pause()
                return True
            self.data.multivariate_results = result
            print(f"\nCreating 3D plot for {analysis['analysis_type']} of {analysis['data_type']} data...")
        else:
            # Multiple analyses exist, let user choose
            print("\nAvailable multivariate analyses:")
            valid_analyses = []
            for i, analysis in enumerate(self.data.multivariate_history):
                if analysis['result']['coordinates'].shape[1] >= 3:
                    valid_analyses.append((i, analysis))
                    print(f"{len(valid_analyses)}. {analysis['analysis_type']} of {analysis['data_type']} data")
                    print(f"   Performed at: {analysis['timestamp'][:19]}")
                    print(f"   Components: {analysis['result']['coordinates'].shape[1]}")
            
            if not valid_analyses:
                print("\nNo analyses with 3 or more components found.")
                print("3D visualization requires at least 3 components.")
                self._pause()
                return True
            
            choice = input("\nWhich analysis would you like to visualize in 3D? (number): ").strip()
            try:
                idx = int(choice) - 1
                if 0 <= idx < len(valid_analyses):
                    orig_idx, analysis = valid_analyses[idx]
                    self.data.multivariate_results = analysis['result']
                    print(f"\nCreating 3D plot for {analysis['analysis_type']} of {analysis['data_type']} data...")
                else:
                    print("Invalid choice.")
                    self._pause()
                    return True
            except ValueError:
                print("Invalid input.")
                self._pause()
                return True

        # Create and display the 3D plot
        output_path = self.file_manager.get_output_path(
            self.data, None, "multivariate_3d", f".{self.config.visualization_format}")
        self.data.visualize('multivariate_3d', output_path=output_path)
        self._pause()
        return True

    def _custom_scatter_plot(self):
        """Create a custom scatter plot by allowing the user to select metrics to plot."""
        print("\nCustom Scatter Plot Creation")
        print("============================")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

        # First, ensure we have calculated all possible metrics
        self._ensure_all_metrics_calculated()

        # Get available metrics to plot
        available_metrics = self._get_available_metrics()

        if not available_metrics:
            print("\nNo metrics available for plotting. Please perform analysis first.")
            self._pause()
            return True

        # Display available metrics
        print("\nAvailable metrics for plotting:")
        for i, (metric_name, metric_info) in enumerate(
                available_metrics.items(), 1):
            print(f"{i}. {metric_name} - {metric_info['description']}")

        # Get user selection for X-axis
        try:
            x_selection = int(
                input("\nSelect metric for X-axis (number): ").strip())
            if x_selection < 1 or x_selection > len(available_metrics):
                print("\nInvalid selection. Please try again.")
                self._pause()
                return True
        except ValueError:
            print("\nInvalid input. Please enter a number.")
            self._pause()
            return True

        # Get user selection for Y-axis
        try:
            y_selection = int(
                input("\nSelect metric for Y-axis (number): ").strip())
            if y_selection < 1 or y_selection > len(available_metrics):
                print("\nInvalid selection. Please try again.")
                self._pause()
                return True
        except ValueError:
            print("\nInvalid input. Please enter a number.")
            self._pause()
            return True

        # Convert selections to metric names
        x_metric = list(available_metrics.keys())[x_selection - 1]
        y_metric = list(available_metrics.keys())[y_selection - 1]

        # Create the scatter plot
        return self._create_custom_scatter_plot(
            x_metric, y_metric, available_metrics)

    def _ensure_all_metrics_calculated(self):
        """Calculate only missing metrics for the loaded sequences."""
        if not self._check_data():
            return True

        # Check what metrics are already calculated
        needs_calculation = []
        already_calculated = []
        
        # Check multivariate analysis
        if not hasattr(self.data, 'multivariate_results') or self.data.multivariate_results is None:
            needs_calculation.append("multivariate analysis (CA on RSCU)")
        else:
            already_calculated.append("multivariate analysis")
            
        # Check ENC values
        if not hasattr(self.data, 'enc_values') or self.data.enc_values is None:
            needs_calculation.append("ENC values")
        else:
            already_calculated.append("ENC values")
            
        # Check optimal codons (needed for CAI/Fop)
        if not hasattr(self.data, 'optimal_codons') or self.data.optimal_codons is None:
            needs_calculation.append("optimal codons")
        else:
            already_calculated.append("optimal codons")
            
        # Check CAI values
        if not hasattr(self.data, 'cai_values') or self.data.cai_values is None:
            needs_calculation.append("CAI values")
        else:
            already_calculated.append("CAI values")
            
        # Check Fop values
        if not hasattr(self.data, 'fop_values') or self.data.fop_values is None:
            needs_calculation.append("Fop values")
        else:
            already_calculated.append("Fop values")
        
        # If everything is already calculated, ask if user wants to recalculate
        if not needs_calculation:
            print("\nAll metrics already calculated!")
            print(f"Cached metrics: {', '.join(already_calculated)}")
            recalc = input("\nDo you want to recalculate all metrics? (y/n): ").strip().lower()
            if recalc != 'y':
                # Ask about export without recalculating
                export = input("\nDo you want to export all metrics to a file? (y/n): ").strip().lower()
                if export == 'y':
                    self._export_comprehensive_metrics()
                self._pause()
                return True
            # If user wants to recalculate, calculate everything
            needs_calculation = ["multivariate analysis (CA on RSCU)", "ENC values", 
                               "optimal codons", "CAI values", "Fop values"]
        else:
            # Show what's already calculated and what needs calculation
            if already_calculated:
                print(f"\nUsing cached: {', '.join(already_calculated)}")
            print(f"Calculating: {', '.join(needs_calculation)}")

        # Calculate only what's needed
        if "multivariate analysis (CA on RSCU)" in needs_calculation:
            print("\nPerforming multivariate analysis...")
            self.data.perform_multivariate_analysis()

        if "ENC values" in needs_calculation:
            print("Calculating ENC values...")
            self.data.calculate_enc()

        if "optimal codons" in needs_calculation:
            print("Calculating optimal codons...")
            self.data.calculate_optimal_codons()

        if "Fop values" in needs_calculation:
            print("Calculating Fop values...")
            self.data.calculate_fop()

        if "CAI values" in needs_calculation:
            print("Calculating CAI values...")
            self.data.calculate_cai()

        # Always calculate SCUO as it's fast and not stored
        print("Calculating SCUO values...")
        scuo_values = self.data.calculate_scuo()

        print("\nAll metrics calculated successfully.")

        # Ask if user wants to export to file
        export = input(
            "\nDo you want to export all metrics to a file? (y/n): ").strip().lower()
        if export == 'y':
            self._export_comprehensive_metrics()

        self._pause()
        return True

    def _get_available_metrics(self):
        """Get all available metrics that can be plotted."""
        metrics = {}

        # Base composition metrics
        for col in self.data.base_composition.columns:
            if col not in [
                'Length',
                    'AA_count']:  # Skip non-percentage metrics
                metrics[f"base_comp_{col}"] = {
                    "description": f"Base Composition: {col}",
                    "values": self.data.base_composition[col],
                    "type": "Base Composition"
                }

        # ENC values
        if hasattr(self.data, 'enc_values') and self.data.enc_values:
            metrics["ENC"] = {
                "description": "Effective Number of Codons",
                "values": pd.Series(self.data.enc_values),
                "type": "Codon Bias"
            }

        # CAI values
        if hasattr(self.data, 'cai_values') and self.data.cai_values:
            metrics["CAI"] = {
                "description": "Codon Adaptation Index",
                "values": pd.Series(self.data.cai_values),
                "type": "Codon Bias"
            }

        # Fop values
        if hasattr(self.data, 'fop_values') and self.data.fop_values:
            metrics["Fop"] = {
                "description": "Frequency of Optimal Codons",
                "values": pd.Series(self.data.fop_values),
                "type": "Codon Bias"
            }

        # SCUO values
        scuo_values = self.data.calculate_scuo()
        metrics["SCUO"] = {
            "description": "Synonymous Codon Usage Order",
            "values": pd.Series(scuo_values),
            "type": "Codon Bias"
        }

        # Multivariate analysis results
        if hasattr(
                self.data,
                'multivariate_results') and self.data.multivariate_results:
            coords = self.data.multivariate_results['coordinates']
            for col in coords.columns:
                metrics[f"multivar_{col}"] = {
                    "description": f"Multivariate: {col}",
                    "values": coords[col],
                    "type": "Multivariate"
                }

        return metrics

    def _create_custom_scatter_plot(self,x_metric,y_metric,available_metrics):
        """Create a custom scatter plot based on selected metrics."""
        # Get data for both axes
        x_data = available_metrics[x_metric]["values"]
        y_data = available_metrics[y_metric]["values"]
        x_desc = available_metrics[x_metric]["description"]
        y_desc = available_metrics[y_metric]["description"]

        # Ensure we have matching indices
        common_indices = x_data.index.intersection(y_data.index)
        if len(common_indices) == 0:
            print("\nNo common data points between selected metrics.")
            self._pause()
            return True

        x_data = x_data.loc[common_indices]
        y_data = y_data.loc[common_indices]

        # Create a DataFrame for the scatter plot
        plot_data = pd.DataFrame({
            'x': x_data,
            'y': y_data,
            'gene': common_indices
        })

        # Generate plot title and axis labels
        title = f"{y_desc} vs {x_desc}"

        # Create the output path
        output_path = self.file_manager.get_output_path(
            self.data, None, f"custom_plot_{x_metric}_vs_{y_metric}"
        )

        # Create visualization
        vis = CustomScatterPlotVisualization()
        result_path = vis.create(plot_data, output_path)

        # Open in browser if path is available
        if result_path:
            vis.open_in_browser(result_path)

        print(
            f"\nScatter plot of {y_desc} vs {x_desc} created and opened in browser.")
        print(f"Number of data points: {len(plot_data)}")

        # Ask if user wants additional customization
        add_customization = input(
            "\nWould you like to perform outlier analysis? (y/n): ").strip().lower()
        if add_customization == 'y':
            self._analyze_outliers(
                plot_data, x_metric, y_metric, x_desc, y_desc
            )

        self._pause()
        return True

    def _analyze_outliers(self, plot_data, x_metric, y_metric, x_desc, y_desc):
        """Perform outlier analysis on the scatter plot data."""
        print("\nPerforming outlier analysis...")

        # Calculate z-scores for both dimensions
        plot_data['x_zscore'] = stats.zscore(plot_data['x'])
        plot_data['y_zscore'] = stats.zscore(plot_data['y'])

        # Calculate Euclidean distance from the center
        plot_data['distance'] = np.sqrt(
            plot_data['x_zscore']**2 +
            plot_data['y_zscore']**2)

        # Sort by distance to find outliers
        sorted_data = plot_data.sort_values('distance', ascending=False)

        # Define threshold for outliers (e.g., top 5% or z-score > 2)
        threshold = 2.0  # Adjust as needed
        outliers = sorted_data[sorted_data['distance'] > threshold]

        # Display outliers
        if not outliers.empty:
            print(
                f"\nFound {len(outliers)} potential outliers (distance > {threshold}):")
            print("\nTop outliers:")
            print(f"{'Gene':20} {x_desc:20} {y_desc:20} Distance")
            print("-" * 80)

            for i, (idx, row) in enumerate(outliers.iterrows(), 1):
                if i <= 10:  # Show top 10 outliers
                    print(
                        f"{row['gene']:20} {row['x']:20.4f} {row['y']:20.4f} {row['distance']:.4f}")

            if len(outliers) > 10:
                print(f"\n... and {len(outliers) - 10} more outliers.")

            # Ask if user wants to save outlier list
            save_outliers = input(
                "\nDo you want to save the complete list of outliers? (y/n): ").strip().lower()
            if save_outliers == 'y':
                output_path = self.file_manager.get_output_path(
                    self.data, None, f"outliers_{x_metric}_vs_{y_metric}") + '.tsv'

                with open(output_path, 'w') as f:
                    f.write(f"Gene\t{x_desc}\t{y_desc}\tDistance\n")
                    for idx, row in outliers.iterrows():
                        f.write(
                            f"{row['gene']}\t{row['x']:.4f}\t{row['y']:.4f}\t{row['distance']:.4f}\n")

                print(f"\nOutliers saved to {output_path}")
        else:
            print("\nNo significant outliers found.")

    def _codon_bias_overview(self):
        """Display an overview of codon usage bias."""
        print("\nGenerating codon usage bias overview...")

        # Make sure all metrics are calculated
        if not hasattr(self.data, 'enc_values') or not self.data.enc_values:
            print("Calculating ENC values...")
            self.data.calculate_enc()

        if not hasattr(self.data, 'cai_values') or not self.data.cai_values:
            print("Calculating CAI values...")
            self.data.calculate_cai()

        if not hasattr(self.data, 'fop_values') or not self.data.fop_values:
            print("Calculating Fop values...")
            self.data.calculate_fop()

        scuo_values = self.data.calculate_scuo()

        # Create overview table
        print("\nCodon Usage Bias Overview (first 10 genes):")
        print("Gene\tENC\tCAI\tFop\tSCUO\tGC\tGC3s")

        # Sort genes by CAI (highest first)
        sorted_genes = sorted(
            self.data.gene_names,
            key=lambda g: self.data.cai_values[g],
            reverse=True)

        for gene in sorted_genes[:min(10, len(sorted_genes))]:
            enc = self.data.enc_values[gene]
            cai = self.data.cai_values[gene]
            fop = self.data.fop_values[gene]
            scuo = scuo_values[gene]
            gc = self.data.base_composition.at[gene, 'GC']
            gc3s = self.data.base_composition.at[gene, 'GC3s']

            print(
                f"{gene}\t{enc:.2f}\t{cai:.4f}\t{fop:.4f}\t{scuo:.4f}\t{gc:.2f}\t{gc3s:.2f}")

        # Calculate and show correlations
        print("\nCorrelations between metrics:")

        # Create a DataFrame with all metrics
        metrics_df = pd.DataFrame({
            'ENC': pd.Series(self.data.enc_values),
            'CAI': pd.Series(self.data.cai_values),
            'Fop': pd.Series(self.data.fop_values),
            'SCUO': pd.Series(scuo_values),
            'GC': self.data.base_composition['GC'],
            'GC3s': self.data.base_composition['GC3s']
        })

        # Calculate correlation matrix
        corr_matrix = metrics_df.corr()

        # Display correlation matrix
        print("\nCorrelation Matrix:")
        print(corr_matrix.round(3))

        # Add genetic code information
        print(
            f"\nGenetic code used: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

        # Ask if user wants to export
        export = input(
            "\nDo you want to export the overview to a file? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, "codon_bias_overview") + '.tsv'

            with open(output_path, 'w') as f:
                # Add metadata
                f.write(f"# GCUA Codon Bias Overview\n")
                f.write(f"# Version: {VERSION}\n")
                f.write(
                    f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
                f.write(f"# Total Genes: {len(self.data.gene_names)}\n")
                f.write(
                    f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"#{'=' * 50}\n\n")

                f.write("Gene\tENC\tCAI\tFop\tSCUO\tGC\tGC3s\n")

                for gene in self.data.gene_names:
                    enc = self.data.enc_values[gene]
                    cai = self.data.cai_values[gene]
                    fop = self.data.fop_values[gene]
                    scuo = scuo_values[gene]
                    gc = self.data.base_composition.at[gene, 'GC']
                    gc3s = self.data.base_composition.at[gene, 'GC3s']

                    f.write(
                        f"{gene}\t{enc:.4f}\t{cai:.4f}\t{fop:.4f}\t{scuo:.4f}\t{gc:.4f}\t{gc3s:.4f}\n")

                # Add correlation matrix
                f.write("\nCorrelation Matrix:\n")
                f.write("\t" + "\t".join(corr_matrix.columns) + "\n")

                for idx, row in corr_matrix.iterrows():
                    f.write(
                        f"{idx}\t" +
                        "\t".join(
                            f"{val:.4f}" for val in row) +
                        "\n")

            print(f"\nCodon bias overview saved to {output_path}")

        self._pause()
        return True

    def _analysis_menu(self):
        """Handle analysis options with organized groups."""
        if not self._check_data():
            return True

        print("\nCustom Analysis Menu")
        print("=" * 40)
        print("\n--- Individual Analyses ---")
        print("1. Codon usage details")
        print("2. Amino acid usage details")
        print("3. Base composition details")
        print("4. RSCU calculation")
        
        print("\n--- Advanced Analysis ---")
        print("5. Multivariate analysis (CA/PCA)")
        print("6. Codon bias metrics (ENC, CAI, Fop, SCUO)")
        print("7. Identify optimal codons")
        print("8. Compare axis cohorts")
        
        print("\n--- Comprehensive Options ---")
        print("9. Calculate ALL metrics (basic + advanced)")
        print("M. Quick multivariate workflow")
        
        print("\n--- Utilities ---")
        print("C. Clear cached metrics")
        print("S. Show analysis status")
        
        print("\nR. Return to main menu")
        print("H. Show help for all menu items")
        
        # Show help hint on first analysis menu display
        if self.show_help_hint:
            print("\n💡 Tip: Type a number followed by '?' (e.g., '1?') for help on that option")
        
        choice = input("\nSelect an option: ").strip().upper()
        
        # Check for help requests
        if choice == 'H':
            self._show_analysis_menu_help()
            self._pause()
            return True
        elif choice.endswith('?'):
            option_key = choice[:-1]
            self._show_quick_help('analysis_menu', option_key)
            self._pause()
            return True
        
        # Handle choices
        if choice == "1":
            return self._codon_usage_menu()
        elif choice == "2":
            return self._aa_usage_menu()
        elif choice == "3":
            return self._base_composition_menu()
        elif choice == "4":
            return self._rscu_menu()
        elif choice == "5":
            return self._multivariate_menu()
        elif choice == "6":
            return self._codon_bias_submenu()
        elif choice == "7":
            return self._optimal_codons_menu()
        elif choice == "8":
            return self._compare_codon_usage_axis_cohorts()
        elif choice == "9":
            return self._calculate_all_metrics()
        elif choice == "M":
            return self._quick_multivariate_workflow()
        elif choice == "C":
            return self._clear_cached_metrics()
        elif choice == "S":
            return self._show_analysis_status_detail()
        elif choice == "R":
            return False
        else:
            print("\nInvalid option. Please try again.")
            self._pause()
            return True
    
    def _rscu_menu(self):
        """Handle RSCU calculation and display."""
        if not self._check_data():
            return True
            
        print("\nRSCU (Relative Synonymous Codon Usage) Analysis")
        print("=" * 50)
        
        # Check if RSCU is already calculated
        if hasattr(self.data, 'rscu_values') and not self.data.rscu_values.empty:
            print("\nRSCU values already calculated.")
            recalc = input("\nRecalculate? (y/n): ").strip().lower()
            if recalc != 'y':
                # Show options for existing RSCU
                print("\nOptions:")
                print("1. Display RSCU values")
                print("2. Export RSCU to file")
                print("3. Return")
                
                choice = input("\nSelect option: ").strip()
                if choice == "1":
                    return self._display_rscu_values()
                elif choice == "2":
                    return self._export_rscu()
                else:
                    return True
        
        # Calculate RSCU
        print("\nCalculating RSCU values...")
        self.data.get_rscu_values()  # This will calculate if needed
        self.data.analysis_status['rscu_calculated'] = True
        
        print("\nRSCU values calculated successfully!")
        
        # Ask what to do next
        print("\nOptions:")
        print("1. Display RSCU values")
        print("2. Export RSCU to file")
        print("3. Return")
        
        choice = input("\nSelect option: ").strip()
        if choice == "1":
            return self._display_rscu_values()
        elif choice == "2":
            return self._export_rscu()
        
        return True
    
    def _show_analysis_status_detail(self):
        """Show detailed analysis status."""
        print("\nDetailed Analysis Status")
        print("=" * 50)
        
        # Basic metrics
        print("\n--- Basic Metrics ---")
        print(f"  Codon usage: {'✓ Calculated' if self.data.analysis_status.get('codon_usage', False) else '✗ Not calculated'}")
        print(f"  Amino acid usage: {'✓ Calculated' if self.data.analysis_status.get('amino_acid_usage', False) else '✗ Not calculated'}")
        print(f"  Base composition: {'✓ Calculated' if self.data.analysis_status.get('base_composition', False) else '✗ Not calculated'}")
        print(f"  RSCU values: {'✓ Calculated' if self.data.analysis_status.get('rscu_calculated', False) else '✗ Not calculated'}")
        
        # Advanced metrics
        print("\n--- Advanced Metrics ---")
        print(f"  ENC values: {'✓ Calculated' if self.data.analysis_status.get('enc_calculated', False) else '✗ Not calculated'}")
        print(f"  CAI values: {'✓ Calculated' if self.data.analysis_status.get('cai_calculated', False) else '✗ Not calculated'}")
        print(f"  Fop values: {'✓ Calculated' if self.data.analysis_status.get('fop_calculated', False) else '✗ Not calculated'}")
        print(f"  SCUO values: {'✓ Calculated' if self.data.analysis_status.get('scuo_calculated', False) else '✗ Not calculated'}")
        
        # Analysis results
        print("\n--- Analysis Results ---")
        print(f"  Optimal codons: {'✓ Identified' if self.data.analysis_status.get('optimal_codons', False) else '✗ Not identified'}")
        print(f"  Multivariate analysis: {'✓ Performed' if self.data.analysis_status.get('multivariate_performed', False) else '✗ Not performed'}")
        if hasattr(self.data, 'multivariate_history') and self.data.multivariate_history:
            print(f"    - Total analyses: {len(self.data.multivariate_history)}")
            for i, analysis in enumerate(self.data.multivariate_history[-3:], 1):
                print(f"    - {analysis['analysis_type']} of {analysis['data_type']} data")
        
        print(f"  Clusters detected: {'✓ Yes' if self.data.analysis_status.get('clusters_detected', False) else '✗ No'}")
        print(f"  Axis cohorts compared: {'✓ Yes' if self.data.analysis_status.get('axis_cohorts_compared', False) else '✗ No'}")
        
        # Memory usage info
        if hasattr(self.data, 'memory_mode'):
            print(f"\n--- Memory Status ---")
            print(f"  Mode: {self.data.memory_mode}")
            if hasattr(self.data, 'codon_counts'):
                memory_mb = self.data.codon_counts.memory_usage(deep=True).sum() / 1024 / 1024
                print(f"  Codon counts memory: {memory_mb:.1f} MB")
        
        self._pause()
        return True
    
    def _display_rscu_values(self):
        """Display RSCU values for genes."""
        if not hasattr(self.data, 'rscu_values') or self.data.rscu_values.empty:
            print("\nRSCU values not calculated yet.")
            self._pause()
            return True
            
        print("\nRSCU Values")
        print("=" * 50)
        
        # Show first few genes as example
        print("\nShowing RSCU values for first 5 genes:")
        print("\nFormat: Codon (Amino Acid): RSCU value")
        
        aa_map = self.config.get_aa_map()
        
        for i, gene in enumerate(self.data.rscu_values.index[:5]):
            print(f"\n{gene}:")
            print("-" * 40)
            
            # Group by amino acid
            aa_groups = {}
            for codon in CODONS:
                aa = aa_map[codon]
                if aa != 'STOP':
                    if aa not in aa_groups:
                        aa_groups[aa] = []
                    aa_groups[aa].append(codon)
            
            # Display by amino acid
            for aa in sorted(aa_groups.keys()):
                print(f"\n  {aa}:")
                for codon in sorted(aa_groups[aa]):
                    rscu = self.data.rscu_values.at[gene, codon]
                    if rscu > 0:  # Only show codons that are used
                        print(f"    {codon}: {rscu:.3f}")
            
            if i < 4:  # Don't ask after the last gene
                cont = input("\nPress Enter to continue or 'q' to quit: ").strip().lower()
                if cont == 'q':
                    break
        
        # Summary statistics
        print("\n\nSummary Statistics:")
        print("-" * 40)
        
        # Find most and least preferred codons overall
        mean_rscu = self.data.rscu_values.mean()
        
        # Exclude stop codons
        non_stop_codons = [c for c in CODONS if aa_map[c] != 'STOP']
        mean_rscu_filtered = mean_rscu[non_stop_codons]
        
        top_codons = mean_rscu_filtered.nlargest(10)
        bottom_codons = mean_rscu_filtered.nsmallest(10)
        
        print("\nMost preferred codons (mean RSCU):")
        for codon, rscu in top_codons.items():
            print(f"  {codon} ({aa_map[codon]}): {rscu:.3f}")
            
        print("\nLeast preferred codons (mean RSCU):")
        for codon, rscu in bottom_codons.items():
            print(f"  {codon} ({aa_map[codon]}): {rscu:.3f}")
        
        self._pause()
        return True
    
    def _clear_cached_metrics(self):
        """Clear all cached metric calculations."""
        if not self._check_data():
            return True
            
        print("\nCurrent cached metrics:")
        cached = []
        
        if hasattr(self.data, 'multivariate_results') and self.data.multivariate_results is not None:
            cached.append("Multivariate analysis (CA on RSCU)")
        if hasattr(self.data, 'multivariate_history') and self.data.multivariate_history:
            cached.append(f"Multivariate history ({len(self.data.multivariate_history)} analyses)")
        if hasattr(self.data, 'enc_values') and self.data.enc_values is not None:
            cached.append("ENC values")
        if hasattr(self.data, 'cai_values') and self.data.cai_values is not None:
            cached.append("CAI values")
        if hasattr(self.data, 'fop_values') and self.data.fop_values is not None:
            cached.append("Fop values")
        if hasattr(self.data, 'optimal_codons') and self.data.optimal_codons is not None:
            cached.append("Optimal codons")
            
        if not cached:
            print("No cached metrics found.")
            self._pause()
            return True
            
        print("\n".join(f"  • {metric}" for metric in cached))
        
        confirm = input("\nAre you sure you want to clear all cached metrics? (y/n): ").strip().lower()
        if confirm == 'y':
            # Clear all cached values
            if hasattr(self.data, 'multivariate_results'):
                self.data.multivariate_results = None
            if hasattr(self.data, 'multivariate_history'):
                self.data.multivariate_history = []
            if hasattr(self.data, 'enc_values'):
                self.data.enc_values = None
            if hasattr(self.data, 'cai_values'):
                self.data.cai_values = None
            if hasattr(self.data, 'fop_values'):
                self.data.fop_values = None
            if hasattr(self.data, 'optimal_codons'):
                self.data.optimal_codons = None
            if hasattr(self.data, 'reference_genes'):
                self.data.reference_genes = None
                
            print("\nAll cached metrics have been cleared.")
            print("Metrics will be recalculated when needed.")
        else:
            print("\nOperation cancelled.")
            
        self._pause()
        return True

    def _codon_usage_menu(self):
        """Handle codon usage output options."""
        if not self._check_data():
            return True

        options = [
            ("Display codon usage for each gene",
             "1", self._display_gene_codon_usage),
            ("Display cumulative codon usage", "2",
             self._display_cumulative_codon_usage),
            ("Export codon usage to file", "3", self._export_codon_usage),
            ("Return to analysis menu", "R", lambda: False)
        ]

        return self._display_menu("Codon Usage Analysis", options)
    
    def _display_gene_codon_usage(self):
        """Display codon usage for each gene."""
        print("\nCodon Usage by Gene:")
        print("(Press 'q' to quit at any time)")
        
        # Display one gene at a time
        for gene in self.data.gene_names:
            print(f"\n{'='*60}")
            print(f"Gene: {gene}")
            print(f"{'='*60}")
            
            # Get codon counts for this gene
            gene_codons = self.data.codon_counts.loc[gene]
            
            # Display by amino acid
            aa_map = self.config.get_aa_map()
            aa_to_codons = {}
            for codon in CODONS:
                aa = aa_map[codon]
                if aa not in aa_to_codons:
                    aa_to_codons[aa] = []
                aa_to_codons[aa].append(codon)
            
            # Sort amino acids
            for aa in sorted(aa_to_codons.keys()):
                if aa == 'STOP':
                    continue
                print(f"\n{aa}:")
                for codon in sorted(aa_to_codons[aa]):
                    count = gene_codons.get(codon, 0)
                    print(f"  {codon}: {int(count)}")
            
            # Ask to continue
            more = input("\nPress Enter for next gene, or 'q' to quit: ").strip().lower()
            if more == 'q':
                break
        
        self._pause()
        return True
    
    def _display_cumulative_codon_usage(self):
        """Display cumulative codon usage for all genes."""
        print("\nCumulative Codon Usage (All Genes):")
        
        # Calculate totals
        total_codons = self.data.codon_counts.sum()
        total_count = total_codons.sum()
        
        print(f"\nTotal codons analyzed: {int(total_count)}")
        print(f"Total genes: {len(self.data.gene_names)}")
        
        # Display by amino acid
        aa_map = self.config.get_aa_map()
        aa_to_codons = {}
        for codon in CODONS:
            aa = aa_map[codon]
            if aa not in aa_to_codons:
                aa_to_codons[aa] = []
            aa_to_codons[aa].append(codon)
        
        print("\nCodon\tCount\tFrequency")
        print("-" * 30)
        
        for aa in sorted(aa_to_codons.keys()):
            if aa == 'STOP':
                continue
            print(f"\n{aa}:")
            for codon in sorted(aa_to_codons[aa]):
                count = total_codons.get(codon, 0)
                freq = count / total_count if total_count > 0 else 0
                print(f"  {codon}\t{int(count)}\t{freq:.4f}")
        
        self._pause()
        return True
    
    def _export_codon_usage(self):
        """Export codon usage to file."""
        output_path = self.file_manager.get_output_path(
            self.data, None, "codon_usage", ".tsv")
        
        print(f"\nExporting codon usage to {output_path}...")
        
        # Export as matrix (genes x codons)
        self.data.codon_counts.to_csv(output_path, sep='\t')
        
        print(f"Codon usage exported to: {output_path}")
        self._pause()
        return True

    def _aa_usage_menu(self):
        """Handle amino acid usage output options."""
        if not self._check_data():
            return True

        options = [
            ("Display amino acid usage for each gene",
             "1", self._display_gene_aa_usage),
            ("Display cumulative amino acid usage",
             "2", self._display_cumulative_aa_usage),
            ("Export amino acid usage to file", "3", self._export_aa_usage),
            ("Return to analysis menu", "R", lambda: False)
        ]

        return self._display_menu("Amino Acid Usage Analysis", options)
    
    def _display_gene_aa_usage(self):
        """Display amino acid usage for each gene."""
        print("\nAmino Acid Usage by Gene:")
        
        # Ensure amino acid usage is calculated
        aa_usage = self.data.get_amino_acid_usage()
        
        # Display header
        amino_acids = sorted([col for col in aa_usage.columns if col != 'STOP'])
        header = "Gene\t" + "\t".join(amino_acids)
        print(header)
        print("-" * (len(header) + 10))
        
        # Display in batches of 20
        count = 0
        for gene in self.data.gene_names:
            row = [gene]
            for aa in amino_acids:
                count_val = aa_usage.loc[gene, aa]
                row.append(str(int(count_val)))
            print("\t".join(row))
            
            count += 1
            if count % 20 == 0:
                more = input("\nPress Enter to continue, or 'q' to quit: ").strip().lower()
                if more == 'q':
                    break
        
        self._pause()
        return True
    
    def _display_cumulative_aa_usage(self):
        """Display cumulative amino acid usage for all genes."""
        print("\nCumulative Amino Acid Usage (All Genes):")
        
        # Ensure amino acid usage is calculated
        aa_usage = self.data.get_amino_acid_usage()
        
        # Calculate totals
        total_aa = aa_usage.sum()
        total_count = total_aa.sum()
        
        print(f"\nTotal amino acids: {int(total_count)}")
        print(f"Total genes: {len(self.data.gene_names)}")
        
        print("\nAmino Acid\tCount\tFrequency\tPer Gene")
        print("-" * 50)
        
        for aa in sorted(total_aa.index):
            if aa == 'STOP':
                continue
            count = total_aa[aa]
            freq = count / total_count if total_count > 0 else 0
            per_gene = count / len(self.data.gene_names)
            print(f"{aa}\t\t{int(count)}\t{freq:.4f}\t{per_gene:.2f}")
        
        self._pause()
        return True
    
    def _export_aa_usage(self):
        """Export amino acid usage to file."""
        output_path = self.file_manager.get_output_path(
            self.data, None, "amino_acid_usage", ".tsv")
        
        print(f"\nExporting amino acid usage to {output_path}...")
        
        # Ensure amino acid usage is calculated
        aa_usage = self.data.get_amino_acid_usage()
        
        # Export as matrix (genes x amino acids)
        aa_usage.to_csv(output_path, sep='\t')
        
        print(f"Amino acid usage exported to: {output_path}")
        self._pause()
        return True

    def _base_composition_menu(self):
        """Handle base composition output options."""
        if not self._check_data():
            return True

        options = [
            ("Display base composition for each gene",
             "1", self._display_gene_base_comp),
            ("Display cumulative base composition",
             "2", self._display_cumulative_base_comp),
            ("Export base composition to file", "3", self._export_base_comp),
            ("Return to analysis menu", "R", lambda: False)
        ]

        return self._display_menu("Base Composition Analysis", options)
    
    def _display_gene_base_comp(self):
        """Display base composition for each gene."""
        print("\nBase Composition by Gene:")
        print("Gene\tGC%\tGC1%\tGC2%\tGC3%\tGC3s%")
        
        # Display in batches of 20
        count = 0
        for gene in self.data.gene_names:
            bc = self.data.base_composition.loc[gene]
            print(f"{gene}\t{bc['GC']:.1f}\t{bc['GC1']:.1f}\t{bc['GC2']:.1f}\t{bc['GC3']:.1f}\t{bc.get('GC3s', 0):.1f}")
            
            count += 1
            if count % 20 == 0:
                more = input("\nPress Enter to continue, or 'q' to quit: ").strip().lower()
                if more == 'q':
                    break
        
        self._pause()
        return True
    
    def _display_cumulative_base_comp(self):
        """Display cumulative base composition for all genes."""
        print("\nCumulative Base Composition:")
        
        # Calculate cumulative statistics
        bc = self.data.base_composition
        print(f"\nTotal genes analyzed: {len(self.data.gene_names)}")
        print(f"Average GC content: {bc['GC'].mean():.2f}%")
        print(f"GC range: {bc['GC'].min():.2f}% - {bc['GC'].max():.2f}%")
        print(f"Standard deviation: {bc['GC'].std():.2f}%")
        
        print(f"\nAverage GC1: {bc['GC1'].mean():.2f}%")
        print(f"Average GC2: {bc['GC2'].mean():.2f}%")
        print(f"Average GC3: {bc['GC3'].mean():.2f}%")
        if 'GC3s' in bc.columns:
            print(f"Average GC3s: {bc['GC3s'].mean():.2f}%")
        
        self._pause()
        return True
    
    def _export_base_comp(self):
        """Export base composition to file."""
        output_path = self.file_manager.get_output_path(
            self.data, None, "base_composition", ".tsv")
        
        print(f"\nExporting base composition to {output_path}...")
        
        with open(output_path, 'w') as f:
            f.write("Gene\tGC%\tGC1%\tGC2%\tGC3%\tGC3s%\n")
            for gene in self.data.gene_names:
                bc = self.data.base_composition.loc[gene]
                f.write(f"{gene}\t{bc['GC']:.2f}\t{bc['GC1']:.2f}\t{bc['GC2']:.2f}\t{bc['GC3']:.2f}\t{bc.get('GC3s', 0):.2f}\n")
        
        print(f"Base composition exported to: {output_path}")
        self._pause()
        return True

    def _multivariate_menu(self):
        """Handle multivariate analysis options."""
        if not self._check_data() or len(self.data.gene_names) < 2:
            print("\nThere must be more than a single sequence in memory")
            self._pause()
            return True
        
        # Show status of performed analyses
        if hasattr(self.data, 'multivariate_history') and self.data.multivariate_history:
            print("\nPerformed analyses:")
            for analysis in self.data.multivariate_history:
                print(f"  • {analysis['analysis_type']} of {analysis['data_type']} data")
        
        options = [
            ("Perform Correspondence Analysis (CA) on RSCU values",
             "1", lambda: self._perform_multivariate('CA', 'RSCU')),
            ("Perform Principal Component Analysis (PCA) on RSCU values",
             "2", lambda: self._perform_multivariate('PCA', 'RSCU')),
            ("Perform CA on amino acid usage", "3",
             lambda: self._perform_multivariate('CA', 'AA')),
            ("Perform PCA on amino acid usage", "4",
             lambda: self._perform_multivariate('PCA', 'AA')),
            ("View/Visualize multivariate analysis results", "5", self._multivariate_plot),
            ("Export multivariate analysis results",
             "6", self._export_multivariate),
            ("Return to analysis menu", "R", lambda: False)
        ]

        return self._display_menu("Multivariate Analysis", options)

    def _perform_multivariate(self, analysis_type, data_type):
        """Perform multivariate analysis and show summary."""
        print(f"\nPerforming {analysis_type} on {data_type} values...")
        result = self.data.perform_multivariate_analysis(
            analysis_type, data_type)
        self._show_multivariate_summary(result, analysis_type, data_type)
        return True

    def _show_multivariate_summary(self, result, analysis_type, data_type):
        """Show a summary of multivariate analysis results."""
        if not result:
            print("Analysis failed. No results to display.")
            self._pause()
            return True

        print(f"\n{analysis_type} of {data_type} data - Summary")
        print("\nExplained variance:")
        for i, var in enumerate(result['explained_variance']):
            print(
                f"  {result['coordinates'].columns[i]}: {var:.4f} ({var*100:.2f}%)")

        # Show top 10 genes on first two axes
        print("\nTop 10 genes on the first axis:")
        sorted_genes = result['coordinates'].sort_values(
            by=result['coordinates'].columns[0], ascending=False)
        for i, (gene, row) in enumerate(sorted_genes.iterrows()):
            if i < 5:
                print(f"  {gene}: {row.iloc[0]:.4f}")
            if i >= len(sorted_genes) - 5:
                print(f"  {gene}: {row.iloc[0]:.4f}")
            if i == 4 and len(sorted_genes) > 10:
                print("  ...")

        # Add genetic code information
        print(
            f"\nGenetic code used: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        
        print(f"\nAnalysis completed successfully!")
        print("You can visualize or export these results from the multivariate analysis menu.")
        
        self._pause()
        return True

    def _visualize_multivariate(self, result, analysis_type, data_type):
        """Visualize multivariate analysis results."""
        # Add analysis type and data type to result
        if 'analysis_type' not in result:
            result['analysis_type'] = analysis_type
        if 'data_type' not in result:
            result['data_type'] = data_type

        # Create visualization using the data's visualization method
        output_path = self.file_manager.get_output_path(
            self.data, None, f"{analysis_type}_{data_type}_plot")
        self.data.visualize('multivariate', result, output_path)

    def _enc_menu(self):
        """Handle ENC calculation options."""
        if not self._check_data():
            return True

        options = [
            ("Calculate and display ENC values", "1", self._display_enc_values),
            ("Create ENC vs GC3s plot (Wright's plot)", "2", self._enc_plot),
            ("Export ENC values to file", "3", self._export_enc_values),
            ("Return to analysis menu", "R", lambda: False)
        ]

        return self._display_menu(
            "Effective Number of Codons (ENC) Analysis", options)
    
    def _display_enc_values(self):
        """Calculate and display ENC values."""
        print("\nCalculating ENC (Effective Number of Codons)...")
        enc_values = self.data.calculate_enc()
        
        print("\nENC Values:")
        print("Gene\tENC")
        print("-" * 30)
        
        # Sort genes by ENC value
        sorted_genes = sorted(enc_values.keys(), key=lambda g: enc_values[g])
        
        # Display in batches of 20
        count = 0
        for gene in sorted_genes:
            print(f"{gene}\t{enc_values[gene]:.2f}")
            
            count += 1
            if count % 20 == 0:
                more = input("\nPress Enter to continue, or 'q' to quit: ").strip().lower()
                if more == 'q':
                    break
        
        # Show summary statistics
        if len(enc_values) > 1:
            enc_vals = list(enc_values.values())
            print(f"\nSummary Statistics:")
            print(f"  Mean ENC: {np.mean(enc_vals):.2f}")
            print(f"  Median ENC: {np.median(enc_vals):.2f}")
            print(f"  Min ENC: {min(enc_vals):.2f}")
            print(f"  Max ENC: {max(enc_vals):.2f}")
            print(f"  Std Dev: {np.std(enc_vals):.2f}")
        
        self._pause()
        return True
    
    def _export_enc_values(self):
        """Export ENC values to file."""
        print("\nCalculating ENC values...")
        enc_values = self.data.calculate_enc()
        
        output_path = self.file_manager.get_output_path(
            self.data, None, "enc_values", ".tsv")
        
        print(f"Exporting ENC values to {output_path}...")
        
        with open(output_path, 'w') as f:
            f.write("Gene\tENC\n")
            for gene in self.data.gene_names:
                enc = enc_values.get(gene, 'NA')
                f.write(f"{gene}\t{enc}\n")
        
        print(f"ENC values exported to: {output_path}")
        self._pause()
        return True

    def _fop_cai_menu(self):
        """Handle Fop and CAI calculation options."""
        if not self._check_data():
            return True

        options = [
            ("Use all genes as reference", "1",
             lambda: self._calculate_and_display_fop_cai(None)),
            ("Select reference genes manually", "2",
             self._calculate_fop_cai_with_manual_ref),
            ("Use multivariate analysis to identify reference genes",
             "3", self._calculate_and_display_fop_cai_multivariate),
            ("Display Fop and CAI values (using current reference)",
             "4", lambda: self._display_fop_cai()),
            ("Create CAI distribution plot", "5", self._cai_distribution_plot),
            ("Export Fop and CAI values", "6", self._export_fop_cai),
            ("Return to analysis menu", "R", lambda: False)
        ]

        return self._display_menu("Fop and CAI Analysis", options)

    def _calculate_fop_cai_with_manual_ref(self):
        """Calculate Fop and CAI using manually selected reference genes."""
        reference_genes = self._select_reference_genes()
        if reference_genes:
            self._calculate_and_display_fop_cai(reference_genes)
        return True

    def _display_fop_cai(self):
        """Display Fop and CAI values."""
        # Check if values are calculated
        if not hasattr(self.data, 'fop_values') or not self.data.fop_values:
            print("\nFop values not calculated. Please calculate Fop values first.")
            self._pause()
            return True

        if not hasattr(self.data, 'cai_values') or not self.data.cai_values:
            print("\nCAI values not calculated. Please calculate CAI values first.")
            self._pause()
            return True

        # Display values sorted by CAI
        print("\nFrequency of Optimal Codons (Fop) and Codon Adaptation Index (CAI):")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print("\nGene\tFop\tCAI")

        # Sort genes by CAI (highest first)
        sorted_genes = sorted(
            self.data.gene_names,
            key=lambda g: self.data.cai_values[g],
            reverse=True)

        for gene in sorted_genes:
            print(
                f"{gene}\t{self.data.fop_values[gene]:.4f}\t{self.data.cai_values[gene]:.4f}")

            # Display in batches of 20
            if sorted_genes.index(gene) > 0 and (
                    sorted_genes.index(gene) + 1) % 20 == 0:
                more = input(
                    "\nPress Enter to continue, or 'q' to quit: ").strip().lower()
                if more == 'q':
                    break

        self._pause()
        return True

    def _calculate_and_display_fop_cai(self, reference_genes=None):
        """Calculate and display Fop and CAI values."""
        print("\nCalculating Frequency of Optimal Codons (Fop)...")
        fop_values = self.data.calculate_fop(reference_genes)

        print("Calculating Codon Adaptation Index (CAI)...")
        cai_values = self.data.calculate_cai(reference_genes)

        # Display values sorted by CAI
        print("\nFrequency of Optimal Codons (Fop) and Codon Adaptation Index (CAI):")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print("Gene\tFop\tCAI")

        # Sort genes by CAI (highest first)
        sorted_genes = sorted(
            self.data.gene_names,
            key=lambda g: cai_values[g],
            reverse=True)

        for gene in sorted_genes:
            print(f"{gene}\t{fop_values[gene]:.4f}\t{cai_values[gene]:.4f}")

            # Display in batches of 20
            if sorted_genes.index(gene) > 0 and (
                    sorted_genes.index(gene) + 1) % 20 == 0:
                more = input(
                    "\nPress Enter to continue, or 'q' to quit: ").strip().lower()
                if more == 'q':
                    break

        # Show summary statistics
        if len(fop_values) > 1:
            fop_list = list(fop_values.values())
            cai_list = list(cai_values.values())

            avg_fop = sum(fop_list) / len(fop_list)
            avg_cai = sum(cai_list) / len(cai_list)

            print("\nSummary statistics:")
            print(f"Average Fop: {avg_fop:.4f}")
            print(f"Average CAI: {avg_cai:.4f}")
            print(f"Maximum Fop: {max(fop_list):.4f}")
            print(f"Maximum CAI: {max(cai_list):.4f}")
            print(f"Minimum Fop: {min(fop_list):.4f}")
            print(f"Minimum CAI: {min(cai_list):.4f}")

        # Ask if user wants to export
        export = input(
            "\nDo you want to export these values to a file? (y/n): ").strip().lower()
        if export == 'y':
            self._export_fop_cai(fop_values, cai_values, reference_genes)

        self._pause()
        return True

    def _calculate_and_display_fop_cai_multivariate(self):
        """Calculate and display Fop and CAI values using multivariate analysis."""
        print("\nUsing multivariate analysis to select reference genes...")

        # Ask for percentage of genes to use
        try:
            percentage = float(input(
                "\nEnter percentage of genes to use from multivariate analysis (5-20 recommended): ").strip())
            if percentage < 1 or percentage > 50:
                print("Invalid percentage. Using default of 10%.")
                percentage = 10
        except ValueError:
            print("Invalid input. Using default of 10%.")
            percentage = 10

        print(
            "\nCalculating Frequency of Optimal Codons (Fop) using multivariate analysis...")
        fop_values = self.data.calculate_fop_multivariate(percentage)

        print("Calculating Codon Adaptation Index (CAI) using multivariate analysis...")
        cai_values = self.data.calculate_cai_multivariate(percentage)

        # Display values sorted by CAI
        print("\nFrequency of Optimal Codons (Fop) and Codon Adaptation Index (CAI):")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print("Gene\tFop\tCAI")

        # Sort genes by CAI (highest first)
        sorted_genes = sorted(
            self.data.gene_names,
            key=lambda g: cai_values[g],
            reverse=True)

        for gene in sorted_genes:
            print(f"{gene}\t{fop_values[gene]:.4f}\t{cai_values[gene]:.4f}")

            # Display in batches of 20
            if sorted_genes.index(gene) > 0 and (
                    sorted_genes.index(gene) + 1) % 20 == 0:
                more = input(
                    "\nPress Enter to continue, or 'q' to quit: ").strip().lower()
                if more == 'q':
                    break

        # Show summary statistics
        if len(fop_values) > 1:
            fop_list = list(fop_values.values())
            cai_list = list(cai_values.values())

            avg_fop = sum(fop_list) / len(fop_list)
            avg_cai = sum(cai_list) / len(cai_list)

            print("\nSummary statistics:")
            print(f"Average Fop: {avg_fop:.4f}")
            print(f"Average CAI: {avg_cai:.4f}")
            print(f"Maximum Fop: {max(fop_list):.4f}")
            print(f"Maximum CAI: {max(cai_list):.4f}")
            print(f"Minimum Fop: {min(fop_list):.4f}")
            print(f"Minimum CAI: {min(cai_list):.4f}")

        # Ask if user wants to export
        export = input(
            "\nDo you want to export these values to a file? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, "fop_cai_multivariate") + '.tsv'

            with open(output_path, 'w') as f:
                # Add metadata
                f.write("# GCUA Fop and CAI Values (Multivariate Analysis)\n")
                f.write(f"# Version: {VERSION}\n")
                f.write(
                    f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
                f.write(f"# Multivariate Percentage: {percentage}%\n")
                f.write(f"# Total Genes: {len(self.data.gene_names)}\n")
                f.write(
                    f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"#{'=' * 50}\n\n")

                f.write("Gene\tFop\tCAI\n")
                for gene in sorted(self.data.gene_names):
                    f.write(
                        f"{gene}\t{fop_values[gene]:.4f}\t{cai_values[gene]:.4f}\n")

            print(f"\nFop and CAI values saved to {output_path}")

        self._pause()
        return True

    def _scuo_menu(self):
        """Handle SCUO calculation options."""
        if not self._check_data():
            return True

        options = [
            ("Calculate and display SCUO values", "1", self._display_scuo_values),
            ("Export SCUO values to file", "2", self._export_scuo_values),
            ("Return to analysis menu", "R", lambda: False)
        ]

        return self._display_menu(
            "Synonymous Codon Usage Order (SCUO) Analysis", options)

    def _display_scuo_values(self):
        """Calculate and display SCUO values."""
        print("\nCalculating Synonymous Codon Usage Order (SCUO)...")
        scuo_values = self.data.calculate_scuo()

        # Display values
        print("\nSynonymous Codon Usage Order (SCUO):")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print("Gene\tSCUO")

        # Sort genes by SCUO (highest first)
        sorted_genes = sorted(
            scuo_values.keys(),
            key=lambda g: scuo_values[g],
            reverse=True)

        for gene in sorted_genes:
            print(f"{gene}\t{scuo_values[gene]:.4f}")

            # Display in batches of 20
            if sorted_genes.index(gene) > 0 and (
                    sorted_genes.index(gene) + 1) % 20 == 0:
                more = input(
                    "\nPress Enter to continue, or 'q' to quit: ").strip().lower()
                if more == 'q':
                    break

        # Show summary statistics
        if len(scuo_values) > 1:
            scuo_list = list(scuo_values.values())
            avg_scuo = sum(scuo_list) / len(scuo_list)
            min_scuo = min(scuo_list)
            max_scuo = max(scuo_list)

            print("\nSummary statistics:")
            print(f"Average SCUO: {avg_scuo:.4f}")
            print(f"Minimum SCUO: {min_scuo:.4f}")
            print(f"Maximum SCUO: {max_scuo:.4f}")

        self._pause()
        return True

    def _export_scuo_values(self):
        """Export SCUO values to a file."""
        print("\nCalculating Synonymous Codon Usage Order (SCUO)...")
        scuo_values = self.data.calculate_scuo()

        output_path = self.file_manager.get_output_path(
            self.data, None, 'scuo') + '.tsv'

        with open(output_path, 'w') as f:
            # Add metadata
            f.write("# GCUA SCUO Values\n")
            f.write(f"# Version: {VERSION}\n")
            f.write(
                f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
            f.write(f"# Total Genes: {len(self.data.gene_names)}\n")
            f.write(
                f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"#{'=' * 50}\n\n")

            f.write("Gene\tSCUO\n")
            for gene in sorted(self.data.gene_names):
                f.write(f"{gene}\t{scuo_values[gene]:.4f}\n")

        print(f"\nSCUO values saved to {output_path}")
        self._pause()
        return True
    
    def _export_extreme_genes_menu(self):
        """Handle export of extreme genes from various metrics."""
        if not self._check_data():
            return True
        
        # Build list of available metrics
        available_metrics = []
        option_num = 1
        
        print("\n" + "="*60)
        print("BASE METRICS:")
        print("="*60)
        
        # Base metrics always available
        available_metrics.append(("GC content", str(option_num), lambda: self._export_extreme_by_metric("GC", "Base Composition")))
        option_num += 1
        available_metrics.append(("GC3s content", str(option_num), lambda: self._export_extreme_by_metric("GC3s", "Base Composition")))
        option_num += 1
        
        # Check for calculated metrics
        if hasattr(self.data, 'enc_values') and self.data.enc_values:
            available_metrics.append(("ENC values", str(option_num), lambda: self._export_extreme_by_metric("ENC", "Codon Bias", self.data.enc_values)))
            option_num += 1
        
        if hasattr(self.data, 'cai_values') and self.data.cai_values:
            available_metrics.append(("CAI values", str(option_num), lambda: self._export_extreme_by_metric("CAI", "Codon Adaptation", self.data.cai_values)))
            option_num += 1
        
        if hasattr(self.data, 'fop_values') and self.data.fop_values:
            available_metrics.append(("Fop values", str(option_num), lambda: self._export_extreme_by_metric("Fop", "Optimal Codon Frequency", self.data.fop_values)))
            option_num += 1
        
        # Check for SCUO
        available_metrics.append(("SCUO values", str(option_num), lambda: self._export_extreme_by_metric("SCUO", "Codon Usage Order", self.data.calculate_scuo())))
        option_num += 1
        
        # Check for multivariate results in history
        if hasattr(self.data, 'multivariate_history') and self.data.multivariate_history:
            print("\nMULTIVARIATE ANALYSES:")
            print("="*60)
            
            for analysis_idx, analysis in enumerate(self.data.multivariate_history):
                analysis_type = analysis['analysis_type']
                data_type = analysis['data_type']
                result = analysis['result']
                
                if 'explained_variance' in result and 'coordinates' in result:
                    # Add axes from this analysis
                    n_axes = min(3, len(result['explained_variance']))  # Show up to 3 axes per analysis
                    for i in range(n_axes):
                        axis_name = f'Axis{i+1}' if analysis_type == 'CA' else f'PC{i+1}'
                        if axis_name in result['coordinates'].columns:
                            variance_pct = result['explained_variance'][i] * 100
                            label = f"{analysis_type} on {data_type} - {axis_name} ({variance_pct:.1f}%)"
                            
                            # Create closure with proper variable capture
                            def make_export_func(anal_type, data_t, ax_name, coords):
                                return lambda: self._export_extreme_by_metric(
                                    f"{anal_type}_{data_t}_{ax_name}", 
                                    f"{anal_type} on {data_t}",
                                    coords[ax_name]
                                )
                            
                            available_metrics.append((
                                label, 
                                str(option_num),
                                make_export_func(analysis_type, data_type, axis_name, result['coordinates'])
                            ))
                            option_num += 1
        
        available_metrics.append(("Return to analysis menu", "R", lambda: False))
        
        print("\nNote: Percentages show explained variance for each axis.")
        return self._display_menu("Export Extreme Genes", available_metrics)
    
    def _export_extreme_by_metric(self, metric_name, metric_type, metric_values=None):
        """Export genes from extremes of a given metric distribution."""
        # Create a display name for the metric
        if '_' in metric_name and metric_name.count('_') >= 2:
            # It's a multivariate metric like "CA_RSCU_Axis1"
            parts = metric_name.split('_')
            display_name = f"{parts[0]} on {parts[1]} - {parts[2]}"
        else:
            display_name = metric_name
            
        print(f"\nExporting extreme genes by {display_name}")
        
        # Get metric values if not provided
        if metric_values is None:
            if metric_name in ["GC", "GC3s"]:
                metric_values = self.data.base_composition[metric_name]
            else:
                print(f"Error: {metric_name} values not available")
                self._pause()
                return True
        
        # Convert dict to Series if needed
        if isinstance(metric_values, dict):
            metric_values = pd.Series(metric_values)
        
        # Ensure index matches gene names
        if isinstance(metric_values, pd.Series):
            metric_values = metric_values.reindex(self.data.gene_names)
        
        # Remove any NaN values
        valid_mask = ~metric_values.isna()
        valid_values = metric_values[valid_mask]
        valid_genes = [g for g, v in zip(self.data.gene_names, valid_mask) if v]
        
        if len(valid_values) == 0:
            print(f"No valid {metric_name} values found.")
            self._pause()
            return True
        
        # Get percentage from user
        print(f"\nTotal genes with valid {display_name} values: {len(valid_values)}")
        print(f"Range: {valid_values.min():.4f} to {valid_values.max():.4f}")
        print(f"Mean: {valid_values.mean():.4f}, Median: {valid_values.median():.4f}")
        
        while True:
            try:
                percentage = float(input("\nEnter percentage of genes to export from each extreme (e.g., 10): "))
                if 0 < percentage <= 50:
                    break
                else:
                    print("Please enter a value between 0 and 50.")
            except ValueError:
                print("Invalid input. Please enter a number.")
        
        # Calculate number of genes
        n_genes = max(1, int(len(valid_values) * percentage / 100))
        
        # Get export preference
        print(f"\nWill export {n_genes} genes from extremes.")
        print("Export options:")
        print("1. Top genes only (highest values)")
        print("2. Bottom genes only (lowest values)")
        print("3. Both extremes (separate files)")
        
        choice = input("\nSelect option (1-3): ").strip()
        
        # Sort genes by metric
        sorted_indices = valid_values.sort_values(ascending=True).index
        
        # Prepare export based on choice
        exports = []
        if choice in ['1', '3']:
            top_genes = sorted_indices[-n_genes:].tolist()
            exports.append(('top', top_genes, valid_values[top_genes]))
        
        if choice in ['2', '3']:
            bottom_genes = sorted_indices[:n_genes].tolist()
            exports.append(('bottom', bottom_genes, valid_values[bottom_genes]))
        
        # Export files
        for extreme_type, gene_list, gene_values in exports:
            output_path = self.file_manager.get_output_path(
                self.data, None, f"extreme_genes_{metric_name}_{extreme_type}_{int(percentage)}pct", ".tsv")
            
            with open(output_path, 'w') as f:
                # Write header with metadata
                f.write(f"# Extreme genes export\n")
                f.write(f"# Metric: {metric_name} ({metric_type})\n")
                f.write(f"# Selection: {extreme_type} {percentage}%\n")
                f.write(f"# Total genes analyzed: {len(valid_values)}\n")
                f.write(f"# Number of genes exported: {len(gene_list)}\n")
                f.write(f"# Threshold: {extreme_type} = {'>' if extreme_type == 'top' else '<'} {gene_values.iloc[0 if extreme_type == 'bottom' else -1]:.4f}\n")
                f.write("#\n")
                
                # Write column headers
                headers = ["Gene", display_name, "Rank", "Percentile"]
                
                # Add other available metrics
                if hasattr(self.data, 'enc_values') and self.data.enc_values:
                    headers.append("ENC")
                if hasattr(self.data, 'cai_values') and self.data.cai_values:
                    headers.append("CAI")
                if hasattr(self.data, 'fop_values') and self.data.fop_values:
                    headers.append("Fop")
                headers.extend(["GC", "GC3s"])
                
                f.write("\t".join(headers) + "\n")
                
                # Write gene data
                for i, (gene, value) in enumerate(zip(gene_list, gene_values)):
                    rank = i + 1 if extreme_type == 'bottom' else len(valid_values) - len(gene_list) + i + 1
                    percentile = (rank / len(valid_values)) * 100
                    
                    row = [gene, f"{value:.4f}", str(rank), f"{percentile:.1f}"]
                    
                    # Add other metrics
                    if hasattr(self.data, 'enc_values') and self.data.enc_values:
                        enc_val = self.data.enc_values.get(gene, 'NA')
                        row.append(f"{enc_val:.2f}" if enc_val != 'NA' else 'NA')
                    
                    if hasattr(self.data, 'cai_values') and self.data.cai_values:
                        cai_val = self.data.cai_values.get(gene, 'NA')
                        row.append(f"{cai_val:.4f}" if cai_val != 'NA' else 'NA')
                    
                    if hasattr(self.data, 'fop_values') and self.data.fop_values:
                        fop_val = self.data.fop_values.get(gene, 'NA')
                        row.append(f"{fop_val:.4f}" if fop_val != 'NA' else 'NA')
                    
                    # Add base composition
                    gc = self.data.base_composition.loc[gene, 'GC']
                    gc3s = self.data.base_composition.loc[gene, 'GC3s']
                    row.extend([f"{gc:.2f}", f"{gc3s:.2f}"])
                    
                    f.write("\t".join(row) + "\n")
            
            print(f"\nExported {len(gene_list)} {extreme_type} genes to: {output_path}")
        
        # Summary statistics
        print(f"\nSummary for exported genes:")
        for extreme_type, gene_list, gene_values in exports:
            print(f"\n{extreme_type.capitalize()} {percentage}%:")
            print(f"  Mean {display_name}: {gene_values.mean():.4f}")
            print(f"  Range: {gene_values.min():.4f} to {gene_values.max():.4f}")
            
            # Show additional metric averages if available
            if hasattr(self.data, 'enc_values') and self.data.enc_values:
                enc_vals = [self.data.enc_values.get(g, np.nan) for g in gene_list]
                enc_mean = np.nanmean(enc_vals)
                if not np.isnan(enc_mean):
                    print(f"  Mean ENC: {enc_mean:.2f}")
            
            if hasattr(self.data, 'cai_values') and self.data.cai_values:
                cai_vals = [self.data.cai_values.get(g, np.nan) for g in gene_list]
                cai_mean = np.nanmean(cai_vals)
                if not np.isnan(cai_mean):
                    print(f"  Mean CAI: {cai_mean:.4f}")
        
        self._pause()
        return True

    def _optimal_codons_menu(self):
        """Handle optimal codons identification options."""
        if not self._check_data():
            return True

        options = [
            ("Identify optimal codons using all genes (frequency-based)",
             "1", self._identify_and_display_optimal_codons),
            ("Identify optimal codons using manually selected reference genes",
             "2", self._identify_optimal_codons_manual_ref),
            ("Identify optimal codons using multivariate analysis", "3",
             self._identify_and_display_optimal_codons_multivariate),
            ("Identify optimal codons using highest RSCU values",
             "4", self._identify_and_display_optimal_codons_rscu),
            ("Identify optimal codons using most common codons", "5",
             self._identify_and_display_optimal_codons_most_common),
            ("Load reference genes from file", "6",
             self._load_reference_genes_from_file),
            ("Load optimal codons from file", "7",
             self._load_optimal_codons_from_file),
            ("Compare codon usage between axis cohorts",
             "8", self._compare_codon_usage_axis_cohorts),
            ("Export optimal codons to file", "9", self._export_optimal_codons),
            ("Return to analysis menu", "R", lambda: False)
        ]

        return self._display_menu("Optimal Codons Identification", options)

    def _identify_and_display_optimal_codons(self, reference_genes=None):
        """Identify and display optimal codons."""
        print("\nIdentifying optimal codons...")
        optimal_codons = self.data.calculate_optimal_codons(reference_genes)

        print("\nOptimal codons:")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print("Amino Acid\tOptimal Codon\tCodon DNA")

        # Group optimal codons by amino acid
        for aa in sorted(optimal_codons.keys()):
            rna_codon = optimal_codons[aa]
            dna_codon = rna_codon.replace('U', 'T')
            print(f"{aa}\t{rna_codon}\t{dna_codon}")

        # Ask if user wants to export
        export = input(
            "\nDo you want to export optimal codons to a file? (y/n): ").strip().lower()
        if export == 'y':
            self._export_optimal_codons(optimal_codons, reference_genes)

        self._pause()
        return True

    def _identify_optimal_codons_manual_ref(self):
        """Identify optimal codons using manually selected reference genes."""
        reference_genes = self._select_reference_genes()
        if reference_genes:
            self._identify_and_display_optimal_codons(reference_genes)
        return True

    def _identify_and_display_optimal_codons_multivariate(self):
        """Identify and display optimal codons using multivariate analysis."""
        # Ask for percentage of genes to use
        try:
            percentage = float(input(
                "\nEnter percentage of genes to use from multivariate analysis (5-20 recommended): ").strip())
            if percentage < 1 or percentage > 50:
                print("Invalid percentage. Using default of 10%.")
                percentage = 10
        except ValueError:
            print("Invalid input. Using default of 10%.")
            percentage = 10

        print("\nIdentifying optimal codons using multivariate analysis...")
        optimal_codons = self.data.calculate_optimal_codons_multivariate(
            percentage)

        print("\nOptimal codons (determined by multivariate analysis):")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print("Amino Acid\tOptimal Codon\tCodon DNA")

        # Group optimal codons by amino acid
        for aa in sorted(optimal_codons.keys()):
            rna_codon = optimal_codons[aa]
            dna_codon = rna_codon.replace('U', 'T')
            print(f"{aa}\t{rna_codon}\t{dna_codon}")

        # Ask if user wants to export
        export = input(
            "\nDo you want to export optimal codons to a file? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, "optimal_codons_multivariate") + '.tsv'

            with open(output_path, 'w') as f:
                # Add metadata
                f.write("# GCUA Optimal Codons (Multivariate Analysis)\n")
                f.write(f"# Version: {VERSION}\n")
                f.write(
                    f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
                f.write(f"# Multivariate Percentage: {percentage}%\n")
                f.write(f"# Total Genes: {len(self.data.gene_names)}\n")
                f.write(
                    f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"#{'=' * 50}\n\n")

                f.write("AA\tRNA_Codon\tDNA_Codon\n")
                for aa in sorted(optimal_codons.keys()):
                    rna_codon = optimal_codons[aa]
                    dna_codon = rna_codon.replace('U', 'T')
                    f.write(f"{aa}\t{rna_codon}\t{dna_codon}\n")

            print(f"\nOptimal codons saved to {output_path}")

        self._pause()
        return True

    def _identify_and_display_optimal_codons_rscu(self):
        """Identify and display optimal codons using highest RSCU values."""
        print("\nIdentifying optimal codons using highest RSCU values...")
        optimal_codons = self.data.calculate_optimal_codons(None, "rscu")

        print("\nOptimal codons (determined by highest RSCU values):")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print("Amino Acid\tOptimal Codon\tCodon DNA")

        # Group optimal codons by amino acid
        for aa in sorted(optimal_codons.keys()):
            rna_codon = optimal_codons[aa]
            dna_codon = rna_codon.replace('U', 'T')
            print(f"{aa}\t{rna_codon}\t{dna_codon}")

        # Ask if user wants to export
        export = input(
            "\nDo you want to export optimal codons to a file? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, "optimal_codons_rscu") + '.tsv'

            with open(output_path, 'w') as f:
                # Add metadata
                f.write("# GCUA Optimal Codons (Highest RSCU)\n")
                f.write(f"# Version: {VERSION}\n")
                f.write(
                    f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
                f.write(f"# Total Genes: {len(self.data.gene_names)}\n")
                f.write(
                    f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"#{'=' * 50}\n\n")

                f.write("AA\tRNA_Codon\tDNA_Codon\n")
                for aa in sorted(optimal_codons.keys()):
                    rna_codon = optimal_codons[aa]
                    dna_codon = rna_codon.replace('U', 'T')
                    f.write(f"{aa}\t{rna_codon}\t{dna_codon}\n")

            print(f"\nOptimal codons saved to {output_path}")

        self._pause()
        return True

    def _identify_and_display_optimal_codons_most_common(self):
        """Identify and display optimal codons using most common codons."""
        print("\nIdentifying optimal codons using most common codons...")
        optimal_codons = self.data.calculate_optimal_codons(None, "raw_count")

        print("\nOptimal codons (determined by most common usage):")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print("Amino Acid\tOptimal Codon\tCodon DNA")

        # Group optimal codons by amino acid
        for aa in sorted(optimal_codons.keys()):
            rna_codon = optimal_codons[aa]
            dna_codon = rna_codon.replace('U', 'T')
            print(f"{aa}\t{rna_codon}\t{dna_codon}")

        # Ask if user wants to export
        export = input(
            "\nDo you want to export optimal codons to a file? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, "optimal_codons_most_common") + '.tsv'

            with open(output_path, 'w') as f:
                # Add metadata
                f.write("# GCUA Optimal Codons (Most Common)\n")
                f.write(f"# Version: {VERSION}\n")
                f.write(
                    f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
                f.write(f"# Total Genes: {len(self.data.gene_names)}\n")
                f.write(
                    f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"#{'=' * 50}\n\n")

                f.write("AA\tRNA_Codon\tDNA_Codon\n")
                for aa in sorted(optimal_codons.keys()):
                    rna_codon = optimal_codons[aa]
                    dna_codon = rna_codon.replace('U', 'T')
                    f.write(f"{aa}\t{rna_codon}\t{dna_codon}\n")

            print(f"\nOptimal codons saved to {output_path}")

        self._pause()
        return True

    def _load_reference_genes_from_file(self):
        """Load reference genes from file and use them to identify optimal codons."""
        file_path = input(
            "\nEnter the path to the file containing reference gene names (one per line): ").strip()

        reference_genes = self.data.load_reference_genes_from_file(file_path)

        if reference_genes:
            # Use these genes to identify optimal codons
            print(
                f"\nUsing {len(reference_genes)} reference genes from file...")
            self._identify_and_display_optimal_codons(reference_genes)
        else:
            print("\nFailed to load reference genes from file.")
            self._pause()
        return True

    def _load_optimal_codons_from_file(self):
        """Load optimal codons from file with clear format instructions."""
        print("\nLoad Optimal Codons from File")
        print("============================")

        print("\nSupported file formats:")
        print("1. Tab-separated values (TSV) format:")
        print("   - Header row with columns: AA, RNA_Codon, DNA_Codon")
        print("   - One row per amino acid with its optimal codon")
        print("   - Example:")
        print("     AA\tRNA_Codon\tDNA_Codon")
        print("     Phe\tUUU\tTTT")
        print("     Leu\tCUG\tCTG")
        print("     ...\t...\t...")

        print("\n2. JSON format:")
        print("   - Structure with 'optimal_codons' dictionary mapping AAs to codons")
        print("   - Optional metadata section")
        print("   - Example:")
        print('     {"optimal_codons": {')
        print('       "Phe": {"rna_codon": "UUU", "dna_codon": "TTT"},')
        print('       "Leu": {"rna_codon": "CUG", "dna_codon": "CTG"},')
        print('       ...}')
        print('     }')

        print(
            "\nNote: The file can be generated by GCUA's 'Export optimal codons' function")

        file_path = input(
            "\nEnter the path to the file containing optimal codons: ").strip()

        if not file_path:
            print("No file path provided.")
            self._pause()
            return True

        optimal_codons = self.data.load_optimal_codons_from_file(file_path)

        if optimal_codons:
            print("\nLoaded optimal codons:")
            print(
                f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
            print("Amino Acid\tOptimal Codon\tCodon DNA")
            print("-" * 40)

            for aa in sorted(optimal_codons.keys()):
                rna_codon = optimal_codons[aa]
                dna_codon = rna_codon.replace('U', 'T')
                print(f"{aa}\t\t{rna_codon}\t\t{dna_codon}")

            print(
                f"\nSuccessfully loaded optimal codons for {len(optimal_codons)} amino acids.")

            # Ask if the user wants to use these codons for sequence
            # optimization
            optimize = input(
                "\nWould you like to optimize sequences using these codons? (y/n): ").strip().lower()
            if optimize == 'y':
                self._optimize_with_external_optimal_codons()
        else:
            print("\nFailed to load optimal codons from the file.")
            print(
                "Please check that the file exists and is in one of the supported formats.")

        self._pause()
        return True

    def _export_optimal_codons(self,optimal_codons=None,reference_genes=None):
        """Export optimal codons to a file."""
        # Calculate if not provided
        if not optimal_codons:
            if not hasattr(
                    self.data,
                    'optimal_codons') or not self.data.optimal_codons:
                print("\nIdentifying optimal codons...")
                optimal_codons = self.data.calculate_optimal_codons(
                    reference_genes)
            else:
                optimal_codons = self.data.optimal_codons

        ref_info = "_with_ref" if reference_genes else ""
        output_path = self.file_manager.get_output_path(
            self.data, None, f"optimal_codons{ref_info}") + '.tsv'

        with open(output_path, 'w') as f:
            # Add metadata
            f.write("# GCUA Optimal Codons\n")
            f.write(f"# Version: {VERSION}\n")
            f.write(
                f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
            f.write(f"# Total Genes: {len(self.data.gene_names)}\n")
            if reference_genes:
                f.write(f"# Reference Genes: {len(reference_genes)}\n")
            f.write(
                f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"#{'=' * 50}\n\n")

            f.write("AA\tRNA_Codon\tDNA_Codon\n")
            for aa in sorted(optimal_codons.keys()):
                rna_codon = optimal_codons[aa]
                dna_codon = rna_codon.replace('U', 'T')
                f.write(f"{aa}\t{rna_codon}\t{dna_codon}\n")

        print(f"\nOptimal codons saved to {output_path}")

        # If reference genes were used, also save them
        if reference_genes:
            ref_path = output_path.replace('.tsv', '_reference_genes.txt')
            with open(ref_path, 'w') as f:
                for gene in reference_genes:
                    f.write(f"{gene}\n")
            print(f"Reference genes list saved to {ref_path}")

        self._pause()
        return output_path

    def _compare_codon_usage_axis_cohorts(self):
        """Compare codon usage between cohorts from opposite ends of Axis 1."""
        if not self._check_data():
            return True

        # Check if we have enough sequences
        if len(self.data.gene_names) < 10:
            print("\nThis analysis requires at least 10 sequences to be meaningful.")
            print(
                f"Currently only {len(self.data.gene_names)} sequences are loaded.")
            self._pause()
            return True

        # Ask for percentage to use for cohorts
        try:
            default_percentage = 10
            percentage_input = input(
                f"\nEnter percentage of genes to use at each axis extreme (default {default_percentage}%): ").strip()
            percentage = float(
                percentage_input) if percentage_input else default_percentage

            if percentage < 1 or percentage > 50:
                print(
                    f"Invalid percentage. Using default of {default_percentage}%.")
                percentage = default_percentage
        except ValueError:
            print(f"Invalid input. Using default of {default_percentage}%.")
            percentage = default_percentage

        # Ask for significance threshold
        try:
            default_threshold = 0.05
            threshold_input = input(
                f"\nEnter p-value threshold for significance (default {default_threshold}): ").strip()
            threshold = float(
                threshold_input) if threshold_input else default_threshold

            if threshold <= 0 or threshold >= 1:
                print(
                    f"Invalid threshold. Using default of {default_threshold}.")
                threshold = default_threshold
        except ValueError:
            print(f"Invalid input. Using default of {default_threshold}.")
            threshold = default_threshold

        # Perform the comparison
        print(
            f"\nComparing codon usage between {percentage}% of genes at each end of Axis 1...")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        results = self.data.compare_codon_usage_axis_cohorts(
            percentage, threshold)

        if not results:
            print("Comparison failed.")
            self._pause()
            return True

        # Ask if user wants to export results
        export = input(
            "\nDo you want to export these results to files? (y/n): ").strip().lower()
        if export == 'y':
            output_path = self.file_manager.get_output_path(
                self.data, None, "codon_comparison")
            self.data.save_codon_comparison_results(results, output_path)

        self._pause()
        return True

    def _export_codon_comparison_results(self):
        """Export results from codon usage comparison between axis cohorts."""
        print("\nThis will perform a comparison of codon usage between cohorts at opposite ends of Axis 1.")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

        # Ask if user wants to perform the comparison now
        perform = input(
            "Do you want to perform the comparison now? (y/n): ").strip().lower()

        if perform == 'y':
            # Execute the comparison
            results = self.data.compare_codon_usage_axis_cohorts()

            if results:
                # Export the results
                output_path = self.file_manager.get_output_path(
                    self.data, None, "codon_comparison")
                self.data.save_codon_comparison_results(results, output_path)
        else:
            print("\nOperation cancelled.")

        self._pause()
        return True

    def _export_menu(self):
        """Handle data export options."""
        if not self._check_data():
            return True

        options = [
            ("Export comprehensive metrics", "1",
             self._export_comprehensive_metrics),
            ("Export codon usage data", "2", self._export_codon_usage),
            ("Export amino acid usage data", "3", self._export_aa_usage),
            ("Export multivariate analysis results",
             "4", self._export_multivariate),
            ("Export optimized sequences", "5", self._optimize_all_genes),
            ("Export optimal codons", "6", self._export_optimal_codons_to_file_menu),
            ("Export codon usage comparison results",
             "7", self._export_codon_comparison_results),
            ("Return to main menu", "R", lambda: False)
        ]

        return self._display_menu("Export Data Menu", options)

    def _export_optimal_codons_to_file_menu(self):
        """Handle exporting optimal codons to file in different formats."""
        if not self._check_data():
            return True

        # Check if optimal codons have been calculated
        if not hasattr(
                self.data,
                'optimal_codons') or not self.data.optimal_codons:
            print("\nNo optimal codons have been identified yet.")
            print("Please identify optimal codons first.")
            self._pause()
            return True

        # Ask for file format
        print("\nChoose output format:")
        print("1. Tab-separated values (TSV)")
        print("2. JSON (with metadata)")

        choice = input("\nEnter your choice: ").strip()

        file_format = "tsv" if choice != '2' else "json"
        format_suffix = ".tsv" if file_format == "tsv" else ".json"

        # Get output path
        output_path = self.file_manager.get_output_path(
            self.data, None, "optimal_codons") + format_suffix

        # Save to file
        self.data.save_optimal_codons_to_file(output_path, format=file_format)

        self._pause()
        return True
    
    def _export_rscu(self):
        """Export RSCU values to file."""
        if not hasattr(self.data, 'rscu_values') or self.data.rscu_values.empty:
            print("\nRSCU values not calculated yet.")
            self._pause()
            return True
            
        output_path = self.file_manager.get_output_path(
            self.data, None, "rscu_values") + '.tsv'
            
        print(f"\nExporting RSCU values to {output_path}...")
        
        # Write with metadata
        with open(output_path, 'w') as f:
            # Add metadata
            f.write(f"# GCUA RSCU Values\n")
            f.write(f"# Version: {VERSION}\n")
            f.write(f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
            f.write(f"# Total Genes: {len(self.data.rscu_values)}\n")
            f.write(f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"#{'=' * 50}\n\n")
            
        # Append the data
        self.data.rscu_values.to_csv(output_path, sep='\t', mode='a')
        
        print(f"RSCU values exported to: {output_path}")
        self._pause()
        return True

    def _export_comprehensive_metrics(self):
        """Export comprehensive metrics to a file."""
        print("\nGenerating comprehensive metrics file...")

        # Get all metrics
        metrics = self.data.get_comprehensive_metrics()

        # Export to file
        output_path = self.file_manager.get_output_path(
            self.data, None, "comprehensive_metrics") + '.tsv'

        # Add metadata about the genetic code used
        metadata = {
            "Genetic Code ID": self.config.genetic_code,
            "Genetic Code Name": self.config.get_genetic_code_name(),
            "Total Sequences": len(self.data.gene_names),
            "Generated": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
            "GCUA Version": VERSION
        }

        with open(output_path, 'w') as f:
            # Write metadata as comments
            for key, value in metadata.items():
                f.write(f"# {key}: {value}\n")
            f.write("#" + "-" * 50 + "\n\n")

            # Write data
            metrics.to_csv(f, sep='\t', float_format='%.4f')

        print(f"\nComprehensive metrics saved to {output_path}")

        # Display summary
        print(f"\nSummary of exported metrics:")
        print(f"Number of genes: {len(metrics)}")
        print(f"Metrics included: {', '.join(metrics.columns)}")
        print(
            f"Using genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

        self._pause()
        return True

    def _detect_clusters_menu(self):
        """Handle cluster detection in multivariate results."""
        if not self._check_data():
            return True
        
        # Check if multivariate analysis has been performed
        if not hasattr(self.data, 'multivariate_results') or not self.data.multivariate_results:
            print("\nNo multivariate analysis results found.")
            print("Please perform a multivariate analysis first.")
            self._pause()
            return True
        
        # Show current multivariate analysis
        analysis_type = self.data.multivariate_results.get('analysis_type', 'Unknown')
        data_type = self.data.multivariate_results.get('data_type', 'Unknown')
        print(f"\nCurrent multivariate analysis: {analysis_type} on {data_type} data")
        
        # Show clustering options
        print("\nCluster Detection Methods:")
        print("1. K-means clustering (spherical clusters)")
        print("2. DBSCAN (density-based, finds arbitrary shapes)")
        print("3. Auto-select method based on dataset size")
        print("R. Return to analysis menu")
        
        choice = input("\nSelect clustering method (1-3, R): ").strip().upper()
        
        if choice == 'R':
            return True
        
        if choice == '1':
            method = 'kmeans'
            # Ask for number of clusters
            k_choice = input("\nNumber of clusters (press Enter for auto-detection): ").strip()
            if k_choice:
                try:
                    n_clusters = int(k_choice)
                except ValueError:
                    print("Invalid number. Using auto-detection.")
                    n_clusters = None
            else:
                n_clusters = None
                print("\nAuto-detecting optimal number of clusters...")
            
            # Perform clustering
            results = self.data.detect_clusters(method=method, n_clusters=n_clusters)
            
        elif choice == '2':
            method = 'dbscan'
            # Ask for parameters
            print("\nDBSCAN Parameters (press Enter for defaults):")
            min_samples = input("Minimum samples per cluster (default=5): ").strip()
            min_samples = int(min_samples) if min_samples else 5
            
            eps = input("Epsilon (neighborhood size, default=auto): ").strip()
            eps = float(eps) if eps else None
            
            # Perform clustering
            results = self.data.detect_clusters(method=method, min_samples=min_samples, eps=eps)
            
        elif choice == '3':
            method = 'auto'
            print("\nAuto-selecting clustering method based on dataset size...")
            results = self.data.detect_clusters(method=method)
            
        else:
            print("Invalid choice.")
            self._pause()
            return True
        
        if results:
            # Show results
            print(f"\nClustering completed using {results['method']}")
            print(f"Number of clusters found: {results['n_clusters']}")
            
            # Show cluster sizes
            print("\nCluster sizes:")
            for cluster_name, stats in results['statistics'].items():
                print(f"  {cluster_name}: {stats['size']} genes")
            
            # Ask if user wants to characterize clusters
            char_choice = input("\nCharacterize clusters by their patterns? (y/n): ").strip().lower()
            if char_choice == 'y':
                characterization = self.data.characterize_clusters()
                
                if characterization:
                    # Display characterization
                    for cluster_name, char_data in characterization.items():
                        print(f"\n{cluster_name} Characteristics:")
                        print(f"  Size: {char_data['size']} genes")
                        if char_data['mean_GC'] is not None:
                            print(f"  Mean GC: {char_data['mean_GC']:.1f}%")
                        if char_data['mean_GC3s'] is not None:
                            print(f"  Mean GC3s: {char_data['mean_GC3s']:.1f}%")
                        if char_data['mean_ENC'] is not None:
                            print(f"  Mean ENC: {char_data['mean_ENC']:.1f}")
                        
                        print(f"\n  Top distinctive {data_type} features:")
                        for i, feature in enumerate(char_data['distinctive_features'][:5]):
                            if data_type == 'RSCU':
                                print(f"    {feature['codon']}: {feature['cluster_rscu']:.3f} vs {feature['other_rscu']:.3f} (enrichment: {feature['enrichment']:.2f}x)")
                            else:
                                print(f"    {feature['amino_acid']}: {feature['cluster_freq']:.2f}% vs {feature['other_freq']:.2f}% (enrichment: {feature['enrichment']:.2f}x)")
            
            # Ask if user wants to visualize clusters
            vis_choice = input("\nVisualize clusters? (y/n): ").strip().lower()
            if vis_choice == 'y':
                self._visualize_clusters()
            
            # Ask if user wants to export cluster results
            export_choice = input("\nExport cluster results? (y/n): ").strip().lower()
            if export_choice == 'y':
                self._export_cluster_results()
        
        self._pause()
        return True
    
    def _visualize_clusters(self):
        """Create a visualization of the detected clusters."""
        if not hasattr(self.data, 'cluster_results') or not self.data.cluster_results:
            print("\nNo cluster results to visualize.")
            return
        
        results = self.data.cluster_results
        coords = self.data.multivariate_results['coordinates']
        labels = results['labels']
        
        # Check how many dimensions were used for clustering
        n_dims = len(results['coordinates'])
        
        # Ask for 2D vs 3D if applicable
        use_3d = False
        if n_dims >= 3 and coords.shape[1] >= 3:
            viz_choice = input("\nVisualization type:\n1. 2D plot (Axis 1 vs Axis 2)\n2. 3D plot (Axes 1, 2, 3)\n\nSelect (1-2): ").strip()
            use_3d = (viz_choice == '2')
        
        # Create figure
        fig = go.Figure()
        
        # Plot each cluster with a different color
        unique_labels = sorted(set(labels))
        colors = px.colors.qualitative.Plotly  # Use Plotly's color palette
        
        for i, label in enumerate(unique_labels):
            if label == -1:  # Noise points from DBSCAN
                color = 'gray'
                name = 'Noise'
            else:
                color = colors[i % len(colors)]
                name = f'Cluster {label}'
            
            mask = labels == label
            
            if use_3d:
                fig.add_trace(go.Scatter3d(
                    x=coords.iloc[mask, 0],
                    y=coords.iloc[mask, 1],
                    z=coords.iloc[mask, 2],
                    mode='markers',
                    marker=dict(
                        color=color,
                        size=6 if label != -1 else 4,
                        opacity=0.7 if label != -1 else 0.3
                    ),
                    text=coords.index[mask],
                    hovertemplate='<b>%{text}</b><br>%{xaxis.title.text}: %{x:.3f}<br>%{yaxis.title.text}: %{y:.3f}<br>%{zaxis.title.text}: %{z:.3f}<extra></extra>',
                    name=name
                ))
            else:
                fig.add_trace(go.Scatter(
                    x=coords.iloc[mask, 0],
                    y=coords.iloc[mask, 1],
                    mode='markers',
                    marker=dict(
                        color=color,
                        size=8 if label != -1 else 5,
                        opacity=0.7 if label != -1 else 0.3
                    ),
                    text=coords.index[mask],
                    hovertemplate='<b>%{text}</b><br>%{xaxis.title.text}: %{x:.3f}<br>%{yaxis.title.text}: %{y:.3f}<extra></extra>',
                    name=name
                ))
        
        # Update layout
        explained_var = self.data.multivariate_results['explained_variance']
        
        if use_3d:
            fig.update_layout(
                title=f"Cluster Analysis - {results['method']} (3D)",
                scene=dict(
                    xaxis=dict(
                        title=f"{coords.columns[0]} ({explained_var[0]:.2%})",
                        showgrid=True,
                        zeroline=True
                    ),
                    yaxis=dict(
                        title=f"{coords.columns[1]} ({explained_var[1]:.2%})",
                        showgrid=True,
                        zeroline=True
                    ),
                    zaxis=dict(
                        title=f"{coords.columns[2]} ({explained_var[2]:.2%})",
                        showgrid=True,
                        zeroline=True
                    ),
                    bgcolor='white'
                ),
                width=1000,
                height=800,
                margin=dict(l=80, r=80, t=100, b=80)
            )
        else:
            fig.update_layout(
                title=f"Cluster Analysis - {results['method']}",
                xaxis=dict(
                    title=f"{coords.columns[0]} ({explained_var[0]:.2%})",
                    showgrid=True,
                    zeroline=True
                ),
                yaxis=dict(
                    title=f"{coords.columns[1]} ({explained_var[1]:.2%})",
                    showgrid=True,
                    zeroline=True
                ),
                plot_bgcolor='white',
                width=1000,
                height=800,
                margin=dict(l=80, r=80, t=100, b=80)
            )
        
        # Save and display
        output_path = self.file_manager.get_output_path(
            self.data, None, "cluster_plot", f".{self.config.visualization_format}")
        
        # Export the figure
        if self.config.visualization_format == 'html':
            fig.write_html(output_path)
            # Try to open in browser
            try:
                webbrowser.open(f"file:///{os.path.abspath(output_path)}")
                print(f"\nCluster visualization opened in browser.")
            except:
                print(f"\nCluster visualization saved to {output_path}")
        else:
            # For other formats, use plotly's export functions
            if self.config.visualization_format == 'svg':
                fig.write_image(output_path, format='svg')
            elif self.config.visualization_format == 'png':
                fig.write_image(output_path, format='png')
            elif self.config.visualization_format == 'pdf':
                fig.write_image(output_path, format='pdf')
            print(f"\nCluster visualization saved to {output_path}")
    
    def _export_cluster_results(self):
        """Export cluster detection results."""
        if not hasattr(self.data, 'cluster_results') or not self.data.cluster_results:
            print("\nNo cluster results to export.")
            return
        
        results = self.data.cluster_results
        labels = results['labels']
        
        # Create main results file
        output_path = self.file_manager.get_output_path(
            self.data, None, "cluster_results") + '.tsv'
        
        # Create DataFrame with gene names and cluster assignments
        cluster_df = pd.DataFrame({
            'Gene': self.data.gene_names,
            'Cluster': labels
        })
        
        # Add multivariate coordinates
        coords = self.data.multivariate_results['coordinates']
        for col in coords.columns[:3]:  # Include up to 3 dimensions
            cluster_df[col] = coords[col].values
        
        # Write main results
        with open(output_path, 'w') as f:
            # Metadata
            f.write(f"# GCUA Cluster Analysis Results\n")
            f.write(f"# Version: {VERSION}\n")
            f.write(f"# Method: {results['method']}\n")
            f.write(f"# Analysis: {results['analysis_type']} on {results['data_type']} data\n")
            f.write(f"# Number of clusters: {results['n_clusters']}\n")
            f.write(f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("#" + "-" * 50 + "\n\n")
            
            cluster_df.to_csv(f, sep='\t', index=False)
        
        print(f"\nCluster assignments saved to {output_path}")
        
        # Export cluster statistics
        stats_path = output_path.replace('.tsv', '_statistics.tsv')
        
        with open(stats_path, 'w') as f:
            f.write("# Cluster Statistics\n")
            f.write("#" + "-" * 50 + "\n\n")
            
            # Create statistics DataFrame
            stats_data = []
            for cluster_name, stats in results['statistics'].items():
                stats_row = {
                    'Cluster': cluster_name,
                    'Size': stats['size'],
                    'Mean_Distance': stats['mean_distance'],
                    'Std_Distance': stats['std_distance'],
                    'Compactness': stats['compactness']
                }
                # Add centroid coordinates
                for i, coord in enumerate(stats['centroid']):
                    stats_row[f'Centroid_{results["coordinates"][i]}'] = coord
                stats_data.append(stats_row)
            
            stats_df = pd.DataFrame(stats_data)
            stats_df.to_csv(f, sep='\t', index=False, float_format='%.4f')
        
        print(f"Cluster statistics saved to {stats_path}")
        
        # Export characterization if available
        if hasattr(self.data, 'cluster_characterization'):
            char_path = output_path.replace('.tsv', '_characterization.tsv')
            
            with open(char_path, 'w') as f:
                f.write("# Cluster Characterization\n")
                f.write("#" + "-" * 50 + "\n\n")
                
                for cluster_name, char_data in self.data.cluster_characterization.items():
                    f.write(f"\n## {cluster_name}\n")
                    f.write(f"Size: {char_data['size']}\n")
                    f.write(f"Mean GC: {char_data['mean_GC']:.2f}%\n")
                    f.write(f"Mean GC3s: {char_data['mean_GC3s']:.2f}%\n")
                    if char_data['mean_ENC']:
                        f.write(f"Mean ENC: {char_data['mean_ENC']:.2f}\n")
                    
                    f.write("\nDistinctive features:\n")
                    if self.data.cluster_results['data_type'] == 'RSCU':
                        f.write("Codon\tCluster_RSCU\tOther_RSCU\tDifference\tEnrichment\n")
                        for feature in char_data['distinctive_features']:
                            f.write(f"{feature['codon']}\t{feature['cluster_rscu']:.3f}\t"
                                   f"{feature['other_rscu']:.3f}\t{feature['difference']:.3f}\t"
                                   f"{feature['enrichment']:.2f}\n")
                    else:
                        f.write("Amino_Acid\tCluster_Freq%\tOther_Freq%\tDifference%\tEnrichment\n")
                        for feature in char_data['distinctive_features']:
                            f.write(f"{feature['amino_acid']}\t{feature['cluster_freq']:.2f}\t"
                                   f"{feature['other_freq']:.2f}\t{feature['difference']:.2f}\t"
                                   f"{feature['enrichment']:.2f}\n")
            
            print(f"Cluster characterization saved to {char_path}")
        
        # Export gene lists for each cluster
        for cluster_name, stats in results['statistics'].items():
            if 'genes' in stats and stats['genes']:
                gene_list_path = output_path.replace('.tsv', f'_{cluster_name}_genes.txt')
                with open(gene_list_path, 'w') as f:
                    f.write(f"# Gene list for {cluster_name}\n")
                    f.write(f"# Size: {stats['size']} genes\n")
                    f.write("#" + "-" * 30 + "\n\n")
                    for gene in stats['genes']:
                        f.write(f"{gene}\n")
                print(f"Gene list for {cluster_name} saved to {gene_list_path}")

        self._pause()
        return True

    def _export_fop_cai(self):
        """Export Fop and CAI values to a file."""
        if not self._check_data():
            return True

        # Check if Fop and CAI values have been calculated
        if not hasattr(self.data, 'fop_values') or not self.data.fop_values:
            print("\nCalculating Fop values...")
            self.data.calculate_fop()

        if not hasattr(self.data, 'cai_values') or not self.data.cai_values:
            print("\nCalculating CAI values...")
            self.data.calculate_cai()

        output_path = self.file_manager.get_output_path(
            self.data, None, "fop_cai") + '.tsv'

        with open(output_path, 'w') as f:
            # Add metadata
            f.write("# GCUA Fop and CAI Values\n")
            f.write(f"# Version: {VERSION}\n")
            f.write(
                f"# Genetic Code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}\n")
            f.write(f"# Total Genes: {len(self.data.gene_names)}\n")
            f.write(
                f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"#{'=' * 50}\n\n")

            f.write("Gene\tFop\tCAI\n")
            for gene in sorted(self.data.gene_names):
                f.write(
                    f"{gene}\t{self.data.fop_values.get(gene, 0):.4f}\t{self.data.cai_values.get(gene, 0):.4f}\n")

        print(f"\nFop and CAI values saved to {output_path}")
        self._pause()
        return True

    def _optimization_menu(self):
        """Handle sequence optimization options."""
        if not self._check_data():
            return True

        options = [
            ("Optimize a single gene using all genes as reference", "1", self._optimize_single_gene),
            ("Optimize a single gene using multivariate analysis", "2", self._optimize_single_gene_multivariate),
            ("Optimize a single gene using manually selected reference genes", "3", self._optimize_single_gene_manual_ref),
            ("Optimize all genes", "4", self._optimize_all_genes),
            ("Load optimal codons from file for optimization", "5", self._optimize_with_external_optimal_codons),
            ("Return to main menu", "R", lambda: False)
        ]

        return self._display_menu("Sequence Optimization Menu", options)

    def _optimize_single_gene(self, reference_genes=None, method="frequency"):
        """Optimize a single gene sequence."""
        if not self._check_data():
            return True

        # Display list of genes
        print("\nAvailable genes:")
        for i, gene in enumerate(self.data.gene_names, 1):
            print(f"{i}. {gene}")
            if i % 20 == 0 and i < len(self.data.gene_names):
                more = input("\nPress Enter to see more genes, or enter a number to select: ").strip()
                if more.isdigit() and 1 <= int(more) <= len(self.data.gene_names):
                    gene_idx = int(more) - 1
                    gene_id = self.data.gene_names[gene_idx]
                    return self._perform_gene_optimization(gene_id, reference_genes, method)

        # Get gene selection
        try:
            selection = input("\nEnter gene number to optimize: ").strip()
            if selection.isdigit() and 1 <= int(selection) <= len(self.data.gene_names):
                gene_idx = int(selection) - 1
                gene_id = self.data.gene_names[gene_idx]
                return self._perform_gene_optimization(gene_id, reference_genes, method)
            else:
                print("Invalid selection.")
        except ValueError:
            print("Invalid input.")

        self._pause()
        return True

    def _perform_gene_optimization(self, gene_id, reference_genes=None, method="frequency"):
        """Perform optimization on the selected gene."""
        print(f"\nOptimizing gene: {gene_id}")
        print(f"Using genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print(f"Method: {method}")

        # Get original and optimized sequences
        result = self.data.optimize_gene_sequence(gene_id, method, reference_genes, output_format='both')

        if not result:
            print(f"Failed to optimize gene {gene_id}.")
            self._pause()
            return True

        original_seq, optimized_seq = result

        # Display comparison
        print("\nOriginal sequence (first 60 bases):")
        print(str(original_seq.seq)[:60] + "...")
        print("\nOptimized sequence (first 60 bases):")
        print(optimized_seq[:60] + "...")

        print(f"\nSequence length: {len(original_seq.seq)} bp")

        # Count differences
        diff_count = sum(1 for a, b in zip(str(original_seq.seq), optimized_seq) if a != b)
        diff_percentage = (diff_count / len(original_seq.seq)) * 100
        print(f"Changes made: {diff_count} ({diff_percentage:.2f}%)")

        # Ask if user wants to save the optimized sequence
        save = input("\nDo you want to save the optimized sequence? (y/n): ").strip().lower()
        if save == 'y':
            output_path = self.file_manager.get_output_path(self.data, None, f"{gene_id}_optimized") + '.fasta'
            with open(output_path, 'w') as f:
                f.write(f">{gene_id}_optimized\n{optimized_seq}")
            print(f"\nOptimized sequence saved to {output_path}")

        self._pause()
        return True

    def _optimize_single_gene_multivariate(self):
        """Optimize a single gene using multivariate analysis for reference genes."""
        try:
            percentage = float(input("\nEnter percentage of genes to use from multivariate analysis (5-20 recommended): ").strip())
            if percentage < 1 or percentage > 50:
                print("Invalid percentage. Using default of 10%.")
                percentage = 10
        except ValueError:
            print("Invalid input. Using default of 10%.")
            percentage = 10

        return self._optimize_single_gene(None, "multivariate")

    def _optimize_single_gene_manual_ref(self):
        """Optimize a single gene using manually selected reference genes."""
        reference_genes = self._select_reference_genes()
        if reference_genes:
            return self._optimize_single_gene(reference_genes)
        return True

    def _optimize_all_genes(self):
        """Optimize all genes and save to a file."""
        print("\nOptimizing all genes...")
        print(f"Using genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

        # Check if optimal codons have been calculated
        if not hasattr(self.data, 'optimal_codons') or not self.data.optimal_codons:
            print("Identifying optimal codons...")
            self.data.calculate_optimal_codons()

        # Get output path
        output_path = self.file_manager.get_output_path(self.data, None, "all_optimized") + '.fasta'

        # Optimize and save
        result = self.data.optimize_all_genes(output_path)

        if result:
            print(f"\nAll genes optimized and saved to {output_path}")
        else:
            print("\nFailed to optimize genes.")

        self._pause()
        return True

    def _optimize_with_external_optimal_codons(self):
        """Optimize sequences using optimal codons loaded from an external file."""
        file_path = input("\nEnter the path to the file containing optimal codons: ").strip()

        if not file_path:
            print("No file path provided.")
            self._pause()
            return True

        # Load the optimal codons
        result = self.data.load_optimal_codons_from_file(file_path)

        if not result:
            print("Failed to load optimal codons.")
            self._pause()
            return True

        # Ask if user wants to optimize a single gene or all genes
        choice = input("\nOptimize (1) a single gene or (2) all genes? Enter 1 or 2: ").strip()

        if choice == '1':
            return self._optimize_single_gene()
        elif choice == '2':
            return self._optimize_all_genes()
        else:
            print("Invalid choice.")
            self._pause()
            return True

    def _check_data(self):
        """Check if data is loaded and return True if it is, False otherwise."""
        if not self.data.gene_names:
            print("\nNo data loaded. Please load a FASTA file first.")
            self._pause()
            return False
        return True

    def _pause(self):
        """Pause execution until user presses Enter."""
        input("\nPress Enter to continue...")
        return True

    def _clear_screen(self):
        """Clear the terminal screen."""
        # This is a simplified version that works on most terminals
        print("\033[H\033[J", end="")
        return True
    
    def _show_quick_help(self, menu_key, option_key):
        """Display quick help for a specific menu option."""
        if not menu_key or menu_key not in self.help_content:
            print(f"\nNo help available for this menu.")
            return
            
        menu_help = self.help_content[menu_key]
        if option_key not in menu_help:
            print(f"\nNo help available for option '{option_key}'.")
            return
            
        help_info = menu_help[option_key]
        print(f"\nHelp for option {option_key}:")
        print("-" * 40)
        print(f"\n{help_info['brief']}")
        print(f"\n{help_info['detailed']}")
        
    def _show_menu_help(self, menu_key, options):
        """Display comprehensive help for all menu items."""
        print("\nMenu Help")
        print("=" * 60)
        
        if menu_key and menu_key in self.help_content:
            menu_help = self.help_content[menu_key]
            for option_text, option_key, _ in options:
                if option_key in menu_help:
                    help_info = menu_help[option_key]
                    print(f"\n{option_key}. {option_text}")
                    print("-" * 40)
                    print(help_info['brief'])
                    print()
        else:
            print("No detailed help available for this menu.")
            print("\nGeneral navigation:")
            print("- Enter the number/letter to select an option")
            print("- Type 'R' to return to the previous menu")
            print("- Type a number followed by '?' for help on that option")
            print("- Type 'H' to show this help screen")
    
    def _show_analysis_menu_help(self):
        """Display comprehensive help for the analysis menu."""
        print("\nCustom Analysis Menu Help")
        print("=" * 60)
        
        if 'analysis_menu' in self.help_content:
            menu_help = self.help_content['analysis_menu']
            
            print("\n--- Individual Analyses ---")
            for key in ['1', '2', '3', '4']:
                if key in menu_help:
                    help_info = menu_help[key]
                    print(f"\n{key}. {help_info['brief']}")
                    print("-" * 40)
                    print(help_info['detailed'].strip())
            
            print("\n\n--- Advanced Analysis ---")
            for key in ['5', '6', '7', '8']:
                if key in menu_help:
                    help_info = menu_help[key]
                    print(f"\n{key}. {help_info['brief']}")
                    print("-" * 40)
                    print(help_info['detailed'].strip())
            
            print("\n\n--- Comprehensive Options ---")
            if '9' in menu_help:
                help_info = menu_help['9']
                print(f"\n9. {help_info['brief']}")
                print("-" * 40)
                print(help_info['detailed'].strip())
        
        print("\n\nNavigation:")
        print("- Enter the number/letter to select an option")
        print("- Type 'R' to return to the main menu")
        print("- Type a number followed by '?' for quick help on that option")

    def _select_reference_genes(self):
        """Allow user to select reference genes manually."""
        if not self._check_data():
            return None

        print("\nSelect reference genes:")
        print("1. Select by numbers")
        print("2. Select by names")
        print("3. Load from file")
        print("R. Return")

        choice = input("\nEnter your choice: ").strip().upper()

        if choice == '1':
            return self._select_genes_by_numbers()
        elif choice == '2':
            return self._select_genes_by_names()
        elif choice == '3':
            return self._load_reference_genes()
        else:
            return None

    def _select_genes_by_numbers(self):
        """Select genes by their numbers."""
        print("\nAvailable genes:")
        for i, gene in enumerate(self.data.gene_names, 1):
            print(f"{i}. {gene}")
            if i % 20 == 0 and i < len(self.data.gene_names):
                more = input("\nPress Enter to see more genes, or 'q' to stop listing: ").strip().lower()
                if more == 'q':
                    break

        selected = input("\nEnter gene numbers separated by commas: ").strip()
        try:
            indices = [int(x.strip()) - 1 for x in selected.split(',') if x.strip()]
            selected_genes = [self.data.gene_names[i] for i in indices if 0 <= i < len(self.data.gene_names)]
            return selected_genes if selected_genes else None
        except ValueError:
            print("Invalid input.")
            return None

    def _select_genes_by_names(self):
        """Select genes by their names."""
        print("\nEnter gene names one by one. Enter empty line when done.")
        selected_genes = []
        while True:
            gene_name = input("Gene name (or empty to finish): ").strip()
            if not gene_name:
                break
            if gene_name in self.data.gene_names:
                selected_genes.append(gene_name)
            else:
                print(f"Gene '{gene_name}' not found.")

        return selected_genes if selected_genes else None

    def _load_reference_genes(self):
        """Load reference genes from a file."""
        file_path = input("\nEnter the path to the file containing reference gene names: ").strip()
        return self.data.load_reference_genes_from_file(file_path)

    def _help(self):
        """Display help information."""
        self._clear_screen()
        self.display_banner()

        help_text = """
    GCUA (General Codon Usage Analysis) Help
    ========================================

    GCUA is a comprehensive tool for analyzing codon usage patterns in DNA sequences.

    NEW MENU STRUCTURE:
    ==================
    1. Data Management - Load, view, filter, and manage sequence data
    2. Quick Analysis - Guided workflow for standard analyses
    3. Custom Analysis - Detailed analysis options organized by type
    4. Visualization & Export - Create plots and export results
    5. Advanced Tools - Sequence optimization, clustering, outlier detection
    6. Settings & Preferences - Configure genetic codes and analysis options

    QUICK START WORKFLOW:
    ====================
    1. Select "Data Management" → Load FASTA file
    2. Select "Quick Analysis" for guided analysis workflow
       - Automatically calculates all basic metrics
       - Performs multivariate analysis
       - Identifies optimal codons
       - Creates standard visualizations
       - Exports results

    CUSTOM ANALYSIS TYPES:
    =====================
    Basic Metrics:
    - Codon Usage: Frequencies and RSCU values
    - Amino Acid Usage: Composition analysis
    - Base Composition: GC content and positional analysis

    Advanced Analysis:
    - Multivariate: CA/PCA with preprocessing options
    - Codon Bias: ENC, CAI, Fop, SCUO metrics
    - Comparative: Optimal codons, axis cohort comparison

    VISUALIZATION OPTIONS:
    =====================
    - Multivariate plots with biplots
    - Scree plots for variance explained
    - GC content analysis plots
    - Wright's plot (ENC vs GC3s)
    - RSCU heatmaps (adaptive for large datasets)
    - CAI distribution plots
    - Custom scatter plots of any metrics

    ADVANCED FEATURES:
    ==================
    - Sequence Optimization: Multiple strategies available
    - Cluster Detection: K-means and DBSCAN methods
    - Outlier Analysis: Based on multivariate results
    - Export Extreme Genes: Top/bottom percentiles
    - Cohort Comparison: Statistical analysis of axis extremes

    PERFORMANCE FEATURES:
    ====================
    - Parallel processing for large datasets
    - Automatic checkpointing for recovery
    - Smart memory management
    - Handles millions of sequences efficiently

    For detailed documentation, see the manual.md file
    """

        print(help_text)

        # Add information about the current genetic code
        print(f"\nCurrent genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

        # Add citation information
        print("\nCitation:")
        print("McInerney JO. GCUA: general codon usage analysis.")
        print("Bioinformatics. 1998;14(4):372-3.")
        print("doi: 10.1093/bioinformatics/14.4.372. PMID: 9632833.")

        self._pause()
        return True

    def _calculate_all_metrics(self):
        """Calculate all available metrics for the loaded sequences."""
        if not self._check_data():
            return True

        print("\nCalculating all metrics for loaded sequences...")
        print(f"Using genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        print("This may take some time for large datasets.")

        # Calculate base composition (already done during loading)
        print("\n1. Base composition: Already calculated")

        # Calculate codon usage and RSCU
        print("2. Analyzing codon usage...")
        # This was already calculated during sequence loading

        # Calculate amino acid usage
        print("3. Analyzing amino acid usage...")
        # This was already calculated during sequence loading

        # Perform multivariate analysis
        print("4. Performing multivariate analysis...")
        self.data.perform_multivariate_analysis()

        # Calculate ENC values
        print("5. Calculating Effective Number of Codons (ENC)...")
        enc_values = self.data.calculate_enc()

        # Calculate optimal codons
        print("6. Identifying optimal codons...")
        self.data.calculate_optimal_codons()

        # Calculate Fop values
        print("7. Calculating Frequency of Optimal Codons (Fop)...")
        fop_values = self.data.calculate_fop()

        # Calculate CAI values
        print("8. Calculating Codon Adaptation Index (CAI)...")
        cai_values = self.data.calculate_cai()

        # Calculate SCUO values
        print("9. Calculating Synonymous Codon Usage Order (SCUO)...")
        scuo_values = self.data.calculate_scuo()

        print("\nAll metrics calculated successfully.")

        # Summary statistics
        print("\nSummary Statistics:")
        print("-" * 30)

        # ENC stats
        enc_list = list(enc_values.values())
        if enc_list:
            print(f"ENC - Average: {sum(enc_list)/len(enc_list):.2f}, Min: {min(enc_list):.2f}, Max: {max(enc_list):.2f}")

        # CAI stats
        cai_list = list(cai_values.values())
        if cai_list:
            print(f"CAI - Average: {sum(cai_list)/len(cai_list):.4f}, Min: {min(cai_list):.4f}, Max: {max(cai_list):.4f}")

        # Fop stats
        fop_list = list(fop_values.values())
        if fop_list:
            print(f"Fop - Average: {sum(fop_list)/len(fop_list):.4f}, Min: {min(fop_list):.4f}, Max: {max(fop_list):.4f}")

        # SCUO stats
        scuo_list = list(scuo_values.values())
        if scuo_list:
            print(f"SCUO - Average: {sum(scuo_list)/len(scuo_list):.4f}, Min: {min(scuo_list):.4f}, Max: {max(scuo_list):.4f}")

        # GC content stats
        gc_list = self.data.base_composition['GC'].values
        if len(gc_list) > 0:
            print(f"GC% - Average: {gc_list.mean():.2f}, Min: {gc_list.min():.2f}, Max: {gc_list.max():.2f}")

        # GC3 content stats
        gc3_list = self.data.base_composition['GC3'].values
        if len(gc3_list) > 0:
            print(f"GC3% - Average: {gc3_list.mean():.2f}, Min: {gc3_list.min():.2f}, Max: {gc3_list.max():.2f}")

        # Ask if user wants to export comprehensive metrics
        export = input("\nDo you want to export all metrics to a file? (y/n): ").strip().lower()
        if export == 'y':
            self._export_comprehensive_metrics()

        self._pause()
        return True


# Main execution block
if __name__ == "__main__":
    print("Starting GCUA...")

    # Create configuration
    config = Config()
    interface = GCUAInterface(config)

    # Display banner and start the program
    print("Launching interface...")
    interface.main_menu()

    print("GCUA exited.")
