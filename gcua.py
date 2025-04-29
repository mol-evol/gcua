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
from Bio.SeqUtils import GC
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import webbrowser
from scipy.spatial.distance import cdist
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Optional, Union, Any, Callable
from scipy.stats import chi2_contingency
import json

# Program constants
VERSION = "2.3.0"
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
    'AUA': 'Met', 'UGA': 'Trp',
    'AUN': 'Met',  # AUN -> Met
    'AGA': 'STOP', 'AGG': 'STOP'
})

# Yeast Mitochondrial code (NCBI translation table 3)
AA_MAP_3 = AA_MAP_1.copy()
AA_MAP_3.update({
    'AUA': 'Met', 'UGA': 'Trp',
    'CUN': 'Thr',  # CUN -> Thr
    'CUU': 'Thr', 'CUC': 'Thr', 'CUA': 'Thr', 'CUG': 'Thr',
    'AUA': 'Met'
})

# Mold, Protozoan, Coelenterate Mitochondrial & Mycoplasma/Spiroplasma
# (NCBI translation table 4)
AA_MAP_4 = AA_MAP_1.copy()
AA_MAP_4.update({'UGA': 'Trp'})

# Invertebrate Mitochondrial code (NCBI translation table 5)
AA_MAP_5 = AA_MAP_1.copy()
AA_MAP_5.update({
    'AUA': 'Met', 'UGA': 'Trp',
    'AUN': 'Met',  # AUN -> Met
    'AAA': 'Lys', 'AGA': 'Ser', 'AGG': 'Ser'
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
AA_MAP_12.update({'CUG': 'Ala'})

# Ascidian Mitochondrial (NCBI translation table 13)
AA_MAP_13 = AA_MAP_1.copy()
AA_MAP_13.update({
    'AUA': 'Met', 'UGA': 'Trp',
    'AAA': 'Asn', 'AGA': 'Gly', 'AGG': 'Gly'
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
    'UAA': 'Gln', 'UAG': 'Gln'
})

# Condylostoma Nuclear (NCBI translation table 28)
AA_MAP_28 = AA_MAP_1.copy()
AA_MAP_28.update({
    'UAA': 'Gln', 'UAG': 'Gln'
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
AA_MAP_33.update({'UAA': 'Tyr', 'UGA': 'Trp'})

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


class Config:
    def __init__(self):
        self.genetic_code = 1  # Default: Standard (Universal)
        self.detailed_output = 1  # 1=basic, 2=some details, 3=all details
        self.beep_on = False
        self.visualization_method = "plotly"  # "plotly" or "text"

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
            html_path = os.path.splitext(output_path)[0] + '.html'
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

        # Calculate outliers based on distance from center
        centroid = coords.iloc[:, 0:2].mean().values.reshape(1, -1)
        distances = cdist(coords.iloc[:, 0:2].values,
                          centroid, 'euclidean').flatten()
        threshold = np.mean(distances) + 2 * np.std(distances)
        is_outlier = distances > threshold

        # Create the figure
        fig = go.Figure()

        # Add main cluster points
        fig.add_trace(go.Scatter(
            x=coords.iloc[~is_outlier, 0],
            y=coords.iloc[~is_outlier, 1],
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
            x=coords.iloc[is_outlier, 0],
            y=coords.iloc[is_outlier, 1],
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

                # Add to the button list
                button = dict(
                    method='update',
                    label=f"{dim1_name} vs {dim2_name}",
                    args=[
                        {'x': [main_cluster_trace.x, outlier_trace.x],
                         'y': [main_cluster_trace.y, outlier_trace.y]},
                        {'xaxis.title': f"{dim1_name} ({explained_var[dim1]:.2%})",
                         'yaxis.title': f"{dim2_name} ({explained_var[dim2]:.2%})"}
                    ]
                )
                buttons.append(button)

        # Improved layout
        fig.update_layout(
            title=f"{analysis_type} of {data_type} data",
            xaxis=dict(
                title=f"{coords.columns[0]} ({explained_var[0]:.2%})",
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
                title=f"{coords.columns[1]} ({explained_var[1]:.2%})",
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
        fig.add_annotation(
            x=0.98,
            y=0.02,
            xref="paper",
            yref="paper",
            text=f"Total genes: {len(coords)}<br>Outliers: {sum(is_outlier)}",
            showarrow=False,
            font=dict(size=12),
            align="right",
            bordercolor="black",
            borderwidth=1,
            bgcolor="white",
            opacity=0.8
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
        rscu_values = data.rscu_values

        # Filter out stop codons
        non_stop_codons = [c for c in CODONS if data.config.get_aa_map()[
            c] != 'STOP']
        rscu_filtered = rscu_values[non_stop_codons]

        # Create the heatmap
        fig = go.Figure(data=go.Heatmap(
            z=rscu_filtered.values.T,
            x=rscu_filtered.index,
            y=non_stop_codons,
            colorscale='Viridis',
            colorbar=dict(title='RSCU Value'),
            hovertemplate='Gene: %{x}<br>Codon: %{y}<br>RSCU: %{z:.2f}<extra></extra>'
        ))

        # Update layout
        fig.update_layout(
            title='RSCU Values Heatmap',
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
        else:  # text visualizations
            # Implement text visualization classes in the next version
            pass

        raise ValueError(
            f"Unsupported visualization type: {vis_type} with method: {vis_method}")


class FileManager:
    @staticmethod
    def get_output_path(data, base_filename, output_type):
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

        return os.path.join(
            output_dir,
            f"{base_filename}_{output_type}_{timestamp}")

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

    def load_fasta(self, file_path):
        """
        Load sequences from a FASTA file with optimized performance.

        Key optimizations:
        - Use Bio.SeqIO more efficiently
        - Minimize repeated computations
        - Use generators and iterators
        - Reduce memory overhead
        """
        try:
            # Use a generator to reduce memory usage
            def sequence_generator():
                for record in SeqIO.parse(file_path, "fasta"):
                    yield record

            # First pass: quick validation and count
            print(f"\nLoading sequences from {file_path}...")
            sequence_count = sum(1 for _ in SeqIO.parse(file_path, "fasta"))

            if sequence_count == 0:
                print(f"\nNo sequences found in {file_path}")
                return False

            # Reset file pointer
            print(f"Total sequences detected: {sequence_count}")

            # Direct loading with progress tracking
            self.sequences = []
            self.gene_names = []

            print("Loading sequences...", end=" ", flush=True)

            # Use a more memory-efficient loading approach
            progress_interval = max(1, sequence_count // 10)  # Update every 10%

            for i, record in enumerate(SeqIO.parse(file_path, "fasta"), 1):
                self.sequences.append(record)
                self.gene_names.append(record.id)

                # Progress indicator
                if i % progress_interval == 0:
                    print(f"{i/sequence_count*100:.0f}%", end=" ", flush=True)

            print("100% Done!")

            # Immediately process sequences to avoid multiple passes
            self._process_sequences()

            print(f"\nSuccessfully loaded {len(self.sequences)} sequences from {file_path}")
            self.file_path = file_path

            return True

        except Exception as e:
            print(f"\nError loading FASTA file: {e}")
            return False

    def _process_sequences(self):
        """Process all loaded sequences to generate codon counts and other statistics."""
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

        print(f"Processing {len(self.sequences)} sequences...")
        print("Progress: ", end="", flush=True)

        last_percent = -1
        total_seqs = len(self.sequences)

        for i, seq in enumerate(self.sequences):
            seq_id = seq.id
            dna_seq = str(seq.seq).upper()

            # Update progress indicator
            current_percent = int((i / total_seqs) * 100)
            if current_percent > last_percent:
                if current_percent % 10 == 0 and current_percent != 0:
                    print(f"{current_percent}%", end="", flush=True)
                else:
                    print(".", end="", flush=True)
                last_percent = current_percent

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

        print(" Done!")

        # Check total codon counts
        total_codons = self.codon_counts.sum().sum()
        print(
            f"Total codons counted across all sequences: {int(total_codons):,}")
        # Calculate RSCUs and amino acid usage
        self._calculate_rscu()
        self._calculate_aa_usage()
        self._calculate_gc3s()

    def _calculate_rscu(self):
        """Calculate RSCU (Relative Synonymous Codon Usage) values."""
        aa_map = self.config.get_aa_map()
        # Get the correct AA_TO_CODONS mapping for current genetic code
        aa_to_codons = AA_TO_CODONS_MAP.get(
            self.config.genetic_code, AA_TO_CODONS_MAP[1])

        self.rscu_values = pd.DataFrame(
            0.0, index=self.gene_names, columns=CODONS)

        for gene in self.gene_names:
            # Group codons by amino acid
            aa_counts = defaultdict(int)
            for codon in CODONS:
                aa = aa_map[codon]
                if aa != 'STOP':
                    aa_counts[aa] += self.codon_counts.at[gene, codon]

            # Calculate RSCU for each codon
            for codon in CODONS:
                aa = aa_map[codon]
                if aa != 'STOP' and aa_counts[aa] > 0:
                    codon_count = self.codon_counts.at[gene, codon]
                    synonymous_codons = len(aa_to_codons[aa])
                    self.rscu_values.at[gene, codon] = codon_count * \
                        synonymous_codons / aa_counts[aa]

    def _calculate_aa_usage(self):
        """Calculate amino acid usage."""
        aa_map = self.config.get_aa_map()
        # Get unique amino acids excluding stop codons
        unique_aas = sorted(
            set([aa for aa in aa_map.values() if aa != 'STOP']) | {'STOP'})

        self.amino_acid_usage = pd.DataFrame(
            0, index=self.gene_names, columns=unique_aas)

        for gene in self.gene_names:
            for codon in CODONS:
                aa = aa_map[codon]
                codon_count = self.codon_counts.at[gene, codon]
                if aa == 'STOP':
                    self.amino_acid_usage.at[gene, 'STOP'] += codon_count
                else:
                    self.amino_acid_usage.at[gene, aa] += codon_count

    def _calculate_gc3s(self):
        """Calculate GC3s (GC content at synonymous third positions)."""
        aa_map = self.config.get_aa_map()
        # Get the correct AA_TO_CODONS mapping for current genetic code
        aa_to_codons = AA_TO_CODONS_MAP.get(
            self.config.genetic_code, AA_TO_CODONS_MAP[1])

        for gene in self.gene_names:
            synonymous_third_positions = 0
            gc3s_count = 0

            for codon in CODONS:
                aa = aa_map[codon]
                if aa != 'STOP' and len(
                        aa_to_codons[aa]) > 1:  # Only count synonymous codons
                    codon_count = self.codon_counts.at[gene, codon]
                    if codon_count > 0:
                        synonymous_third_positions += codon_count
                        if codon[-1] in 'GC':
                            gc3s_count += codon_count

            if synonymous_third_positions > 0:
                self.base_composition.at[gene, 'GC3s'] = gc3s_count / \
                    synonymous_third_positions * 100
            else:
                self.base_composition.at[gene, 'GC3s'] = 0

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
            mean_rscu = self.rscu_values.loc[reference_genes].mean()

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

        # Prepare data matrix
        print("Preparing data matrix...", end="", flush=True)
        if data_type == 'RSCU':
            data_matrix = self.rscu_values.copy()
            # Filter out stop codons
            aa_map = self.config.get_aa_map()
            stop_codons = [codon for codon in CODONS if aa_map[codon] == 'STOP']
            for codon in stop_codons:
                if codon in data_matrix.columns:
                    data_matrix = data_matrix.drop(codon, axis=1)
        elif data_type == 'AA':
            data_matrix = self.amino_acid_usage.copy()
            if 'STOP' in data_matrix.columns:
                data_matrix = data_matrix.drop('STOP', axis=1)
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
        print("Running analysis calculations...")
        print("Progress: ", end="", flush=True)

        if analysis_type == 'PCA':
            print("10%", end="", flush=True)
            scaled_data = stats.zscore(data_matrix, axis=0)
            print("...30%", end="", flush=True)

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

            Dr_inv_sqrt = np.diag(1.0 / np.sqrt(r.flatten()))
            Dc_inv_sqrt = np.diag(1.0 / np.sqrt(c.flatten()))
            print("...60%", end="", flush=True)

            S = Dr_inv_sqrt @ (P - np.outer(r, c)) @ Dc_inv_sqrt
            print("...70%", end="", flush=True)

            u, s, vh = np.linalg.svd(S, full_matrices=False)
            print("...80%", end="", flush=True)

            n_components = min(5, len(s))

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

        self.multivariate_results = result
        return result

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

        # Find the gene in our sequences
        gene_seq = None
        for seq in self.sequences:
            if seq.id == gene_id:
                gene_seq = seq
                break

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

            # Open the visualization if saved to a path
            if result_path:
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
        left_data._calculate_rscu()

        right_data = GCUAData(self.config)
        right_data.gene_names = ["Right_Cohort"]
        right_data.codon_counts = pd.DataFrame([right_codon_usage], index=["Right_Cohort"])
        right_data._calculate_rscu()

        aa_map = self.config.get_aa_map()
        aa_to_codons = AA_TO_CODONS_MAP.get(self.config.genetic_code, AA_TO_CODONS_MAP[1])

        results = {
            "left_cohort": left_extreme_genes,
            "right_cohort": right_extreme_genes,
            "left_cohort_usage": {codon: left_codon_usage[codon] for codon in CODONS},
            "right_cohort_usage": {codon: right_codon_usage[codon] for codon in CODONS},
            "left_cohort_rscu": {codon: left_data.rscu_values.at["Left_Cohort", codon] for codon in CODONS},
            "right_cohort_rscu": {codon: right_data.rscu_values.at["Right_Cohort", codon] for codon in CODONS},
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
        menu_options = [
            ("Load FASTA file", "1", self._input_fasta),
            ("Analysis", "2", self._analysis_menu),
            ("Visualization", "3", self._visualization_menu),
            ("Sequence Optimization", "4", self._optimization_menu),
            ("Export Data", "5", self._export_menu),
            ("Preferences", "6", self._preferences_menu),
            ("Help", "7", self._help),
            ("Quit program", "Q", lambda: False)
        ]

        while True:
            self._clear_screen()
            self.display_banner()

            # Show data status
            if self.data.sequences:
                print(
                    f"\nCurrent data: {len(self.data.sequences)} sequences from '{self.data.file_path}'")
                print(
                    f"Using genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
            else:
                print("\nNo data loaded.")
                print(
                    f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

            # Display menu with improved formatting
            print("\n+----------------------+")
            print("|      MAIN MENU       |")
            print("+----------------------+")

            # Calculate padding for alignment
            max_option_length = max(len(option)
                                    for option, _, _ in menu_options)

            for option, key, _ in menu_options:
                # Right-align the key and add padding to align all options
                padding = " " * (max_option_length - len(option))
                print(f"  {key}. {option}{padding}")

            choice = input("\nEnter your choice: ").strip().upper()

            for option, key, handler in menu_options:
                if choice == key:
                    if key == "Q":  # only exit program on Q
                        return
                    handler()
                    break
            else:
                if choice != 'Q':  # 'Q' is already handled in the loop
                    print("\nUnrecognized command.")
                    self._pause()

    def _display_menu(self, title, options):
        """
        Enhanced menu display with a cleaner, more minimalist layout.
        """
        while True:
            self._clear_screen()

            # Title with underline
            print(title)
            print("=" * len(title))

            # Show data status if applicable
            if self.data.sequences:
                print(f"\nCurrent Data: {len(self.data.sequences)} sequences")
                print(f"Source: {self.data.file_path or 'N/A'}")

            print("\nOptions:")
            # Calculate padding for alignment
            max_option_length = max(len(f"{option_key}. {option_text}") for option_text, option_key, _ in options)

            for option_text, option_key, _ in options:
                # Right-align the key and add padding to align all options
                full_option = f"{option_key}. {option_text}"
                print(full_option)

            print("\nR. Return to previous menu")

            choice = input("\nSelect an option: ").strip().upper()

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

    def _input_fasta(self):
        """Handle FASTA file input."""
        file_path = input("\nName of the FASTA-formatted file: ").strip()
        if file_path:
            self.data.load_fasta(file_path)
            self._pause()
        return True

    def _visualization_menu(self):
        """Handle visualization options with enhanced plotting capabilities."""
        if not self._check_data():
            return True

        options = [
            ("Create multivariate analysis plot", "1", self._multivariate_plot),
            ("Create GC content vs GC3 plot", "2", self._gc_content_plot),
            ("Create ENC vs GC3s plot (Wright's plot)", "3", self._enc_plot),
            ("Create RSCU heatmap", "4", self._rscu_heatmap),
            ("Create CAI distribution plot", "5", self._cai_distribution_plot),
            ("Show codon bias overview", "6", self._codon_bias_overview),
            ("Create custom scatter plot", "7", self._custom_scatter_plot),
            ("Return to main menu", "R", lambda: False)
        ]

        return self._display_menu("Visualization Menu", options)

    def _preferences_menu(self):
        """Handle program preferences menu."""
        options = [
            ("Genetic Code", "1", self._genetic_code_menu),
            ("Output Detail Level", "2", self._detail_level_menu),
            ("Visualization Method", "3", self._visualization_method_menu),
            ("Toggle Beep", "4", self._toggle_beep),
            ("Return to main menu", "R", lambda: False)
        ]

        return self._display_menu("Preferences Menu", options)

    def _toggle_beep(self):
        """Toggle beep setting."""
        self.config.beep_on = not self.config.beep_on
        print(f"\nBeep {'enabled' if self.config.beep_on else 'disabled'}.")
        self._pause()
        return True

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
            if self.data.sequences:
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

    def _multivariate_plot(self):
        """Create and display a multivariate analysis plot."""
        if not hasattr(
                self.data,
                'multivariate_results') or not self.data.multivariate_results:
            # Perform analysis if not done yet
            print("\nNo multivariate analysis results found. Performing analysis...")
            result = self.data.perform_multivariate_analysis()
        else:
            result = self.data.multivariate_results

        # Create and display the plot
        output_path = self.file_manager.get_output_path(
            self.data, None, "multivariate_plot")
        self.data.visualize('multivariate', result, output_path)
        self._pause()
        return True

    def _gc_content_plot(self):
        """Create and display GC content plot."""
        print("\nCreating GC content vs GC3 plot...")
        output_path = self.file_manager.get_output_path(
            self.data, None, "gc_content_plot")
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
            self.data, None, "enc_plot")
        self.data.visualize('enc', output_path=output_path)
        self._pause()
        return True

    def _rscu_heatmap(self):
        """Create and display RSCU heatmap."""
        print("\nCreating RSCU heatmap...")
        print(
            f"Current genetic code: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")
        output_path = self.file_manager.get_output_path(
            self.data, None, "rscu_heatmap")
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
            self.data, None, "cai_distribution")
        self.data.visualize('cai_distribution', cai_values, output_path)
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
        """Calculate all metrics for the loaded sequences."""
        if not self._check_data():
            return True

        print("\nCalculating all metrics...")

        # Perform all analyses
        print("Analyzing base composition...")
        # This was done during sequence loading

        print("Performing multivariate analysis...")
        self.data.perform_multivariate_analysis()

        print("Calculating ENC values...")
        self.data.calculate_enc()

        print("Calculating optimal codons...")
        self.data.calculate_optimal_codons()

        print("Calculating Fop values...")
        self.data.calculate_fop()

        print("Calculating CAI values...")
        self.data.calculate_cai()

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
        """Handle analysis options."""
        if not self._check_data():
            return True

        options = [
            ("Calculate all metrics", "1", self._calculate_all_metrics),
            ("Calculate codon usage", "2", self._codon_usage_menu),
            ("Calculate amino acid usage", "3", self._aa_usage_menu),
            ("Calculate base composition", "4", self._base_composition_menu),
            ("Perform multivariate analysis", "5", self._multivariate_menu),
            ("Calculate ENC values", "6", self._enc_menu),
            ("Calculate Fop and CAI values", "7", self._fop_cai_menu),
            ("Calculate SCUO values", "8", self._scuo_menu),
            ("Identify optimal codons", "9", self._optimal_codons_menu),
            ("Compare codon usage between axis cohorts",
             "0", self._compare_codon_usage_axis_cohorts),
            ("Return to main menu", "R", lambda: False)
        ]

        return self._display_menu("Analysis Menu", options)

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

    def _multivariate_menu(self):
        """Handle multivariate analysis options."""
        if not self._check_data() or len(self.data.gene_names) < 2:
            print("\nThere must be more than a single sequence in memory")
            self._pause()
            return True

        options = [
            ("Perform Correspondence Analysis (CA) on RSCU values",
             "1", lambda: self._perform_multivariate('CA', 'RSCU')),
            ("Perform Principal Component Analysis (PCA) on RSCU values",
             "2", lambda: self._perform_multivariate('PCA', 'RSCU')),
            ("Perform CA on amino acid usage", "3",
             lambda: self._perform_multivariate('CA', 'AA')),
            ("Perform PCA on amino acid usage", "4",
             lambda: self._perform_multivariate('PCA', 'AA')),
            ("Visualize multivariate analysis results", "5", self._multivariate_plot),
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
                print(f"  {gene}: {row[0]:.4f}")
            if i >= len(sorted_genes) - 5:
                print(f"  {gene}: {row[0]:.4f}")
            if i == 4 and len(sorted_genes) > 10:
                print("  ...")

        # Add genetic code information
        print(
            f"\nGenetic code used: [{self.config.genetic_code}] {self.config.get_genetic_code_name()}")

        # Ask if user wants to visualize or export
        print("\nOptions:")
        print("1. Visualize these results")
        print("2. Export these results to a file")
        print("R. Return to multivariate menu")

        choice = input("\nEnter your choice: ").strip().upper()

        if choice == '1':
            # Visualize the results
            self._visualize_multivariate(result, analysis_type, data_type)
        elif choice == '2':
            self._export_multivariate(result, analysis_type, data_type)

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
        if not self.data.sequences:
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

    GCUA is a tool for analyzing codon usage patterns in DNA sequences.

    Basic Usage:
    1. Load a FASTA file containing DNA sequences
    2. Perform analysis on the loaded sequences
    3. Visualize the results
    4. Export data for further analysis

    Available Analysis Types:
    - Codon Usage: Calculates codon frequencies and RSCU values
    - Amino Acid Usage: Analyzes amino acid composition
    - Base Composition: Examines nucleotide content and GC bias
    - Multivariate Analysis: Performs correspondence analysis or PCA
    - ENC Values: Calculates Effective Number of Codons
    - Fop and CAI: Measures codon optimization and adaptation
    - SCUO: Analyzes synonymous codon usage order
    - Optimal Codons: Identifies preferred codons

    Visualization Options:
    - Multivariate plots: View correspondence analysis results
    - GC content plots: Compare GC content with GC3 composition
    - ENC plots: Create Wright's plot of ENC vs GC3s
    - RSCU heatmaps: Visualize RSCU values across genes
    - CAI distribution: Examine distribution of CAI values
    - Custom scatter plots: Create plots of any two metrics

    Sequence Optimization:
    - Optimize gene sequences using various reference strategies
    - Export optimized sequences for expression experiments

    Genetic Codes:
    - Multiple genetic code tables are supported (NCBI translation tables)
    - Change the genetic code in the Preferences menu

    For detailed information, visit https://github.com/mol-evol/gcua
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
