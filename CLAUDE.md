# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GCUA (General Codon Usage Analysis) is a comprehensive Python bioinformatics tool for analyzing codon usage patterns in DNA sequences. It's a single-file application (`gcua.py`) that provides:
- Codon usage analysis and bias metrics (ENC, CAI, Fop, SCUO)
- Multivariate analysis (CA and PCA)
- Interactive visualizations using Plotly
- Sequence optimization based on codon usage patterns
- Support for 33 different genetic codes

## Recent Updates (2025 Menu Reorganization)

### Menu System Overhaul
The menu system has been completely reorganized for better workflow:

1. **Data Management** - Load, view, filter, and manage sequence data
2. **Quick Analysis** - Guided workflow for standard analyses
3. **Custom Analysis** - Detailed analysis options organized by type
4. **Visualization & Export** - Unified menu for plots and data export  
5. **Advanced Tools** - Sequence optimization, clustering, outlier detection
6. **Settings & Preferences** - Configure genetic codes and analysis options

### Key Improvements
- **Workflow-based organization** replaces function-based menus
- **Status indicators** show which analyses have been completed
- **Context-aware options** gray out unavailable choices
- **Quick Actions** for common workflows (e.g., "Calculate All Basic Metrics")
- **Grouped analysis options** (Basic Metrics, Advanced Analysis, Comparative)
- **Unified export menu** streamlines data output

## Development Commands

### Setup and Installation
```bash
# Install all dependencies
pip install -r requirements.txt

# Verify installation
python gcua.py --version
```

### Running the Application
```bash
# Launch the interactive menu-driven interface
python gcua.py
```

Note: There are no command-line arguments. All functionality is accessed through the interactive menu system.

### Development Workflow
- **No formal testing framework** - Testing is done manually through the interactive interface
- **No linting or type checking commands** are currently set up
- **No build process** - This is a standalone Python script
- All outputs are saved to `gcua_outputs/` directory with timestamped filenames

## Architecture and Code Structure

### Single-File Architecture
The entire application is contained in `gcua.py` (~5500+ lines). Key components:

1. **GCUAProcessor Class**: Main application class containing all functionality
   - Data loading and validation
   - Metric calculations
   - Multivariate analysis
   - Visualization generation
   - Sequence optimization

2. **Menu-Driven Interface**: Hierarchical menu system for user interaction
   - Main Menu → Sub-menus → Functions
   - Input validation and error handling
   - Clear screen functionality for better UX

3. **Data Flow**:
   - FASTA sequences → Internal DataFrame storage
   - Calculations stored in class attributes
   - Results exported to `gcua_outputs/` directory

### Key Technical Considerations

1. **Genetic Code Handling**:
   - Uses BioPython's genetic code tables
   - CRITICAL: Genetic code must be set correctly before analysis
   - Changing genetic code requires reloading FASTA files
   - Supports 33 different genetic codes (NCBI translation tables)

2. **Memory Management (Enhanced in v2.5.0)**:
   - **Hybrid approach**: Automatic mode selection based on file size
   - **Small files (<100MB)**: Loaded entirely into memory for speed
   - **Large files (>100MB)**: Sequences loaded on-demand to conserve memory
   - **Sequence indexing**: Maps gene IDs to file positions for fast random access
   - DataFrames store all calculated metrics
   - Memory usage scales with number of genes, not sequence length

3. **Performance Optimizations (New in v2.5.0)**:
   - **Parallel Processing**: 
     - Uses multiprocessing for codon counting
     - Configurable number of processes
     - Batch processing for optimal throughput
   - **Checkpointing**:
     - Automatic progress saving at configurable intervals
     - Resume from last checkpoint after interruption
     - Checkpoint files: `[filename]_checkpoint.pkl`
   - **Batch Processing**:
     - Sequences processed in configurable batches
     - Default: 1000 sequences per batch

4. **File I/O Patterns**:
   - Input: FASTA format only (must be coding sequences divisible by 3)
   - Output: TSV (data), HTML (plots), FASTA (optimized sequences), JSON (complex data)
   - All outputs include timestamps in filenames
   - Checkpoint files for progress recovery

5. **Visualization Strategy**:
   - Plotly generates standalone HTML files
   - Files auto-open in default browser
   - Interactive features: zoom, pan, hover tooltips

6. **Error Handling**:
   - Extensive input validation
   - Try-except blocks around file operations
   - User-friendly error messages in the CLI
   - Graceful handling of interrupted processing

### Common Development Tasks

When modifying the code:
1. Follow existing code style (no formal style guide)
2. Maintain the monolithic structure (don't split into modules)
3. Add new menu options by extending the menu dictionaries
4. Use existing visualization patterns for new plots
5. Store new metrics as DataFrame columns in the class

### Working with Performance Features

When extending or modifying performance features:

1. **Adding New Parallel Operations**:
   - Create a module-level function that accepts serializable arguments
   - Use the `process_sequence_batch` pattern for batch processing
   - Return dictionaries that can be merged into DataFrames

2. **Extending Checkpointing**:
   - Add new data to the `checkpoint_data` dictionary in `_save_checkpoint()`
   - Update `_load_checkpoint()` to restore the new data
   - Ensure all data is pickleable

3. **Memory Management**:
   - Use `get_sequence()` method for unified sequence access
   - Check `self.memory_mode` before accessing `self.sequences`
   - Update sequence index when adding new sequence access patterns

4. **Configuration**:
   - Add new settings to the `Config` class `__init__` method
   - Create corresponding menu options in preferences
   - Document default values and ranges

### Dependencies and Compatibility
- Python 3.8+ required
- Cross-platform (Windows, macOS, Linux)
- No external configuration files
- No environment variables needed
- Optional: scikit-learn for PCA (gracefully handled if missing)
- Optional: kaleido for static image export (SVG, PNG, PDF)

## Recent Improvements and Features

### Smart Visualizations (2025)
- **RSCU Heatmap**: Automatically adapts to dataset size
  - Small datasets (<1000 genes): Traditional gene × codon heatmap
  - Large datasets (≥1000 genes): Summary statistics heatmap with mean, std dev, percentiles
  - Prevents browser crashes with 500K+ genes
  - Shows meaningful patterns instead of overwhelming detail
  - Codons grouped by amino acid in summary mode

### Export Formats (2025)
- Added support for SVG, PNG, PDF export in addition to HTML
- Configurable via Preferences > Visualization Format
- Requires kaleido package for non-HTML formats
- Graceful fallback to HTML if kaleido not installed

### Menu System Fixes (2025)
- Implemented all missing menu methods for Analysis submenus
- Fixed AttributeError issues in base composition, codon usage, amino acid usage, and ENC menus
- Added proper pagination for large result displays

### Extreme Genes Export (2025)
- Added comprehensive extreme genes export functionality (Analysis Menu → Option E)
- Features:
  - Export genes from extremes of any calculated metric (GC, GC3s, ENC, CAI, Fop, SCUO, multivariate axes)
  - Flexible selection: top X%, bottom X%, or both extremes
  - Customizable percentage thresholds (1-50%)
  - Comprehensive output includes all available metrics for exported genes
  - Metadata headers with selection criteria and thresholds
  - Automatic ranking and percentile calculation
- Use cases:
  - Identify potential highly/lowly expressed genes
  - Find genes with extreme codon bias
  - Extract outliers for experimental validation
  - Create gene sets for downstream analysis

### CA Biplot Visualization (2025)
- Enhanced Correspondence Analysis visualization to show true biplots
- Features:
  - Displays both genes (rows) and columns in the same coordinate space
  - Columns are either codons (for RSCU analysis) or amino acids (for AA analysis)
  - Visual distinction: genes as circles, codons/amino acids as green diamonds with labels
  - Interactive: hover to see exact coordinates
  - Labels can be toggled on/off via legend for cleaner visualization
  - Works for all axis combinations (Axis1 vs Axis2, etc.)
  - Correctly labels items as "Codons" or "Amino acids" based on analysis type
- Benefits:
  - See which codons/amino acids are associated with gene clusters
  - Identify variables driving separation along each axis
  - Understand preferences of different gene groups
  - Standard CA biplot interpretation

### Enhanced Export Extreme Genes Menu (2025)
- Improved to handle multiple multivariate analyses properly
- Features:
  - Now accesses ALL performed multivariate analyses from history
  - Clear labeling: "CA on RSCU - Axis1 (45.2%)" shows analysis type, data type, axis, and variance
  - Organized menu with BASE METRICS and MULTIVARIATE ANALYSES sections
  - Sequential numbering for easy selection
  - Shows up to 3 axes per multivariate analysis
- File naming:
  - Multivariate exports include analysis details: "extreme_genes_CA_RSCU_Axis1_top_10pct.tsv"
  - Metadata headers clearly indicate which analysis was used
- Benefits:
  - Access extremes from any performed analysis, not just the last one
  - Clear understanding of what each option represents
  - Better organization for datasets with multiple analyses

### Efficient Metric Caching (2025)
- Added intelligent caching to avoid redundant calculations
- Features:
  - Checks which metrics are already calculated before recalculating
  - Shows "Using cached" vs "Calculating" messages
  - Option to force recalculation if needed
  - "Clear cached metrics" option in Analysis menu (option C)
- Benefits:
  - Second call to custom scatter plot uses cached values instantly
  - No redundant expensive calculations (CA, ENC, CAI, etc.)
  - Significant time savings for large datasets (200K+ genes)
  - Memory overhead is minimal (~15-20 MB for 200K genes)
- Clear cache when:
  - Changing genetic code
  - Modifying reference genes
  - Switching between different analyses

### Proper CA Dimensions for RSCU (2025)
- Fixed multivariate analysis to exclude single-codon amino acids
- Changes:
  - Automatically removes Met (ATG) and Trp (TGG) for standard genetic code
  - These codons have no synonymous alternatives (RSCU always = 1.0)
  - Reduces dimensions from 61 to 59 for standard code
  - Shows excluded codons in output: "excluded 2 single-codon amino acids: ATG (Met), TGG (Trp)"
- Benefits:
  - Focuses analysis on actual codon usage choice, not amino acid presence
  - Cleaner biplots without uninformative points
  - Matches standard practice in codon usage analysis
  - More meaningful principal axes

### CLR Default for Amino Acid PCA (2025)
- Made CLR transformation the default for amino acid compositional data analysis
- Changes:
  - When performing PCA on amino acid data, CLR is now default (option 2)
  - User prompt shows "default=2 CLR for amino acids"
  - Added note explaining CLR is mathematically appropriate for compositional data
  - Still allows users to choose standard frequencies or log transformation
- Benefits:
  - Proper handling of compositional constraints (amino acids sum to 100%)
  - Avoids spurious correlations from closure effect
  - First PC no longer dominated by protein length variation
  - More meaningful biological patterns revealed
  - Follows best practices for compositional data analysis

### Cluster Detection Fix (2025)
- Fixed ValueError when characterizing clusters: "Item wrong length X instead of Y"
- Issue: Cluster labels were assigned to genes in multivariate results (subset) but characterization tried to apply to full dataset
- Solution:
  - Filter data matrices to only include genes that were actually clustered
  - Added gene_indices to cluster_results to track which genes were included
  - Updated GC content calculations to handle missing genes gracefully
  - Added None checks for display of mean values
- The fix ensures cluster characterization works correctly with filtered multivariate data

### Quick Analysis Workflow Fix (2025)
- Fixed AttributeError: 'GCUAData' object has no attribute 'calculate_all'
- Issue: Quick Analysis workflow was calling non-existent calculate_all() method
- Solution:
  - Updated workflow to use existing calculation methods
  - Basic metrics (codon usage, amino acid usage, base composition) are calculated during loading
  - Added RSCU calculation as the additional basic metric
  - Added proper status flag updates throughout the workflow
  - Created non-interactive versions of plot creation for workflow
  - Updated _create_standard_plots to generate plots without user prompts
- The workflow now completes successfully without errors

### Custom Analysis Menu Reorganization (2025)
- Fixed confusing duplicate options ("Calculate all metrics" vs "Calculate All Basic Metrics")
- Reorganized menu into clearer sections:
  - **Individual Analyses**: Specific metric details (codon, AA, base comp, RSCU)
  - **Advanced Analysis**: Multivariate, codon bias metrics, optimal codons
  - **Comprehensive Options**: Calculate ALL metrics, quick workflows
  - **Utilities**: Clear cache, show detailed status
- Added new features:
  - RSCU menu with display and export options
  - Detailed analysis status display (option S)
  - Clearer distinction between individual vs comprehensive calculations
- Menu is now more intuitive and avoids redundancy

### Cluster Visualization Fix (2025)
- Fixed AttributeError: 'FigWrapper' object has no attribute 'create'
- Issue: Cluster visualization was using incorrect method to export plots
- Solution:
  - Replaced custom FigWrapper approach with direct plotly export methods
  - Added support for both 2D and 3D cluster visualizations
  - When 3+ dimensions used in clustering, user can choose 2D or 3D view
  - Properly handles different export formats (HTML, SVG, PNG, PDF)
  - Auto-opens HTML plots in browser
- Enhanced visualization with proper axis labels showing explained variance

### Outlier Analysis Implementation (2025)
- Implemented comprehensive outlier detection functionality (was placeholder)
- **Four detection methods**:
  1. **Statistical outliers**: Z-score based detection on any metric (GC, GC3s, length, ENC, CAI)
  2. **Multivariate outliers**: Mahalanobis distance in PCA/CA space
  3. **Domain-specific outliers**: Extreme codon usage patterns (RSCU >2 or <0.5)
  4. **Combined approach**: Consensus outliers flagged by multiple methods
- **Key features**:
  - Interactive metric selection with availability checking
  - Customizable thresholds (Z-score, percentiles)
  - Visualization of outliers on multivariate plots (2D/3D)
  - Detection of stop codon readthrough (potential pseudogenes)
  - Export functionality for all outlier lists
  - Consensus analysis showing genes flagged by multiple methods
- **Use cases**:
  - Quality control (identify problematic sequences)
  - Find genes with unusual expression patterns
  - Detect potential annotation errors
  - Identify horizontally transferred genes

### Text-Based Help System Implementation (2025)
- Added comprehensive help system to improve usability for unfamiliar menu items
- **Help System Features**:
  - Quick help: Type option number followed by '?' (e.g., '1?') for specific help
  - Full menu help: Type 'H' to see help for all menu items  
  - Context-sensitive help content stored in `self.help_content` dictionary
  - Help hint shown on first menu display with 💡 tip (can be disabled)
  - Both brief (2-3 lines) and detailed explanations for each option
- **Implementation Details**:
  - Modified `_display_menu()` to accept `menu_key` parameter and handle help requests
  - Added `_show_quick_help()` method for single-option help display
  - Added `_show_menu_help()` method for comprehensive menu help
  - Added `_show_analysis_menu_help()` for custom analysis menu (which doesn't use _display_menu)
  - Created `_init_help_content()` method with help text for all major menus:
    - Main menu (explaining each major section)
    - Custom Analysis menu (detailed metric descriptions)
    - Advanced Tools menu (clustering, optimization, outliers)
    - Data Management menu (loading, filtering, saving)
    - Visualization & Export menu (plot types and formats)
    - Sequence Optimization menu (methods and usage)
- **User Experience**:
  - Non-intrusive - help only shown when requested
  - Clear navigation instructions in help screens
  - Tip shown once per session to inform users about help availability
  - Works with both standard menus (using _display_menu) and custom menus
  - Help content includes usage tips, requirements, and applications
- **Example help content**:
  - Brief: "Calculate codon bias metrics (ENC, CAI, Fop, SCUO)"
  - Detailed: Explains what each metric measures, value ranges, and interpretation

### Multivariate Plot Axis Selection (2025)
- Added interactive axis selection when creating multivariate plots
- **Changes**:
  - When visualizing multivariate results, users are now prompted to select which axes to plot
  - Shows available axes with their explained variance percentages
  - Allows selection of any valid axis combination (e.g., Axis1 vs Axis3)
  - Validates input to prevent plotting same axis against itself
  - Works for all output formats (HTML, SVG, PNG, PDF), not just HTML
- **Implementation**:
  - Modified `_multivariate_plot()` to prompt for axis selection
  - Updated `MultivariateVisualization` class to use selected axes
  - Axis information included in filename (e.g., "multivariate_plot_axis1_vs_axis3.svg")
  - Biplot loadings also use the selected axes
  - Axis labels show correct names and variance percentages
- **Benefits**:
  - Users can explore different dimension combinations
  - Especially useful for static formats (SVG, PNG, PDF) that don't have interactive buttons
  - Better understanding of data structure beyond just first two dimensions
  - Consistent experience across all visualization formats

### Quick Multivariate Workflow Fix (2025)
- Fixed AttributeError in Quick Multivariate Workflow: 'GCUAData' object has no attribute 'calculate_rscu'
- **Issue**: Code was calling `self.data.calculate_rscu()` which doesn't exist as a public method
- **Solution**: 
  - Replaced all calls to `calculate_rscu()` with `get_rscu_values()`
  - The `get_rscu_values()` method automatically calculates RSCU if needed using the private `_calculate_rscu()` method
  - This follows the existing pattern where RSCU is calculated on-demand
- **Affected areas**:
  - Quick Analysis Workflow
  - Quick Multivariate Workflow 
  - RSCU menu
  - Various other locations that calculate RSCU
- The fix ensures all RSCU calculations use the proper public interface

### Genetic Code Bug Fixes (2025)
- Fixed critical bugs in genetic code implementations where some codes had 65 codons instead of 64
- **Issues found and fixed**:
  - **Genetic Code 2 (Vertebrate Mitochondrial)**: Had invalid codon 'AUN' with ambiguous nucleotide
  - **Genetic Code 3 (Yeast Mitochondrial)**: Had invalid codon 'CUN' with ambiguous nucleotide
  - **Genetic Code 5 (Invertebrate Mitochondrial)**: Had invalid codon 'AUN' with ambiguous nucleotide
- **Solution**:
  - Removed all invalid codons containing 'N' (ambiguous nucleotide)
  - Valid codons only use A, U, G, C (or A, T, G, C for DNA)
  - All genetic codes now have exactly 64 codons as required
- **Verification**:
  - Tested all 25 genetic codes - each now has exactly 64 codons
  - Verified specific changes against NCBI Genetic Code standards
  - Confirmed amino acid assignments match NCBI specifications
- **Impact**: This was a critical fix as incorrect genetic codes would produce wrong amino acid translations and RSCU calculations

### Complete NCBI Genetic Code Verification (2025)
- Performed comprehensive verification of all 25 genetic codes against NCBI standards
- **Additional translation errors found and fixed**:
  - **Code 12 (Alternative Yeast Nuclear)**: CTG was incorrectly Ala, should be Ser
  - **Code 13 (Ascidian Mitochondrial)**: AAA was incorrectly Asn, should be Lys (standard)
  - **Code 27 (Karyorelict Nuclear)**: Missing TGA → Trp translation
  - **Code 28 (Condylostoma Nuclear)**: Missing TGA → Trp translation  
  - **Code 33 (Cephalodiscidae Mitochondrial)**: Missing AGA → Ser and AGG → Lys translations
- **Verification method**:
  - Compared all 64 codons for each genetic code against NCBI's official translation tables
  - Source: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
  - Tested both RNA (U) and DNA (T) versions of codons
- **Results**: 
  - All 25 genetic codes now correctly implement NCBI standards
  - Every codon translates to the correct amino acid
  - Proper handling of alternative start codons and stop codon reassignments
- **Critical importance**: Ensures accurate protein translation, codon usage analysis, and sequence optimization across all supported organisms