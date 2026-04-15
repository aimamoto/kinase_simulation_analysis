# Dual-Engine Kinase Allostery & Structural Analysis Pipeline (`analyzeR`)

This directory contains the master statistical, machine-learning, and structural evaluation suite for high-throughput AlphaFold 3 kinase ensembles. Upgraded to a **Dual-Engine Architecture**, it seamlessly handles both standard monomeric/homodimeric kinase systems and complex, asymmetric heterodimers (such as the ERBB family).

## 🛠️ System Requirements

- **Python 3.8+** (Standard libraries only: `subprocess`, `csv`, `glob`, `itertools`, `re`)
- **R 4.0+**
- **Required R Packages**: 
  ```R
  install.packages(c("tidyverse", "scales", "rstatix", "RVAideMemoire", "patchwork", 
                     "corrplot", "factoextra", "RColorBrewer", "cluster", "mclust", 
                     "plotly", "ggrepel", "pals"))
  ```
  *(Note for Ubuntu 22.04/24.04 users: Ensure system-level dependencies are installed via `sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev` prior to installing the `tidyverse` and `plotly` R packages).*

## 📂 Directory Structure

Ensure your `analyzeR` folder is structured as follows before execution. The Python wrapper expects the R scripts to reside within a `modules/` subdirectory.

```text
analyzeR/
│
├── run_allostery_discovery.py        # Master interactive Python wrapper
├── experiment.csv                    # Required ONLY for the Standard Engine (Ligand QC)
├── METRICS_CHEATSHEET.md             # Guide to biological mapping of CSV variables
│
└── modules/
    ├── multimer_core_engine.R        # Standard Kinase Engine (Phases 1-7)
    ├── erbb_asymmetric_engine.R      # ERBB/Asymmetric Engine (Phases 1-7)
    └── posthoc_differentiate.R       # Phase 8: Differential Allosteric Drivers
```

## 📥 Expected Inputs

1. **HMM Metadata CSVs**: The Python wrapper will automatically perform a recursive global search for any files ending in `_results_v6.csv` or `*_results.csv` in your active directories. It automatically ignores folders named `archive`, `old`, or `temp_chimerax_chunks`.
2. **`experiment.csv` (Standard Engine Only)**: If running Engine 1, you **must** have an `experiment.csv` in the root directory. The pipeline uses this to map the original experimental design (chains, mutations, ligands) and perform automated Truth-Overwrite for AF3 ligand QC (detecting when AF3 physically overrides a designed Apo/Holo configuration).

---

## 🚀 Execution Guide

The entire pipeline is driven by an interactive command-line interface. 

### Step 1: Launch the Pipeline
Navigate to the `analyzeR` directory and execute the Python wrapper:
```bash
./run_allostery_discovery.py
# or
python3 run_allostery_discovery.py
```

### Step 2: Choose Your Starting Point
You will be prompted to either:
1. **Run the Full Core Pipeline (Phases 1-7)**: Aggregates CSVs, runs the biological engine, and generates global clustering models.
2. **Skip to Post-Hoc Analysis (Phase 8)**: Bypasses the core engine if you already have a `Phase7_Complete_Structural_Metadata.csv` generated from a previous run and simply want to extract new pairwise differential drivers.

### Step 3: Select the Biological Engine
If running the Full Core Pipeline, you must select the appropriate thermodynamic logic:

* **Engine 1: Standard Kinase Engine**
  * **Best for**: Monomers, homodimers, or simple heterocomplexes (e.g., SRC, CSK).
  * **Key Features**: Aggressive AF3 Ligand QC (Truth Overwrite via `experiment.csv`), Global Macro-State categorical mapping, Dunbrack 2D Phase Space, and continuous Network Coupling (PCA & MAC).
* **Engine 2: ERBB Asymmetric Engine**
  * **Best for**: Complex asymmetric dimers evaluating Activator vs. Receiver roles.
  * **Key Features**: Computes global thermodynamic assembly probabilities (ATP-agnostic receiver probability), strictly maps Catalytic Spine (C-Spine) integrity, evaluates competitive binding preferences (1ATP logic), and measures catalytic clamping (2ATP logic).

### Step 4: Clustering Methodology
Choose how the pipeline should discover Meta-Stable states in Phase 6:
* **Gaussian Mixture Models (GMM) [Default]**: Best for continuous, dynamic structural transitions. Allows for non-spherical, probabilistic energetic basins.
* **K-Means + Gap Statistic**: Stricter, spherical basin mapping with mathematically rigorous K-optimization.

---

## 🔬 Pipeline Architecture (Phases 1-8)

Regardless of the engine selected, the pipeline strictly adheres to an 8-Phase analytical architecture:

* **Phase 1: Data Aggregation & QC** – Parses thousands of simulated models. The Standard engine applies Truth Overwrite for ligand mismatches; the ERBB engine maps Activator/Receiver directionality.
* **Phase 2: Macro-States & Global Assemblies** – Bar charts and Fisher's Exact tests defining the coarse global conformation (e.g., DFGin/out, Receiver %).
* **Phase 3: Structural Integrity & Phase Space** – Dunbrack Dihedral mapping (Standard) or C-Spine categorical matrices vs. ATP availability (ERBB).
* **Phase 4: Local Allosteric Micro-Metrics** – High-resolution violin plots of specific catalytic dials (e.g., Catalytic Clamping `HRD_ATP_Dist`, N/C-lobe anchors `aCb4_aE_Dist`).
* **Phase 5: Network Density & Correlation** – Generates Spearman correlation heatmaps and computes the Mean Absolute Correlation (MAC) to quantify the global network rigidity of the kinase.
* **Phase 6: Unsupervised 3D Discovery** – PCA and clustering (GMM/K-Means) to identify meta-stable states. Features an automated noise-reduction filter to merge extraneous data points. Outputs interactive HTML 3D scatter plots and static convex hull mappings.
* **Phase 7: KinCore Biological Signatures** – Maps the discovered mathematical clusters back to known biological categories (BLAminus, R-Spine Intact) using Fisher's exact tests. Outputs the master `Phase7_Complete_Structural_Metadata.csv`.
* **Phase 8: Post-Hoc Differential Drivers (Dual-Mode)**
  * **Trigger**: Prompts automatically after Phase 7, or can be run independently.
  * **Function**: Executes `posthoc_differentiate.R`. Allows for manual, sequential, or "lazy" pairwise comparisons between any two groups (e.g., State 1 vs State 2, or WT vs L858R).
  * **Output**: Volcano plots mapping Cohen's *d* Effect Size against FDR-adjusted p-values to isolate the specific structural drivers responsible for phenotypic shifts.

---

## 🔀 Under the Hood: Orchestration & Routing Mechanics

The `run_allostery_discovery.py` wrapper acts as a robust traffic controller, executing specific modules via `subprocess` calls and managing data hand-offs. 

### 1. Data Discovery & Engine Routing
To ensure a "Single Source of Truth," the Python wrapper handles all directory parsing. It identifies valid CSVs and writes their absolute paths to a temporary file (`.temp_active_csv_list.txt`). 
Based on your Engine selection, the wrapper triggers either `multimer_core_engine.R` or `erbb_asymmetric_engine.R`, passing three critical arguments:
1. The target protein name.
2. The clustering method string (`gmm` or `kmeans`).
3. The path to the temporary CSV list.
*(The R script reads this list, loads the data, and the Python wrapper subsequently deletes the temporary file to keep the environment clean).*

### 2. Meta-Stable State Discovery & Noise Reduction (Phase 6)
The pipeline maps complex, multi-dimensional distance metrics into discrete states using one of two selectable methodologies:
* **Gaussian Mixture Models (GMM) [Default]**: Utilizes the `mclust` package to identify non-spherical, probabilistic energetic basins. Optimal for continuous structural transitions and "breathing" motions.
* **K-Means + Gap Statistic**: Utilizes the `cluster` package for strict, spherical basin mapping. Computes the Gap Statistic across a range of $K$ values to mathematically prove the optimal number of clusters. **Note:** The pipeline explicitly defaults to the `globalSEmax` method after Dudoit and Fridlyand (2002) for optimal $K$ selection, abandoning the standard `Tibs2001SEmax` (Tibshirani et al., 2001) as internal testing proved `globalSEmax` to be more robust for evaluating these specific conformational ensembles.

**Micro-Cluster Consolidation (The Noise Filter):**
To prevent overfitting and statistical noise, both routing paths include an automated extraneous state filter. If a discovered state captures less than 5% of the total simulation models (with a hard minimum threshold of 10 models), it is flagged as a statistically insignificant "micro-cluster." The algorithm then calculates the Euclidean distance between cluster centroids and merges the micro-cluster into its nearest logical macroscopic neighbor.

### 3. Post-Hoc Differential Analysis Orchestration (Phase 8)
When Phase 8 is triggered, the Python wrapper shifts from a single-execution model to an iterative loop. It reads the resulting `Phase7_Complete_Structural_Metadata.csv` to identify all available comparison groups (e.g., distinct Macro-States or Experimental Conditions).

The user selects a routing mode:
* **Lazy Mode**: Automatically computes `itertools.combinations()` to run every possible unique pairwise comparison.
* **Sequential Mode**: Compares neighboring states linearly (State 1 vs 2, 2 vs 3, etc.).
* **Manual Mode**: Allows explicit target selection.

For each pair, the Python wrapper sanitizes the string inputs (converting newlines to `___` to prevent OS command-line breaks) and executes `posthoc_differentiate.R`. The R script decodes these arguments, filters the dataset exclusively to the two targets, calculates Cohen's $d$ and FDR-adjusted Wilcoxon statistics, and plots the differential volcano metrics.

---

## 📊 Locating Your Outputs

Upon completion, all data, PDFs, and interactive HTML files are routed to a dynamically generated folder based on your target, clustering method, and engine:

* **Standard Engine**: `./plots_and_stats_TARGET_GMM/`
* **ERBB Engine**: `./plots_and_stats_TARGET_GMM_ERBB/`
* **Phase 8 Results**: Found within a nested `./Phase8_Volcanos/` subdirectory containing specific `Stats_A_vs_B.csv` and high-resolution Volcano `.pdf` graphics.
