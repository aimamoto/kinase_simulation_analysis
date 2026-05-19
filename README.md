# KINOME Structural Bioinformatic Pipeline (Part 1 v6r6)

This pipeline automates the extraction, sequence alignment, and structural analysis of kinase protein simulations (AlphaFold 3 or MD). It is optimized for batch processing large sets of simulations, such as combinatorial mutant matrices (e.g., 4x8 EGFR/ERBB variants) and single-protein runs (e.g., ABL1 wild-type).

## Project Goal
This pipeline is developed as a sister project of another pipeline specific for the ERBB family dimers (github.com/aimamoto/erbb_dimer_complex). While the sister project relies on the landmark residues in the ERBB family, this project aims to extend the analysis schemes to a broad range of KINOME-wide kinases using hmmer. So far, tested in CDK1 (with or without Cyclin B1 and key phosphorylation), SRC-CSK heterodimer, the ERBB family members, ABL, FAK, JAK2, and PKA catalytic subunit. Part 1 is orchestrated by the wrapper 'run_kinase_pipeline_v6r6.py', while Part 2 is orchestrated by another wrapper 'run_allostery_discovery_v6r9.py' (not detailed in this README).

## Pipeline Capabilities & Evolution

**Core Analytical Features (Established in v4r1):**
1. **Advanced Structural Metrics:** The pipeline evaluates both the dynamic hydrophobic core and the αC-β4 loop allosteric network. 
   * *Hydrophobic Core:* Tracks the cohesiveness of the Hydrophobic Shell and the integrity of critical Core Bridges (V104 to R-Spine, I150 to αE helix).
   * *Allosteric Networks:* Quantifies the K105 toggle switch (sensing active vs. apo states), the αE helix anchor (Y156-N99), and the rigid deep αF-helix scaffold (D220).
2. **Visual Alignments (1D):** Automatically generates publication-ready Multiple Sequence Alignment (MSA) text panels and abstract 1D topological schematics mapping structural motifs perfectly to sequence coordinates.

**Recent Upgrades (v6r6):** 
1. **Universal Steric Co-Factor Integration:** Seamless integration and structural tracking of non-kinase co-factors (e.g., CDK-Cyclin, RAF-14-3-3). The pipeline now preserves non-kinase chains and maps their steric influence on the kinase αC-helix and Activation Loop.
2. **AlphaFold 3 Confidence Metrics:** Parallel extraction of AF3 predicted align errors (PAE) and template modeling scores (ipTM/pTM), automatically merged with the geometric data into a singular master CSV dataset.
3. **Parallel Execution Engine:** Massive speed improvements for high-throughput batch runs. Both sequence extraction (`extract_fasta.py`) and ChimeraX structural analysis now utilize an orchestrator/worker multiprocessing architecture. You can control this via the `-c` or `--cores` flag.
4. **Smart Filtering:** The pipeline automatically detects AF3 nested seeds (e.g., `seed-m_sample-n/model.cif`) and safely excludes redundant top-level `*_model.cif` summary files, preventing duplicate processing and skewed datasets.
5. **Expanded Catalytic & Ligand Tracking:** Explicit distance tracking for ATP and Magnesium coordination. Calculates critical distances for HRD-Asp to ATP, DFG-Asp to Mg/ATP, and P-Loop to ATP.

## Getting Started

### 1. System Requirements
Ensure your OS (Ubuntu, macOS, Windows) has the following installed and accessible via your system PATH:
* **UCSF ChimeraX**: Typing `chimerax` in the terminal must successfully launch the application (or run it headlessly).
* **HMMER**: Specifically the `hmmalign` command.
  * *Ubuntu:* `sudo apt install hmmer`
  * *macOS:* `brew install hmmer`
* **Pfam Kinase Profile**: Download `Pkinase.hmm` and place it at `~/pfam/Pkinase.hmm`. *(If stored elsewhere, modify the wrapper script to pass the `--hmm /your/path/` flag to the landmark extractor).*

### 2. Python Dependencies
Install the required Python packages (requires Python 3.x):
```bash
pip install -r requirements.txt
```
*(This installs `biopython`, `pandas`, `pyyaml`, and `matplotlib`. The multiprocessing features utilize standard built-in Python libraries.)*

### 3. The Two-Stage Workflow

#### Stage A: Discovery & Audit
Scan your current directory for simulation folders:
```bash
python3 run_kinase_pipeline_v6r6.py
```
* The script safely ignores non-simulation folders (like `modules/` or `archives/`).
* It strictly reads your folder names (e.g., `a-erbb2cattail_b-egfrcattail-t790m`) and generates a `generated_matrix.csv` manifest.
* **Select Option 2** when prompted to pause the pipeline. Open the CSV to verify your simulation matrix. You can easily drop unwanted runs by deleting the "x" in the corresponding cell.

#### Stage B: Parallel Analysis & Visualization
Resume the pipeline after verification, specifying the number of CPU cores to utilize:
```bash
python3 run_kinase_pipeline_v6r6.py --resume -c 16
```
The pipeline automatically handles:
1. **Archiving**: Moves any previous results to a timestamped folder in `archives/` to prevent data mixing.
2. **Setup**: Generates `proteins.yaml` and extracts `sequences.fasta` via parallel workers.
3. **Landmark Mapping**: Aligns sequences against the HMM profile to generate `hmm_landmarks.json`.
4. **1D Sequence Visualizations**: Generates `MSA_Annotated_Panel.pdf` and individual topological schematics (`kinome_VISalign.py`).
5. **3D Structural Execution**: Spawns parallel headless ChimeraX instances to extract 3D measurements.
6. **AF3 Metric Extraction**: Parses AF3 JSON files to compile PAE and ipTM confidence scores.
7. **Final Consolidation**: Merges geometric output and AF3 metrics into `master_kinase_analysis_results_v6r6.csv`.

### 4. Advanced: Custom Entry Points
If you already have a curated configuration or sequence file, you can bypass the discovery stages.

* **Start from YAML (--use-yaml)**:
```bash
python3 run_kinase_pipeline_v6r6.py --use-yaml -c 8
```
* **Start from FASTA (--use-fasta)**:
```bash
python3 run_kinase_pipeline_v6r6.py --use-fasta -c 8
```

## Utilities & Troubleshooting
* **`debug_landmarks.py`**: Run `python3 modules/debug_landmarks.py sequences.fasta` to print a terminal table verifying the exact amino acid letters matching your extracted HMM coordinate indices.
* **`make_bingo.py`**: Generate conditional association matrices and isolate single-variable changes across complex simulation datasets. 

## Outputs
* `master_kinase_analysis_results_v6r6.csv`: Comprehensive master dataset for all analyzed chains. Includes spatial/dihedral states, AlphaFold 3 confidences (`ipTM`, `PAE`), Co-factor proximity mapping, Allosteric nodes (`Y156_N99_Dist`, `K105_E107_Dist`), and Ligand metrics (`HRD_ATP_Dist`, `DFG_Mg_Dist`).
* `MSA_Annotated_Panel.pdf`: Wrapped, monospaced sequence alignment visually mapped to canonical landmarks and structural topology tracks.
* `[sequence_name]_schematic.pdf`: Individual 1D abstract topologies showing N/C-lobes, helices, sheets, and catalytic loop positioning.
* **Visualization Directories:**
  * `cx_viz_core/`: `.cxc` ChimeraX macros isolating the catalytic core, R/C-Spines, Salt Bridge, and Hydrophobic Shell.
  * `cx_viz_allosteric/`: `.cxc` ChimeraX macros isolating the αC-β4 toggle switch, αE anchor, and αF scaffold.

## Project Structure
```text
. (Working Directory)
├── run_kinase_pipeline_v6r6.py        # Main entry point orchestrator
├── requirements.txt                   # Python package dependencies
├── modules/                           # Core logic directory
│   ├── generate_config.py
│   ├── extract_fasta.py               # Parallel FASTA extraction
│   ├── extract_landmarks.py           # HMM-based structural mapping
│   ├── kinome_VISalign.py             # 1D Schematics and MSA plotting
│   ├── debug_landmarks.py             # Verification utility for HMM mapping
│   ├── 1_run_parallel_chimerax_hmm_v6r6.py
│   ├── chimerax_hmm_worker_v6r6.py    # Headless 3D geometric calculations
│   ├── extract_af3_metrics.py         # Extracts ipTM and PAE metrics
│   └── make_bingo.py                  # Statistical crossover matrix generator
├── a-erbb2_b-egfr.../                 # Input simulation directories (.cif/.pdb)
│
│   # --- Pipeline Generated Outputs (Created during execution) ---
├── archives/                          # Historical pipeline runs (auto-archived)
├── generated_matrix.csv               # Matrix manifest of discovered simulations
├── proteins.yaml                      # Strict naming configuration map
├── sequences.fasta                    # Extracted sequence FASTA
├── hmm_landmarks.json                 # HMM mapped canonical/allosteric residues
├── temp_af3_metrics.csv               # Temporary extracted confidence metrics
├── hmm_kinase_analysis_results_v6r6.csv # Temporary geometric measurements
├── master_kinase_analysis_results_v6r6.csv # Final merged dataset (Geometry + AF3)
├── MSA_Annotated_Panel.pdf            # Unified MSA structural annotation graphic
├── [kinase]_schematic.pdf             # Individual 1D topology maps
├── cx_viz_core/                       # 3D macros for the catalytic core (.cxc)
└── cx_viz_allosteric/                 # 3D macros for the allosteric network (.cxc)
```

---

## References
1. **Kinase Classification:** Modi, V., & Dunbrack, R. L., Jr (2019). "Defining a new nomenclature for the structures of active and inactive kinases." *PNAS*, 116(14), 6818-6827.
2. **Hydrophobic Core:** Kim, J. et al. (2017) "A dynamic hydrophobic core orchestrates allostery in protein kinases," *Science Advances*, 3(4). doi: 10.1126/sciadv.1600663.
3. **DFG-in/out Conformational Coupling:** Levinson, N. M. et al. (2006). "A Src-like inactive conformation in the abl tyrosine kinase domain." *PLoS Biology*, 4(5), e144.
4. **αC-β4 Loop Allostery:** Wu, J., Jonniya, N. A., et al. (2024). "Role of the aC-b4 loop in protein kinase structure and dynamics." *eLife*. doi: 10.7554/eLife.91980
