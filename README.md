# KINOME Structural Bioinformatic Pipeline (v6)

This pipeline automates the extraction, sequence alignment, and structural analysis of kinase protein simulations (AlphaFold3 or MD). It is optimized for batch processing large sets of simulations, such as combinatorial mutant matrices (e.g., 4x8 EGFR/ERBB variants) and single-protein runs (e.g., ABL1 wild-type).

## Project Goal
This pipeline is developed as a sister project of another pipeline specific for the ERBB family dimers (github.com/aimamoto/erbb_dimer_complex). While the sister project relies on the landmark residues in the ERBB family, this project aims to extend the analysis schemes to a broad range of KINOME-wide kinases using hmmer. So far, tested in the ERBB family members, SRC, CSK, ABL, FAK, JAK2, and PKA catalytic subunit.

## Pipeline Capabilities & Evolution

**Core Analytical Features (Established in v4r1):**
1. **Advanced Structural Metrics:** The pipeline evaluates both the dynamic hydrophobic core and the αC-β4 loop allosteric network. 
   * *Hydrophobic Core:* Tracks the cohesiveness of the Hydrophobic Shell and the integrity of critical Core Bridges (V104 to R-Spine, I150 to αE helix).
   * *Allosteric Networks:* Quantifies the K105 toggle switch (sensing active vs. apo states), the αE helix anchor (Y156-N99), and the rigid deep αF-helix scaffold (D220).
2. **Visual Alignments (1D):** Automatically generates publication-ready Multiple Sequence Alignment (MSA) text panels and abstract 1D topological schematics mapping structural motifs perfectly to sequence coordinates.

**Recent Upgrades (v6):** 
1. **Parallel Execution Engine:** Massive speed improvements for high-throughput batch runs. Both sequence extraction (`extract_fasta.py`) and ChimeraX structural analysis now utilize an orchestrator/worker multiprocessing architecture. You can control this via the `-c` or `--cores` flag.
2. **AlphaFold 3 Smart Filtering:** The pipeline now automatically detects AF3 nested seeds (e.g., `seed-m_sample-n/model.cif`) and safely excludes redundant top-level `*_model.cif` summary files, preventing duplicate processing and skewed datasets.
3. **Expanded Catalytic & Ligand Tracking:** Explicit distance tracking for ATP and Magnesium coordination. The pipeline now calculates critical distances for HRD-Asp to ATP, DFG-Asp to Mg/ATP, and P-Loop to ATP, alongside the existing hydrophobic core and αC-β4 loop metrics.
4. **Universal Header Normalization:** Improved regex parsing ensures bullet-proof mapping between highly variable 3D CIF folder paths and 1D FASTA sequence headers.

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
*(This installs `biopython`, `pandas`, `pyyaml`, and `matplotlib`. The new multiprocessing features utilize standard built-in Python libraries.)*

### 3. The Two-Stage Workflow

#### Stage A: Discovery & Audit
Scan your current directory for simulation folders:
```bash
python3 run_kinase_pipeline_v6.py
```
* The script safely ignores non-simulation folders (like `modules/` or `archives/`).
* It strictly reads your folder names (e.g., `a-erbb2cattail_b-egfrcattail-t790m`) and generates a `generated_matrix.csv` manifest.
* **Select Option 2** when prompted to pause the pipeline. Open the CSV to verify your simulation matrix. You can easily drop unwanted runs by deleting the "x" in the corresponding cell.

#### Stage B: Parallel Analysis & Visualization
Resume the pipeline after verification, specifying the number of CPU cores to utilize:
```bash
python3 run_kinase_pipeline_v6.py --resume -c 16
```
The pipeline automatically handles:
1. **Archiving**: Moves any previous results to a timestamped folder in `archives/` to prevent data mixing (safely bypasses active input files).
2. **Setup**: Generates `proteins.yaml` and extracts `sequences.fasta` via parallel workers.
3. **Landmark Mapping**: Aligns sequences against the HMM profile to generate `hmm_landmarks.json` (extracting canonical motifs, shell residues, and allosteric loop anchors).
4. **1D Sequence Visualizations**: Generates `MSA_Annotated_Panel.pdf` and individual topological schematics (`kinome_VISalign.py`).
5. **3D Structural Execution**: Spawns parallel headless ChimeraX instances, compiles a merged statistical dataset, and dynamically generates split visualization macros.

### 4. Advanced: Custom Entry Points
If you already have a curated configuration or sequence file, you can bypass the discovery stages. The pipeline will automatically protect your provided files from being archived.

* **Start from YAML (--use-yaml)**:
```bash
python3 run_kinase_pipeline_v6.py --use-yaml -c 8
```
*(Skips directory scanning; uses your existing `proteins.yaml` to extract FASTA sequences in parallel, then proceeds to analysis.)*

* **Start from FASTA (--use-fasta)**:
```bash
python3 run_kinase_pipeline_v6.py --use-fasta -c 8
```
*(Skips scanning and extraction entirely; uses your existing `sequences.fasta` for HMM alignment, PDF schematic generation, and ChimeraX structural analysis.)*

## Utilities & Troubleshooting
* **`debug_landmarks.py`**: If alignment visual markers or structural distance calculations look incorrect, you can run `python3 modules/debug_landmarks.py sequences.fasta` to print a terminal table verifying the exact amino acid letters matching your extracted HMM coordinate indices.

## Outputs
* `hmm_kinase_analysis_results_v6.csv`: Comprehensive structural measurement dataset for all analyzed chains. Includes spatial/dihedral states, Spine/Shell integrity, Allosteric nodes (`Y156_N99_Dist`, `K105_E107_Dist`, etc.), and Ligand metrics (`HRD_ATP_Dist`, `DFG_Mg_Dist`, `PLoop_ATP_Dist`). This file contains the source data for the *analyzeR* pipeline (find the pipeline scripts and README files there for details).
* `MSA_Annotated_Panel.pdf`: Wrapped, monospaced sequence alignment visually mapped to canonical landmarks and structural topology tracks.
* `[sequence_name]_schematic.pdf`: Individual 1D abstract topologies showing N/C-lobes, helices, sheets, and catalytic loop positioning.
* **Visualization Directories:**
  * `cx_viz_core/`: `.cxc` ChimeraX macros isolating the catalytic core, R/C-Spines, Salt Bridge, and Hydrophobic Shell.
  * `cx_viz_allosteric/`: `.cxc` ChimeraX macros isolating the αC-β4 toggle switch, αE anchor, and αF scaffold.
* `archives/`: Historical record of previous runs.

## Project Structure
To maintain a clean working environment, all core logic is isolated in the `modules/` directory. Your working directory should look like this before and after running the pipeline:

```text
. (Working Directory)
├── run_kinase_pipeline_v6.py        # The main entry point orchestrator
├── requirements.txt                 # Python package dependencies
├── modules/                         # Core logic directory
│   ├── generate_config.py
│   ├── extract_fasta.py             # Parallel FASTA extraction
│   ├── extract_landmarks.py
│   ├── kinome_VISalign.py           # 1D Schematics and MSA plotting
│   ├── debug_landmarks.py           # Verification utility for HMM mapping
│   ├── 1_run_parallel_chimerax_hmm_v6.py
│   └── chimerax_hmm_worker_v6.py    # Headless calculation & logic
├── a-erbb2cattail_b-egfrcattail.../ # Your input simulation directories (.cif/.pdb)
├── abl1-wtcat_py161_apo/            # Your input simulation directories (.cif/.pdb)
│
│   # --- Pipeline Generated Outputs (Created during execution) ---
├── archives/                        # Historical pipeline runs (auto-archived)
│   └── run_YYYYMMDD_HHMMSS/
├── generated_matrix.csv             # Matrix manifest of discovered simulations
├── proteins.yaml                    # Strict naming configuration map
├── sequences.fasta                  # Extracted sequence FASTA
├── hmm_landmarks.json               # HMM mapped canonical and allosteric residues
├── hmm_kinase_analysis_results_v6.csv # Final structural measurements dataset
├── MSA_Annotated_Panel.pdf          # Unified MSA structural annotation graphic
├── [kinase]_schematic.pdf           # Individual 1D topology maps
├── cx_viz_core/                     # 3D macros for the catalytic core (.cxc)
└── cx_viz_allosteric/               # 3D macros for the allosteric network (.cxc)
```

---

## References
1. **Kinase Classification:** Modi, V., & Dunbrack, R. L., Jr (2019). "Defining a new nomenclature for the structures of active and inactive kinases." *PNAS*, 116(14), 6818-6827.
2. **Hydrophobic Core:** Kim, J. et al. (2017) "A dynamic hydrophobic core orchestrates allostery in protein kinases," *Science Advances*, 3(4). doi: 10.1126/sciadv.1600663.
3. **DFG-in/out Conformational Coupling:** Levinson, N. M. et al. (2006). "A Src-like inactive conformation in the abl tyrosine kinase domain." *PLoS Biology*, 4(5), e144.
4. **αC-β4 Loop Allostery:** Wu, J., Jonniya, N. A., et al. (2024). "Role of the aC-b4 loop in protein kinase structure and dynamics." *eLife*. doi: 10.7554/eLife.91980
