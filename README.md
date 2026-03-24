# Kinase Structural Bioinformatic Pipeline

This pipeline automates the extraction, sequence alignment, and structural analysis of kinase protein simulations (AlphaFold3 or MD). It is optimized for batch processing large sets of simulations, such as combinatorial mutant matrices (e.g., 4x8 EGFR/ERBB variants) and single-protein runs (e.g., ABL1 wild-type).

## Project Goal
This pipeline is developed as a sister project of another pipeline specific for erbb family dimers (github.com/aimamoto/erbb_dimer_complex). While the sister project relies on the landmark residues in the ERBB family, this project aims to extend the analysis schemes to broad range of tyrosine kinases using hmmer. Sofar, tested in the ERBB family members, SRC, CSK, and ABL.

## Getting Started

### 1. System Requirements
Ensure your OS (Ubuntu, macOS, etc.) has the following installed and accessible via your system PATH:
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
*(This installs `biopython`, `pandas`, and `pyyaml`)*

### 3. The Two-Stage Workflow

#### Stage A: Discovery & Audit
Scan your current directory for simulation folders:
```bash
python3 run_kinase_pipeline.py
```
* The script safely ignores non-simulation folders (like `modules/` or `archives/`).
* It strictly reads your folder names (e.g., `a-erbb2cattail_b-egfrcattail-t790m`) and generates a `generated_matrix.csv` manifest.
* **Select Option 2** when prompted to pause the pipeline. Open the CSV to verify your simulation matrix. You can easily drop unwanted runs by deleting the "x" in the corresponding cell.

#### Stage B: Analysis & Visualization
Resume the pipeline after verification:
```bash
python3 run_kinase_pipeline.py --resume
```
The pipeline automatically handles:
1. **Archiving**: Moves any previous results to a timestamped folder in `archives/` to prevent data mixing.
2. **Setup**: Generates `proteins.yaml` and strictly extracts `sequences.fasta` without forcefully altering your folder names.
3. **Mapping**: Aligns sequences against the HMM profile to generate `hmm_landmarks.json`.
4. **Execution**: Spawns parallel ChimeraX workers for high-throughput structural analysis of salt bridges, R- and C-spines, and A-loop dynamics.

5. **Advanced**: Custom Entry Points
If you already have a curated configuration or sequence file, you can bypass the discovery stages. The pipeline will automatically protect your provided files from being archived.

* **Start from YAML (--use-yaml)**:

```bash
python3 run_kinase_pipeline.py --use-yaml
```
(Skips directory scanning; uses your existing proteins.yaml to extract FASTA sequences, then proceeds to analysis.)

* **Start from FASTA (--use-fasta)*:

```bash
python3 run_kinase_pipeline.py --use-fasta
```
(Skips scanning and extraction entirely; uses your existing sequences.fasta for HMM alignment and ChimeraX structural analysis.)

## Outputs
* `hmm_aloop_analysis_results.csv`: Comprehensive structural measurement dataset for all analyzed chains.
* `cx_pocket_visualization_hmm/`: `.cxc` ChimeraX macros for 3D visualization of each simulation.
* `archives/`: Historical record of previous runs.

## Project Structure
To maintain a clean working environment, all core logic is isolated in the `modules/` directory. Your working directory should look like this before and after running the pipeline:

```
. (Working Directory)
├── run_kinase_pipeline.py           # The main entry point orchestrator
├── requirements.txt                 # Python package dependencies
├── modules/                         # Core logic directory
│   ├── generate_config.py
│   ├── extract_fasta.py
│   ├── extract_landmarks.py
│   ├── 1_run_parallel_chimerax_hmm_v2.py
│   └── chimerax_hmm_worker_v2upgrade.py
├── a-erbb2cattail_b-egfrcattail.../ # Your input simulation directories (.cif/.pdb)
├── abl1-wtcat_py161_apo/            # Your input simulation directories (.cif/.pdb)
│
│   # --- Pipeline Generated Outputs (Created during execution) ---
├── archives/                        # Historical pipeline runs (auto-archived)
│   └── run_YYYYMMDD_HHMMSS/
├── generated_matrix.csv             # Matrix manifest of discovered simulations
├── proteins.yaml                    # Strict naming configuration map
├── sequences.fasta                  # Extracted sequence FASTA
├── hmm_landmarks.json               # HMM mapped canonical residues
├── hmm_aloop_analysis_results.csv   # Final structural measurements dataset
└── cx_pocket_visualization_hmm/     # 3D visualization macros (.cxc)

```
