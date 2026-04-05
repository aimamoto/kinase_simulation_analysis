# Kinase Structural Bioinformatic Pipeline

This pipeline automates the extraction, sequence alignment, and structural analysis of kinase protein simulations (AlphaFold3 or MD). It is optimized for batch processing large sets of simulations, such as combinatorial mutant matrices (e.g., 4x8 EGFR/ERBB variants) and single-protein runs (e.g., ABL1 wild-type).

## Project Goal
This pipeline is developed as a sister project of another pipeline specific for erbb family dimers (github.com/aimamoto/erbb_dimer_complex). While the sister project relies on the landmark residues in the ERBB family, this project aims to extend the analysis schemes to a broad range of tyrosine kinases using hmmer. So far, tested in the ERBB family members, SRC, CSK, ABL, FAK, JAK2, and PKA catalytic subunit.

**Recent Upgrades (v4r1):** 
1. **Advanced Structural Metrics:** The pipeline now evaluates both the dynamic hydrophobic core and the <span class="math inline">&alpha;</span>C-<span class="math inline">&beta;</span>4 loop allosteric network. 
   * *Hydrophobic Core:* Tracks the cohesiveness of the Hydrophobic Shell and the integrity of critical Core Bridges (V104 to R-Spine, I150 to <span class="math inline">&alpha;</span>E helix).
   * *Allosteric Networks:* Quantifies the K105 toggle switch (sensing active vs. apo states), the <span class="math inline">&alpha;</span>E helix anchor (Y156-N99), and the rigid deep <span class="math inline">&alpha;</span>F-helix scaffold (D220).
2. **Visual Alignments (1D):** Automatically generates publication-ready Multiple Sequence Alignment (MSA) text panels and abstract 1D topological schematics mapping structural motifs perfectly to sequence coordinates.

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
*(This installs `biopython`, `pandas`, `pyyaml`, and `matplotlib`)*

### 3. The Two-Stage Workflow

#### Stage A: Discovery & Audit
Scan your current directory for simulation folders:
```bash
python3 run_kinase_pipeline_v4r1.py
```
* The script safely ignores non-simulation folders (like `modules/` or `archives/`).
* It strictly reads your folder names (e.g., `a-erbb2cattail_b-egfrcattail-t790m`) and generates a `generated_matrix.csv` manifest.
* **Select Option 2** when prompted to pause the pipeline. Open the CSV to verify your simulation matrix. You can easily drop unwanted runs by deleting the "x" in the corresponding cell.

#### Stage B: Analysis & Visualization
Resume the pipeline after verification:
```bash
python3 run_kinase_pipeline_v4r1.py --resume
```
The pipeline automatically handles:
1. **Archiving**: Moves any previous results to a timestamped folder in `archives/` to prevent data mixing (safely bypasses active input files).
2. **Setup**: Generates `proteins.yaml` and strictly extracts `sequences.fasta` without forcefully altering your folder names.
3. **Landmark Mapping**: Aligns sequences against the HMM profile to generate `hmm_landmarks.json` (extracting canonical motifs, shell residues, and allosteric loop anchors).
4. **1D Sequence Visualizations**: Generates `MSA_Annotated_Panel.pdf` and individual topological schematics (`kinome_VISalign.py`).
5. **3D Structural Execution**: Spawns parallel ChimeraX workers for high-throughput structural analysis, dynamically generating split visualization macros.

### 4. Advanced: Custom Entry Points
If you already have a curated configuration or sequence file, you can bypass the discovery stages. The pipeline will automatically protect your provided files from being archived.

* **Start from YAML (--use-yaml)**:
```bash
python3 run_kinase_pipeline_v4r1.py --use-yaml
```
*(Skips directory scanning; uses your existing proteins.yaml to extract FASTA sequences, then proceeds to analysis.)*

* **Start from FASTA (--use-fasta)**:
```bash
python3 run_kinase_pipeline_v4r1.py --use-fasta
```
*(Skips scanning and extraction entirely; uses your existing sequences.fasta for HMM alignment, PDF schematic generation, and ChimeraX structural analysis.)*

## Utilities & Troubleshooting
* **`debug_landmarks.py`**: If alignment visual markers or structural distance calculations look incorrect, you can run `python3 modules/debug_landmarks.py sequences.fasta` to print a terminal table verifying the exact amino acid letters matching your extracted HMM coordinate indices.

## Outputs
* `hmm_aloop_analysis_results.csv`: Comprehensive structural measurement dataset for all analyzed chains. Includes spatial/dihedral states, Spine/Shell integrity, and explicit distance metrics for allosteric nodes (`Y156_N99_Dist`, `K105_E107_Dist`, `D220_HRD_Dist`, etc.).
* `MSA_Annotated_Panel.pdf`: Wrapped, monospaced sequence alignment visually mapped to canonical landmarks and structural topology tracks.
* `[sequence_name]_schematic.pdf`: Individual 1D abstract topologies showing N/C-lobes, helices, sheets, and catalytic loop positioning.
* **Split Visualization Directories:**
  * `cx_viz_core/`: `.cxc` ChimeraX macros isolating the catalytic core, R/C-Spines, Salt Bridge, and Hydrophobic Shell.
  * `cx_viz_allosteric/`: `.cxc` ChimeraX macros isolating the <span class="math inline">&alpha;</span>C-<span class="math inline">&beta;</span>4 toggle switch, <span class="math inline">&alpha;</span>E anchor, and <span class="math inline">&alpha;</span>F scaffold.
* `archives/`: Historical record of previous runs.

## Project Structure
To maintain a clean working environment, all core logic is isolated in the `modules/` directory. Your working directory should look like this before and after running the pipeline:

```text
. (Working Directory)
├── run_kinase_pipeline_v4r1.py      # The main entry point orchestrator
├── requirements.txt                 # Python package dependencies
├── modules/                         # Core logic directory
│   ├── generate_config.py
│   ├── extract_fasta.py
│   ├── extract_landmarks.py
│   ├── kinome_VISalign.py           # 1D Schematics and MSA plotting
│   ├── debug_landmarks.py           # Verification utility for HMM mapping
│   ├── 1_run_parallel_chimerax_hmm_v4r1.py
│   └── chimerax_hmm_worker_v4r1.py  # Calculates core domains and allosteric networks
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
├── hmm_aloop_analysis_results.csv   # Final structural measurements dataset
├── MSA_Annotated_Panel.pdf          # Unified MSA structural annotation graphic
├── [kinase]_schematic.pdf           # Individual 1D topology maps
├── cx_viz_core/                     # 3D visualization macros for the catalytic core (.cxc)
└── cx_viz_allosteric/               # 3D visualization macros for the allosteric network (.cxc)
```
---

## References
1. **Kinase Classification:** Modi, V., & Dunbrack, R. L., Jr (2019). "Defining a new nomenclature for the structures of active and inactive kinases." *PNAS*, 116(14), 6818-6827.
2. **Hydrophobic Core:** Kim, J. et al. (2017) "A dynamic hydrophobic core orchestrates allostery in protein kinases," *Science Advances*, 3(4). doi: 10.1126/sciadv.1600663.
3. **DFG-in/out Conformational Coupling:** Levinson, N. M. et al. (2006). "A Src-like inactive conformation in the abl tyrosine kinase domain." *PLoS Biology*, 4(5), e144.
4. **<span class="math inline">&alpha;</span>C-<span class="math inline">&beta;</span>4 Loop Allostery:** Wu, J., Jonniya, N. A., et al. (2024). "Role of the aC-b4 loop in protein kinase structure and dynamics." *eLife*. doi: 10.7554/eLife.91980
