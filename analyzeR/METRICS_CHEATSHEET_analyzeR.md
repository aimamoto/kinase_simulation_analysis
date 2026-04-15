# Kinase Structural Metrics Cheat Sheet

This guide maps the raw measurements and derived classifications in the structural metadata CSVs to their biological meaning, their visual representation in the R statistics pipeline (Dual-Engine architecture), and where to find them in your 3D ChimeraX sessions. 

## 1. The Catalytic Core & Hydrophobic Shell
**3D Visual Location:** `cx_viz_core/` macros.  
These metrics define the fundamental "engine" of the kinase, focusing on the hydrophobic architecture that dynamically compresses to link the N-lobe and C-lobe during activation *(Kim, 2017)*.

* **`SB_Dist` (The Salt Bridge)**
  * **Biology:** Distance between the β3-Lysine of the VAIK motif (or the unique VCIK motif in pseudokinases like ERBB3) and the αC-Glutamate. Short distances (< 4.0 Å) indicate the αC-helix is swung "In", bridging ATP to the catalytic machinery. 
  * **Plot Reference:** Phase 4 (Violin Plots), Phase 5 (PCA/Heatmaps).
  * **Citation:** (Huse & Kuriyan, 2002).
* **`V104_RS2_Dist` (Core Compression)**
  * **Biology:** Distance between the deep Hydrophobic Shell and the Regulatory Spine (RS2). Measures the "tightness" of the core. In active kinases, the shell compresses tightly against the R-spine.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kim, 2017).
* **`Shell_M118_M120_Dist` (Hinge Breathing)**
  * **Biology:** Measures the local flexibility/compaction of the hinge region connecting the two lobes.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kim, 2017).
* **`I150_HRD_Dist` (Base Rigidity)**
  * **Biology:** Distance between the αE helix (e.g., I150, or residue 147 in ERBB3) and the catalytic HRD motif. The αE helix acts as the stable "floor" of the kinase. This ensures the catalytic loop is firmly mounted.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kim, 2017).

## 2. The Allosteric Network & Regulatory Dials
**3D Visual Location:** `cx_viz_allosteric/` macros.  
These networks sense external signals (ligand binding, mutations) and transmit them to the catalytic core via the highly dynamic αC-β4 loop *(Wu, 2024)*.

* **`K105_E107_Dist` & `K105_E121_Dist` (The Toggle Switch)**
  * **Biology:** K105 sits on the flexible αC-β4 loop and acts as a mechanical switch. When active, it flips inward to pair with E107. When inactive, it flips outward to pair with E121 or face the solvent.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Wu, 2024).
* **`K105_N99_Dist` (Loop Cohesiveness)**
  * **Biology:** Distance between the K105 switch and the top of the αC-β4 loop (N99). Defines how tightly folded the regulatory loop is.
  * **Plot Reference:** Phase 4 (Violin Plots), Phase 5 (Heatmaps).
  * **Citation:** (Wu, 2024).
* **`Y156_N99_Dist` (The αE Anchor)**
  * **Biology:** A massive cross-lobe bridge. It anchors the flexible regulatory loop of the N-lobe (N99) down to the rigid C-lobe (Y156). Breaking this distance functionally decouples the two halves of the kinase.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Wu, 2024).
* **`D220_HRD_Dist` (The αF Scaffold)**
  * **Biology:** Distance from the highly conserved Aspartate on the deep αF helix to the backbone of the catalytic loop. The αF helix acts as the ultimate rigid "hub" of the kinase, anchoring the catalytic loop in place for structural stability.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kornev et al., 2008; Taylor & Kornev, 2011).

## 3. DFG Conformation & Phase Space
**3D Visual Location:** Present in both Core and Allosteric macros (focus on the Activation Loop).  
These metrics strictly define the spatial positioning and backbone torsion of the DFG Motif (Asp-Phe-Gly), dictating ATP-Mg2+ coordination *(Modi, 2019; Levinson, 2006)*.

* **`D1_Dist` & `D2_Dist` (Dunbrack Spatial Coordinates)**
  * **Biology:** `D1` (DFG-Asp to αC-Glu) separates DFG-*in* (short) from DFG-*out* (long). `D2` (DFG-Asp to HRD-His/Tyr) separates active states from intermediates.
  * **Plot Reference:** Phase 3 (2D Phase Space, lower panels).
  * **Citation:** (Modi, 2019).
* **`Phi_D` & `Psi_D` (DFG Dihedrals)**
  * **Biology:** The Ramachandran backbone torsion angles of the DFG-Aspartate. Maps the exact rotation the backbone undergoes when transitioning between states.
  * **Plot Reference:** Phase 3 (2D Phase Space, upper panels).
  * **Citation:** (Modi, 2019).

## 4. Asymmetric ERBB & Dimer Architecture (Dual-Engine Expansion)
**3D Visual Location:** Dimer interface mapping and comparative structural overlays.
Metrics specific to the evaluation of asymmetric heterodimers and competitive multi-state complexes.

* **`Role` (Activator vs. Receiver)**
  * **Biology:** Evaluates the conformational identity of subunits in asymmetric interfaces. Calculates receiver probabilities based on global thermodynamic assemblies across the simulation ensemble.
  * **Plot Reference:** Phase 2 (Global Receiver Probability), Phase 6 (State-Role Composition).
* **`HRD_ATP_Dist` (Catalytic Clamping)**
  * **Biology:** Measures the physical clamping distance between the catalytic HRD aspartate (e.g., Asp 813 in mature EGFR / 837 full-length) and the captured ATP ligand. Serves as a primary indicator of 2ATP complex viability.
  * **Plot Reference:** ERBB Engine Phase 4b (2ATP Clamping Violins).
* **`aCb4_aE_Dist` (N/C-Lobe Anchor Shift)**
  * **Biology:** Tracks the long-range structural displacement defining the global shift from an inactive monomeric state into a tightly coupled asymmetric receiver state. 
  * **Plot Reference:** ERBB Engine Phase 4c.

## 5. AlphaFold 3 Ligand QC & Thermodynamics
**3D Visual Location:** Pocket rendering and small-molecule coordination views.
Metrics derived from the pipeline’s built-in "Truth Overwrite" algorithm, which detects instances where AlphaFold's thermodynamic behavior supersedes human-designed experimental parameters.

* **`Sim_Has_Mismatch` (AF3 Truth Overwrite)**
  * **Biology:** A boolean flag triggered when the pipeline detects that AF3 preferentially bound a ligand to an unintended chain (e.g., placing ATP in a designed Apo pocket). Maps competitive binding affinities across dimer pairs.
  * **Plot Reference:** Phase 1 (Ligand Preference Rates).
* **`C_Spine` (Catalytic Spine Integrity)**
  * **Biology:** Now dynamically categorized as "Intact", "Ligand Distant", or "No Ligand". Assesses whether the ligand is properly coordinated by the hydrophobic spine residues.
  * **Plot Reference:** Phase 3 (C-Spine Integrity vs ATP Availability).

## 6. Global Network Density & Differential Drivers
**Pipeline Location:** Phase 5 (Coupling) and Phase 8 (Post-Hoc Differential Analysis).
Quantitative readouts of how allosteric networks rewire in response to mutations or binding events.

* **`Global_Coupling_Score` (MAC)**
  * **Calculation:** The Mean Absolute Correlation (MAC) of all internal dynamic features. 
  * **Biology:** A measure of "Global Network Rigidity." Higher MAC scores indicate a tightly coupled, mechanically rigid kinase (often highly active or trapped in a specific state), while lower scores indicate a more flexible, decoupled state.
  * **Plot Reference:** Phase 5 (Global Network Density Bar Charts).
* **`Signed_EffSize` (Cohen's d Directionality)**
  * **Calculation:** A post-hoc extraction identifying structural drivers between two defined populations.
  * **Biology:** Represents the magnitude and direction (Expansion vs Compression) of specific allosteric distances when comparing variant A to variant B, isolating the structural components responsible for phenotypic shifts.
  * **Plot Reference:** Phase 8 (Differential Volcano Plots).

## 7. Categorical Macro-States
**3D Visual Location:** Global structure visualization.  
Discrete, rule-based classifications derived from the continuous metrics, offering a shorthand for the global state of the kinase *(Modi, 2019)*.

* **`State` (Global Conformation)**
  * **Plot Reference:** Phase 2, Phase 7.
* **`Macro_State` (Unsupervised 3D Discovery)**
  * **Biology:** The cluster assignment (e.g., GMM or K-Means) mapping the simulation to a specific meta-stable energetic basin. Consolidates multi-dimensional shifts into distinct structural phenotypes.
  * **Plot Reference:** Phase 6 (PCA/3D Phase Space).

## 🛑 Important Disclaimer: CSV vs. 3D Visuals
If you compare the distances in the `hmm_kinase_analysis_results.csv` (or the Phase 7 Master Metadata) to the dashed lines drawn in the ChimeraX `.cxc` 3D visuals, you will notice different values. **This is intentional.**

1. **The CSV (Quantitative Data):** Calculates the **minimum heavy-atom sidechain distance** (e.g., Nitrogen to Oxygen). This represents the true stereochemical contact/bond distance used for all rigorous statistical plots.
2. **The 3D Visuals (Qualitative Guide):** For structural features spanning large distances (Allosteric networks, R-Spine, Hydrophobic Shell), the macros draw lines between **Alpha-Carbons (Cα)**. This standard visualization technique is used to track macroscopic backbone shifts without drawing lines that clip through the protein ribbons. (Precise catalytic interactions, such as ATP-Mg2+ coordination, are still shown via exact heavy-atom distances in the 3D view).

---

## References

1. **Huse, M., & Kuriyan, J. (2002).** "The conformational plasticity of protein kinases." *Cell*, 109(3), 275-282. 
2. **Kim, J. et al. (2017).** "A dynamic hydrophobic core orchestrates allostery in protein kinases." *Science Advances*, 3(4), e1600663. 
3. **Wu, J., Jonniya, N. A., et al. (2024).** "Role of the aC-b4 loop in protein kinase structure and dynamics." *eLife*, 13, e91980.
4. **Kornev, A. P., Taylor, S. S., & Ten Eyck, L. F. (2008).** "A helix scaffold for the assembly of active protein kinases." *PNAS*, 105(38), 14377-14382.
5. **Taylor, S. S., & Kornev, A. P. (2011).** "Protein kinases: evolution of dynamic regulatory proteins." *Trends in Biochemical Sciences*, 36(2), 65-77.
6. **Modi, V., & Dunbrack, R. L., Jr. (2019).** "Defining a new nomenclature for the structures of active and inactive kinases." *PNAS*, 116(14), 6818-6827.
7. **Levinson, N. M. et al. (2006).** "A Src-like inactive conformation in the abl tyrosine kinase domain." *PLoS Biology*, 4(5), e144.
