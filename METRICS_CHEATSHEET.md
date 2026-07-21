[METRICS_CHEATSHEET.md](https://github.com/user-attachments/files/26516007/METRICS_CHEATSHEET.md)
# Kinase Structural Metrics Cheat Sheet

This guide maps the raw measurements in `master_kinase_analysis_results_v7r2.csv` to their biological meaning, their visual representation in the R statistics pipeline (`scripts/modules/`), and where to find them in your 3D ChimeraX sessions. 

## 1. The Catalytic Core & Hydrophobic Shell
**3D Visual Location:** `cx_viz_core/` macros.  
These metrics define the fundamental "engine" of the kinase, focusing on the hydrophobic architecture that dynamically compresses to link the N-lobe and C-lobe during activation *(Kim, 2017)*.

* **`SB_Dist` (The Salt Bridge)**
  * **Biology:** Distance between the β3-Lysine (VAIK) and the αC-Glutamate. Short distances (< 4.0 Å) indicate the αC-helix is swung "In", bridging ATP to the catalytic machinery. 
  * **Plot Reference:** Phase 4 (Violin Plots), Phase 5 (PCA/Heatmaps).
  * **Citation:** (Huse & Kuriyan, 2002) [1].
* **`V104_RS2_Dist` (Core Compression)**
  * **Biology:** Distance between the deep Hydrophobic Shell (V104) and the Regulatory Spine (RS2). Measures the "tightness" of the core. In active kinases, the shell compresses tightly against the R-spine.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kim, 2017) [2].
* **`Shell_M118_M120_Dist` (Hinge Breathing)**
  * **Biology:** Measures the local flexibility/compaction of the hinge region connecting the two lobes.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kim, 2017) [2].
* **`I150_HRD_Dist` (Base Rigidity)**
  * **Biology:** Distance between the αE helix (I150) and the catalytic HRD motif. The αE helix acts as the stable "floor" of the kinase. This ensures the catalytic loop is firmly mounted to the floor.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kim, 2017) [2].

## 2. The Allosteric Network & Regulatory Dials
**3D Visual Location:** `cx_viz_allosteric/` macros.  
These networks sense external signals (ligand binding, mutations) and transmit them to the catalytic core via the highly dynamic αC-β4 loop *(Wu, 2024)*.

* **`K105_E107_Dist` & `K105_E121_Dist` (The Toggle Switch)**
  * **Biology:** K105 sits on the flexible αC-β4 loop and acts as a mechanical switch. When active, it flips inward to pair with E107. When inactive, it flips outward to pair with E121 or face the solvent.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Wu, 2024) [3].
* **`K105_N99_Dist` (Loop Cohesiveness)**
  * **Biology:** Distance between the K105 switch and the top of the αC-β4 loop (N99). Defines how tightly folded the regulatory loop is.
  * **Plot Reference:** Phase 4 (Violin Plots), Phase 5 (Heatmaps).
  * **Citation:** (Wu, 2024) [3].
* **`Y156_N99_Dist` (The αE Anchor)**
  * **Biology:** A massive cross-lobe bridge. It anchors the flexible regulatory loop of the N-lobe (N99) down to the rigid C-lobe (Y156). Breaking this distance functionally decouples the two halves of the kinase.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Wu, 2024) [3].
* **`D220_HRD_Dist` (The αF Scaffold)**
  * **Biology:** Distance from the highly conserved Aspartate on the deep αF helix (D220) to the backbone of the catalytic loop (HRD). The αF helix acts as the ultimate rigid "hub" of the kinase, anchoring the catalytic loop in place for structural stability.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kornev et al., 2008; Taylor & Kornev, 2011)[4, 5].

## 3. DFG Conformation & Phase Space
**3D Visual Location:** Present in both Core and Allosteric macros (focus on the Activation Loop).  
These metrics strictly define the spatial positioning and backbone torsion of the DFG Motif (Asp-Phe-Gly), dictating ATP-Mg2+ coordination *(Modi, 2019; Levinson, 2006)*.

* **`D1_Dist` & `D2_Dist` (Dunbrack Spatial Coordinates)**
  * **Biology:** `D1` (DFG-Asp to αC-Glu) separates DFG-*in* (short) from DFG-*out* (long). `D2` (DFG-Asp to HRD-His/Tyr) separates active states from intermediates.
  * **Plot Reference:** Phase 3 (2D Phase Space, lower panels).
  * **Citation:** (Modi, 2019) [6].
* **`Phi_D` & `Psi_D` (DFG Dihedrals)**
  * **Biology:** The Ramachandran backbone torsion angles of the DFG-Aspartate. Maps the exact rotation the backbone undergoes when transitioning between states.
  * **Plot Reference:** Phase 3 (2D Phase Space, upper panels).
  * **Citation:** (Modi, 2019) [6].

## 4. Categorical Macro-States
**3D Visual Location:** Global structure visualization.  
These are discrete, rule-based classifications derived from the continuous metrics above, providing a shorthand for the global state of the kinase *(Modi, 2019)*.

* **`State` (Global Conformation)**
  * **Biology:** Classifies the overall architecture (e.g., "Active (BLAminus)", "Inactive (BLAplus)", "DFGout").
  * **Plot Reference:** Phase 2 (Macro-State Bar Charts).
  * **Citation:** (Modi, 2019 [6]; Levinson, 2006 [7]).
* **`C_Helix` (In vs Out)**
  * **Biology:** Classified strictly based on the integrity of the `SB_Dist` (Salt Bridge).
  * **Plot Reference:** Phase 2 (Macro-State Bar Charts).
  * **Citation:** (Huse & Kuriyan, 2002)[1].
* **`R_Spine` (Intact vs Broken)**
  * **Biology:** Evaluates if the Regulatory Spine residues are stacked contiguously.
  * **Plot Reference:** Phase 2 (Macro-State Bar Charts).
  * **Citation:** (Kim, 2017) [2].
* **`Spatial` (DFGin vs DFGout)**
  * **Biology:** Classified strictly based on `D1_Dist` and `D2_Dist` boundaries.
  * **Plot Reference:** Phase 2 (Macro-State Bar Charts).
  * **Citation:** (Modi, 2019) [6].

## 🛑 Important Disclaimer: CSV vs. 3D Visuals
If you compare the distances in the `hmm_kinase_analysis_results.csv` to the dashed lines drawn in the ChimeraX `.cxc` 3D visuals, you will notice different values. **This is intentional.**

1. **The CSV (Quantitative Data):** Calculates the **minimum heavy-atom sidechain distance** (e.g., Nitrogen to Oxygen). This represents the true stereochemical contact/bond distance used for all rigorous statistical plots.
2. **The 3D Visuals (Qualitative Guide):** For structural features spanning large distances (Allosteric networks, R-Spine, Hydrophobic Shell), the macros draw lines between **Alpha-Carbons (Cα)**. This is a standard visualization technique used to track macroscopic backbone shifts without drawing messy lines that clip through the protein ribbons. (Note: Precise catalytic interactions, such as ATP-Mg2+ coordination, are still shown via precise heavy-atom distances in the 3D view).

---

## References

1. **Huse, M., & Kuriyan, J. (2002).** "The conformational plasticity of protein kinases." *Cell*, 109(3), 275-282. 
2. **Kim, J. et al. (2017).** "A dynamic hydrophobic core orchestrates allostery in protein kinases." *Science Advances*, 3(4), e1600663. 
3. **Wu, J., Jonniya, N. A., et al. (2024).** "Role of the aC-b4 loop in protein kinase structure and dynamics." *eLife*, 13, e91980.
4. **Kornev, A. P., Taylor, S. S., & Ten Eyck, L. F. (2008).** "A helix scaffold for the assembly of active protein kinases." *PNAS*, 105(38), 14377-14382.
5. **Taylor, S. S., & Kornev, A. P. (2011).** "Protein kinases: evolution of dynamic regulatory proteins." *Trends in Biochemical Sciences*, 36(2), 65-77.
6. **Modi, V., & Dunbrack, R. L., Jr. (2019).** "Defining a new nomenclature for the structures of active and inactive kinases." *PNAS*, 116(14), 6818-6827.
7. **Levinson, N. M. et al. (2006).** "A Src-like inactive conformation in the abl tyrosine kinase domain." *PLoS Biology*, 4(5), e144.
