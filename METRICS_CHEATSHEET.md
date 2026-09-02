[METRICS_CHEATSHEET.md](https://github.com/user-attachments/files/26516007/METRICS_CHEATSHEET.md)
# Kinase Structural Metrics Cheat Sheet

This guide maps the raw measurements in `master_kinase_analysis_results_v7r3.csv` to their biological meaning, their visual representation in the R statistics pipeline (`scripts/modules/`), and where to find them in your 3D ChimeraX sessions. 

## 0. Read this first: what the landmark names mean

Every numbered landmark in the column names (`N99`, `V104`, `K105`, `E107`, `M118`, `M120`, `E121`, `I150`, `Y156`, `D220`) is a **position label, not a residue identity**. Each is a fixed node of the Pfam Pkinase HMM profile, named after the residue that occupies it in **PKA**, the canonical reference kinase. The pipeline aligns every target to that profile (`modules/extract_landmarks.py`, `TARGET_NODES`), so one column means the same structural position across the kinome regardless of what residue sits there or how the target is numbered.

The consequence: **the letter in the column name usually does not describe the residue actually measured.** For the three model kinases in this study, six of the ten numbered landmarks are occupied by a residue of a different type in all three:

| Column label | HMM node | CSK | SRC | CDK1 | Letter holds? |
|---|---|---|---|---|---|
| `N99`  | 56  | **R**59 | **R**61 | **R**59 | No — Arg, never Asn |
| `V104` | 61  | V64 | V66 | V64 | Yes |
| `K105` | 62  | **Q**65 | **Q**67 | **S**65 | No — Gln/Ser, never Lys |
| `E107` | 64  | **L**67 | **Y**69 | **Q**67 | No — Leu/Tyr/Gln, never Glu |
| `M118` | 75  | **I**79 | **I**79 | **L**78 | No — Ile/Leu, never Met |
| `M120` | 77  | **T**81 | **T**81 | **F**80 | No — Thr/Phe, never Met |
| `E121` | 78  | E82 | E82 | E81 | Yes — Glu in all three |
| `I150` | 113 | **Y**119 | **Y**119 | **F**118 | No — Tyr/Phe, never Ile |
| `Y156` | 119 | **F**125 | Y125 | **V**124 | Partial — Tyr in SRC only |
| `D220` | 179 | D183 | D187 | D186 | Yes |

Positions are given in the numbering of the deposited construct (`data/*/inputs/sequences.fasta`). CDK1 is full-length, so its numbers are UniProt numbering directly. For the CSK and SRC constructs, add **+185** and **+257** respectively to recover canonical full-length numbering — e.g. CSK `K` → K222, SRC `K` → K295, SRC DFG-Phe → F405. Those offsets are constant across four independent landmarks and agree with Kincore's own residue calls on 3D7T (`data/validation_3d7t/`).

**So when reading the entries below, take the biology as the description of a structural position and its role, not as a claim about side-chain chemistry.** Where the PKA-derived name implies an interaction that cannot occur in these kinases, the entry says so.

## 1. The Catalytic Core & Hydrophobic Shell
**3D Visual Location:** `cx_viz_core/` macros.  
These metrics define the fundamental "engine" of the kinase, focusing on the hydrophobic architecture that dynamically compresses to link the N-lobe and C-lobe during activation *(Kim, 2017)*.

* **`SB_Dist` (The Salt Bridge)**
  * **Biology:** Minimum distance from the β3-Lys (VAIK) **Nζ** to the nearest polar (O\*/N\*) side-chain atom of the αC-Glu. Short values (~2.5-3 Å intact; < 4 Å broadly bridged) mean the αC is swung "In" and bridging ATP to the catalytic machinery. Both partners are as named in every kinase here — this is one of the few genuinely conserved positions (β3-Lys and αC-Glu). Note the categorical `C_Helix` column is **not** derived from this value; see its entry below. 
  * **Plot Reference:** Phase 4 (Violin Plots), Phase 5 (PCA/Heatmaps).
  * **Citation:** (Huse & Kuriyan, 2002) [1].
* **`V104_RS2_Dist` (Core Compression)**
  * **Biology:** Minimum side-chain distance between the deep Hydrophobic Shell (`V104`) and the Regulatory Spine (`RS2`, the β4 spine residue). Measures the "tightness" of the core. In active kinases, the shell compresses tightly against the R-spine. This is one of the labels whose letter does hold: the `V104` position is a valine in all three model kinases (Val64/Val66/Val64).
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kim, 2017) [2].
* **`Shell_M118_M120_Dist` (Hinge Breathing)**
  * **Biology:** Measures the local flexibility/compaction of the hinge region connecting the two lobes, as the minimum side-chain distance between the two shell positions (drives `Shell_State`: `Packed` if < 5 Å, else `Loose`). Despite the name, **neither position is a methionine in the model kinases** — they are Ile79/Thr81 in both CSK and SRC, and Leu78/Phe80 in CDK1. The packing readout is a hydrophobic-contact distance, not a Met-Met interaction.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kim, 2017) [2].
* **`I150_HRD_Dist` (Base Rigidity)**
  * **Biology:** Minimum side-chain distance between the αE helix position (`I150`) and the catalytic **HRD-His**. The αE helix acts as the stable "floor" of the kinase, and this ensures the catalytic loop is firmly mounted to it. The `I150` position is an aromatic, not an isoleucine, in the model kinases — Tyr119 in both CSK and SRC, Phe118 in CDK1 — so the contact is aromatic packing against the catalytic histidine.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kim, 2017) [2].

## 2. The Allosteric Network & Regulatory Dials
**3D Visual Location:** `cx_viz_allosteric/` macros.  
These networks sense external signals (ligand binding, mutations) and transmit them to the catalytic core via the highly dynamic αC-β4 loop *(Wu, 2024)*.

* **`K105_E107_Dist` & `K105_E121_Dist` (The Toggle Switch)**
  * **Biology:** The `K105` position sits on the flexible αC-β4 loop and acts as a mechanical switch, reorienting between an inward partner (`E107`) and an outward one (`E121`) as the loop repacks. Both columns are minimum side-chain distances. **The PKA-derived Lys-Glu chemistry does not transfer to these kinases:** the `K105` position is Gln65 (CSK), Gln67 (SRC) and Ser65 (CDK1), and the `E107` position is Leu67 (CSK), Tyr69 (SRC) and Gln67 (CDK1) — so `K105_E107_Dist` is not a salt bridge in any model kinase, and in CSK it is a polar-to-aliphatic contact with no possible hydrogen bond. `K105_E121_Dist` is the better-behaved of the two, since the `E121` position genuinely is a glutamate in all three (Glu82/Glu82/Glu81), making a Gln/Ser-Glu hydrogen bond chemically available. Read both as loop-conformation distances that report the switch geometry, not as bond-formation readouts. Neither is the β3-Lys/αC-Glu regulatory salt bridge — that is `SB_Dist`.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Wu, 2024) [3].
* **`K105_N99_Dist` (Loop Cohesiveness)**
  * **Biology:** Minimum side-chain distance between the `K105` switch position and the top of the αC-β4 loop (`N99`). Defines how tightly folded the regulatory loop is. The `N99` position is an **arginine** in all three model kinases (Arg59 in CSK and CDK1, Arg61 in SRC), not the asparagine its label implies.
  * **Plot Reference:** Phase 4 (Violin Plots), Phase 5 (Heatmaps).
  * **Citation:** (Wu, 2024) [3].
* **`Y156_N99_Dist` (The αE Anchor)**
  * **Biology:** A massive cross-lobe bridge. It anchors the flexible regulatory loop of the N-lobe (`N99`) down to the rigid C-lobe (`Y156`). Breaking this distance functionally decouples the two halves of the kinase. Note this is the one **all-atom** ("residue-min") distance in this group — every other entry here is side-chain-only — so backbone atoms can set the value. In the model kinases the pair is Phe125-Arg59 (CSK), Tyr125-Arg61 (SRC) and Val124-Arg59 (CDK1): a cation-π or hydrophobic contact against an arginine, and only in SRC is the `Y156` position the tyrosine its label implies.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Wu, 2024) [3].
* **`D220_HRD_Dist` (The αF Scaffold)**
  * **Biology:** Minimum **side-chain** distance from the highly conserved Aspartate on the deep αF helix (`D220`) to the catalytic **HRD-His** (not to the catalytic loop backbone). The αF helix acts as the ultimate rigid "hub" of the kinase, anchoring the catalytic loop in place for structural stability. Both partners are as named in the model kinases — Asp183 (CSK), Asp187 (SRC), Asp186 (CDK1) against the HRD histidine — so this entry's chemistry transfers directly.
  * **Plot Reference:** Phase 4 (Violin Plots).
  * **Citation:** (Kornev et al., 2008; Taylor & Kornev, 2011)[4, 5].

## 3. DFG Conformation & Phase Space
**3D Visual Location:** Present in both Core and Allosteric macros (focus on the Activation Loop).  
These metrics strictly define the spatial positioning and backbone torsion of the DFG Motif (Asp-Phe-Gly), dictating ATP-Mg2+ coordination *(Modi, 2019; Levinson, 2006)*.

* **`D1_Dist` & `D2_Dist` (Dunbrack Spatial Coordinates)**
  * **Biology:** Both coordinates measure to the same moving atom — the DFG-Phe Cζ (or the farthest side-chain heavy atom from Cα, if Cζ is absent as in a mutated Phe) — and differ only in the fixed anchor they measure it from. `D1` (αC-Glu(+4) Cα → DFG-Phe Cζ) separates DFG-*in* (short, ≤ 11 Å) from DFG-*out* (long, > 11 Å). `D2` (β3-Lys/VAIK-Lys Cα → DFG-Phe Cζ) then distinguishes DFG-*in* (≥ 11 Å) from the DFG-*inter* intermediate (≤ 11 Å). The D1 anchor is the fourth residue *after* the αC-Glu, not the αC-Glu itself: the published 11/14 Å cutoffs are calibrated on the (+4) anchor, and measuring from the αC-Glu was the v7r2 error corrected at v7r3.
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
  * **Biology:** An αC rotamer call: `In` if the β3-Lys **Cβ** to αC-Glu **Cβ** distance is ≤ 10 Å, else `Out`. This is a *separate* measurement from `SB_Dist`, **not a threshold applied to it** — it asks whether the αC is positioned so that a salt bridge is possible ("salt-bridge competent"), which is a looser condition than the bridge being formed. In practice the two agree closely (across the 2,200 deposited chains, `C_Helix` matches a naive `SB_Dist ≤ 4 Å` cut on 99.3% of rows), so the distinction rarely changes a call — but they are not the same criterion, and the 15 rows where they disagree do so in both directions.
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
If you compare the distances in the `master_kinase_analysis_results_v7r3.csv` to the dashed lines drawn in the ChimeraX `.cxc` 3D visuals, you will notice that some values differ. **This is intentional.**

1. **The CSV (Quantitative Data):** Calculates the **minimum heavy-atom sidechain distance** (e.g., Nitrogen to Oxygen). This represents the true stereochemical contact/bond distance used for all rigorous statistical plots.
2. **The 3D Visuals (Qualitative Guide):** The macros use two drawing conventions depending on the interaction. For elements tracked as macroscopic backbone shifts (the R-Spine, the Hydrophobic Shell, and Cα landmarks such as I150–HRD and Y156–N99), lines are drawn between **Alpha-Carbons (Cα)** to track the backbone without drawing messy lines that clip through the protein ribbons; these Cα values differ from the CSV. For specific chemical contacts, the macros instead draw the **exact heavy-atom (atom-to-atom) shortest distance**, which reproduces the CSV value: this includes the β3-Lys–αC-Glu salt bridge (`SB_Dist`, drawn Nζ→nearest carboxylate oxygen), the αC-β4/αE polar contact, and catalytic coordination such as HRD–ATP and ATP-Mg2+.

---

## References

1. **Huse, M., & Kuriyan, J. (2002).** "The conformational plasticity of protein kinases." *Cell*, 109(3), 275-282. 
2. **Kim, J. et al. (2017).** "A dynamic hydrophobic core orchestrates allostery in protein kinases." *Science Advances*, 3(4), e1600663. 
3. **Wu, J., Jonniya, N. A., et al. (2024).** "Role of the aC-b4 loop in protein kinase structure and dynamics." *eLife*, 13, e91980.
4. **Kornev, A. P., Taylor, S. S., & Ten Eyck, L. F. (2008).** "A helix scaffold for the assembly of active protein kinases." *PNAS*, 105(38), 14377-14382.
5. **Taylor, S. S., & Kornev, A. P. (2011).** "Protein kinases: evolution of dynamic regulatory proteins." *Trends in Biochemical Sciences*, 36(2), 65-77.
6. **Modi, V., & Dunbrack, R. L., Jr. (2019).** "Defining a new nomenclature for the structures of active and inactive kinases." *PNAS*, 116(14), 6818-6827.
7. **Levinson, N. M. et al. (2006).** "A Src-like inactive conformation in the abl tyrosine kinase domain." *PLoS Biology*, 4(5), e144.
