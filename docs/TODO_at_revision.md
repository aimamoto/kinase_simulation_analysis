# Deferred documentation corrections — action at manuscript revision

**Created:** 2026-09-02. **Updated: 2026-09-19 — Items 1, 2 and 5 are ACTIONED in the repo
(staged for the editor, not sent); Item 4, cosmetic and repo-only, is cleared. 2026-09-26 — Item 6
(the `Type` definition) and Item 7 (`README.md` has no v7r4 section) added, both OPEN. 2026-09-27 —
Item 8 (Module 2: `Spine_Bridge_Dist` in the state discovery) added, OPEN.** **Item 8 is not
documentation: it changes reported state-level results. Read it before the others.**

> **Status as of 2026-09-18.** Items 1 and 2 have been applied to
> `docs/AlloQuant_master_CSV_data_dictionary_v7r3.xlsx`/`.pdf`, together with a third
> correction to the `ActLoop_CT` entry found on 2026-09-18 (see below). They were bundled so a
> single corrected S1 Dataset can be supplied at revision instead of piecemeal changes.
>
> **The deposited/submitted files have NOT been changed and no update has been sent to the
> editor.** `plos_submission/` still holds the 2026-09-01 files exactly as submitted. The
> corrected file, an unmodified copy of the submitted one, and a full change note are staged in
> `plos_submission/updates/2026-09-19/` for use at revision.
>
> Verified: exactly three cells differ from the as-submitted workbook (`A9`, `I38`, `I55`).

The manuscript, including **S1 Dataset** (`AlloQuant_master_CSV_data_dictionary_v7r3.xlsx`/`.pdf`),
was submitted to *PLOS Computational Biology* on **2026-09-01**. The items below were found on
2026-09-02, one day after submission; Item 5 was found on 2026-09-18. They were held so the
deposited supplementary files matched what the reviewers received, and have now been applied to
the repository copy only — see the status block above. Supply the corrected S1 Dataset at
revision, alongside the response to reviewers.

## Bottom line first

**Items 1–7: no reported number changes, and nothing there is an erratum.** Each is a
prose/description defect in documentation.

**Item 8 is different.** It is a Module 2 analysis defect that changes the published CSK and SRC
metastable states (and so Fig 4B/C labels, Fig 6 state names, Fig S7 and Suppl. Note 10/Table 7).
The Fig 4B/4C conclusions survive the fix, and CDK1 is largely unaffected. See Item 8 for the
evidence and the per-element impact. Everything below about Items 1–7 still holds.

For Items 1–7: the pipeline code is correct, and every value in
`master_kinase_analysis_results_v7r3.csv`, the `Phase*` tables, the figures and the paper is
unaffected. Items 1, 2 and 5 sit inside a submitted supplementary file (S1 Dataset); each is a
single sentence, and Items 1 and 2 contradict the same document's own per-column text. Item 5
additionally has a framing consequence for two rows of S1 Table — see that item.

---

## Item 1 — S1 Dataset: the "CSK numbering" sentence is wrong  ✅ ACTIONED 2026-09-18

**File:** `docs/AlloQuant_master_CSV_data_dictionary_v7r3.xlsx` (and the `.pdf`), general notes
section — the row beginning "Distances in Angstrom (A), angles in degrees…".

**It currently says:**

> numbered landmarks (n99, v104, k105, e107, m118, m120, e121, i150, y156, d220) are labelled in
> **CSK numbering** and map to the equivalent target position via the HMM.

**Why it is wrong:** the labels are **PKA** numbering, not CSK. The document's own per-column
entries say so correctly — e.g. `K105_E121_Dist` reads "Pfam Pkinase HMM node 62 (**PKA canonical
K105 position**; GLN65 in CSK domain numbering)". CSK's residue at that node is Gln65, so the
label `K105` cannot be CSK numbering. Confirmed by regenerating the mapping (see the reference
table below): none of `N99`, `K105`, `E107`, `M118`, `M120`, `I150` is the named residue type in
CSK.

**Fix:** change "labelled in CSK numbering" to "labelled in **PKA** numbering (the canonical
reference kinase), and map to the equivalent target position via the HMM regardless of the residue
type actually present". Leave every per-column description alone — they are already correct.

## Item 2 — S1 Dataset: "the two shell methionines"  ✅ ACTIONED 2026-09-18

**File:** same, the `Shell_M118_M120_Dist` row.

**It currently says:** "Min side-chain distance between the two shell **methionines** M118 and
M120 (drives Shell_State)."

**Why it is wrong:** neither position is a methionine in any model kinase — Ile79/Thr81 in both CSK
and SRC, Leu78/Phe80 in CDK1. "Methionine" is the PKA identity, inherited from the label.

**Fix:** "…between the two hydrophobic-shell positions M118 and M120 (PKA numbering; Ile79/Thr81 in
CSK and SRC, Leu78/Phe80 in CDK1)".

## Item 3 — `data/MANIFEST_md5.txt` needs **no** change for Items 1-2 (nor for Item 5)

Checked 2026-09-02: the manifest's paths are relative to `data/` and cover only
`cdk1_ccnb1_260517/`, `csk_monomer_260806_fullmsa/`, `csk_src_dimer_260718/`,
`interface_analyses/`, `validation_3d7t/` and `data/README.md`. **No `docs/` file is listed**, so
editing the dictionary does not invalidate it and it must not be rebuilt on account of Items 1-2.

Re-verify only if a file under `data/` changes:

```bash
cd data && md5sum -c MANIFEST_md5.txt
```

As of 2026-09-02 all 203 entries verify clean.

## Item 4 — cosmetic, repo only (no deposit impact)  ✅ ACTIONED 2026-09-19

* `METRICS_CHEATSHEET.md` line 1 is a stray pasted GitHub attachment URL sitting above the title.
  — removed; the file now opens on its `# Kinase Structural Metrics Cheat Sheet` title.
* `scripts/modules/placeholder.txt` is leftover cruft, absent from the documented repo layout.
  — deleted. Nothing referenced it (the only mention in the repo was this list), and
  `scripts/modules/` holds ten real files, so the directory is not at risk of disappearing.

Neither file is part of the submission; no deposited or staged file changed.

## Item 5 — S1 Dataset: the `ActLoop_CT` Kincore comparison  ✅ ACTIONED 2026-09-18

**Found:** 2026-09-18. **File:** same dictionary, the `ActLoop_CT` row (cell `I38`).

**It said:** "The 5.5 A cutoff is AlloQuant's own sensitivity choice, deliberately looser than
Kincore's 6.0 A."

**Why it was wrong:** 5.5 Å is *stricter* than 6.0 Å, so the cutoff alone makes the criterion
less permissive, not more; and the sentence implies the two tools measure the same contact
differing only in cutoff. They do not — Kincore uses two named atoms of one residue
(APE9 Cα ↔ HRD-Arg backbone O), this uses an all-atom minimum over a four-residue window. The
window, not the cutoff, is what makes AlloQuant's criterion the more permissive of the two. Git
history shows the window, all-atom minimum and 5.5 Å cutoff are inherited unchanged from
`chimerax_hmm_worker_v6r6.py`; v7r3 changed only the anchor, so the framing predates v7r3.

**⚠️ Knock-on for S1 Table:** S1 Table compares "A-loop C-terminal" against Kincore and marks
the SRC row **agree** and the CSK row **divergent**. If the two quantities are not comparable,
that agreement is between different measurements. The same caveat applies to the "A-loop
N-terminal" row (AlloQuant scans `f+3…f+6` against `hrd-1`; Kincore uses the single DFG6–XHRD
pair). **This is a framing issue in a submitted table and is the item most likely to need
raising with the editor.** No measured value in S1 Table changes. Full analysis:
`plos_submission/updates/2026-09-19/CHANGES_2026-09-19.md` §3, and
`S1_Table_note_2026-09-19.md` in the same directory for proposed footnote wording.

Also corrected in `README.md` (v7r3 release note #2) and the worker source comment — commits
`346eefb` and `b1fc720` on branch `docs-cheatsheet-corrections`.

## Item 6 — S1 Dataset: what `Type` is  ⏳ OPEN (found 2026-09-26)

**File:** same dictionary, the `Type` row (cell `I22`).

**It says:** "Kinase identity assigned by best sequence-identity match to the FASTA landmark set
(get_best_landmark_for_chain)."

**Why it is imprecise:** sequence identity picks *which landmark entry* measures the chain, but
the value written is that entry's **name**. The name comes from `proteins.yaml` or the folder
token, via `extract_fasta.py`, so `Type` encodes whatever the operator's name encodes: family
(`CDK1`), allele (`CSK-WT`), or construct. It never encodes modification. The row's own second
sentence (SRC pY159 reads `SRC-WT`) already shows this. Every consumer treats it as the
operator's name: the Module 2 `TARGET_TYPE` filter, the `Condition_reviewed` join, and the addons'
exact matches (`Type == 'CSK-WT'`). Found while fixing v7r4 (`6d638ee`, `91858a6` on
`v7r4-bugfix`). There, one name covering several sequences (crystal copies, isoforms, constructs)
produced dedup variants `NAME_1`, `NAME_2` …, and v7r4 now reports these as `NAME`.

**Proposed wording:** "Protein name (from proteins.yaml / the folder token) of the FASTA landmark
entry that best matches the chain's sequence (get_best_landmark_for_chain); its landmark positions
are the ones used to measure this chain. Carries what the name carries (family, allele, construct),
not modification: the SRC pY159 variant reads SRC-WT, and the pY159 vs WT distinction lives in
Simulation_ID/Directory."

**No reported number changes.** The published datasets have no dedup variants, so their `Type`
values and every measured column are what v7r4 gives. This clarifies a definition; bundle it with
Items 1, 2 and 5 only if the S1 Dataset is being re-supplied anyway.

**Related Module 2 code issues (not documentation; not fixed; published outputs unaffected):**
* `multimer_core_engine.R:276` selects the target by *prefix* (`^(?:[A-Z]-)?TARGET`), so `CDK1`
  would also take `CDK10`/`CDK12`. No published dataset has two such names.
* `multimer_core_engine.R:211` rebuilds `Condition_reviewed` by finding `tolower(Type)` inside the
  condition text, so a chain whose name is not spelled out in its condition string keeps its
  designed apo/holo state without any warning.

## Item 7 — `README.md` has no v7r4 section  ⏳ OPEN (found 2026-09-26)

**Repo only, not part of the submission.** Merging `v7r4-bugfix` into `main` (after this branch)
brings in the three v7r4 files but not a README that describes them. `README.md` on
`v7r4-bugfix` never mentions v7r4. Its title, Deployment and Two-Stage Workflow commands, Outputs,
Documentation, Repository Layout and Runtime Working Directory all name only the v7r3 files
(lines 1, 6, 52–91, 99–110, 121–155 at `756e478`). Right now the v7r4 changes are recorded only in
the UPDATE LOG header of `scripts/modules/chimerax_hmm_worker_v7r4.py`. The README of the local
deployment copy (`../alloquant_module1_v7r4/`, outside every repository) is written for that
copy and does **not** replace this one, although parts of it can be reused.

**What the section needs:**
* **Two versions, and which is which.** Line 30 says "v7r3 is the only version shipped here and
  is the version of record for every published result". The second half stays true. The first
  must change: v7r3 remains frozen and byte-identical, and v7r4 is the recommended version for
  new runs.
* **What v7r4 fixes**, so a user can tell whether their v7r3 output is affected:
  (1) whole arms silently dropped at exit 0 when the landmark name gate matches nothing
  (`c91474a`); (2) worker failures propagated, model coverage asserted, and the chunk directory
  kept on failure (`c91474a`, `dfe90c7`); (3) `NAME_n` landmark entries, where one name covers
  several sequences, now measured with their own map instead of the bare `NAME`'s (`6d638ee`),
  with `Type` still reading `NAME` (`91858a6`).
* **The value statement, stated exactly.** v7r4 changes no reported value **except** for runs
  whose `sequences.fasta` has `NAME_n` kinase entries. For those runs v7r4 is correct and v7r3
  is not. None of the published datasets has such entries (checked 2026-09-26), so every
  published number stands. The label stays v7r4 on purpose (Akira, 2026-09-26), because the
  difference depends on the run, not on the version.
* **New outputs:** `landmark_refs_v7r4.csv` (only when a chain used a `NAME_n` entry) and the
  printed variant table, whose substitution counts are lower bounds (`756e478`). The v7r4
  output filenames (`hmm_…_v7r4.csv`, `master_…_v7r4.csv`) belong in Outputs, Repository Layout
  and Runtime Working Directory.
* **Naming rule:** isoforms, constructs or mutants that should stay apart in Part 2 need
  distinct `NAME-variant` names in `proteins.yaml`. `Type` carries the operator's name (see Item 6).
* **Still the operator's job in v7r4:** `run_step` still calls `python3` directly, so activate the
  environment first. Count rows against the design, because v7r4 checks model coverage, not the
  experimental design.

**Related:** `docs/AlloQuant_output_file_manifest_v7r3` and the data dictionary are v7r3
documents. Decide whether v7r4 gets its own versions or a note on the v7r3 ones. The master CSV
schema is unchanged, so the dictionary needs only Item 6. The manifest lacks
`landmark_refs_v7r4.csv` and the v7r4 filenames.


## Item 8 — Module 2: `Spine_Bridge_Dist` enters state discovery with two definitions  ⏳ OPEN (found 2026-09-27)

**Not documentation. It changes reported results at the state level.** Found while running
Module 2 on a new AF3 condition set (Aur A ± TPX2 ± CEP192 ± pT ± ATP·Mg), where every GMM state
was either all-apo or all-holo.

**The defect.** Module 1 computes `Spine_Bridge_Dist` as the αC–β4 loop to the **ligand** when
one is bound, and to the **k−3 (β3) residue** when none is. This is documented in S1 (row 35:
"fallback for apo"), so the *metric* is disclosed correctly. The two branches measure different
things on different scales (medians CDK1 3.91 vs 7.92 Å, CSK 4.15 vs 7.54, SRC 4.12 vs 7.72).
`multimer_core_engine.R` (v7r1) blanks ligand columns in apo chains only when their names match
`ATP|Mg` (line 304 on `main`), so `Spine_Bridge_Dist`:
1. enters the PCA/GMM feature set (`universal_dist_cols`), which separates apo from holo chains
   by construction;
2. is tested across apo-vs-holo groups as one quantity. It ranks 1st–2nd in most published
   CDK1 apo-vs-holo Phase 8 contrasts (p.adj ≈ 1e-33); that ranking comes from the definition,
   not biology;
3. enters per-state MAC wherever a state mixes apo and holo chains.

**Evidence** (outside this repo; nothing here changed):
- `../alloquant_spine_bridge_check_2026-09-27/NOTE_spine_bridge_published_runs.md`. All three
  published Module 2 runs were reproduced exactly (state labels identical: CDK1 600/600, CSK
  795/795, SRC 795/795), then repeated with the column dropped.
- `../SRC/af3_output/260927_spine_bridge_masked_engine_test/README.md`. Engine copy with the
  one-line mask fix, run on the deposited CSK–SRC 260718 data. Per-figure comparison in
  `compare_report.txt`.

| run | states pub → fixed | ARI | Cramér's V, state ~ apo/holo |
|---|---|---|---|
| CDK1 | 7 → 8 | 0.81 | 0.96 → 0.93 (the split is real conformation) |
| CSK | 8 → 5 | 0.62 | 0.73 → 0.21 |
| SRC | 7 → 5 | 0.58 | 0.93 → 0.38 |

**Impact by manuscript element** (masked engine vs published):

| element | effect |
|---|---|
| Fig 4A, pooled MAC per condition | holo conditions unchanged; apo conditions ≈ −0.007 (CSK apo/apo 0.134 → 0.127). "CSK rises on priming, SRC flat" unchanged |
| Fig 4B, "pY419 flips SRC state" | **conclusion holds**: both-ATP and primed SRC still occupy disjoint state sets. Labels change; the main primed state (old S8) now shares a state with apo SRC |
| Fig 4C, state handshake | **conclusion holds, stronger**: primed V 0.44 → 0.58 (p 4e-6 → 3e-16); unprimed n.s. (0.26 → 0.19). Axis labels change |
| Fig 6 | the representative models and their SB distances are unchanged; only their state names change (CSK 3→5, 9→3, 5→4; SRC 8→4) |
| Fig S7 (SRC State 4 vs 7/8) | **must be redone**: old S8 merges with apo SRC and old S4 splits, so the contrasts no longer exist as defined |
| Suppl. Note 10 / Table 7 (per-state MAC) | recompute on the new states |
| Fig 2C (CDK1 state composition) | light check: States 4/5/6/8 map 1:1, State 1 mostly; States 7 and 9 re-split |
| Any text naming `Spine_Bridge` as an apo/holo discriminator | remove or qualify |

**Not affected:** every Module 1 value and the S1 Dataset; within-condition analyses (the engine
appends the physical ligand state to each condition, so no published condition mixes apo and
holo chains, 0/7 for CSK and SRC), including the paired-coupling figure and per-condition MAC
among holo conditions.

**Fix: repo side.**
1. Engine: add `Spine_Bridge` to the apo mask:
   `ligand_cols <- grep("ATP|Mg|Spine_Bridge", all_numeric_cols, ignore.case = TRUE, value = TRUE)`.
   Ship it as a **new** engine version so v7r1 stays byte-identical for the published outputs
   (the same policy as v7r3/v7r4), and have the driver call it.
2. Longer term, in a future Module 1 version: split the column into `Spine_Bridge_Ligand_Dist`
   and `Spine_Bridge_B3_Dist` (a schema change: S1 Dataset and data dictionary).
3. Regenerate `data/csk_src_dimer_260718/plots_and_stats_{CSK,SRC}_GMM` (and CDK1 if the paper
   adopts the fixed engine everywhere) with the new engine, including Phase 8 for the new SRC
   state pairs. **`data/MANIFEST_md5.txt` must then be rebuilt** (Item 3 applies: files under
   `data/` change).
4. Figure scripts hard-code state names and must be updated: `figures/fig4_model2.py`
   (`src_states`, `SRC_ROWS`, `CSK_COLS`, the V values in the panel titles), `figures/compose_fig6.py`
   (docstring state names), `figures/figS7_src_state_volcanos.py` (`PAIRS`).
5. README: record the fix and the value statement (which outputs change and why) in the section
   for the new engine version, next to Item 7.

**Fix: manuscript side.**
1. Methods (Module 2): state that ligand-dependent distances, including `Spine_Bridge_Dist`, are
   excluded from state discovery when apo and holo chains are pooled.
2. Replace Fig 4B/C (labels, V values, legend footnote), Fig S7, and Suppl. Note 10/Table 7 with
   the re-run outputs; rename the states in Fig 6 and its legend; check Fig 2C States 7/9.
3. Response to reviewers: disclose it as a correction found by the authors after submission. The
   conclusions of Fig 4B/C are unchanged (4C strengthens); the state labels and the S7 contrasts
   change.

**Related, not fixed (decide at revision):** because the ATP/Mg columns are masked in apo chains,
per-condition MAC is already computed on different column sets for apo and holo conditions (Fig
4A compares them). Report that or equalise the feature set.

---

## Already fixed on 2026-09-02 — no action needed, listed for the record

These were corrected in `METRICS_CHEATSHEET.md`, which is **repo documentation only and not part of
any submitted supplementary file**, so fixing it immediately carries no risk to the submission. All
were verified against `scripts/modules/chimerax_hmm_worker_v7r3.py`.

| Column | Cheatsheet had said | Code actually does |
|---|---|---|
| `D1_Dist` | "DFG-Asp to αC-Glu" | αC-Glu**(+4)** Cα → DFG-**Phe** Cζ |
| `D2_Dist` | "DFG-Asp to HRD-His/Tyr" | β3-Lys Cα → DFG-**Phe** Cζ |
| `C_Helix` | "based on the integrity of `SB_Dist`" | β3-Lys **Cβ** → αC-Glu **Cβ** ≤ 10 Å, an independent measurement |
| `D220_HRD_Dist` | "to the **backbone** of the catalytic loop" | min **side-chain** distance to the **HRD-His** |

Only the `D1_Dist` item was a missed v7r3 sweep. The other three were unchanged at v7r3 and so fell
outside a sweep correctly scoped to what v7r3 altered; they appear never to have been right.

A new "Section 0" was also added to the cheatsheet explaining that numbered landmarks are PKA-named
HMM *positions* rather than residue identities, carrying the reference table below, and the
`K105`/`N99`/`Y156`/`I150`/`M118`/`M120` entries were reworded so they no longer assert PKA
chemistry that does not exist in the model kinases. The most substantive of those: the old
`K105_E107_Dist` text described a Lys→Glu salt bridge that occurs in **none** of the three model
kinases (Gln65-Leu67 in CSK, Gln67-Tyr69 in SRC, Ser65-Gln67 in CDK1 — and Gln against Leu cannot
hydrogen bond at all).

---

## Reference: HMM landmark → actual residue, by model kinase

Positions are in the numbering of the deposited construct (`data/*/inputs/sequences.fasta`).
**CDK1 is full-length, so its numbers are UniProt numbering directly.** For CSK and SRC add
**+185** and **+257** respectively for canonical full-length numbering.

| Column label | HMM node | CSK | SRC | CDK1 | Label letter holds? |
|---|---|---|---|---|---|
| `N99`  | 56  | R59 | R61 | R59 | No — Arg, never Asn |
| `V104` | 61  | V64 | V66 | V64 | Yes |
| `K105` | 62  | Q65 | Q67 | S65 | No — Gln/Ser, never Lys |
| `E107` | 64  | L67 | Y69 | Q67 | No — Leu/Tyr/Gln, never Glu |
| `M118` | 75  | I79 | I79 | L78 | No — Ile/Leu, never Met |
| `M120` | 77  | T81 | T81 | F80 | No — Thr/Phe, never Met |
| `E121` | 78  | E82 | E82 | E81 | Yes — Glu in all three |
| `I150` | 113 | Y119 | Y119 | F118 | No — Tyr/Phe, never Ile |
| `Y156` | 119 | F125 | Y125 | V124 | Partial — Tyr in SRC only |
| `D220` | 179 | D183 | D187 | D186 | Yes |

Non-numbered landmarks are genuinely conserved and need no caveat: `k` (β3-Lys) = K37/K38/K33,
`c` (αC-Glu) = E51/E53/E51, `hrd` (HRD-His) = H127/H127/H126, `f` (DFG-Phe) = F148/F148/F147,
`ape` = E171/E175/E173.

### How to re-derive this table

Two independent checks confirmed it, and both should be repeated if the table is ever regenerated:

1. It reproduces the dictionary's own stated identities, `k105 → GLN65 (CSK)` and
   `e121 → GLU82 (CSK)`.
2. The construct→canonical offsets are **constant across four independent landmarks**
   (CSK +185, SRC +257) and agree with Kincore's own residue calls on 3D7T in
   `data/validation_3d7t/kincore_3d7t_output.txt` — yielding SRC K295/E310/H384/F405, the
   canonical numbers. CDK1's numbering is corroborated by the condition names themselves
   (`pt14`, `py15`, `pt161` → CDK1 T14/Y15/T161).

Procedure: run `hmmalign --outformat A2M` with `~/pfam/Pkinase.hmm` against the deposited
`inputs/sequences.fasta`, then walk the A2M string exactly as
`scripts/modules/extract_landmarks.py::extract_nodes_from_aligned_seq` does (uppercase = match
state, advancing both the HMM node counter and the sequence index; lowercase = insert, advancing
only the sequence index; `-` = delete, advancing only the node counter), using the `TARGET_NODES`
map from that same file. Do not reimplement the walk from scratch — mirror the shipped function, or
the node numbering will drift.
