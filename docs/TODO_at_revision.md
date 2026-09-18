# Deferred documentation corrections — action at manuscript revision

**Created:** 2026-09-02. **Updated: 2026-09-18 — Items 1 and 2 are now ACTIONED in the repo.**

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

**No reported number changes. Nothing here is an erratum.** Every item is a prose/description
defect in documentation. The pipeline code is correct, and every value in
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

## Item 4 — cosmetic, repo only (no deposit impact)

* `METRICS_CHEATSHEET.md` line 1 is a stray pasted GitHub attachment URL sitting above the title.
* `scripts/modules/placeholder.txt` is leftover cruft, absent from the documented repo layout.

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
