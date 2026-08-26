# Pipeline add-ons

Complementary, **post-hoc** analyses that run on the outputs of the core kinase/allostery
pipeline. These are intentionally **kept separate from `modules/` and the wrappers** — they do
not modify or depend on the core pipeline internals, and the core pipeline runs identically
whether or not they are present.

The two general-purpose add-ons (`paired_chain_coupling.py`, `mac_confidence_control.py`) are
argparse-driven and non-interactive. The four manuscript-specific controls added below them read
their inputs from the `ALLOQUANT_*` environment variables documented in `figures/README.md`.
Paths below are written relative to a pipeline output directory (the one holding
`plots_and_stats_*/`).

## `paired_chain_coupling.py` — cross-chain conformational coupling

Tests whether the two kinase domains of a heterodimer are conformationally coupled at the
**single-structure** level, exploiting the fact that each multimeric AlphaFold3 model contains
both chains at once (their per-model geometry is 1:1 paired).

**Input:** the two per-chain `Phase7_Complete_Structural_Metadata.csv` tables produced by
Module 2 (e.g. `plots_and_stats_CSK_GMM/` and `plots_and_stats_SRC_GMM/`).

**Method:** join on `Simulation_ID`; for every chainA-metric × chainB-metric pair, Spearman ρ
computed **within each experimental condition** and combined by Fisher-z meta-analysis (removes
the shared-design confound); optionally recomputed as **partial** correlation controlling for
per-model AF3 confidence (`ipTM`, `PAE_Mean_AB`) to rule out a model-quality artifact;
Benjamini–Hochberg FDR across the grid.

Partialling uses the closed-form partial-correlation identity **applied sequentially, one
covariate at a time**, taking the marginal ρ(x,z)/ρ(y,z) at each step — not the fully recursive
form. Anything that needs to reproduce these numbers must mirror that rather than "improve" it.

**Run (conda env `main`):**
```bash
conda run -n main python addons/paired_chain_coupling.py \
    --chainA plots_and_stats_CSK_GMM/Phase7_Complete_Structural_Metadata.csv \
    --chainB plots_and_stats_SRC_GMM/Phase7_Complete_Structural_Metadata.csv \
    --out-dir addons/coupling_out
```
Add `--confidence ""` to skip the confidence control if `ipTM`/`PAE_*` are all-NaN (runs
predating the AF3-assessment pipeline upgrade — re-run stage 1 with the confidence JSONs to
populate them; distances are unaffected). `--min-n` (default 25) is the minimum number of models
a condition needs to enter the meta-analysis.

**Outputs:** `paired_coupling_full_grid.csv` (all pairs, raw + partial ρ/p, FDR) and
`paired_coupling_homologous.csv` (same-named "conformational mirroring" pairs).

**Interpretation notes.** Use `fdr_partial`, not `p_partial` — the uncorrected p-values make
weakly coupled metrics look significant. And the coupling is *concentrated in, but not strictly
confined to,* the regulatory-spine module: in the CSK–SRC analysis two non-spine active-site
readouts (the P-loop–ATP contact and the DFG-Asp ψ backbone dihedral) also survive FDR.

**Provenance:** used for the CSK–SRC manuscript **Figure 5**, rendered by
`figures/fig_coupling.py` in this repository.

## `mac_confidence_control.py` — is MAC a model-confidence artefact?

Recomputes Global Network Rigidity (**MAC** = Mean Absolute Correlation) while partialling
per-model AlphaFold3 confidence out of **every edge** of the Spearman distance-correlation
matrix, then re-runs the same between-condition contrasts on the partialled edges. This answers
the reasonable objection that condition-to-condition MAC differences reflect differences in
prediction quality — better-predicted ensembles being more internally self-consistent — rather
than genuine mechanical coupling.

**Method:** the same sequential closed-form partial-correlation identity used above. To stay
comparable with the published MAC values it mirrors the R engine
(`modules/multimer_core_engine.R`, Phase 5 / Phase 6b): only `*_Dist` columns enter the network
(angles break on periodicity), a column is used within a group only if it has enough non-NA
values, uncomputable correlations are set to 0, MAC is the mean of |ρ| over the upper triangle,
and conditions are contrasted by Mann–Whitney on the edge distributions.

**Run (conda env `main`):**
```bash
conda run -n main python addons/mac_confidence_control.py \
    --metadata plots_and_stats_CSK_GMM/Phase7_Complete_Structural_Metadata.csv \
    --confidence master_kinase_analysis_results_v7r3.csv \
    --label CSK --out-dir addons/mac_confidence_out \
    --phase5 plots_and_stats_CSK_GMM/Phase5_Global_Network_Density.csv
```
* `--group-col` selects which MAC is tested: `Condition` (default, Phase 5 *global* MAC) or
  `Macro_State` (Phase 6b *intrinsic* MAC).
* `--covars` defaults to `ipTM PAE_Mean_AB`. For a **monomeric** series pass `--covars pTM`:
  a monomer has no interface, so it has no ipTM or interface PAE.
* `--covar-mode common` (the default) partials out only the covariates populated in *every*
  group, so all groups lose the same variance. **`per-group` makes MAC_partial non-comparable
  across groups and must not be used for contrasts** — this matters whenever a monomeric
  condition is compared against a complex.
* `--phase5` is optional but recommended: it checks that MAC_raw reproduces the pipeline's own
  Phase 5 values before anything is concluded from MAC_partial.

**Outputs** (each prefixed with `--label`): `{label}_mac_raw_vs_partial.csv` (per-group MAC before
and after partialling, with Δ), `{label}_mac_pairwise_tests.csv` (every contrast, raw and
partialled p/FDR, and whether the direction is preserved), `{label}_confidence_by_group.csv` (the
covariate distributions themselves), and `{label}_mac_raw_vs_partial.png` (diagnostic plot).

**Interpretation note — the important one.** Partialling removes variance, so |ρ| and therefore
MAC fall for **every** group. That is arithmetic, not evidence of weakening. The claim under test
is that the *between-condition contrast* survives: direction preserved and still significant.

**Provenance:** used for the manuscript's confidence-control supplementary note, tables and
figure; its output directory is what `figures/figSI_confidence.py` reads.

## The composition control — `csk_mac_by_condition.py`, `mac_nmatched_control.py`, `cdk1_mac_decompose.py`, `cdk1_matched_panel.py`

Global (condition-level) MAC is one correlation network built across all of a condition's models,
so it responds to the condition's **state composition** as well as to the coupling *within* its
states: pooling models from structurally distinct metastable states induces correlation among
every metric that separates them. These three scripts decompose a global-MAC difference into
those two contributions. They are the code behind the manuscript's composition control.

* **`csk_mac_by_condition.py`** — for each of the four biological CSK conditions, fix that
  condition's own global-MAC metric panel, then recompute the same statistic separately within
  each occupied metastable state **on that fixed panel**, so the per-state values are
  panel-matched to the condition's own bar. Reproduces the pipeline's
  `Phase5_Global_Network_Density.csv` values as a check, and writes
  `csk_mac_by_condition.csv`. These are the per-state values tabulated in the supplement.
  States with fewer than 8 models are reported but not emitted. (Up to 2026-08-24 this CSV also
  fed a per-state overlay on Figure 4 panel A; that overlay was removed, see below.)
* **`cdk1_mac_decompose.py`** — the same decomposition applied to CDK1, testing whether the
  CDK1 "thaw" survives it. It does: the thaw is overwhelmingly a within-state coupling change.
* **`cdk1_matched_panel.py`** — checks the CDK1 result is not an artefact of comparing different
  metric panels, by recomputing each contrast on the metrics viable in **both** of its
  conditions. Replicates the R engine's semantics exactly, including the fact that R *retains*
  zero-variance columns (their correlations become NA, then 0, dragging MAC down) — do not
  "fix" that by dropping them.

* **`mac_nmatched_control.py`** — tests whether a state is *less coupled than the ensemble it
  sits in*, which the floor below cannot answer. Each state is compared with 2000 random subsets
  of its own condition **at that state's own sample size**, state labels ignored, and reports the
  departure in null standard deviations. In the primed CSK condition every state falls short of
  its size-matched expectation and the occupancy-weighted value falls 0.091 below it; the
  both-ATP baseline is null. So conditioning on state removes the whole of the pooled rise in
  the primed condition and none of it in the baseline.

Note the small-sample floor: the expected mean |ρ| between independent variables is
`sqrt(2 / (pi * (n - 1)))` — 0.080 at n = 100, 0.125 at 42, 0.216 at 15 — so every per-state
value is inflated in absolute terms and is only meaningful as a matched comparison.

**That floor is the wrong reference for the composition question.** It is the expectation for
*independent* variables, which establishes whether a cell carries information at all. A random
subset of a *mixed* condition inherits the mixture's own induced correlation and sits far above
it: 0.257 against a floor of 0.125 at n = 42 in the primed condition. Use
`mac_nmatched_control.py` for that comparison, not the floor.

The same *n* dependence is why Figure 4 panel A no longer overlays per-state values on the
pooled bars: values at n = 15-45 and n = 100 do not belong on one linear axis.

**Run:**
```bash
ALLOQUANT_SRC=/path/to/csk-src_output  python addons/csk_mac_by_condition.py
ALLOQUANT_SRC=/path/to/csk-src_output  python addons/mac_nmatched_control.py
ALLOQUANT_CDK1=/path/to/cdk1_output    python addons/cdk1_mac_decompose.py
ALLOQUANT_CDK1=/path/to/cdk1_output ALLOQUANT_SRC=/path/to/csk-src_output \
    python addons/cdk1_matched_panel.py
```
`ALLOQUANT_CSK_MAC_CSV` and `ALLOQUANT_CSK_NMATCHED_CSV` set where `csk_mac_by_condition.py`
and `mac_nmatched_control.py` write their CSVs (defaults: `csk_mac_by_condition.csv` and
`csk_mac_nmatched.csv` in the current directory). `mac_nmatched_control.py` is seeded, so it
reproduces exactly. The two CDK1 scripts print their tables
to stdout and write nothing.

**Provenance:** this control is independent of `mac_confidence_control.py` above, which
addresses prediction confidence and says nothing about composition.

## `iface_py419.py` — does SRC pY419 change the docking interface?

A direct, assumption-free measurement on the AF3 ensembles rather than on pipeline summaries:
reads each `model.cif` with `gemmi` and computes inter-chain contacts, buried surface area
(numerical Shrake–Rupley with a 92-point sphere), polar bridges, and per-residue-pair contact
frequency, for four conditions spanning the unprimed/primed and CSK-apo/CSK-holo contrasts.

**Run:**
```bash
ALLOQUANT_SRC=/path/to/csk-src_output python addons/iface_py419.py
```
Writes `analysis_iface_py419.json` (override with `ALLOQUANT_IFACE_JSON`), which
`figures/figSI_interface.py` reads. It walks every model of four conditions, so expect it to
take a while.

**Interpretation note.** The displacements it finds are small enough to need a resolution check
before they can be called real: compare them against AlphaFold3's own predicted aligned error at
the same residue pairs. In the manuscript they come to a few percent of the per-pair PAE, which
is why the conclusion is "below resolution" rather than "invariant".

**Provenance:** Supplementary Note 8.

## `chain_local_rmsf.py` — chain-local per-residue Cα RMSF

Measures how much each chain's own fold deforms across an ensemble, without letting the two
subunits' motion relative to each other leak into the number.

**Input:** an AlphaFold3 output directory holding one subdirectory per condition, each with
`seed-*_sample-*/model.cif`. It reads the coordinates directly and does not touch any
`plots_and_stats_*` output, so it carries no AlloQuant version tag and its result is unaffected
by the v7r3 metric corrections.

**Method:** each chain is fitted on its own conformation in the **first model of that
condition**, using that chain's Cα atoms only (Kabsch). Per-residue RMSF is then taken about
that chain's mean position across the condition's models, over residues present in every model.
Reports median, 90th percentile, the percentage of residues above a mobility threshold
(0.5 Å by default), and mean pLDDT, plus the full per-residue profile on request.

**Run:**
```bash
python3 addons/chain_local_rmsf.py --root "$ALLOQUANT_CDK1" \
    --conditions a-cdk1_b-ccnb1-166_0atp=+CCNB1 \
                 a-cdk1_b-ccnb1-166_1atp=+CCNB1+ATP \
                 "a-cdk1-pt161_b-ccnb1-166_1atp=pT161 (active)" \
                 a-cdk1-pt14-py15-pt161_b-ccnb1-166_1atp=triple-P \
    --chains A=CDK1 B=CCNB1 \
    --out chain_local_rmsf.csv --profiles-out chain_local_rmsf_profiles.csv
```

**Two traps.** Chain-local fitting is not interchangeable with the shared-frame superposition
used for a whole complex: a shared frame retains inter-chain motion, so the two sets of values
must not be compared or placed on a common axis. And polymer membership must be decided by
`label_seq_id`, not by `group_PDB`: AlphaFold3 writes modified residues (TPO, PTR, SEP) as
HETATM even though they belong to the chain, while true ligands (ATP, MG) carry no
`label_seq_id`. Filtering on `group_PDB` silently shortens every phosphorylated chain and shifts
its statistics.

**Provenance:** S9 Table and Supplementary Note 15.
