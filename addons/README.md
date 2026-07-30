# Pipeline add-ons

Complementary, **post-hoc** analyses that run on the outputs of the core kinase/allostery
pipeline. These are intentionally **kept separate from `modules/` and the wrappers** — they do
not modify or depend on the core pipeline internals, and the core pipeline runs identically
whether or not they are present.

Both add-ons are argparse-driven and non-interactive. Paths below are written relative to a
pipeline output directory (the one holding `plots_and_stats_*/`).

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

**Provenance:** used for the manuscript's confidence-control supplementary note, tables and figure.
