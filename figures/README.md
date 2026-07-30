# Manuscript figure-generation scripts

Scripts that render the main-text and supplementary figures for the
AlloQuant / trans-allosteric kinase manuscript. Figures are regenerated
for presentation from AlloQuant's computed outputs using a shared color
palette and type scale (`palette.py`); the underlying values are those
produced directly by the pipeline.

These are the **release-adapted** copies: they resolve their input
directories from environment variables (below) and import `palette.py`
relative to their own location, so they run from a clean checkout. The
working copies used while drafting the manuscript carry hardcoded
absolute paths instead.

## Script → figure

| Script | Figure |
|---|---|
| `fig2_model1.py` | Figure 2 — CDK1–CCNB1 thawing / re-tensioning |
| `fig3_cdk1_states.py` | Figure 3 — CDK1 activation structural series |
| `fig4_model2.py` | Figure 4 — CSK–SRC handshake |
| `fig_coupling.py` | Figure 5 — CSK/SRC spine coupling |
| `figSI_cage.py` | Figure S1 — P-loop pT14 phosphate clash |

Shared: `palette.py` (colors + type scale, imported by every script);
`PALETTE.md` (palette documentation).

### Not in this directory yet

Figure 6 (the CSK–SRC structural composite) and Figure S4 (the MAC
confidence control) are rendered by scripts that have not been ported
here yet. Figure 6 also depends on pre-rendered ChimeraX panels and
their `.cxc` macros rather than on pipeline data alone.

### Retired

`figSI_substates.py` rendered a CDK1 rigid/floppy substate figure that
was **removed from the manuscript at v7r3**: the substates it showed
were an artefact of the v7r2 `D1_Dist` anchor, and under v7r3 the two
populations merge into a single state. The script is kept for the
record; it does not correspond to a current figure.

## Version status

The science in these scripts is on **v7r2** values and terminology. The
v7r3 rebase (corrected `D1_Dist` anchor, corrected `ActLoop_CT` anchor,
"metastable" rather than "macro-state") has been applied to the drafting
copies but **not yet ported here**, so re-running these will not
reproduce the current published figures. Port the science changes while
keeping the environment-variable adaptation above — do not overwrite
these files with the hardcoded-path working copies.

Specifically outstanding: `fig_coupling.py` here partials **`ipTM` only**,
whereas the analysis it depicts — and the published Figure 5 legend —
partial both `ipTM` and `PAE_Mean_AB`, matching
`addons/paired_chain_coupling.py`. The ported version must partial both
and should gate its values against that add-on's published
`rho_partial` at render time.

## Dependencies

Python 3 with `matplotlib`, `numpy`, `pandas`, `scipy`, and `Pillow`
(the last only for `figSI_cage.py`).

## Data locations (environment variables)

Scripts read AlloQuant pipeline outputs (`Phase*` /
`plots_and_stats_*` directories). Set these to your local output dirs
(defaults are relative placeholders):

| Variable | Points to |
|---|---|
| `ALLOQUANT_CDK1` | CDK1–CCNB1 pipeline output directory |
| `ALLOQUANT_SRC` | CSK–SRC pipeline output directory |
| `ALLOQUANT_STRUCT` | folder of ChimeraX-rendered structural panels (`struct/`) |
| `ALLOQUANT_FIGDIR` | output directory for rendered figures (default: current dir) |

Run each script from this directory (so `palette.py` is importable), e.g.:

```bash
ALLOQUANT_SRC=/path/to/csk-src_output python fig4_model2.py
```

## Note on the structural panels

`fig3_cdk1_states.py` and `figSI_cage.py` compose pre-rendered ChimeraX
images (from `ALLOQUANT_STRUCT`) rather than generating them from data.
Those renders are not produced by this pipeline; the final composited
figures are the published figure files.
