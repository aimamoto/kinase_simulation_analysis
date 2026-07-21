# Manuscript figure-generation scripts

Scripts that render the main-text and supplementary figures for the
AlloQuant / trans-allosteric kinase manuscript. Figures are regenerated
for presentation from AlloQuant's computed outputs using a shared color
palette and type scale (`palette.py`); the underlying values are those
produced directly by the pipeline.

## Script → figure

| Script | Figure |
|---|---|
| `fig2_model1.py` | Figure 2 — CDK1–CCNB1 thawing / re-tensioning |
| `fig3_cdk1_states.py` | Figure 3 — CDK1 activation structural series |
| `fig4_model2.py` | Figure 4 — CSK–SRC handshake |
| `fig_coupling.py` | Figure 5 — CSK/SRC spine coupling |
| `figSI_cage.py` | Figure S2 — P-loop pT14 phosphate clash |
| `figSI_substates.py` | Figure S1 — CDK1 rigid/floppy substates |

Shared: `palette.py` (colors + type scale, imported by every script);
`PALETTE.md` (palette documentation).

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

`fig3_cdk1_states.py`, `figSI_cage.py`, and `figSI_substates.py` compose
pre-rendered ChimeraX images (from `ALLOQUANT_STRUCT`) rather than
generating them from data. Those renders are not produced by this
pipeline; the final composited figures are the published figure files.
