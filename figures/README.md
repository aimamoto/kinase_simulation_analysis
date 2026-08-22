# Manuscript figure-generation scripts

Scripts that render the main-text and supplementary figures for the
AlloQuant / trans-allosteric kinase manuscript. Figures are regenerated
for presentation from AlloQuant's computed outputs using a shared color
palette and type scale (`palette.py`); the underlying values are those
produced directly by the pipeline.

These are the **release-adapted** copies: they resolve their input
directories from environment variables (below), import `palette.py`
relative to their own location, and write to `$ALLOQUANT_FIGDIR`, so
they run from a clean checkout without writing into it. The working
copies used while drafting the manuscript carry hardcoded absolute paths
instead.

## Script → figure

| Script | Figure |
|---|---|
| `fig2_model1.py` | Figure 2 — CDK1–CCNB1 thawing / re-tensioning |
| `fig3_cdk1_states.py` | Figure 3 — CDK1 activation structural series |
| `fig4_model2.py` | Figure 4 — CSK–SRC handshake |
| `fig_coupling.py` | Figure 5 — CSK/SRC spine coupling |
| `compose_fig6.py` | Figure 6 — CSK–SRC structural composite |
| `figSI_cage.py` | Figure S1 — P-loop pT14 phosphate clash |
| `figs2alt.py` | Figure S2 — monomeric CSK, apo vs holo |
| `figSI_activation.py` | Figure S4 — categorical activation status |
| `compose_si_model2.py` | Figure S5 — Model 2 dynamics (RMSF, C-tail) |
| `figSI_confidence.py` | Figure S6 — MAC confidence control |
| `figS7_src_state_volcanos.py` | Figure S7 — SRC state volcano panels |

Shared: `palette.py` (colors + type scale, imported by every script);
`PALETTE.md` (palette documentation).

Two output filenames still carry an earlier figure numbering, because
the SI was renumbered after they were written: `figSI_activation.py`
writes `figS3_activation.*` for Figure S4, and `figSI_confidence.py`
writes `FigureS5_confidence_control.*` for Figure S6. The filenames were
left alone so they match the deposited figure files.

### Analysis scripts with no figure in the paper

| Script | Backs |
|---|---|
| `figSI_interface.py` | Supplementary Note 8 — the docking interface is unchanged on SRC priming |
| `figSI_monomer.py` | the monomeric-CSK αC-In occupancies quoted in Results and Discussion |

Both render diagnostic panels that were not submitted; they are included
because the numbers they produce are quoted in the text.
`figSI_interface.py` reads the contact/BSA measurement written by
`addons/iface_py419.py` — run that first, or point
`ALLOQUANT_IFACE_JSON` at an existing copy.

### Figures with no generator

Figure 1 (the framework schematic) and Figure S3 (the AlphaFold3
thermodynamic-override diagram) were made by hand and have no script
here or anywhere else.

Figure S5 is a partial case. `compose_si_model2.py` is its generator,
but the deposited PDF was **manually adjusted** afterwards (panel-label
positions) and that adjusted file is the master, as the script's own
in-line note records. Re-running the script reproduces the panels and
approximates the manual label placement; it does not reproduce the
deposited file byte-for-byte.

### Removed

`figSI_substates.py` was deleted at this revision. It rendered a CDK1
rigid/floppy substate figure that was **removed from the manuscript at
v7r3**: the substates it showed were an artefact of the v7r2 `D1_Dist`
anchor, and under v7r3 the two populations merge into a single state. It
corresponded to no current figure; the code remains in git history.

## Version status

The science in these scripts is **v7r3**, matching the manuscript. Every
script here was re-run from a clean output directory against the
deposited pipeline outputs, and **12 of the 13 rendered PNGs came back
byte-identical to the figures in the paper** — the exception being
Figure S5, for the manual-adjustment reason above. PDFs differ only in
their embedded creation timestamp.

## Dependencies

Python 3 with `matplotlib`, `numpy`, `pandas`, `scipy`, `Pillow`, and
`gemmi` (the last only for `figSI_cage.py`, which measures phosphate
distances directly from `model.cif`).

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
| `ALLOQUANT_AF3_ROOT` | parent directory holding the AF3 run dirs (`figSI_monomer.py` only, for the 40-seed monomer runs) |
| `ALLOQUANT_CSK_MAC_CSV` | `csk_mac_by_condition.csv` from `addons/csk_mac_by_condition.py` (`fig4_model2.py` panel A) |
| `ALLOQUANT_MAC_CONFIDENCE_OUT` | `mac_confidence_out_v7r3/` from `addons/mac_confidence_control.py` (`figSI_confidence.py`; defaults to `$ALLOQUANT_SRC/addons/mac_confidence_out_v7r3`) |
| `ALLOQUANT_IFACE_JSON` | `analysis_iface_py419.json` from `addons/iface_py419.py` (`figSI_interface.py`) |

Run each script from this directory (so `palette.py` is importable), e.g.:

```bash
ALLOQUANT_SRC=/path/to/csk-src_output ALLOQUANT_FIGDIR=/tmp/figs python fig4_model2.py
```

Two figures depend on an add-on being run first: Figure 4 needs
`addons/csk_mac_by_condition.py`, and Figure S6 needs
`addons/mac_confidence_control.py`.

## Note on the structural panels

`fig3_cdk1_states.py`, `figSI_cage.py`, `figs2alt.py`,
`compose_fig6.py` and `compose_si_model2.py` compose pre-rendered
ChimeraX images (from `ALLOQUANT_STRUCT`) rather than generating them
from data. Those renders are not produced by this pipeline; the final
composited figures are the published figure files.
