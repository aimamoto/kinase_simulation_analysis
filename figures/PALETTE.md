# AlloQuant manuscript — LOCKED figure design system

Extracted from `Fig1.odg` (embedded SVG, Flat-UI palette) so Figs 2–N read as one system.

## Colours
| Role | Hex | Use |
|---|---|---|
| Ink (primary) | `#2C3E50` | text, axes, borders, box outlines |
| Ink-2 | `#34495E` | secondary text / darker fills |
| Muted grey | `#7F8C8D` | annotations, tertiary text |
| Gridline | `#BDC3C7` | recessive grid |
| Light fill | `#ECF0F1` | panel/box backgrounds |
| Blue / dark | `#3498DB` / `#2980B9` | structure, N-lobe, primary data series |
| Red / dark | `#E74C3C` / `#C0392B` | tense / rigid pole, C-lobe features |
| Green / dark | `#27AE60` / `#1E8449` | active state / output / positive |
| Amber | `#F39C12` | nucleotide / intermediate |
| Purple | `#8E44AD` | special landmark / 4th category |

**Sequential "rigidity" ramp** (magnitude of MAC, fluid→rigid): `#EAF2FB → #AED4F0 → #5DAEE3 → #2E86C1 → #21618C → #1B2F45`.

**Diverging** (for correlation/Δ heatmaps): red `#C0392B` ↔ neutral `#ECF0F1` ↔ blue `#2980B9`.

## Type
Sans-serif, Liberation Sans / Arial / Noto Sans (Fig 1 uses these). Ink `#2C3E50`.

## CVD rule (validated, OKLab ΔE under deuteranopia/protanopia)
The full Flat-UI categorical set passes normal vision (min ΔE 18.2) but **red–green fails
deuteranopia (ΔE 4.8)** and green–amber fails protanopia (6.6). Therefore:
- **Magnitude → sequential blue ramp** (never a categorical rainbow).
- **Red = tense/rigid, Green = active** are *semantic poles*, always **direct-labelled**; never
  distinguish two comparable series by red-vs-green alone.
- Categorical series ≤ ~4, drawn from {blue, amber, purple, navy, grey} + red/green-as-semantic,
  and **always labelled** (identity never by colour alone).

## Fonts/tooling
Build with `conda run -n main python` (matplotlib in env `main`). Structure renders via ChimeraX (`/usr/bin/chimerax`). Export vector PDF/SVG for submission + PNG for review.
