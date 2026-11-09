# Validation and small benchmark

Validated on 2026-11-09 with GCC/G++ on Windows. Original potential tables and
upstream C source were kept unchanged. The extension was built with `-O3` and
floating-point contraction disabled; no fast-math flags were used.

- The regression suite checks native energies and all 16 alternatives for
  adjacent, distant, reversed and inter-chain pairs against the original C
  program. It also checks weighted pair/component sums, identity controls,
  fixed coordinates, model selection, invalid inputs and the Python adapter.
- A separate end-to-end integration check used 127 native contacts of 2K96
  supplied by the existing BRiQ reference annotation. All 2032 variant energies
  were compared with independent full-structure rescoring by upstream cgRNASP.
- The largest absolute difference was **4.996824714e-7 kBT**. Upstream writes six
  digits after the decimal point; this agrees within its rounding precision.
- Reference total energy for that 2K96 preparation was -291.0030368413926 kBT.
- The existing contact map, contact density and contact fraction plots were
  generated successfully using the Python integration example. The BRiQ-specific
  component plot was intentionally not used with cgRNASP component columns.

## Timing scope

The following measurements include loading/preparing a reference with Gemmi,
starting one native scoring process, its cache and outputs, reading all variant
energies and constructing the existing analysis results. They exclude building
the program, Python imports, determining the contact list, and plotting.

| Input | Contacts | Variant energies | Measured pipeline wall time |
| --- | ---: | ---: | ---: |
| 2K96, model 0 | 127 | 2032 | 0.193, 0.199, 0.200 s |
| Historical local 1RNA input, model 0, canonical RNA | 60 | 960 | 0.245 s |

These are local small-RNA measurements, not a guarantee for other molecules or
hardware. The separate full upstream-oracle validation of all 2032 2K96 variants,
including writing its intermediate PDBs, took about 24 s.

Full 1EHZ is not benchmarked: it includes modified residues (e.g. 2MG) that are
rejected by default. No chemically unjustified automatic mapping to A/C/G/U was
introduced just to obtain a timing. A protocol for modifications is required
before interpreting an analysis of that complete structure.

This validates numerical implementation of fixed-geometry cgRNASP scoring,
not biological interpretation or frustration-class thresholds for RNA.
