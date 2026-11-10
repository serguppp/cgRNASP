# Fixed-geometry pair and mutation scoring

This fork adds `cgRNASP_frustration`. The original C programs and all potential
tables remain unchanged. The added implementation follows the three-bead model,
distance bins and energy weights of the upstream code at
`ff8c5a23727b725e636d0b52dead2cecdb624a87`.

## Build

Requires Python 3 and GCC/G++ with C++17 support. From the repository root:

```sh
python build.py --with-upstream
python -m unittest discover -s tests -v
```

On Windows the executable is `build/cgRNASP_frustration.exe`; on Linux it is
`build/cgRNASP_frustration`. `--with-upstream` additionally compiles the unchanged
original program as a numerical regression oracle, not as a runtime dependency.
The extension exits with 0 on success and 1 on error. The original program still
returns its historical exit code 6.

## Reference scoring

```sh
./build/cgRNASP_frustration score input.pdb \
  --data-dir ./cgRNASP/data --output-dir ./test-output/reference
```

The output directory must be new or empty. Files:

| File | Contents |
| --- | --- |
| `reference_cg.pdb` | Selected model reduced to P, C4' and N9/N1 |
| `pairs.csv` | Weighted residue-pair energies and optional supplied annotations |
| `energy_components.csv` | Weighted reference energy components |
| `run.json` | Version, options, table fingerprints, omissions and native energy |

`run.json` is written last as the successful completion marker. A failed write
can leave partial outputs; rerun in a fresh directory.

Pairs outside every applicable distance bin are exactly zero for every identity
and omitted from the cache/export. Explicitly requested contacts are exported
even if their direct pair energy is zero. Therefore the sum of all `E_pair`
entries still reproduces the native `E_total`, up to summation roundoff.
`score --contacts contacts.csv ...` can annotate selected reference pairs without
generating variants. An energetic pair is not automatically a base pair or stack.

## Contact selection (default)

With no `--contacts`, both `score` and `mutate` select the native cgRNASP
geometric pair graph. A residue pair is selected when at least one of its
present bead pairs has a valid distance bin in the applicable upstream table.
There is no external annotation dependency and no filter on native pair energy.
Zero native energy, including cancellation, does not exclude a pair.

This means model interactions, not necessarily base pairs or stacks. The upper
range depends on sequence separation: approximately 5.1 A (1), 9 A (2), 12.9 A
(3-4), and 24 A (>4 or inter-chain). Exact bin tests are shared with scoring.
`contact_type` is `cgrnasp_range`; it is not a biological interaction annotation.
Adjacent residues and longer-range contacts are included. Separate them when
interpreting or visualising results rather than applying a hidden filter.

`run.json` records `contact_selection` and `contact_definition`. Default selection
has `contacts_file: null`. Each selected pair receives all 16 variants in one
process, using the existing cache and unchanged potential tables.

## Contact mutation scan

```sh
./build/cgRNASP_frustration mutate input.pdb \
  --data-dir ./cgRNASP/data \
  --output-dir ./test-output/mutations
```

To deliberately restrict the scan, add `--contacts contacts.csv`. Minimal file:

```csv
chain_i,residue_i,chain_j,residue_j,contact_type
A,1,A,2,example
A,3,A,11,example
```

The `example` rows demonstrate the format, not a biological annotation.
An explicit CSV replaces the default cgRNASP pair selection. Use it only for
a deliberately restricted or matched comparison; no external annotation is
needed for the default scan. The four identifier columns are required;
`contact_type` is optional and defaults to `contact`. If it is provided, rows
are all evaluated regardless of their annotation value. Empty labels default to
`contact`. Other columns are ignored.
Self-contacts, duplicates (including reversed duplicates) and absent identifiers
are rejected. CSV orientation and annotation are preserved in variant output.
Blank chain IDs are supported. Multiline CSV fields are not supported.

In addition to the reference files, `mutate` writes `contact_variants.csv`:

```text
chain_i,residue_i,chain_j,residue_j,contact_type,base_i,base_j,
native_E_total,E_total,delta_E_total,native_E_pair,E_pair
```

Each contact has exactly 16 rows in `ACGU x ACGU` order. `E_total` includes the
effects of both substitutions on all their interactions with the environment.
`E_pair` is only the direct interaction of the two selected positions; it is
provided for diagnostics, not as a substitute for the total-energy protocol.

`energy_components.csv` contains the reference plus all variants, identified by
`state`, the four contact columns and `base_i,base_j`. These components differ
from per-pair energies: they are the four weighted terms of the entire structure.

## Energy and exact incremental calculation

For every residue pair, the energy sums up to nine bead-bead contributions.
The table is selected by residue separation 1, 2, 3-4, or >4/inter-chain.
As in the original program, distance bins have width 0.3 A and the four table
sizes are 17, 30, 43 and 80. The original integer-bin tests, including their
precise upper limits, are retained. Historical table labels are not contiguous
CG type IDs: the original code reads by row order. The extension validates this
layout and uses the same ordering.

The unweighted group sums are combined as:

```text
E_total = E1 + 1.5*E2 + 2.4*E3_4 + 6.5*E_long/f(L)
f(L) = -355/sqrt(L+16) + 72
```

The CSV columns `E_seq1,E_seq2,E_seq3_4,E_long` are ALREADY weighted, so their
sum is `E_total`. The output energy unit follows upstream: kBT.
Only the original inter-residue terms are present; no additional energy terms
are introduced. Short-distance penalties remain in
the original tables.

The engine computes distance bins once and caches the 4x4 identity matrix for
each geometrically contributing pair. A contact mutation updates the union of
pairs incident on its two positions, counting their shared pair only once.
Adding these differences to the native group sums gives the same total energy
as rescoring the complete mutant, apart from floating-point summation roundoff.
Coordinates, bead count, table assignments and L stay fixed.

All three bead types depend on nucleotide identity: e.g. AP/AC4'/AN9 become
UP/UC4'/UN1. N9-to-N1 is a relabeling of one fixed point, not addition of a point.
No atomistic rebuilding or relaxation is performed. This does not validate the
steric feasibility of a subsequent all-atom reconstruction.

For a sparse geometric pair graph with M edges, cache construction costs
O(L^2 + 16*M + M*log(M)) and O(L+16*M) memory (the lookup uses an ordered map).
Updating a contact costs approximately
O(16*(degree(i)+degree(j))), plus CSV output. There is no subprocess per mutant.

## Input choices and limits

- Canonical A/C/G/U only. Modified nucleotides are not automatically mapped.
- Full-atom and three-bead PDB inputs work; only the model's selected beads
  contribute. Both ATOM and canonical HETATM records are accepted.
- Ligands/water/unsupported residues cause an error by default.
  `--skip-unsupported-residues` is an explicit opt-in: exclusions are recorded
  in `run.json`, and contacts referring to excluded residues still fail.
  Prefer preparing the RNA reference explicitly in your pipeline.
- Model 0 is selected by default. `--model-index N` is zero-based; model count
  and selected model are recorded. No concatenation of NMR models occurs.
- C4' and the native glycosidic N are required. A missing P is allowed only
  at the first observed residue of a chain and reported, without invented
  coordinates. Interior gaps require explicit preparation.
- Alternate locations and duplicate selected atoms are rejected.
- `--sequence-separation pdb` is the default and matches upstream number-based
  separation. Numbering gaps are retained. Insertion codes are rejected in
  this mode because upstream does not distinguish them.
- `--sequence-separation ordinal` explicitly uses observed per-chain order
  and supports insertion codes. It can change the model for missing segments;
  it is not automatically equivalent to true biological sequence positions.
- A table fingerprint is FNV-1a of file bytes, labeled as such, not a security
  checksum. LF/CRLF conversion can change it without changing parsed energies.

## Python integration

`examples/python/cgrnasp.py` is a standard-library-only adapter that can be copied
into the Python package. It accepts a prepared reference PDB, optional contacts CSV and
paths to this binary and its data directory. It returns native energy and paths
to the same kinds of output files used by the existing analysis pipeline.

Keep structure preparation in Gemmi. Bypass ModeRNA and the atomistic mutation
cache on this scoring path. Keep the 16-state population-standard-deviation
frustration formula and the same native contact list. Use the full-atom reference
for structural plots; component plots use only E_seq1, E_seq2, E_seq3_4 and E_long.

The native state occurs exactly once among the 16 combinations and retains its
original geometry. Uniform weights are 1/16. This extension supplies energies;
it does not establish biological interpretation or RNA frustration thresholds.

## Tests

Tests compare native energy and every variant of adjacent, distant, reversed
and inter-chain pairs against the unchanged upstream scorer. They check weighted
component and pair sums, native-identity controls, NMR selection, identifier
handling, fixed coordinates and rejected invalid inputs. The compact fixture
is canonical RNA from model 0 of PDB 1NC0, reduced to CG without rebuilding atoms.

The original authors' software, reference and license are retained in the root
README and LICENSE. This extension is distributed under the same GPL license.
