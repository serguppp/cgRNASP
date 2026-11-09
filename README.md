# cgRNASP ? fixed-geometry frustration extension

This fork adds fast residue-pair exports and all 16 A/C/G/U substitutions for
native contacts. It reuses the original three-bead potentials and energy weights.
The upstream programs and data remain unchanged.

```sh
python build.py --with-upstream
python -m unittest discover -s tests -v
```

```sh
./build/cgRNASP_frustration mutate input.pdb --contacts contacts.csv \
  --data-dir ./cgRNASP/data --output-dir ./test-output/mutations
```

On Windows use `build/cgRNASP_frustration.exe`.

The scorer keeps P/C4'/N9-or-N1 coordinates fixed and updates nucleotide types in
memory. `E_total` includes the effects on both residues' interactions with their
environment. No atomistic mutation or relaxation is performed.

See [usage, file formats and input choices](docs/fixed-geometry.md) and the
[Python adapter](examples/python/cgrnasp.py). Outputs are `reference_cg.pdb`,
`pairs.csv`, `energy_components.csv`, `run.json`, and, in mutation mode,
`contact_variants.csv`.

## Upstream project

cgRNASP: coarse-grained statistical potentials with residue separation for RNA structure evaluation.

We developed our residue-separation-based CG statistical potentials (cgRNASP) at several CG levels: 

(i) three CG beads at P, C4’, and N9 atoms for purine (or N1 atom for pyrimidine); 

(ii) two CG beads at P and C4’ atoms;

(iii) one CG bead on C4’ atom. 

Correspondingly, for simplicity, our CG statistical potentials were named as cgRNASP, cgRNASP-PC, and cgRNASP-C, respectively. The first potential was named in such way since the potential of 3-bead representation was regarded as a representative one of the CG potentials.

Please follow the readmes in the corresponding subfolders for using the respective statistical potentials. And the potentials were set to a high value of 30 kBT for atom pairs within 3 Å to exclude the structure candidates with severe steric conflicts; updated on August 8, 2024.


Please feel free to contact us if you have any question!
zjtan@whu.edu.cn

Reference:                                      
[1] Tan YL, Wang X, Yu S, Zhang B, & Tan ZJ. 2023. cgRNASP: coarse-grained statistical potentials with residue separation for RNA structure evaluation. NAR Genom Bioinform. 5(1): lqad016.
