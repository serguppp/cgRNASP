"""Numerical regression against upstream plus CLI/input validation tests."""

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUFFIX = ".exe" if os.name == "nt" else ""
BINARY = ROOT / "build" / f"cgRNASP_frustration{SUFFIX}"
UPSTREAM = ROOT / "build" / f"cgRNASP_upstream{SUFFIX}"
DATA = ROOT / "cgRNASP/data"
FIXTURE = ROOT / "tests/fixtures/1NC0_model0_cg.pdb"
BASES = "ACGU"


def rows(path):
    with path.open(newline="", encoding="utf-8") as file:
        return list(csv.DictReader(file))


def residue_key(line):
    return line[21].strip(), str(int(line[22:26])) + line[26].strip()


def mutated(text, changes):
    lines = []
    for line in text.splitlines():
        if line.startswith("ATOM  ") and residue_key(line) in changes:
            base = changes[residue_key(line)]
            atom = line[12:16].strip()
            if atom in ("N1", "N9"):
                atom = "N9" if base in "AG" else "N1"
            line = line[:12] + f" {atom:<3}" + line[16:17] + f"{base:>3}" + line[20:]
        lines.append(line)
    return "\n".join(lines) + "\n"


class ScoringTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cgrnasp-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.counter = 0

    def invoke(self, mode="score", source=FIXTURE, contacts=None, extra=(), success=True):
        self.counter += 1
        destination = self.root / f"result-{self.counter}"
        command = [str(BINARY), mode, str(source), "--data-dir", str(DATA), "--output-dir", str(destination)]
        if contacts is not None:
            command += ["--contacts", str(contacts)]
        process = subprocess.run([*command, *extra], capture_output=True, text=True)
        if success:
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertTrue((destination / "run.json").is_file())
        else:
            self.assertNotEqual(process.returncode, 0)
            self.assertFalse(destination.exists(), process.stderr)
        return destination, process

    def contact_file(self, contacts):
        path = self.root / "contacts.csv"
        with path.open("w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerow(["chain_i", "residue_i", "chain_j", "residue_j", "contact_type"])
            for record in contacts:
                writer.writerow(record)
        return path

    def oracle(self, pdbs):
        self.counter += 1
        directory = self.root / f"oracle-{self.counter}"
        directory.mkdir()
        for name, text in pdbs.items():
            (directory / name).write_text(text, encoding="ascii")
        output = directory / "energies.txt"
        process = subprocess.run([str(UPSTREAM), str(directory), str(len(pdbs)), str(output)],
                                 cwd=ROOT / "cgRNASP", text=True, capture_output=True)
        # The original, unchanged program deliberately returns 6 after success.
        self.assertEqual(process.returncode, 6, process.stderr)
        self.assertIn("spendtime", process.stdout)
        energies = {line.split()[0]: float(line.split()[1]) for line in output.read_text().splitlines()}
        self.assertEqual(set(energies), set(pdbs))
        return energies

    def test_automatic_contacts_cover_the_scoring_graph_and_match_explicit_scan(self):
        scored, _ = self.invoke()
        native_pairs = rows(scored / "pairs.csv")
        automatic, _ = self.invoke("mutate")
        run = json.loads((automatic / "run.json").read_text())
        self.assertEqual(run["contact_selection"], "cgrnasp_distance_ranges")
        self.assertIsNone(run["contacts_file"])
        self.assertEqual(run["contact_count"], run["cached_pair_count"])
        self.assertEqual(run["contact_count"], len(native_pairs))
        self.assertTrue(all(r["contact_type"] == "cgrnasp_range" for r in rows(automatic / "pairs.csv")))
        self.assertEqual(len(rows(automatic / "contact_variants.csv")), 16 * len(native_pairs))
        explicit, _ = self.invoke("mutate", contacts=scored / "pairs.csv")
        self.assertEqual((automatic / "contact_variants.csv").read_bytes(),
                         (explicit / "contact_variants.csv").read_bytes())

    def test_automatic_contacts_use_geometry_even_with_zero_potential(self):
        data = self.root / "zero-data"
        data.mkdir()
        for source in DATA.iterdir():
            if source.suffix == ".potential":
                lines = [" ".join(line.split()[:-1] + ["0"]) for line in source.read_text().splitlines() if line.strip()]
                (data / source.name).write_text("\n".join(lines) + "\n")
            else:
                (data / source.name).write_bytes(source.read_bytes())
        original, _ = self.invoke()
        zero, _ = self.invoke("mutate", extra=["--data-dir", str(data)])
        original_pairs = rows(original / "pairs.csv")
        zero_pairs = rows(zero / "pairs.csv")
        pair_key = lambda r: (r["chain_i"], r["residue_i"], r["chain_j"], r["residue_j"])
        self.assertEqual({pair_key(r) for r in original_pairs}, {pair_key(r) for r in zero_pairs})
        self.assertTrue(zero_pairs)
        self.assertTrue(all(float(r["E_pair"]) == 0 for r in zero_pairs))
        self.assertEqual(len(rows(zero / "contact_variants.csv")), 16 * len(zero_pairs))

    def test_automatic_contacts_exclude_geometrically_isolated_residue(self):
        source = self.root / "isolated.pdb"
        lines = []
        for line in FIXTURE.read_text().splitlines():
            if line.startswith("ATOM  ") and residue_key(line) == ("A", "24"):
                line = line[:30] + f"{float(line[30:38]) + 200:8.3f}" + line[38:]
            lines.append(line)
        source.write_text("\n".join(lines) + "\n")
        destination, _ = self.invoke("mutate", source=source)
        selected = rows(destination / "pairs.csv")
        self.assertTrue(selected)
        self.assertFalse(any(r["residue_i"] == "24" or r["residue_j"] == "24" for r in selected))

    def test_python_adapter_defaults_to_native_cgrnasp_contacts(self):
        sys.path.insert(0, str(ROOT / "examples/python"))
        try:
            from cgrnasp import run_cgrnasp_mutations
        finally:
            sys.path.pop(0)
        result = run_cgrnasp_mutations(reference_path=FIXTURE, output_dir=self.root / "auto-python",
                                      executable=BINARY, data_dir=DATA)
        run = json.loads(result.run_file.read_text())
        self.assertIsNone(run["contacts_file"])
        self.assertEqual(run["contact_count"], run["cached_pair_count"])
        self.assertEqual(len(rows(result.variants_file)), 16 * run["contact_count"])

    def test_native_matches_original_program(self):
        destination, _ = self.invoke()
        native = json.loads((destination / "run.json").read_text())["E_total"]
        oracle = self.oracle({"reference.pdb": (destination / "reference_cg.pdb").read_text()})
        self.assertAlmostEqual(native, oracle["reference.pdb"], delta=6e-7)

    def test_pair_and_component_sums(self):
        destination, _ = self.invoke()
        run = json.loads((destination / "run.json").read_text())
        pair_rows = rows(destination / "pairs.csv")
        self.assertAlmostEqual(sum(float(row["E_pair"]) for row in pair_rows), run["E_total"], delta=1e-9)
        for row in pair_rows:
            self.assertAlmostEqual(sum(float(row[name]) for name in ["E_seq1", "E_seq2", "E_seq3_4", "E_long"]), float(row["E_pair"]), delta=1e-12)
        self.assertAlmostEqual(sum(run["energy_components"].values()), run["E_total"], delta=1e-12)
        self.assertEqual(len(rows(destination / "energy_components.csv")), 1)

    def test_every_variant_matches_full_upstream_rescoring(self):
        contacts = self.contact_file([("A", "1", "A", "2", "test"),
                                      ("A", "3", "A", "11", "test"),
                                      ("A", "20", "A", "5", "reversed")])
        destination, _ = self.invoke("mutate", contacts=contacts)
        variants = rows(destination / "contact_variants.csv")
        self.assertEqual(len(variants), 48)
        reference = (destination / "reference_cg.pdb").read_text()
        pdbs = {}
        for index, row in enumerate(variants):
            changes = {(row["chain_i"], row["residue_i"]): row["base_i"],
                       (row["chain_j"], row["residue_j"]): row["base_j"]}
            pdbs[f"variant_{index:04d}.pdb"] = mutated(reference, changes)
        oracle = self.oracle(pdbs)
        for index, row in enumerate(variants):
            self.assertAlmostEqual(float(row["E_total"]), oracle[f"variant_{index:04d}.pdb"], delta=6e-7)
            self.assertAlmostEqual(float(row["E_total"]) - float(row["native_E_total"]), float(row["delta_E_total"]), delta=1e-12)
        for row in rows(destination / "energy_components.csv"):
            self.assertAlmostEqual(sum(float(row[name]) for name in ["E_seq1", "E_seq2", "E_seq3_4", "E_long"]), float(row["E_total"]), delta=1e-10)
        self.assertEqual(len(rows(destination / "energy_components.csv")), 49)

    def test_native_combination_has_exactly_zero_delta_and_all_16_states(self):
        contacts = self.contact_file([("A", "1", "A", "2", "test")])
        destination, _ = self.invoke("mutate", contacts=contacts)
        variants = rows(destination / "contact_variants.csv")
        self.assertEqual({(row["base_i"], row["base_j"]) for row in variants}, set(product(BASES, repeat=2)))
        native = {residue_key(line): line[19] for line in FIXTURE.read_text().splitlines() if line.startswith("ATOM  ")}
        identity = next(row for row in variants if row["base_i"] == native[("A", "1")] and row["base_j"] == native[("A", "2")])
        self.assertEqual(float(identity["delta_E_total"]), 0.0)
        self.assertEqual(identity["E_total"], identity["native_E_total"])

    def test_cross_chain_energy_and_variant_match_upstream(self):
        lines = []
        for line in FIXTURE.read_text().splitlines():
            if line.startswith("ATOM  ") and int(line[22:26]) > 12:
                line = line[:21] + "B" + f"{int(line[22:26]) - 12:4d}" + line[26:]
            lines.append(line)
        source = self.root / "two_chains.pdb"
        source.write_text("\n".join(lines) + "\n")
        contacts = self.contact_file([("A", "3", "B", "3", "cross_chain")])
        destination, _ = self.invoke("mutate", source=source, contacts=contacts)
        reference = (destination / "reference_cg.pdb").read_text()
        variants = rows(destination / "contact_variants.csv")
        oracle = self.oracle({f"variant_{i:04d}.pdb": mutated(reference, {("A", "3"): row["base_i"], ("B", "3"): row["base_j"]}) for i, row in enumerate(variants)})
        for i, row in enumerate(variants):
            self.assertAlmostEqual(float(row["E_total"]), oracle[f"variant_{i:04d}.pdb"], delta=6e-7)

    def test_geometry_and_bead_count_stay_fixed_under_purine_pyrimidine_changes(self):
        text = FIXTURE.read_text()
        expected = [line[30:54] for line in text.splitlines() if line.startswith("ATOM  ")]
        for a, b in product(BASES, repeat=2):
            result = mutated(text, {("A", "1"): a, ("A", "2"): b})
            self.assertEqual([line[30:54] for line in result.splitlines() if line.startswith("ATOM  ")], expected)

    def test_models_are_selected_explicitly(self):
        text = FIXTURE.read_text().replace("END\n", "")
        alternative = mutated(text, {("A", "1"): "C"})
        source = self.root / "ensemble.pdb"
        source.write_text("MODEL        1\n" + text + "ENDMDL\nMODEL        2\n" + alternative + "ENDMDL\nEND\n")
        original, _ = self.invoke()
        first, _ = self.invoke(source=source)
        second, _ = self.invoke(source=source, extra=["--model-index", "1"])
        self.assertEqual(json.loads((original / "run.json").read_text())["E_total"], json.loads((first / "run.json").read_text())["E_total"])
        self.assertEqual(json.loads((first / "run.json").read_text())["input_model_count"], 2)
        expected = self.oracle({"reference.pdb": (second / "reference_cg.pdb").read_text()})
        self.assertAlmostEqual(json.loads((second / "run.json").read_text())["E_total"], expected["reference.pdb"], delta=6e-7)
        self.invoke(source=source, extra=["--model-index", "2"], success=False)

    def test_explicit_contact_labels_do_not_control_selection(self):
        source = self.root / "distant.pdb"
        text = FIXTURE.read_text()
        for line in text.splitlines():
            if line.startswith("ATOM  ") and residue_key(line) == ("A", "24"):
                text = text.replace(line, line[:30] + f"{float(line[30:38]) + 100:8.3f}" + line[38:])
        source.write_text(text)
        contacts = self.contact_file([("A", "1", "A", "24", "distant"), ("A", "2", "A", "3", "unannotated")])
        destination, _ = self.invoke("mutate", source=source, contacts=contacts)
        selected = [row for row in rows(destination / "pairs.csv") if row["contact_type"] != "unannotated"]
        self.assertEqual(len(selected), 1)
        self.assertEqual(float(selected[0]["E_pair"]), 0)
        variants = rows(destination / "contact_variants.csv")
        self.assertEqual(len(variants), 32)
        self.assertTrue(any(row["contact_type"] == "unannotated" for row in variants))
        self.assertTrue(any(float(row["delta_E_total"]) != 0 for row in variants))

    def test_insertion_codes_require_explicit_separation(self):
        source = self.root / "insertions.pdb"
        source.write_text("\n".join(line[:26] + "A" + line[27:] if line.startswith("ATOM  ") and residue_key(line) == ("A", "2") else line for line in FIXTURE.read_text().splitlines()) + "\n")
        self.invoke(source=source, success=False)
        destination, _ = self.invoke(source=source, extra=["--sequence-separation", "ordinal"])
        self.assertTrue(any(row["residue_i"] == "2A" or row["residue_j"] == "2A" for row in rows(destination / "pairs.csv")))

    def test_missing_terminal_phosphate_is_reported(self):
        destination, _ = self.invoke()
        self.assertEqual(json.loads((destination / "run.json").read_text())["missing_terminal_phosphates"], ["A:1"])

    def test_missing_internal_phosphate_or_base_atom_is_rejected(self):
        base_atom = next(line[12:16].strip() for line in FIXTURE.read_text().splitlines()
                         if line.startswith("ATOM  ") and residue_key(line) == ("A", "2") and line[12:16].strip() in ("N1", "N9"))
        for atom in ["P", "C4'", base_atom]:
            source = self.root / "missing.pdb"
            source.write_text("\n".join(line for line in FIXTURE.read_text().splitlines() if not (line.startswith("ATOM  ") and line[12:16].strip() == atom and residue_key(line) == ("A", "2"))) + "\n")
            self.invoke(source=source, success=False)

    def test_unsupported_residues_are_not_silently_skipped(self):
        source = self.root / "unsupported.pdb"
        extra = "HETATM 9999  O   HOH Z 999       1.000   2.000   3.000  1.00  0.00           O\n"
        source.write_text(FIXTURE.read_text().replace("END\n", extra + "END\n"))
        self.invoke(source=source, success=False)
        destination, _ = self.invoke(source=source, extra=["--skip-unsupported-residues"])
        self.assertEqual(json.loads((destination / "run.json").read_text())["skipped_residue_names"], ["HOH"])

    def test_duplicate_or_alternate_atoms_are_rejected(self):
        lines = FIXTURE.read_text().splitlines()
        for text in [lines[0] + "\n" + FIXTURE.read_text(), lines[0][:16] + "A" + lines[0][17:] + "\n" + "\n".join(lines[1:])]:
            source = self.root / "ambiguous.pdb"
            source.write_text(text)
            self.invoke(source=source, success=False)

    def test_bad_contacts_are_rejected(self):
        for contacts in [
            [("A", "1", "A", "1", "test")],
            [("A", "1", "A", "9999", "test")],
            [("A", "1", "A", "2", "test"), ("A", "2", "A", "1", "test")],
        ]:
            self.invoke("mutate", contacts=self.contact_file(contacts), success=False)

    def test_nonempty_output_is_never_overwritten(self):
        destination, _ = self.invoke()
        previous = (destination / "run.json").read_bytes()
        process = subprocess.run([str(BINARY), "score", str(FIXTURE), "--data-dir", str(DATA), "--output-dir", str(destination)], capture_output=True)
        self.assertNotEqual(process.returncode, 0)
        self.assertEqual(previous, (destination / "run.json").read_bytes())

    def test_truncated_tables_fail_before_output_is_created(self):
        data = self.root / "bad-data"
        data.mkdir()
        for source in DATA.iterdir():
            (data / source.name).write_bytes(source.read_bytes())
        (data / "long-ranged.potential").write_text("0 0 0 30\n")
        self.invoke(extra=["--data-dir", str(data)], success=False)

    def test_blank_chain_csv_identifiers(self):
        source = self.root / "blank_chain.pdb"
        source.write_text("\n".join(line[:21] + " " + line[22:] if line.startswith("ATOM  ") else line for line in FIXTURE.read_text().splitlines()) + "\n")
        destination, _ = self.invoke("mutate", source=source, contacts=self.contact_file([("", "1", "", "2", "test")]))
        self.assertEqual(len(rows(destination / "contact_variants.csv")), 16)

    def test_full_atom_input_ignores_non_cg_atoms_without_changing_geometry(self):
        text = FIXTURE.read_text()
        first = text.splitlines()[0]
        extra = first[:12] + " C1'" + first[16:]
        source = self.root / "extra_atoms.pdb"
        source.write_text(first + "\n" + extra + "\n" + "\n".join(text.splitlines()[1:]) + "\n")
        original, _ = self.invoke()
        full, _ = self.invoke(source=source)
        self.assertEqual(json.loads((original / "run.json").read_text())["E_total"], json.loads((full / "run.json").read_text())["E_total"])

    def test_python_adapter(self):
        sys.path.insert(0, str(ROOT / "examples/python"))
        try:
            from cgrnasp import run_cgrnasp_mutations
        finally:
            sys.path.pop(0)
        result = run_cgrnasp_mutations(
            reference_path=FIXTURE,
            contacts_file=self.contact_file([("A", "1", "A", "2", "test")]),
            output_dir=self.root / "python output with spaces",
            executable=BINARY,
            data_dir=DATA,
        )
        self.assertEqual(len(rows(result.variants_file)), 16)
        self.assertEqual(result.total_energy, json.loads(result.run_file.read_text())["E_total"])
        self.assertTrue(result.pairs_file.is_file())


if __name__ == "__main__":
    unittest.main()
