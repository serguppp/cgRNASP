"""Copy into your Python package; this module only invokes the compiled scorer."""

import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Literal


@dataclass(frozen=True)
class CgrnaspMutationResult:
    total_energy: float
    pairs_file: Path
    variants_file: Path
    run_file: Path


def run_cgrnasp_mutations(
    reference_path: Path,
    contacts_file: Path,
    output_dir: Path,
    *,
    executable: Path,
    data_dir: Path,
    sequence_separation: Literal["pdb", "ordinal"] = "pdb",
) -> CgrnaspMutationResult:
    """Score all 16 fixed-geometry variants for every supplied native contact.

    Prepare a single-model, canonical RNA reference before calling this function.
    Noncanonical residues and unresolved alternate locations are rejected by C++.
    """
    reference_path = Path(reference_path).resolve()
    contacts_file = Path(contacts_file).resolve()
    output_dir = Path(output_dir).resolve()
    executable = Path(executable).resolve()
    data_dir = Path(data_dir).resolve()

    for path in (reference_path, contacts_file, executable):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not data_dir.is_dir():
        raise NotADirectoryError(data_dir)
    if sequence_separation not in ("pdb", "ordinal"):
        raise ValueError("sequence_separation must be pdb or ordinal")
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise FileExistsError(f"Output directory must be new or empty: {output_dir}")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(executable),
        "mutate",
        str(reference_path),
        "--contacts",
        str(contacts_file),
        "--data-dir",
        str(data_dir),
        "--output-dir",
        str(output_dir),
        "--model-index",
        "0",
        "--sequence-separation",
        sequence_separation,
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"cgRNASP failed ({completed.returncode}): "
            f"{completed.stderr.strip() or completed.stdout.strip()}"
        )

    run_file = output_dir / "run.json"
    metadata = json.loads(run_file.read_text(encoding="utf-8"))
    total_energy = float(metadata["E_total"])
    if not math.isfinite(total_energy):
        raise ValueError("cgRNASP returned a nonfinite native energy")
    if metadata["program"] != "cgRNASP_frustration" or metadata["mode"] != "mutate":
        raise ValueError("Unexpected scorer output")
    if metadata["variant_count"] != 16 * metadata["contact_count"]:
        raise ValueError("Incomplete cgRNASP mutation scan")

    result = CgrnaspMutationResult(
        total_energy=total_energy,
        pairs_file=output_dir / "pairs.csv",
        variants_file=output_dir / "contact_variants.csv",
        run_file=run_file,
    )
    for path in (result.pairs_file, result.variants_file, result.run_file):
        if not path.is_file():
            raise FileNotFoundError(path)
    return result
