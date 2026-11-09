"""Build the extension and, optionally, the unchanged upstream regression oracle."""

import argparse
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "build")
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--cc", default=os.environ.get("CC", "gcc"))
    parser.add_argument("--with-upstream", action="store_true")
    args = parser.parse_args()
    target = args.output_dir.resolve()
    target.mkdir(parents=True, exist_ok=True)
    version = subprocess.run(
        ["git", "describe", "--always", "--dirty"],
        cwd=ROOT, text=True, capture_output=True, check=True,
    ).stdout.strip()
    changed = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=ROOT, text=True, capture_output=True, check=True,
    ).stdout.strip()
    if changed and not version.endswith("-dirty"):
        version += "-dirty"
    suffix = ".exe" if os.name == "nt" else ""
    command = [
        args.cxx, "-std=c++17", "-O3", "-ffp-contract=off", "-Wall", "-Wextra",
        "-Wpedantic", "-Werror", f"-DCGRNASP_SOURCE_VERSION={json.dumps(version)}",
        str(ROOT / "src/cgrnasp_frustration.cpp"),
        "-o", str(target / f"cgRNASP_frustration{suffix}"),
    ]
    subprocess.run(command, check=True)
    if args.with_upstream:
        subprocess.run([
            args.cc, "-O3", "-ffp-contract=off", str(ROOT / "cgRNASP/cgRNASP.c"),
            "-lm", "-o", str(target / f"cgRNASP_upstream{suffix}"),
        ], check=True)
    print(f"Built in {target} (source {version})")


if __name__ == "__main__":
    main()
