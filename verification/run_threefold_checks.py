#!/usr/bin/env python3
"""Run all distributed arithmetic and numerical checks. Standard library only.

These checks concern arithmetic and combinatorial consistency, not analytic
existence, positivity, symplectic comparison, or global Lagrangian gluing.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import shutil


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "build" / "threefold-checks"
CHECKS = (
    ("threefold_v11/verify_threefold.py", (
        "threefold_v11/verification_results.json",
        "threefold_v11/invariants_table.tex",
    )),
    ("threefold_v11/verify_relative_models.py", (
        "threefold_v11/relative_model_results.json",
    )),
    ("local_blocks/verify_local_blocks.py", ("local_blocks/local_block_checks.json",)),
    ("local_blocks/verify_relative_gluing.py", ("local_blocks/relative_gluing_checks.json",)),
    ("local_blocks/verify_corner_models.py", ("local_blocks/corner_model_checks.json",)),
    ("local_blocks/verify_positive_end.py", ("local_blocks/positive_end_checks.json",)),
    ("local_blocks/verify_positive_transport.py", ("local_blocks/positive_transport_checks.json",)),
    ("local_blocks/verify_straightening.py", ("local_blocks/straightening_checks.json",)),
    ("local_blocks/verify_ports.py", ("local_blocks/ports_action_checks.json",)),
    ("local_blocks/verify_cap_extensions.py", ("local_blocks/cap_extension_checks.json",)),
    ("local_blocks/verify_boundary_surgery.py", ("local_blocks/boundary_surgery_checks.json",)),
    ("local_blocks/verify_closed_curve_and_families.py", ("local_blocks/closed_curve_family_checks.json",)),
    ("local_blocks/verify_focus_normal_and_regions.py", ("local_blocks/focus_normal_region_checks.json",)),
)


def main():
    for folder in ("local_blocks", "threefold_v11"):
        shutil.copytree(REPO / "papers/sections/threefold" / folder,
                        ROOT / folder, dirs_exist_ok=True)
    records = []
    for relative, outputs in CHECKS:
        script = ROOT / relative
        completed = subprocess.run(
            [sys.executable, str(script)], cwd=script.parent,
            text=True, encoding="utf-8", capture_output=True,
        )
        if completed.returncode:
            print(completed.stdout)
            print(completed.stderr, file=sys.stderr)
            raise SystemExit(f"FAILED: {relative} ({completed.returncode})")
        generated = []
        for output in outputs:
            data = (ROOT / output).read_bytes()
            generated.append({
                "path": output,
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            })
        records.append({
            "script": relative, "status": "passed", "generated": generated,
        })
        print(f"PASS {relative}")
    report = {
        "scope": (
            "Exact arithmetic, combinatorial, and selected numerical checks; "
            "no claim of analytic or global fibration existence."
        ),
        "checks": records,
        "status": "passed",
    }
    (ROOT / "VALIDATION.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8",
    )
    print("All thirteen checks passed. Report: VALIDATION.json")


if __name__ == "__main__":
    main()
