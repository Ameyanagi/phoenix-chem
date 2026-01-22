#!/usr/bin/env python3
"""Generate llms-full.txt by concatenating all documentation files.

This script creates a single file containing all documentation content,
suitable for LLMs that prefer complete context in one file.

Usage:
    python scripts/generate_llms_full.py
"""

from pathlib import Path


def generate_llms_full() -> None:
    """Generate llms-full.txt from all documentation files."""
    docs_dir = Path("docs")
    output_file = Path("llms-full.txt")

    # Define order to concatenate (matches mkdocs.yml nav structure)
    files_order = [
        "index.md",
        "getting-started/installation.md",
        "getting-started/quickstart.md",
        "user-guide/core-concepts.md",
        "user-guide/compounds.md",
        "user-guide/thermodynamics.md",
        "user-guide/hazard-evaluation.md",
        "user-guide/decomposition.md",
        "user-guide/reactions.md",
        "user-guide/batch-processing.md",
        "user-guide/functional-groups.md",
        "user-guide/error-handling.md",
        "api/index.md",
        "api/compound.md",
        "api/reaction.md",
        "api/hazard.md",
        "api/thermo.md",
        "api/batch.md",
        "api/exceptions.md",
        "examples/index.md",
        "reference/chemistry-background.md",
        "reference/data-sources.md",
        "reference/limitations.md",
    ]

    with open(output_file, "w", encoding="utf-8") as out:
        # Write header
        out.write("# PHOENIX - Full Documentation\n\n")
        out.write("> Auto-generated from all documentation files.\n")
        out.write("> For the latest version, visit https://phoenix-chem.readthedocs.io\n\n")

        # Safety warning
        out.write("## ⚠️ SAFETY NOTICE\n\n")
        out.write("PHOENIX is a **screening tool only**. Results:\n\n")
        out.write("- Must be validated experimentally\n")
        out.write("- Do not replace domain expert review\n")
        out.write("- Supported elements: C, H, N, O, S, P, F, Cl, Br only\n")
        out.write("- Accuracy: ±10 kJ/mol for formation enthalpies\n\n")
        out.write("---\n\n")

        # Concatenate files
        files_found = 0
        files_missing = []

        for file_path in files_order:
            full_path = docs_dir / file_path
            if not full_path.exists():
                files_missing.append(str(full_path))
                continue

            files_found += 1
            with open(full_path, encoding="utf-8") as f:
                content = f.read()

            # Add file header
            out.write(f"# File: {file_path}\n\n")
            out.write(content)
            out.write("\n\n---\n\n")

        # Summary
        print(f"Generated {output_file}")
        print(f"  Files included: {files_found}")
        if files_missing:
            print(f"  Files missing: {len(files_missing)}")
            for f in files_missing:
                print(f"    - {f}")


if __name__ == "__main__":
    generate_llms_full()
