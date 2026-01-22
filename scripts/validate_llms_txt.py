#!/usr/bin/env python3
"""Validate llms.txt file for broken links and required content.

This script checks that:
1. All linked files exist
2. Safety warning is present in the blockquote
3. Format follows llmstxt.org specification

Usage:
    python scripts/validate_llms_txt.py

Exit codes:
    0: All validations passed
    1: Validation failed
"""

import re
import sys
from pathlib import Path


def validate_llms_txt() -> bool:
    """Validate llms.txt file.

    Returns:
        True if all validations pass, False otherwise.
    """
    llms_txt = Path("llms.txt")

    if not llms_txt.exists():
        print("ERROR: llms.txt not found")
        return False

    content = llms_txt.read_text(encoding="utf-8")
    errors = []
    warnings = []

    # Check H1 header
    if not content.startswith("# "):
        errors.append("llms.txt must start with H1 header (# PROJECT_NAME)")

    # Check blockquote with safety warning
    blockquote_match = re.search(r"^> .+$", content, re.MULTILINE)
    if not blockquote_match:
        errors.append("Missing blockquote summary after H1 header")
    else:
        blockquote_text = ""
        for line in content.split("\n"):
            if line.startswith(">"):
                blockquote_text += line[1:].strip() + " "
            elif blockquote_text and not line.startswith(">"):
                break

        safety_keywords = ["screening", "safety", "validation", "experimental"]
        has_safety_warning = any(kw.lower() in blockquote_text.lower() for kw in safety_keywords)
        if not has_safety_warning:
            errors.append(
                "Blockquote must include safety warning "
                "(keywords: screening, safety, validation, experimental)"
            )

    # Check all markdown links exist
    link_pattern = r"\[([^\]]+)\]\(([^)]+)\)"
    links = re.findall(link_pattern, content)

    broken_links = []
    for title, url in links:
        # Skip external URLs
        if url.startswith("http://") or url.startswith("https://"):
            continue

        # Check if file exists
        path = Path(url)
        if not path.exists():
            broken_links.append((title, url))

    if broken_links:
        errors.append(f"Found {len(broken_links)} broken links:")
        for title, url in broken_links:
            errors.append(f"  - [{title}]({url})")

    # Check H2 sections exist
    h2_sections = re.findall(r"^## (.+)$", content, re.MULTILINE)
    if not h2_sections:
        warnings.append("No H2 sections found (expected: Documentation, API Reference, etc.)")

    # Report results
    print("=" * 60)
    print("llms.txt Validation Report")
    print("=" * 60)

    if errors:
        print(f"\n❌ ERRORS ({len(errors)}):")
        for error in errors:
            print(f"  {error}")

    if warnings:
        print(f"\n⚠️  WARNINGS ({len(warnings)}):")
        for warning in warnings:
            print(f"  {warning}")

    if not errors and not warnings:
        print("\n✅ All validations passed!")

    # Summary
    print(f"\nLinks checked: {len(links)}")
    print(f"H2 sections: {len(h2_sections)}")

    return len(errors) == 0


if __name__ == "__main__":
    success = validate_llms_txt()
    sys.exit(0 if success else 1)
