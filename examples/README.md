# PHOENIX Examples

Practical examples for chemical hazard evaluation with PHOENIX.

## Quick Start

```bash
cd examples
python 01_quick_start.py
```

## Examples Overview

| Example | Description | Time |
|---------|-------------|------|
| `01_quick_start.py` | Basic API introduction + NIST data | 2 min |
| `02_energetic_screening.py` | Batch hazard evaluation | 5 min |
| `03_decomposition_analysis.py` | Decomposition products & energy | 5 min |
| `04_reaction_thermodynamics.py` | Reactions + NIST reference table | 5 min |
| `05_solvent_safety.py` | Single-compound safety check | 3 min |

## Data Sources

PHOENIX uses thermodynamic data from:
- **NIST-JANAF Thermochemical Tables** (Chase, 1998)
- **chemicals library** (CalebBell/ChEDL) - includes NIST WebBook data

Examples 01 and 04 demonstrate how to access and verify reference data.

## Running Examples

All examples print results to stdout. No external files required.

```bash
# Run a single example
python examples/01_quick_start.py

# Run all examples
for f in examples/*.py; do python "$f"; done
```

## Safety Warning

These examples use real energetic materials (TNT, RDX, etc.) for
educational purposes. Results are for screening only. Experimental
validation is required before handling energetic materials.
