# PHOENIX

<p align="center">
  <img src="assets/logo.png" alt="PHOENIX logo" width="180">
</p>

[![CI](https://github.com/Ameyanagi/phoenix-chem/actions/workflows/ci.yml/badge.svg)](https://github.com/Ameyanagi/phoenix-chem/actions/workflows/ci.yml)
[![Documentation](https://github.com/Ameyanagi/phoenix-chem/actions/workflows/docs.yml/badge.svg)](https://github.com/Ameyanagi/phoenix-chem/actions/workflows/docs.yml)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Physicochemical Hazard Observation & Energetics Numerical Indexing eXpert**

A Python library for reactive chemical hazard evaluation, successor to ASTM CHETAH.

> **Warning**: PHOENIX is for **screening purposes only**. Results must not be used as the sole basis for safety decisions. Experimental validation is required before handling energetic materials.

## Documentation

**[Read the full documentation](https://Ameyanagi.github.io/phoenix-chem/)**

- [Installation](https://Ameyanagi.github.io/phoenix-chem/getting-started/installation/)
- [Quick Start](https://Ameyanagi.github.io/phoenix-chem/getting-started/quickstart/)
- [User Guide](https://Ameyanagi.github.io/phoenix-chem/user-guide/core-concepts/)
- [API Reference](https://Ameyanagi.github.io/phoenix-chem/api/)

## Features

- **Thermodynamic Estimation**: Calculate ΔHf°, S°, Cp(T) using Benson Group Additivity
- **Hazard Classification**: CHETAH-style criteria for instability screening
- **Decomposition Analysis**: Maximum heat of decomposition (ΔHd) calculation
- **Batch Processing**: Screen large SMILES datasets efficiently
- **Reaction Balancing**: Auto-balance chemical reactions with thermodynamics

## Installation

```bash
pip install phoenix-chem
```

Or with uv:

```bash
uv pip install phoenix-chem
```

## Quick Example

```python
from phoenix import Compound

# Create a compound from SMILES
compound = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")  # TNT

# Get properties
print(f"Formula: {compound.formula}")
print(f"Molecular Weight: {compound.molecular_weight:.2f} g/mol")
print(f"Oxygen Balance: {compound.oxygen_balance:.1f}%")

# Evaluate hazard
result = compound.evaluate_hazard()
print(f"Hazard Class: {result.hazard_class}")
print(f"Max ΔHd: {result.max_decomposition_cal_g:.1f} cal/g")
```

## Batch Screening

```python
from phoenix import screen

smiles_list = ["CCO", "c1ccccc1[N+](=O)[O-]", "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]"]
results = screen(smiles_list)

# Get DataFrame
df = results.dataframe
print(df[["formula", "hazard_class", "delta_hd_cal_g"]])

# Export
results.to_csv("screening_results.csv")
```

## Supported Elements

PHOENIX supports compounds containing: **C, H, N, O, S, P, F, Cl, Br**

## LLM-Friendly Documentation

PHOENIX provides [llms.txt](llms.txt) for LLM consumption following the [llmstxt.org](https://llmstxt.org) specification.

## License

MIT License - see [LICENSE](LICENSE) for details.

## References

- ASTM E659: Standard Test Method for CHETAH
- Benson, S.W. "Thermochemical Kinetics" (Wiley, 1976)
- Meyer et al., "Explosives" 6th Ed. (Wiley-VCH, 2007)
