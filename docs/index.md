# PHOENIX

**Physicochemical Hazard Observation & Energetics Numerical Indexing eXpert**

A Python library for reactive chemical hazard evaluation.

!!! warning "Screening Tool Only"
    PHOENIX is for **screening purposes only**. Results must not be used as the sole basis for safety decisions. Experimental validation is required before handling energetic materials. Consult qualified safety professionals and relevant regulations.

## Features

- **Thermodynamic Estimation**: Calculate ΔHf°, S°, Cp(T) using Benson Group Additivity
- **Hazard Classification**: CHETAH-style criteria for instability screening
- **Decomposition Analysis**: Maximum heat of decomposition (ΔHd) calculation
- **Batch Processing**: Screen large SMILES datasets efficiently
- **NIST Reference Data**: Access thermodynamic data from NIST-JANAF tables

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
print(f"Max ΔHd: {result.max_decomposition_kJ_mol:.1f} kJ/mol")
```

## Supported Elements

PHOENIX supports compounds containing: **C, H, N, O, S, P, F, Cl, Br**

Metals and other elements are not supported. See [Limitations](reference/limitations.md) for details.

## Next Steps

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } **Installation**

    ---

    Install PHOENIX with pip or uv

    [:octicons-arrow-right-24: Getting Started](getting-started/installation.md)

-   :material-rocket-launch:{ .lg .middle } **Quick Start**

    ---

    Get up and running in 5 minutes

    [:octicons-arrow-right-24: Quick Start](getting-started/quickstart.md)

-   :material-book-open-variant:{ .lg .middle } **User Guide**

    ---

    Learn core concepts and workflows

    [:octicons-arrow-right-24: User Guide](user-guide/core-concepts.md)

-   :material-api:{ .lg .middle } **API Reference**

    ---

    Complete API documentation

    [:octicons-arrow-right-24: API Reference](api/index.md)

</div>
