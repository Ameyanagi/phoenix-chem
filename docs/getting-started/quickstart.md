# Quick Start

Get up and running with PHOENIX in 5 minutes.

!!! warning "Safety Notice"
    PHOENIX is for **screening purposes only**. Results indicate potential hazards but must be validated experimentally before making safety decisions.

## Your First Compound

```python
from phoenix import Compound

# Create a compound from SMILES
ethanol = Compound.from_smiles("CCO")

# Basic properties
print(f"Formula: {ethanol.formula}")
print(f"MW: {ethanol.molecular_weight:.2f} g/mol")
print(f"Composition: {ethanol.composition}")
```

**Expected output:**
```
Formula: C2H6O
MW: 46.07 g/mol
Composition: {'C': 2, 'O': 1, 'H': 6}
```

## Thermodynamic Properties

```python
# Enthalpy of formation at 298.15 K (standard state)
hf = ethanol.enthalpy_of_formation
print(f"ΔHf° = {hf.value:.1f} {hf.unit}")

# Standard entropy
s = ethanol.entropy
print(f"S° = {s.value:.1f} {s.unit}")

# Temperature-dependent properties (T in Kelvin)
state = ethanol.thermo_at(T=400)  # at 400 K
print(f"H(400K) = {state.H.value:.1f} kJ/mol")
```

**Expected output:**
```
ΔHf° = -235.1 kJ/mol
S° = 289.9 J/(mol·K)
H(400K) = -227.8 kJ/mol
```

## Hazard Evaluation

```python
# Evaluate hazard for nitrobenzene
nitrobenzene = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")

# Quick hazard indicators
print(f"Oxygen Balance: {nitrobenzene.oxygen_balance:.1f}%")

# Full hazard evaluation
result = nitrobenzene.evaluate_hazard()
print(f"Hazard Class: {result.hazard_class}")
print(f"Max ΔHd: {result.max_decomposition_kJ_mol:.1f} kJ/mol")
print(f"Max ΔHd: {result.max_decomposition_cal_g:.1f} cal/g")
```

**Expected output:**
```
Oxygen Balance: -162.5%
Hazard Class: HIGH
Max ΔHd: -552.2 kJ/mol
Max ΔHd: -1072.0 cal/g
```

### Interpreting Results

| Hazard Class | ΔHd Threshold | Meaning |
|--------------|---------------|---------|
| **HIGH** | < -700 cal/g | Significant decomposition energy; treat as potentially explosive |
| **MEDIUM** | -300 to -700 cal/g | Moderate hazard; requires careful handling |
| **LOW** | > -300 cal/g | Lower hazard but still requires standard safety protocols |

!!! tip "Oxygen Balance (OB%)"
    - **OB% ≈ 0**: Near-optimal oxygen for complete combustion (maximum energy release)
    - **OB% < 0**: Oxygen-deficient (fuel-rich)
    - **OB% > 0**: Oxygen-excess (oxidizer-rich)

    Compounds with OB% between -120% and +80% warrant closer scrutiny.

## Decomposition Analysis

```python
# Analyze decomposition products
decomp = nitrobenzene.max_decomposition()

print(f"ΔHd = {decomp.delta_hd_kJ_mol:.1f} kJ/mol")
print(f"ΔHd = {decomp.delta_hd_cal_g:.1f} cal/g")

print("\nDecomposition Products:")
for product, moles in decomp.products.items():
    if moles > 0.01:
        print(f"  {product}: {moles:.2f} mol")
```

**Expected output:**
```
ΔHd = -552.2 kJ/mol
ΔHd = -1072.0 cal/g

Decomposition Products:
  N2: 0.50 mol
  H2O: 2.00 mol
  C: 6.00 mol
  H2: 0.50 mol
```

## Batch Screening

```python
from phoenix import screen

# Screen multiple compounds
smiles_list = [
    "CCO",                                    # Ethanol
    "CC(=O)C",                                # Acetone
    "c1ccccc1[N+](=O)[O-]",                   # Nitrobenzene
    "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",  # TNT
]

results = screen(smiles_list)

# Access results as DataFrame
df = results.dataframe
print(df[["formula", "hazard_class", "delta_hd_kJ_mol"]])

# Export to CSV
results.to_csv("screening_results.csv")
```

**Expected output:**
```
    formula hazard_class  delta_hd_kJ_mol
0     C2H6O          LOW        -6.712263
1     C3H6O       MEDIUM       -25.122387
2   C6H5NO2         HIGH      -552.160000
3  C7H5N3O6         HIGH     -1294.800000
```

## Chemical Reactions

```python
from phoenix import Reaction, Auto

# Create a reaction with auto-balancing
# SMILES "C" = methane (CH4), "O=O" = O2, "O=C=O" = CO2, "O" = H2O
rxn = Reaction.from_smiles(
    reactants=["C", "O=O"],           # CH4 + O2
    products=["O=C=O", "O"],          # CO2 + H2O
)
rxn.balance()  # Auto-balance the equation

print(f"Reaction: {rxn}")
print(f"ΔH_rxn = {rxn.delta_h:.1f} kJ/mol")
print(f"Exothermic: {rxn.is_exothermic}")
```

**Expected output:**
```
Reaction: CH4 + 2 O2 -> CO2 + 2 H2O
ΔH_rxn = -802.6 kJ/mol
Exothermic: True
```

## Error Handling

PHOENIX raises specific exceptions for invalid inputs:

```python
from phoenix import (
    Compound,
    InvalidSmilesError,
    UnsupportedElementError,
    PhoenixError,
)

# Handle invalid SMILES
try:
    compound = Compound.from_smiles("not-a-smiles")
except InvalidSmilesError as e:
    print(f"Invalid SMILES: {e}")

# Handle unsupported elements (e.g., metals)
try:
    compound = Compound.from_smiles("[Fe]")  # Iron
except UnsupportedElementError as e:
    print(f"Unsupported: {e}")

# Catch any PHOENIX error
try:
    compound = Compound.from_smiles("invalid")
except PhoenixError as e:
    print(f"Error: {e}")
```

**Expected output:**
```
Invalid SMILES: Invalid SMILES string: 'not-a-smiles'
Unsupported: Unsupported elements: Fe
Error: Invalid SMILES string: 'invalid'
```

!!! info "Supported Elements"
    PHOENIX supports: **C, H, N, O, S, P, F, Cl, Br**

    Metals and other elements are not supported and will raise `UnsupportedElementError`.

## Next Steps

- [Core Concepts](../user-guide/core-concepts.md): Understand oxygen balance, ΔHd, and CHETAH criteria
- [Hazard Evaluation](../user-guide/hazard-evaluation.md): Deep dive into hazard screening
- [Error Handling](../user-guide/error-handling.md): Complete exception reference
- [Examples](../examples/index.md): More detailed examples
- [API Reference](../api/index.md): Complete API documentation
