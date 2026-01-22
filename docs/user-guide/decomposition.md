# Decomposition Analysis

PHOENIX calculates the maximum heat of decomposition (ΔHd) to predict energy release potential.

## Overview

The maximum heat of decomposition represents the theoretical maximum energy release when a compound decomposes to its most thermodynamically stable products.

$$
\Delta H_d = \Delta H_f(\text{products}) - \Delta H_f(\text{reactant})
$$

Negative values indicate **exothermic** decomposition (energy release).

## Basic Usage

```python
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")  # Nitrobenzene
decomp = compound.max_decomposition()

print(f"ΔHd = {decomp.delta_hd_kJ_mol:.1f} kJ/mol")
print(f"ΔHd = {decomp.delta_hd_cal_g:.1f} cal/g")
```

## DecompositionResult

The `max_decomposition()` method returns a `DecompositionResult`:

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `delta_hd_kJ_mol` | float | ΔHd in kJ/mol (negative = exothermic) |
| `delta_hd_cal_g` | float | ΔHd in cal/g |
| `products` | dict[str, float] | Decomposition products and moles |
| `reactant_hf_kJ_mol` | float | Reactant ΔHf° |
| `products_hf_kJ_mol` | float | Total products ΔHf° |
| `gas_volume_L_g` | float | Gas volume at temperature |
| `gas_moles` | float | Total moles of gas |
| `gas_composition` | dict[str, float] | Gas mole fractions |
| `gas_temperature_K` | float | Temperature for gas calculation |
| `method` | str | 'hierarchy' or 'lp' |

### Example

```python
decomp = compound.max_decomposition()

# Energy release
print(f"ΔHd = {decomp.delta_hd_cal_g:.1f} cal/g")

# Products
print("\nDecomposition products:")
for product, moles in decomp.products.items():
    if moles > 0.01:
        print(f"  {product}: {moles:.3f} mol")

# Gas generation
print(f"\nGas generation: {decomp.gas_volume_L_g:.4f} L/g")
print(f"Gas moles: {decomp.gas_moles:.2f} mol/mol")
```

## Calculation Methods

### Hierarchy Method (Default)

Uses CHETAH analytical priority rules:

```python
decomp = compound.max_decomposition(method="hierarchy")
```

The hierarchy follows thermodynamic priorities:

| Priority | Reaction | Product |
|----------|----------|---------|
| 1 | F + H → HF | Hydrogen fluoride |
| 2 | N → ½N₂ | Nitrogen gas |
| 3 | P + O → ¼P₄O₁₀ | Phosphorus pentoxide |
| 4 | H + O → ½H₂O | Water |
| 5 | C + O → ½CO₂ | Carbon dioxide |
| 6 | S + O → SO₂ | Sulfur dioxide |
| 7 | C + ½O → CO | Carbon monoxide |
| 8 | Cl + H → HCl | Hydrogen chloride |
| 9 | Br + H → HBr | Hydrogen bromide |
| 10 | C → C(s) | Graphite |
| 11 | H → ½H₂ | Hydrogen gas |
| 12 | S → S(s) | Solid sulfur |
| 13 | Halogens → X₂ | Diatomic halogens |

### LP Method

Uses Linear Programming optimization to find the mathematically optimal solution:

```python
decomp = compound.max_decomposition(method="lp")
```

LP minimizes total product enthalpy subject to atom balance constraints.

### Comparing Methods

Use `method="both"` to compare:

```python
comparison = compound.max_decomposition(method="both")

print(f"Hierarchy: {comparison.hierarchy_result.delta_hd_cal_g:.1f} cal/g")
print(f"LP: {comparison.lp_result.delta_hd_cal_g:.1f} cal/g")
print(f"Deviation: {comparison.deviation_percent:.2f}%")
```

The `DecompositionComparison` object provides:

| Attribute | Type | Description |
|-----------|------|-------------|
| `hierarchy_result` | DecompositionResult | Hierarchy method result |
| `lp_result` | DecompositionResult | LP method result |
| `deviation_percent` | float | % deviation between methods |
| `hierarchy_delta_hd` | float | Shortcut to hierarchy ΔHd |
| `lp_delta_hd` | float | Shortcut to LP ΔHd |

### When Methods Differ

Methods typically agree for:
- Simple molecules
- Complete oxidation scenarios

Methods may differ for:
- Complex molecules with many possible products
- Partially oxidized systems

```python
# Example: Check for significant deviation
comparison = compound.max_decomposition(method="both")

if comparison.deviation_percent > 5:
    print("Warning: Methods disagree significantly")
    print(f"  Hierarchy: {comparison.hierarchy_delta_hd:.1f} kJ/mol")
    print(f"  LP: {comparison.lp_delta_hd:.1f} kJ/mol")
```

## Gas Generation

### Basic Gas Data

```python
decomp = compound.max_decomposition()

print(f"Total gas: {decomp.gas_moles:.2f} mol/mol")
print(f"Volume: {decomp.gas_volume_L_g:.4f} L/g")
```

### Gas Composition

```python
decomp = compound.max_decomposition()

print("Gas composition (mole fractions):")
for gas, fraction in decomp.gas_composition.items():
    print(f"  {gas}: {fraction*100:.1f}%")
```

### Temperature Effects

Gas volume depends on temperature:

```python
# Default (298.15 K)
decomp_298 = compound.max_decomposition()

# Elevated temperature
decomp_500 = compound.max_decomposition(gas_temperature_K=500)

print(f"Gas at 298 K: {decomp_298.gas_volume_L_g:.4f} L/g")
print(f"Gas at 500 K: {decomp_500.gas_volume_L_g:.4f} L/g")
```

The ideal gas law is used:

$$
V = \frac{n \cdot R \cdot T}{P \cdot MW}
$$

## Product Distribution

### Accessing Products

```python
decomp = compound.max_decomposition()

# All products
print("All products:")
for product, moles in decomp.products.items():
    print(f"  {product}: {moles:.4f} mol")

# Filter significant products
print("\nSignificant products (>0.01 mol):")
for product, moles in decomp.products.items():
    if moles > 0.01:
        print(f"  {product}: {moles:.3f} mol")
```

### Product Enthalpies

Reference formation enthalpies used:

| Product | ΔHf° (kJ/mol) | Phase |
|---------|---------------|-------|
| HF | -273.30 | gas |
| N₂ | 0.0 | gas |
| P₄O₁₀ | -2984.0 | solid |
| H₂O | -241.83 | gas |
| CO₂ | -393.52 | gas |
| SO₂ | -296.81 | gas |
| CO | -110.53 | gas |
| HCl | -92.31 | gas |
| HBr | -36.29 | gas |
| C (graphite) | 0.0 | solid |
| H₂ | 0.0 | gas |

## Practical Examples

### Comparing Different Compounds

```python
from phoenix import Compound

smiles_data = [
    ("CCO", "Ethanol"),
    ("c1ccccc1[N+](=O)[O-]", "Nitrobenzene"),
    ("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]", "TNT"),
]

for smiles, name in smiles_data:
    compound = Compound.from_smiles(smiles)
    decomp = compound.max_decomposition()

    print(f"\n{name}:")
    print(f"  ΔHd = {decomp.delta_hd_cal_g:.0f} cal/g")
    print(f"  Gas = {decomp.gas_volume_L_g:.4f} L/g")
```

### Energy Balance Analysis

```python
compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
decomp = compound.max_decomposition()

# Energy balance
print("Energy Balance:")
print(f"  Reactant ΔHf°: {decomp.reactant_hf_kJ_mol:.1f} kJ/mol")
print(f"  Products ΔHf°: {decomp.products_hf_kJ_mol:.1f} kJ/mol")
print(f"  ΔHd: {decomp.delta_hd_kJ_mol:.1f} kJ/mol")

# Verify calculation
calculated = decomp.products_hf_kJ_mol - decomp.reactant_hf_kJ_mol
print(f"\nVerification:")
print(f"  Calculated: {calculated:.1f} kJ/mol")
print(f"  Reported: {decomp.delta_hd_kJ_mol:.1f} kJ/mol")
```

### Oxygen-Limited Decomposition

```python
# Oxygen-deficient compound
compound = Compound.from_smiles("CCCCCC")  # Hexane

decomp = compound.max_decomposition()

print(f"Hexane decomposition:")
print(f"  OB% = {compound.oxygen_balance:.1f}%")
print(f"  ΔHd = {decomp.delta_hd_cal_g:.1f} cal/g")

print("\nProducts (no oxygen available):")
for product, moles in decomp.products.items():
    if moles > 0.01:
        print(f"  {product}: {moles:.2f} mol")
```

## Error Handling

### Decomposition Errors

```python
from phoenix import Compound, DecompositionError

try:
    compound = Compound.from_smiles("CCO")
    decomp = compound.max_decomposition()
except DecompositionError as e:
    print(f"Decomposition failed: {e.reason}")
    if e.formula:
        print(f"  Formula: {e.formula}")
```

### LP Solver Failures

If LP fails, it falls back to hierarchy:

```python
import warnings

# LP failure triggers warning and fallback
decomp = compound.max_decomposition(method="lp")
# Warning: LP solver failed: ... Falling back to hierarchy method.
```

---

## Next Steps

- [Hazard Evaluation](hazard-evaluation.md): Complete hazard assessment
- [Reactions](reactions.md): Reaction thermodynamics
- [API Reference](../api/hazard.md): DecompositionResult API
