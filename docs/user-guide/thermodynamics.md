# Thermodynamic Properties

PHOENIX provides thermodynamic property estimation using Benson Group Additivity (GA) with reference data comparison.

## Overview

### Available Properties

| Property | Symbol | Unit | Description |
|----------|--------|------|-------------|
| Enthalpy of formation | ΔHf° | kJ/mol | Heat of formation from elements |
| Entropy | S° | J/(mol·K) | Standard molar entropy |
| Heat capacity | Cp | J/(mol·K) | Constant pressure heat capacity |
| Gibbs free energy | G° | kJ/mol | Calculated as H - TS |

### Data Sources

1. **Benson Group Additivity** - Primary estimation method
   - Uses pgradd library
   - Covers most organic functional groups
   - Typical accuracy: ±5-10 kJ/mol for ΔHf°

2. **Reference Data** - For validation and comparison
   - NIST-JANAF thermochemical tables
   - chemicals library (CalebBell/ChEDL)

## ThermoProperty Class

The `ThermoProperty` dataclass holds a thermodynamic value with metadata:

```python
from phoenix import Compound

compound = Compound.from_smiles("CCO")
hf = compound.enthalpy_of_formation

# Access value and metadata
print(f"Value: {hf.value} {hf.unit}")
print(f"Source: {hf.source}")
print(f"Uncertainty: {hf.uncertainty}")
print(f"Temperature: {hf.temperature_K} K")
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `value` | float | The property value |
| `unit` | str | Unit of measurement |
| `uncertainty` | float | None | Estimated uncertainty |
| `source` | str | None | Data source description |
| `phase` | str | Phase: 'g', 'l', or 's' |
| `temperature_K` | float | Temperature in Kelvin |
| `breakdown` | tuple | Group contribution details |
| `reference_value` | ThermoValue | None | Reference for comparison |

### Using as Float

ThermoProperty can be used directly as a float:

```python
hf = compound.enthalpy_of_formation

# Arithmetic operations
delta_h = float(hf) - 100.0

# Comparisons
if float(hf) < -200:
    print("Highly exothermic formation")
```

## Temperature-Dependent Calculations

### Single Temperature

```python
compound = Compound.from_smiles("CCO")

# Default (298.15 K)
hf_default = compound.enthalpy_of_formation
print(f"ΔHf°(298K) = {hf_default.value:.1f} kJ/mol")

# Custom temperature
hf_500 = compound.enthalpy_of_formation(T=500)
print(f"ΔHf°(500K) = {hf_500.value:.1f} kJ/mol")
```

### Multiple Temperatures

Use NumPy arrays for vectorized calculations:

```python
import numpy as np
from phoenix import Compound

compound = Compound.from_smiles("CCO")

# Temperature range
temps = np.linspace(300, 1000, 50)

# Get values at all temperatures
hf_values = compound.enthalpy_of_formation(T=temps)

print(f"Type: {type(hf_values)}")  # numpy.ndarray
print(f"Min: {hf_values.min():.1f}, Max: {hf_values.max():.1f} kJ/mol")
```

### Temperature Warnings

PHOENIX warns for extreme temperatures:

```python
import warnings

compound = Compound.from_smiles("CCO")

# Low temperature warning (<200 K)
hf_low = compound.enthalpy_of_formation(T=100)
# Warning: Temperature 100 K is below 200 K. Benson GA extrapolation may be unreliable.

# High temperature warning (>6000 K)
hf_high = compound.enthalpy_of_formation(T=7000)
# Warning: Temperature 7000 K is above 6000 K. NASA polynomial extrapolation may be unreliable.
```

## ThermoState

Get all properties at a single temperature using `thermo_at()`:

```python
compound = Compound.from_smiles("CCO")

# All properties at 500 K
state = compound.thermo_at(T=500)

print(f"Temperature: {state.temperature} K")
print(f"H  = {state.H.value:>8.1f} kJ/mol")
print(f"S  = {state.S.value:>8.1f} J/(mol·K)")
print(f"Cp = {state.Cp.value:>8.1f} J/(mol·K)")
print(f"G  = {state.G.value:>8.1f} kJ/mol")
```

### State Properties

ThermoState provides both short and long names:

```python
state = compound.thermo_at(T=500)

# Short names
print(state.H.value)   # Enthalpy
print(state.S.value)   # Entropy
print(state.Cp.value)  # Heat capacity
print(state.G.value)   # Gibbs energy

# Long names
print(state.enthalpy.value)
print(state.entropy.value)
print(state.heat_capacity.value)
print(state.gibbs_energy.value)
```

### Gibbs Energy Calculation

Gibbs energy is calculated as:

$$
G = H - T \cdot S
$$

```python
state = compound.thermo_at(T=500)

# Verify G calculation
H_kJ = state.H.value
S_kJ_K = state.S.value / 1000  # Convert J to kJ
T = state.temperature

G_calculated = H_kJ - T * S_kJ_K
print(f"G (calculated): {G_calculated:.1f} kJ/mol")
print(f"G (from state): {state.G.value:.1f} kJ/mol")
```

## Group Contribution Breakdown

For Benson GA estimates, view the group breakdown:

```python
compound = Compound.from_smiles("CCC")  # Propane
hf = compound.enthalpy_of_formation

# Check if breakdown is available
if hf.has_breakdown():
    print(hf.format_breakdown())
```

Output example:
```
ENTHALPY OF FORMATION (GAS) CALCULATION
============================================================
Group                         Count    Contribution      Total
------------------------------------------------------------
C-(H)3(C)                        2         -10.00       -20.00
C-(H)2(C)2                       1          -5.00        -5.00
------------------------------------------------------------
ESTIMATED VALUE (Benson GA):                           -25.00 kJ/mol
```

## Reference Value Comparison

When reference data is available:

```python
compound = Compound.from_smiles("CCO")
hf = compound.enthalpy_of_formation

# Check for reference value
if hf.has_reference():
    print(f"Estimated: {hf.value:.1f} kJ/mol")
    print(f"Reference: {hf.reference_value.value:.1f} kJ/mol")
    print(f"Deviation: {hf.deviation:.1f} kJ/mol")
    print(f"Deviation: {hf.deviation_percent:.1f}%")
```

## Backward Compatibility

Legacy float properties are still available:

```python
compound = Compound.from_smiles("CCO")

# Legacy properties (at 298.15 K)
hf_kJ = compound.delta_hf_kJ_mol  # float
s_J = compound.entropy_J_mol_K    # float

print(f"ΔHf° = {hf_kJ:.1f} kJ/mol")
print(f"S° = {s_J:.1f} J/(mol·K)")
```

## Benson GA vs Reference Data

### When to Trust GA

Benson GA is most accurate for:

- Simple hydrocarbons
- Common alcohols, ethers, ketones
- Well-characterized functional groups

### When to Be Cautious

GA may be less accurate for:

- Strained ring systems
- Multiple adjacent functional groups
- Unusual bonding patterns
- Very large molecules (>100 atoms)

### Checking Accuracy

```python
compound = Compound.from_smiles("CCO")
hf = compound.enthalpy_of_formation

# Typical Benson GA uncertainty
typical_uncertainty = 10.0  # kJ/mol

# Check if reference exists for validation
if hf.has_reference():
    actual_deviation = abs(hf.deviation)
    print(f"Deviation from reference: {actual_deviation:.1f} kJ/mol")

    if actual_deviation > 2 * typical_uncertainty:
        print("Warning: Large deviation - validate experimentally")
```

## Temperature Ranges

### Valid Ranges

| Method | Temperature Range |
|--------|------------------|
| Benson GA | 298-1500 K (best) |
| NASA polynomials | 200-6000 K |

### Extrapolation Warnings

```python
# These will generate warnings
state_cold = compound.thermo_at(T=100)   # Below 200 K
state_hot = compound.thermo_at(T=7000)   # Above 6000 K
```

---

## Next Steps

- [Decomposition](decomposition.md): Maximum decomposition calculation
- [Reactions](reactions.md): Reaction thermodynamics
- [API Reference](../api/thermo.md): Complete ThermoProperty API
