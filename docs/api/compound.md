# Compound API

::: phoenix.core.compound.Compound
    options:
      show_root_heading: true
      show_source: false
      members_order: source
      heading_level: 2

---

## Class Overview

The `Compound` class is the primary entry point for PHOENIX. It represents a single chemical compound and provides access to molecular properties, thermodynamic data, and hazard evaluation.

```python
from phoenix import Compound

compound = Compound.from_smiles("CCO")
print(compound)  # C2H6O (MW=46.07)
```

## Factory Methods

### from_smiles

```python
@classmethod
def from_smiles(cls, smiles: str, phase: str = "g") -> Compound
```

Create a Compound from a SMILES string.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `smiles` | str | required | SMILES string |
| `phase` | str | `"g"` | Phase: 'g' (gas), 'l' (liquid), or 's' (solid) |

**Returns:** `Compound`

**Raises:**

- `InvalidSmilesError` - If SMILES cannot be parsed
- `UnsupportedElementError` - If contains unsupported elements
- `UnsupportedStructureError` - If charged or radical

**Example:**

```python
from phoenix import Compound

ethanol = Compound.from_smiles("CCO")
print(f"Formula: {ethanol.formula}")     # C2H6O
print(f"MW: {ethanol.molecular_weight}") # 46.07

# Specify phase for condensed-phase calculations
water_liquid = Compound.from_smiles("O", phase="l")
print(f"Phase: {water_liquid.phase}")    # l
```

## Properties

### Molecular Information

| Property | Type | Description |
|----------|------|-------------|
| `formula` | str | Molecular formula (Hill notation) |
| `molecular_weight` | float | MW in g/mol |
| `composition` | dict[str, int] | Element counts |
| `num_atoms` | int | Total atom count |
| `canonical_smiles` | str | Canonicalized SMILES |
| `original_smiles` | str \| None | Input SMILES |
| `inchikey` | str \| None | InChIKey identifier |
| `phase` | str | Phase ('g', 'l', 's') |
| `warnings` | list[str] | Validation warnings |
| `rdmol` | Chem.Mol | RDKit molecule object |

**Example with Output:**

```python
from phoenix import Compound

tnt = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")

print(f"Formula: {tnt.formula}")            # C7H5N3O6
print(f"MW: {tnt.molecular_weight:.2f}")    # 227.13
print(f"Composition: {tnt.composition}")    # {'C': 7, 'N': 3, 'O': 6, 'H': 5}
print(f"Atom count: {tnt.num_atoms}")       # 21
print(f"InChIKey: {tnt.inchikey}")          # SPSSULHKWOKEEL-UHFFFAOYSA-N
```

### Thermodynamic Properties

| Property | Type | Description |
|----------|------|-------------|
| `enthalpy_of_formation` | ThermoPropertyAccessor | ΔHf° accessor |
| `delta_hf_kJ_mol` | float | ΔHf° at 298.15 K (shortcut) |
| `entropy` | ThermoPropertyAccessor | S° accessor |
| `entropy_J_mol_K` | float | S° at 298.15 K (shortcut) |
| `oxygen_balance` | float | Oxygen balance (OB%) |
| `oxygen_balance_percent` | float | Alias for oxygen_balance |

**Example with Output:**

```python
from phoenix import Compound

ethanol = Compound.from_smiles("CCO")

# Access as property (returns value at 298.15 K)
hf = ethanol.enthalpy_of_formation
print(f"ΔHf° = {hf.value:.1f} {hf.unit}")  # ΔHf° = -235.1 kJ/mol

# Quick access shortcuts
print(f"ΔHf° = {ethanol.delta_hf_kJ_mol:.1f} kJ/mol")  # -235.1 kJ/mol
print(f"S° = {ethanol.entropy_J_mol_K:.1f} J/(mol·K)") # 289.9 J/(mol·K)
print(f"OB% = {ethanol.oxygen_balance:.1f}%")          # -208.6%
```

### ThermoPropertyAccessor

The `enthalpy_of_formation` and `entropy` properties return a `ThermoPropertyAccessor` that supports both property and method access:

```python
from phoenix import Compound
import numpy as np

compound = Compound.from_smiles("CCO")

# As property: value at 298.15 K
hf = compound.enthalpy_of_formation
print(hf.value)  # -235.1

# As method: value at specified temperature
hf_500 = compound.enthalpy_of_formation(T=500)
print(f"ΔHf(500K) = {hf_500.value:.1f} kJ/mol")  # -227.8 kJ/mol

# Vectorized with NumPy arrays
temps = np.array([300, 400, 500, 600])
values = compound.enthalpy_of_formation(T=temps)
print(values)  # Array of ΔHf values at each temperature
```

## Methods

### Thermodynamic Methods

#### heat_capacity

```python
def heat_capacity(self, temperature_K: float = 298.15) -> ThermoProperty
```

Calculate heat capacity at specified temperature.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `temperature_K` | float | 298.15 | Temperature in Kelvin |

**Returns:** `ThermoProperty` with `value` in J/(mol·K)

**Example:**

```python
compound = Compound.from_smiles("CCO")
cp = compound.heat_capacity(temperature_K=400)
print(f"Cp(400K) = {cp.value:.1f} {cp.unit}")  # Cp(400K) = 92.3 J/(mol·K)
```

#### thermo_at

```python
def thermo_at(self, *, T: float) -> ThermoState
```

Get all thermodynamic properties at specified temperature. Returns a `ThermoState` object with cached properties.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `T` | float | Temperature in Kelvin (keyword-only) |

**Returns:** `ThermoState` with properties: `H`, `S`, `Cp`, `G`

**Example:**

```python
compound = Compound.from_smiles("CCO")
state = compound.thermo_at(T=500)

print(f"H(500K)  = {state.H.value:.1f} kJ/mol")      # -227.8 kJ/mol
print(f"S(500K)  = {state.S.value:.1f} J/(mol·K)")   # 336.2 J/(mol·K)
print(f"Cp(500K) = {state.Cp.value:.1f} J/(mol·K)")  # 109.5 J/(mol·K)
print(f"G(500K)  = {state.G.value:.1f} kJ/mol")      # -395.9 kJ/mol

# Alternative property names
print(state.enthalpy.value)      # Same as state.H
print(state.entropy.value)       # Same as state.S
print(state.heat_capacity.value) # Same as state.Cp
print(state.gibbs_energy.value)  # Same as state.G
```

### Hazard Methods

#### max_decomposition

```python
def max_decomposition(
    self,
    *,
    method: str = "hierarchy",
    gas_temperature_K: float = 298.15,
) -> DecompositionResult | DecompositionComparison
```

Calculate maximum heat of decomposition (ΔHd).

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `method` | str | `"hierarchy"` | Method: 'hierarchy', 'lp', or 'both' |
| `gas_temperature_K` | float | 298.15 | Temperature for gas volume calculations |

**Returns:**

- `DecompositionResult` when method is 'hierarchy' or 'lp'
- `DecompositionComparison` when method is 'both'

**Example:**

```python
from phoenix import Compound

nb = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")  # Nitrobenzene

# Default (hierarchy method)
decomp = nb.max_decomposition()
print(f"ΔHd = {decomp.delta_hd_kJ_mol:.1f} kJ/mol")  # -552.2 kJ/mol
print(f"ΔHd = {decomp.delta_hd_cal_g:.1f} cal/g")    # -1072.0 cal/g

# Decomposition products
for product, moles in decomp.products.items():
    if moles > 0.01:
        print(f"  {product}: {moles:.2f} mol")
# Output:
#   N2: 0.50 mol
#   H2O: 2.00 mol
#   C: 6.00 mol
#   H2: 0.50 mol

# Compare both methods
comparison = nb.max_decomposition(method="both")
print(f"Hierarchy: {comparison.hierarchy_delta_hd:.1f} kJ/mol")
print(f"LP: {comparison.lp_delta_hd:.1f} kJ/mol")
print(f"Deviation: {comparison.deviation_percent:.2f}%")
```

#### evaluate_hazard

```python
def evaluate_hazard(self) -> HazardResult
```

Perform full hazard evaluation using CHETAH criteria.

**Returns:** `HazardResult` with classification and decomposition data

**Example:**

```python
from phoenix import Compound

tnt = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")
result = tnt.evaluate_hazard()

print(f"Hazard Class: {result.hazard_class}")        # HIGH
print(f"Max ΔHd: {result.max_decomposition_kJ_mol:.1f} kJ/mol")  # -1294.8 kJ/mol
print(f"Max ΔHd: {result.max_decomposition_cal_g:.1f} cal/g")    # -1362.4 cal/g

# Triggered CHETAH criteria
print(f"Triggered criteria: {result.triggered_criteria}")
# Example: ['Criterion 1: ΔHd < -300 cal/g']

# Functional group alerts
for alert in result.functional_group_alerts:
    print(f"Alert: {alert.group_name} - {alert.severity}")
```

#### detect_functional_groups

```python
def detect_functional_groups(self) -> list[FunctionalGroupAlert]
```

Detect hazardous functional groups in the molecule.

**Returns:** List of `FunctionalGroupAlert` objects

**Example:**

```python
from phoenix import Compound

rdx = Compound.from_smiles("O=[N+]([O-])N1CN([N+](=O)[O-])CN([N+](=O)[O-])C1")
alerts = rdx.detect_functional_groups()

for alert in alerts:
    print(f"{alert.group_name}: {alert.count} occurrence(s)")
    print(f"  Severity: {alert.severity}")
# Output:
#   Nitro (N-NO2): 3 occurrence(s)
#     Severity: HIGH
```

## Magic Methods

| Method | Description | Example Output |
|--------|-------------|----------------|
| `__repr__` | Debug representation | `Compound('CCO')` |
| `__str__` | String representation | `C2H6O (MW=46.07)` |
| `__eq__` | Equality by canonical SMILES | `compound1 == compound2` |
| `__hash__` | Hash by canonical SMILES | Enables use in sets/dicts |

```python
from phoenix import Compound

c1 = Compound.from_smiles("CCO")
c2 = Compound.from_smiles("OCC")  # Same compound, different SMILES

print(repr(c1))    # Compound('CCO')
print(str(c1))     # C2H6O (MW=46.07)
print(c1 == c2)    # True (same canonical SMILES)

# Use in collections
compounds = {c1, c2}  # Set contains only 1 element
print(len(compounds)) # 1
```

## Constants

### Supported Elements

```python
SUPPORTED_ELEMENTS = frozenset({"C", "H", "N", "O", "S", "P", "F", "Cl", "Br"})
```

Compounds containing other elements raise `UnsupportedElementError`.

### Large Molecule Warning

```python
LARGE_MOLECULE_THRESHOLD = 100  # atoms
```

Molecules with more than 100 atoms generate a `UserWarning`:

```python
import warnings
from phoenix import Compound

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    large_mol = Compound.from_smiles("C" * 101)  # Very large alkane
    if w:
        print(w[0].message)  # Warning about large molecule
```
