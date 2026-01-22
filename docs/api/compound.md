# Compound API

::: phoenix.core.compound.Compound
    options:
      show_root_heading: true
      show_source: false
      members_order: source
      heading_level: 2

---

## Class Overview

The `Compound` class is the primary entry point for PHOENIX.

```python
from phoenix import Compound

compound = Compound.from_smiles("CCO")
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
| `phase` | str | `"g"` | Phase: 'g', 'l', or 's' |

**Returns:** `Compound`

**Raises:**

- `InvalidSmilesError` - If SMILES cannot be parsed
- `UnsupportedElementError` - If contains unsupported elements
- `UnsupportedStructureError` - If charged or radical

**Example:**

```python
ethanol = Compound.from_smiles("CCO")
water_liquid = Compound.from_smiles("O", phase="l")
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
| `original_smiles` | str | None | Input SMILES |
| `inchikey` | str | None | InChIKey identifier |
| `phase` | str | Phase ('g', 'l', 's') |
| `warnings` | list[str] | Validation warnings |
| `rdmol` | Chem.Mol | RDKit molecule object |

### Thermodynamic Properties

| Property | Type | Description |
|----------|------|-------------|
| `enthalpy_of_formation` | ThermoPropertyAccessor | ΔHf° accessor |
| `delta_hf_kJ_mol` | float | ΔHf° at 298.15 K |
| `entropy` | ThermoPropertyAccessor | S° accessor |
| `entropy_J_mol_K` | float | S° at 298.15 K |
| `oxygen_balance` | float | OB% |
| `oxygen_balance_percent` | float | Alias for oxygen_balance |

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

**Returns:** `ThermoProperty`

#### thermo_at

```python
def thermo_at(self, *, T: float) -> ThermoState
```

Get all thermodynamic properties at specified temperature.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `T` | float | Temperature in Kelvin (keyword-only) |

**Returns:** `ThermoState`

**Example:**

```python
state = compound.thermo_at(T=500)
print(f"H = {state.H.value} kJ/mol")
print(f"S = {state.S.value} J/(mol·K)")
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

Calculate maximum heat of decomposition.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `method` | str | `"hierarchy"` | 'hierarchy', 'lp', or 'both' |
| `gas_temperature_K` | float | 298.15 | Temperature for gas volume |

**Returns:** `DecompositionResult` or `DecompositionComparison` (when method='both')

#### evaluate_hazard

```python
def evaluate_hazard(self) -> HazardResult
```

Perform full hazard evaluation.

**Returns:** `HazardResult`

## Magic Methods

| Method | Description |
|--------|-------------|
| `__repr__` | `Compound('CCO')` |
| `__str__` | `C2H6O (MW=46.07)` |
| `__eq__` | Equality by canonical SMILES |
| `__hash__` | Hash by canonical SMILES |

## Supported Elements

```python
SUPPORTED_ELEMENTS = frozenset({"C", "H", "N", "O", "S", "P", "F", "Cl", "Br"})
```

## Large Molecule Warning

```python
LARGE_MOLECULE_THRESHOLD = 100  # atoms
```

Molecules with more than 100 atoms generate a `UserWarning`.
