# Thermo API

Thermodynamic property models and data access.

---

## ThermoProperty

::: phoenix.thermo.models.ThermoProperty
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `value` | float | Property value |
| `unit` | str | Unit of measurement |
| `uncertainty` | float | None | Estimated uncertainty |
| `source` | str | None | Data source description |
| `phase` | str | Phase ('g', 'l', 's') |
| `breakdown` | tuple[GroupContribution, ...] | Group contributions |
| `references` | tuple[Reference, ...] | Literature references |
| `reference_value` | ThermoValue | None | Reference for comparison |
| `estimation_method` | str | Estimation method used |
| `temperature_K` | float | Temperature in Kelvin |

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `deviation` | float | None | Estimate - reference |
| `deviation_percent` | float | None | % deviation from reference |

### Methods

```python
def has_breakdown(self) -> bool
```
Check if group contribution breakdown is available.

```python
def has_reference(self) -> bool
```
Check if reference value is available.

```python
def format_breakdown(self, property_name: str = "ENTHALPY OF FORMATION") -> str
```
Format CHETAH-style breakdown table.

### Float Conversion

```python
hf = compound.enthalpy_of_formation
value = float(hf)  # Get raw float
```

---

## ThermoState

::: phoenix.thermo.models.ThermoState
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `temperature` | float | Temperature in Kelvin |

### Properties

| Property | Alias | Type | Unit | Description |
|----------|-------|------|------|-------------|
| `H` | `enthalpy` | ThermoProperty | kJ/mol | Enthalpy of formation |
| `S` | `entropy` | ThermoProperty | J/(mol·K) | Standard entropy |
| `Cp` | `heat_capacity` | ThermoProperty | J/(mol·K) | Heat capacity |
| `G` | `gibbs_energy` | ThermoProperty | kJ/mol | Gibbs free energy |

### Example

```python
from phoenix import Compound

compound = Compound.from_smiles("CCO")
state = compound.thermo_at(T=500)

print(f"T = {state.temperature} K")
print(f"H = {state.H.value:.2f} kJ/mol")
print(f"S = {state.S.value:.2f} J/(mol·K)")
print(f"Cp = {state.Cp.value:.2f} J/(mol·K)")
print(f"G = {state.G.value:.2f} kJ/mol")
```

---

## ThermoPropertyAccessor

::: phoenix.thermo.models.ThermoPropertyAccessor
    options:
      show_root_heading: true
      show_source: false

Enables dual property/method access for temperature-dependent properties.

### Usage

```python
# As property (298.15 K)
hf = compound.enthalpy_of_formation
print(hf.value)  # Value at 298.15 K

# As method with temperature
hf_500 = compound.enthalpy_of_formation(T=500)
print(hf_500.value)  # Value at 500 K

# Vectorized with NumPy
import numpy as np
temps = np.linspace(300, 1000, 100)
values = compound.enthalpy_of_formation(T=temps)  # Returns ndarray
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `value` | float | Value at default temperature |
| `unit` | str | Unit of measurement |
| `uncertainty` | float | None | Uncertainty at default T |
| `source` | str | None | Data source |
| `temperature_K` | float | Default temperature (298.15) |

---

## GroupContribution

::: phoenix.thermo.models.GroupContribution
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `group_name` | str | Group notation |
| `count` | int | Number of occurrences |
| `contribution` | float | Contribution per group |
| `property_type` | str | "Hf", "S", or "Cp" |
| `source` | str | Data source |

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `total` | float | count × contribution |

---

## ThermoValue

::: phoenix.thermo.models.ThermoValue
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `value` | float | Property value |
| `unit` | str | Unit |
| `uncertainty` | float | None | Uncertainty |
| `method` | str | Determination method |
| `references` | tuple[Reference, ...] | Literature references |

---

## Estimation Functions

### estimate_delta_hf

```python
from phoenix.thermo.benson import estimate_delta_hf

def estimate_delta_hf(compound: Compound, T: float = 298.15) -> ThermoProperty
```

Estimate enthalpy of formation using Benson GA.

### estimate_entropy

```python
from phoenix.thermo.benson import estimate_entropy

def estimate_entropy(compound: Compound, T: float = 298.15) -> ThermoProperty
```

Estimate standard entropy using Benson GA.

### estimate_heat_capacity

```python
from phoenix.thermo.benson import estimate_heat_capacity

def estimate_heat_capacity(compound: Compound, temperature_K: float = 298.15) -> ThermoProperty
```

Estimate heat capacity using Benson GA.

---

## Data Access

### get_formation_enthalpy

```python
from phoenix.thermo.data import get_formation_enthalpy

def get_formation_enthalpy(formula: str) -> FormationEnthalpyData
```

Get reference formation enthalpy data.

---

## Temperature Constants

```python
from phoenix.thermo.models import (
    TEMP_MIN_WARN,   # 200.0 K - warn below
    TEMP_MAX_WARN,   # 6000.0 K - warn above
    TEMP_DEFAULT,    # 298.15 K - standard state
)
```

---

## Example: Full Breakdown

```python
from phoenix import Compound

compound = Compound.from_smiles("CCC")  # Propane
hf = compound.enthalpy_of_formation

# Print CHETAH-style breakdown
if hf.has_breakdown():
    print(hf.format_breakdown())

# Access individual contributions
for group in hf.breakdown:
    print(f"{group.group_name}: {group.count} × {group.contribution:.2f} = {group.total:.2f}")
```
