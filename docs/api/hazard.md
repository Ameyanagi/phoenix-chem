# Hazard API

Hazard evaluation and decomposition analysis classes.

---

## HazardResult

::: phoenix.hazard.classification.HazardResult
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | str | Canonical SMILES |
| `formula` | str | Molecular formula |
| `molecular_weight` | float | MW in g/mol |
| `delta_hf_kJ_mol` | float | ΔHf° in kJ/mol |
| `max_decomposition_kJ_mol` | float | ΔHd in kJ/mol |
| `max_decomposition_cal_g` | float | ΔHd in cal/g |
| `oxygen_balance_percent` | float | OB% |
| `hazard_class` | HazardClass | "HIGH", "MEDIUM", or "LOW" |
| `triggered_criteria` | tuple[int, ...] | Criteria numbers (1-4) |
| `functional_group_alerts` | tuple[str, ...] | Alert names |
| `product_breakdown` | dict[str, float] | Decomposition products |
| `gas_volume_L_g` | float | Gas generation in L/g |
| `method` | str | Calculation method |

### Example

```python
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
result = compound.evaluate_hazard()

print(result.hazard_class)
print(result.triggered_criteria)
print(result.functional_group_alerts)
```

---

## DecompositionResult

::: phoenix.hazard.decomposition.DecompositionResult
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `delta_hd_kJ_mol` | float | ΔHd in kJ/mol (negative = exothermic) |
| `delta_hd_cal_g` | float | ΔHd in cal/g |
| `products` | dict[str, float] | Product moles |
| `reactant_hf_kJ_mol` | float | Reactant ΔHf° |
| `products_hf_kJ_mol` | float | Products ΔHf° sum |
| `gas_volume_L_g` | float | Gas volume per gram |
| `gas_moles` | float | Gas moles per mol compound |
| `gas_composition` | dict[str, float] | Gas mole fractions |
| `gas_temperature_K` | float | Temperature for gas calc |
| `method` | str | 'hierarchy' or 'lp' |

### Example

```python
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
decomp = compound.max_decomposition()

print(f"ΔHd = {decomp.delta_hd_cal_g:.1f} cal/g")
for product, moles in decomp.products.items():
    if moles > 0.01:
        print(f"  {product}: {moles:.2f} mol")
```

---

## DecompositionComparison

::: phoenix.hazard.decomposition.DecompositionComparison
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `hierarchy_result` | DecompositionResult | Hierarchy method result |
| `lp_result` | DecompositionResult | LP method result |
| `deviation_percent` | float | % deviation between methods |

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `hierarchy_delta_hd` | float | Hierarchy ΔHd (kJ/mol) |
| `lp_delta_hd` | float | LP ΔHd (kJ/mol) |

### Example

```python
comparison = compound.max_decomposition(method="both")
print(f"Hierarchy: {comparison.hierarchy_delta_hd:.1f} kJ/mol")
print(f"LP: {comparison.lp_delta_hd:.1f} kJ/mol")
print(f"Deviation: {comparison.deviation_percent:.2f}%")
```

---

## FunctionalGroupAlert

::: phoenix.hazard.functional_groups.FunctionalGroupAlert
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Group name |
| `smarts` | str | SMARTS pattern |
| `count` | int | Number of occurrences |
| `description` | str | Hazard description |

### Example

```python
from phoenix.hazard.functional_groups import detect_functional_groups
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
alerts = detect_functional_groups(compound)

for alert in alerts:
    print(f"{alert.name}: {alert.count}x - {alert.description}")
```

---

## Functions

### evaluate_hazard

```python
from phoenix.hazard.classification import evaluate_hazard

def evaluate_hazard(compound: Compound) -> HazardResult
```

Perform complete hazard evaluation.

### calculate_max_decomposition

```python
from phoenix.hazard.decomposition import calculate_max_decomposition

def calculate_max_decomposition(
    compound: Compound,
    *,
    method: Literal["hierarchy", "lp", "both"] = "hierarchy",
    gas_temperature_K: float = 298.15,
) -> DecompositionResult | DecompositionComparison
```

Calculate maximum heat of decomposition.

### detect_functional_groups

```python
from phoenix.hazard.functional_groups import detect_functional_groups

def detect_functional_groups(compound: Compound) -> list[FunctionalGroupAlert]
```

Detect reactive functional groups.

### calculate_oxygen_balance

```python
from phoenix.hazard.oxygen_balance import calculate_oxygen_balance

def calculate_oxygen_balance(
    composition: dict[str, int],
    molecular_weight: float
) -> float
```

Calculate oxygen balance percentage.

---

## CHETAH Criteria Constants

```python
from phoenix.hazard.classification import (
    CRITERION_1_THRESHOLD,  # -300.0 cal/g
    CRITERION_2_THRESHOLD,  # -100.0 cal/g
    CRITERION_3_OB_LOW,     # -200.0 %
    CRITERION_3_OB_HIGH,    # +100.0 %
)
```

---

## Decomposition Products

```python
from phoenix.hazard.decomposition import DECOMPOSITION_PRODUCTS

# Dict[str, Dict] with keys:
# - delta_hf_kJ_mol: Formation enthalpy
# - atoms: Element composition
# - is_gas: Boolean for gas phase
```

Available products:

| Product | ΔHf° (kJ/mol) | Gas |
|---------|---------------|-----|
| HF | -273.30 | Yes |
| N₂ | 0.0 | Yes |
| P₄O₁₀ | -2984.0 | No |
| H₂O | -241.83 | Yes |
| CO₂ | -393.52 | Yes |
| SO₂ | -296.81 | Yes |
| CO | -110.53 | Yes |
| HCl | -92.31 | Yes |
| HBr | -36.29 | Yes |
| C (graphite) | 0.0 | No |
| H₂ | 0.0 | Yes |
| S (rhombic) | 0.0 | No |
| F₂ | 0.0 | Yes |
| Cl₂ | 0.0 | Yes |
| Br₂ | 30.91 | Yes |
| O₂ | 0.0 | Yes |
