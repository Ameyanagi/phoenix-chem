# Exceptions API

All PHOENIX exception classes.

---

## Exception Hierarchy

```
PhoenixError (base)
├── InvalidSmilesError
├── UnsupportedElementError
├── UnsupportedStructureError
├── MissingGroupError
├── DecompositionError
└── BalanceError
    ├── OverconstrainedError
    └── UnderconstrainedError
```

---

## PhoenixError

::: phoenix.exceptions.PhoenixError
    options:
      show_root_heading: true
      show_source: false

Base exception for all PHOENIX errors.

```python
from phoenix import PhoenixError

try:
    # PHOENIX operations
    pass
except PhoenixError as e:
    print(f"PHOENIX error: {e}")
```

---

## InvalidSmilesError

::: phoenix.exceptions.InvalidSmilesError
    options:
      show_root_heading: true
      show_source: false

Raised when a SMILES string cannot be parsed.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | str | The invalid SMILES string |

### Example

```python
from phoenix import Compound, InvalidSmilesError

try:
    compound = Compound.from_smiles("not-valid")
except InvalidSmilesError as e:
    print(f"Invalid SMILES: {e.smiles}")
```

---

## UnsupportedElementError

::: phoenix.exceptions.UnsupportedElementError
    options:
      show_root_heading: true
      show_source: false

Raised when a compound contains unsupported elements.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `elements` | list[str] | List of unsupported element symbols |

### Example

```python
from phoenix import Compound, UnsupportedElementError

try:
    compound = Compound.from_smiles("[Fe]")
except UnsupportedElementError as e:
    print(f"Unsupported elements: {e.elements}")
```

### Supported Elements

`C, H, N, O, S, P, F, Cl, Br`

---

## UnsupportedStructureError

::: phoenix.exceptions.UnsupportedStructureError
    options:
      show_root_heading: true
      show_source: false

Raised for valid SMILES with unsupported molecular structures.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `reason` | str | Why the structure is unsupported |
| `smiles` | str | None | The SMILES if available |

### Example

```python
from phoenix import Compound, UnsupportedStructureError

try:
    compound = Compound.from_smiles("[NH4+]")
except UnsupportedStructureError as e:
    print(f"Reason: {e.reason}")
    print(f"SMILES: {e.smiles}")
```

### Unsupported Structures

- Charged species (ions)
- Radical species
- Coordination complexes

---

## MissingGroupError

::: phoenix.exceptions.MissingGroupError
    options:
      show_root_heading: true
      show_source: false

Raised when Benson GA lacks group contribution data.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `groups` | list[str] | List of missing group names |

### Example

```python
from phoenix import MissingGroupError

try:
    hf = compound.enthalpy_of_formation
except MissingGroupError as e:
    print(f"Missing groups: {e.groups}")
```

---

## DecompositionError

::: phoenix.exceptions.DecompositionError
    options:
      show_root_heading: true
      show_source: false

Raised when decomposition calculation fails.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `reason` | str | Why decomposition failed |
| `formula` | str | None | Molecular formula if available |

### Example

```python
from phoenix import DecompositionError

try:
    decomp = compound.max_decomposition()
except DecompositionError as e:
    print(f"Reason: {e.reason}")
    print(f"Formula: {e.formula}")
```

---

## BalanceError

::: phoenix.exceptions.BalanceError
    options:
      show_root_heading: true
      show_source: false

Base exception for reaction balancing errors.

```python
from phoenix import BalanceError

try:
    rxn.balance()
except BalanceError as e:
    print(f"Balance error: {e}")
```

---

## OverconstrainedError

::: phoenix.exceptions.OverconstrainedError
    options:
      show_root_heading: true
      show_source: false

Raised when reaction constraints are inconsistent.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `imbalances` | dict[str, float] | Element imbalances |

### Imbalance Signs

- Positive: excess on product side
- Negative: excess on reactant side

### Example

```python
from phoenix import Reaction, OverconstrainedError

try:
    rxn = Reaction.from_smiles(
        reactants=[("CH4", 1), ("O=O", 1)],
        products=[("O=C=O", 1), ("O", 2)]
    )
    rxn.balance()
except OverconstrainedError as e:
    print(f"Imbalances: {e.imbalances}")
    for elem, imb in e.imbalances.items():
        side = "products" if imb > 0 else "reactants"
        print(f"  {elem}: {abs(imb):.2f} excess on {side}")
```

---

## UnderconstrainedError

::: phoenix.exceptions.UnderconstrainedError
    options:
      show_root_heading: true
      show_source: false

Raised when multiple valid solutions exist.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `degrees_of_freedom` | int | Additional constraints needed |
| `suggestion` | str | None | How to resolve |

### Example

```python
from phoenix import Reaction, UnderconstrainedError, Auto

try:
    rxn = Reaction.from_smiles(
        reactants=["C", "O=O"],
        products=["O=C=O", "C=O"]  # Both CO2 and CO
    )
    rxn.balance()
except UnderconstrainedError as e:
    print(f"DOF: {e.degrees_of_freedom}")
    print(f"Suggestion: {e.suggestion}")
```

### Resolution

Specify additional coefficients:

```python
rxn = Reaction.from_smiles(
    reactants=[("C", 1), ("O=O", Auto)],
    products=[("O=C=O", 1), ("C=O", 0)]  # Force no CO
)
rxn.balance()
```

---

## Import All Exceptions

```python
from phoenix import (
    PhoenixError,
    InvalidSmilesError,
    UnsupportedElementError,
    UnsupportedStructureError,
    MissingGroupError,
    DecompositionError,
    BalanceError,
    OverconstrainedError,
    UnderconstrainedError,
)
```

---

## Best Practices

### Catch Specific First

```python
try:
    compound = Compound.from_smiles(smiles)
except InvalidSmilesError:
    # Handle invalid SMILES
    pass
except UnsupportedElementError:
    # Handle unsupported elements
    pass
except PhoenixError:
    # Catch-all for other PHOENIX errors
    pass
```

### Batch Processing

```python
from phoenix import Compound, PhoenixError

results = []
errors = []

for smiles in smiles_list:
    try:
        compound = Compound.from_smiles(smiles)
        results.append((smiles, compound))
    except PhoenixError as e:
        errors.append((smiles, type(e).__name__, str(e)))
```
