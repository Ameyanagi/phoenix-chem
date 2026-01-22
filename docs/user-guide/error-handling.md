# Error Handling

PHOENIX provides a hierarchy of exceptions for precise error handling.

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

## PhoenixError

Base exception for all PHOENIX errors:

```python
from phoenix import Compound, PhoenixError

try:
    compound = Compound.from_smiles("invalid")
except PhoenixError as e:
    print(f"PHOENIX error: {e}")
```

Use this for catching all PHOENIX-specific errors while letting other exceptions propagate.

## InvalidSmilesError

Raised when a SMILES string cannot be parsed:

```python
from phoenix import Compound, InvalidSmilesError

try:
    compound = Compound.from_smiles("not-a-valid-smiles")
except InvalidSmilesError as e:
    print(f"Invalid SMILES: {e.smiles}")
    print(f"Message: {e}")
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | str | The invalid SMILES string |

### Common Causes

- Typos in SMILES string
- Unbalanced parentheses or brackets
- Invalid atom symbols
- Malformed ring closures

```python
# Examples of invalid SMILES
invalid_examples = [
    "CCO(",      # Unbalanced parenthesis
    "C1CC",      # Unclosed ring
    "CXC",       # Invalid atom 'X'
    "random",    # Not SMILES at all
]
```

## UnsupportedElementError

Raised when a compound contains elements not supported by PHOENIX:

```python
from phoenix import Compound, UnsupportedElementError

try:
    compound = Compound.from_smiles("[Fe](C)(C)(C)(C)C")  # Ferrocene-like
except UnsupportedElementError as e:
    print(f"Unsupported elements: {e.elements}")
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `elements` | list[str] | List of unsupported element symbols |

### Supported Elements

PHOENIX supports: **C, H, N, O, S, P, F, Cl, Br**

Unsupported include:
- Metals (Fe, Cu, Zn, etc.)
- Iodine (I)
- Silicon (Si)
- Boron (B)
- Selenium (Se)

## UnsupportedStructureError

Raised for valid SMILES with unsupported molecular structures:

```python
from phoenix import Compound, UnsupportedStructureError

# Charged species
try:
    compound = Compound.from_smiles("[NH4+]")
except UnsupportedStructureError as e:
    print(f"Unsupported: {e.reason}")
    print(f"SMILES: {e.smiles}")

# Radical species
try:
    compound = Compound.from_smiles("[CH3]")  # Methyl radical
except UnsupportedStructureError as e:
    print(f"Unsupported: {e.reason}")
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `reason` | str | Why the structure is unsupported |
| `smiles` | str | None | The SMILES if available |

### Unsupported Structures

- Charged species (ions)
- Radical species
- Coordination complexes

## MissingGroupError

Raised when Benson Group Additivity lacks data for a molecule:

```python
from phoenix import Compound, MissingGroupError

try:
    compound = Compound.from_smiles("complex-molecule")
    hf = compound.enthalpy_of_formation
except MissingGroupError as e:
    print(f"Missing groups: {e.groups}")
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `groups` | list[str] | List of missing group names |

### Common Causes

- Unusual bonding patterns
- Rare functional groups
- Strained ring systems not in database

## DecompositionError

Raised when decomposition calculation fails:

```python
from phoenix import Compound, DecompositionError

try:
    compound = Compound.from_smiles("CCO")
    decomp = compound.max_decomposition()
except DecompositionError as e:
    print(f"Decomposition failed: {e.reason}")
    if e.formula:
        print(f"Formula: {e.formula}")
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `reason` | str | Why decomposition failed |
| `formula` | str | None | Molecular formula if available |

## BalanceError

Base exception for reaction balancing errors:

```python
from phoenix import Reaction, BalanceError

try:
    rxn = Reaction.from_smiles(["CH4"], ["CO2"])
    rxn.balance()
except BalanceError as e:
    print(f"Balance error: {e}")
```

## OverconstrainedError

Raised when reaction constraints are inconsistent:

```python
from phoenix import Reaction, OverconstrainedError

try:
    # Impossible: 1 CH4 + 1 O2 cannot give these products
    rxn = Reaction.from_smiles(
        reactants=[("CH4", 1), ("O=O", 1)],
        products=[("O=C=O", 1), ("O", 2)]
    )
    rxn.balance()
except OverconstrainedError as e:
    print(f"Over-constrained: {e}")
    print(f"Imbalances: {e.imbalances}")
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `imbalances` | dict[str, float] | Element imbalances (positive = excess product) |

### Understanding Imbalances

```python
except OverconstrainedError as e:
    for element, imbalance in e.imbalances.items():
        if imbalance > 0:
            print(f"{element}: {imbalance:.2f} excess on product side")
        else:
            print(f"{element}: {-imbalance:.2f} excess on reactant side")
```

## UnderconstrainedError

Raised when multiple valid solutions exist:

```python
from phoenix import Reaction, UnderconstrainedError

try:
    rxn = Reaction.from_smiles(
        reactants=["C", "O=O"],
        products=["O=C=O", "C=O"]  # Both CO2 and CO possible
    )
    rxn.balance()
except UnderconstrainedError as e:
    print(f"Under-constrained: {e}")
    print(f"Degrees of freedom: {e.degrees_of_freedom}")
    print(f"Suggestion: {e.suggestion}")
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `degrees_of_freedom` | int | Additional coefficients needed |
| `suggestion` | str | None | How to resolve the issue |

### Resolving Under-Constrained Reactions

```python
from phoenix import Reaction, Auto

# Fix one more coefficient
rxn = Reaction.from_smiles(
    reactants=[("C", 1), ("O=O", Auto)],
    products=[("O=C=O", 1), ("C=O", 0)]  # Force no CO
)
rxn.balance()
```

## Best Practices

### Catch Specific Exceptions

```python
from phoenix import (
    Compound,
    InvalidSmilesError,
    UnsupportedElementError,
    UnsupportedStructureError,
    PhoenixError,
)

def safe_create_compound(smiles):
    try:
        return Compound.from_smiles(smiles)
    except InvalidSmilesError:
        print(f"Invalid SMILES: {smiles}")
    except UnsupportedElementError as e:
        print(f"Unsupported elements: {e.elements}")
    except UnsupportedStructureError as e:
        print(f"Unsupported structure: {e.reason}")
    except PhoenixError as e:
        print(f"Other PHOENIX error: {e}")
    return None
```

### Batch Processing with Error Recovery

```python
from phoenix import Compound, PhoenixError

smiles_list = ["CCO", "invalid", "[Fe]", "c1ccccc1"]

results = []
errors = []

for smiles in smiles_list:
    try:
        compound = Compound.from_smiles(smiles)
        result = compound.evaluate_hazard()
        results.append((smiles, result))
    except PhoenixError as e:
        errors.append((smiles, type(e).__name__, str(e)))

print(f"Processed: {len(results)}")
print(f"Errors: {len(errors)}")
```

### Logging Errors

```python
import logging
from phoenix import Compound, PhoenixError

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

def process_smiles(smiles):
    try:
        return Compound.from_smiles(smiles)
    except PhoenixError as e:
        logger.warning(f"Failed to process {smiles}: {type(e).__name__}: {e}")
        return None
```

### Using the screen() Function

The `screen()` function handles errors gracefully:

```python
from phoenix import screen

results = screen(["CCO", "invalid", "[Fe]", "c1ccccc1"])

# Check for errors in DataFrame
df = results.dataframe
failed = df[df['error'].notna()]

for _, row in failed.iterrows():
    print(f"SMILES: {row['smiles']}")
    print(f"Error: {row['error']}: {row['error_message']}")
```

## Error Messages

PHOENIX provides informative error messages:

```python
# InvalidSmilesError
"Invalid SMILES string: 'not-valid'"

# UnsupportedElementError
"Unsupported elements: Fe, Cu"
"Iodine (I) is not supported in MVP. Unsupported elements: I"

# UnsupportedStructureError
"Unsupported structure: Charged species not supported (formal charge: 1)"
"Unsupported structure: Radical species not supported"

# OverconstrainedError
"Reaction is over-constrained: atom conservation cannot be satisfied. Imbalances: O: +1.0000"

# UnderconstrainedError
"Reaction is under-constrained: multiple solutions exist (2 degree(s) of freedom). Specify 2 more coefficient(s) to get a unique solution."
```

---

## Next Steps

- [Compounds](compounds.md): Creating compounds
- [Batch Processing](batch-processing.md): Handling errors in batch
- [API Reference](../api/exceptions.md): Complete exceptions reference
