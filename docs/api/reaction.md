# Reaction API

::: phoenix.core.reaction.Reaction
    options:
      show_root_heading: true
      show_source: false

---

## Class Overview

```python
from phoenix import Reaction, Auto

rxn = Reaction.from_smiles(
    reactants=[("CH4", 1), ("O=O", Auto)],
    products=[("O=C=O", Auto), ("O", Auto)]
)
rxn.balance()
```

## Auto Sentinel

```python
from phoenix import Auto
```

`Auto` is a singleton sentinel value for auto-calculated coefficients.

```python
# These are equivalent
rxn1 = Reaction.from_smiles(
    reactants=[("CH4", Auto), ("O2", Auto)],
    products=[("CO2", Auto), ("H2O", Auto)]
)

rxn2 = Reaction.from_smiles(
    reactants=["CH4", "O2"],
    products=["CO2", "H2O"]
)
```

## Factory Methods

### from_smiles

```python
@classmethod
def from_smiles(
    cls,
    reactants: Sequence[str] | Sequence[tuple[str, CoeffSpec]] | None = None,
    products: Sequence[str] | Sequence[tuple[str, CoeffSpec]] | None = None,
    *,
    reactant_smiles: Sequence[...] | None = None,
    product_smiles: Sequence[...] | None = None,
) -> Reaction
```

Create a Reaction from SMILES strings.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `reactants` | Sequence | Reactant SMILES or (SMILES, coeff) tuples |
| `products` | Sequence | Product SMILES or (SMILES, coeff) tuples |
| `reactant_smiles` | Sequence | Legacy keyword argument |
| `product_smiles` | Sequence | Legacy keyword argument |

**Coefficient Spec:**

- `int` or `float`: Explicit coefficient
- `Auto` or `None`: Auto-calculate

**Examples:**

```python
# All auto
rxn = Reaction.from_smiles(["CH4", "O2"], ["CO2", "H2O"])

# Mixed
rxn = Reaction.from_smiles(
    reactants=[("CH4", 1), ("O2", Auto)],
    products=[("CO2", Auto), ("H2O", Auto)]
)
```

### from_reaction_smiles

```python
@classmethod
def from_reaction_smiles(cls, reaction_smiles: str, auto_balance: bool = True) -> Reaction
```

Create from reaction SMILES format.

**Format:** `"coeff SMILES + coeff SMILES >> coeff SMILES + coeff SMILES"`

**Examples:**

```python
rxn = Reaction.from_reaction_smiles("CH4 + 2 O2 >> CO2 + 2 H2O")
rxn = Reaction.from_reaction_smiles("CH4 + O2 >> CO2 + H2O", auto_balance=True)
```

## Properties

### Species Access

| Property | Type | Description |
|----------|------|-------------|
| `reactants` | list[ReactionSpecies] | Reactant species |
| `products` | list[ReactionSpecies] | Product species |
| `all_species` | list[ReactionSpecies] | All species |
| `elements` | set[str] | Elements in reaction |
| `is_balanced` | bool | Whether balanced |
| `coefficients` | dict[str, float | None] | Formula → coefficient |
| `stoichiometry_vector` | np.ndarray | None | Signed coefficients |

### Thermodynamic Properties

| Property | Type | Description |
|----------|------|-------------|
| `enthalpy` / `delta_h` | float | ΔH_rxn in kJ/mol |
| `entropy` / `delta_s` | float | ΔS_rxn in J/(mol·K) |
| `gibbs_free_energy` / `delta_g` | float | ΔG_rxn in kJ/mol |

## Methods

### balance

```python
def balance(
    self,
    *,
    normalize: bool = True,
    prefer_integers: bool = True
) -> Reaction
```

Balance the reaction.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `normalize` | bool | True | Normalize smallest coefficient to 1 |
| `prefer_integers` | bool | True | Convert to integers if possible |

**Returns:** Self (for chaining)

**Raises:**

- `OverconstrainedError` - Inconsistent constraints
- `UnderconstrainedError` - Multiple solutions

### to_equation

```python
def to_equation(self, *, use_names: bool = False) -> str
```

Format as equation string.

**Returns:** `"CH4 + 2 O2 -> CO2 + 2 H2O"`

---

## ReactionSpecies

::: phoenix.core.reaction.ReactionSpecies
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `compound` | Compound | The compound |
| `coefficient` | float | None | Stoichiometric coefficient |
| `is_auto` | bool | Was coefficient auto-calculated |

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `formula` | str | Molecular formula |
| `composition` | dict[str, int] | Element counts |

### Methods

```python
def with_coefficient(self, coeff: float) -> ReactionSpecies
```

Create copy with new coefficient.

---

## Balancing Algorithm

### Atom Conservation

For each element j:

$$
\sum_i \nu_i \cdot a_{ij} = 0
$$

### Null-Space Method

When all coefficients unknown, find null(A) where A is the composition matrix.

### Constrained System

When some coefficients known:

$$
A_{unknown} \cdot \nu_{unknown} = -A_{known} \cdot \nu_{known}
$$

---

## Examples

### Full Auto-Balance

```python
rxn = Reaction.from_smiles(["CH4", "O2"], ["CO2", "H2O"])
rxn.balance()
print(rxn)  # CH4 + 2 O2 -> CO2 + 2 H2O
```

### Partial Constraints

```python
rxn = Reaction.from_smiles(
    reactants=[("OCC(O)CO", 1), ("[H][H]", Auto)],
    products=[("CC(O)CO", 1), ("O", Auto)]
)
rxn.balance()
print(rxn.coefficients)
```

### Thermodynamics

```python
rxn = Reaction.from_smiles(["CH4", "O2"], ["CO2", "H2O"])
rxn.balance()
print(f"ΔH = {rxn.delta_h:.1f} kJ/mol")
print(f"ΔS = {rxn.delta_s:.1f} J/(mol·K)")
print(f"ΔG = {rxn.delta_g:.1f} kJ/mol")
```
