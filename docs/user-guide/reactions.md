# Chemical Reactions

PHOENIX provides reaction balancing and thermodynamic calculations through the `Reaction` class.

## Creating Reactions

### From SMILES with Coefficients

```python
from phoenix import Reaction, Auto

# Explicit coefficients
rxn = Reaction.from_smiles(
    reactants=[("CH4", 1), ("O=O", 2)],
    products=[("O=C=O", 1), ("O", 2)]
)
print(rxn)  # CH4 + 2 O2 -> CO2 + 2 H2O
```

### Auto-Balancing

Let PHOENIX determine coefficients:

```python
# All coefficients auto-determined
rxn = Reaction.from_smiles(
    reactants=["CH4", "O=O"],
    products=["O=C=O", "O"]
)
rxn.balance()
print(rxn)  # CH4 + 2 O2 -> CO2 + 2 H2O
```

### Mixed Explicit/Auto

Fix some coefficients, auto-calculate others:

```python
# Glycerol hydrogenation
rxn = Reaction.from_smiles(
    reactants=[("OCC(O)CO", 1), ("[H][H]", Auto)],  # Fix glycerol = 1
    products=[("CC(O)CO", 1), ("O", Auto)]          # Fix propanediol = 1
)
rxn.balance()
print(rxn.coefficients)
# {'C3H8O3': 1.0, 'H2': 1.0, 'C3H8O2': 1.0, 'H2O': 1.0}
```

### From Reaction SMILES

```python
# Standard reaction SMILES format
rxn = Reaction.from_reaction_smiles("CH4 + 2 O2 >> CO2 + 2 H2O")

# Auto-balance if coefficients not specified
rxn = Reaction.from_reaction_smiles("CH4 + O2 >> CO2 + H2O", auto_balance=True)
```

## The Auto Sentinel

`Auto` marks coefficients for auto-calculation:

```python
from phoenix import Auto

# Auto is a singleton
print(Auto)        # Auto
print(type(Auto))  # <class 'phoenix.core.reaction._AutoType'>

# Auto is falsy for conditional checks
if not Auto:
    print("Auto evaluates to False")
```

Using `Auto` vs `None`:

```python
# These are equivalent
rxn1 = Reaction.from_smiles(
    reactants=[("CH4", Auto), ("O=O", Auto)],
    products=[("O=C=O", Auto), ("O", Auto)]
)

rxn2 = Reaction.from_smiles(
    reactants=["CH4", "O=O"],  # Plain strings = Auto
    products=["O=C=O", "O"]
)
```

## Balancing Algorithm

PHOENIX uses atom conservation to balance reactions.

### Atom Conservation

For each element, atoms must balance:

$$
\sum_i \nu_i \cdot a_{ij} = 0
$$

Where $\nu_i$ is the stoichiometric coefficient and $a_{ij}$ is atoms of element $j$ in species $i$.

### Null-Space Method

When all coefficients are unknown:

1. Build composition matrix $A$
2. Compute null space of $A$
3. Select simplest integer solution

```python
rxn = Reaction.from_smiles(
    reactants=["CH4", "O=O"],
    products=["O=C=O", "O"]
)
rxn.balance()

# Check if balanced
print(rxn.is_balanced)  # True
```

### Constrained System

When some coefficients are known:

1. Partition into known/unknown
2. Solve linear system for unknowns

## Reaction Properties

### Accessing Species

```python
rxn = Reaction.from_smiles(
    reactants=[("CH4", 1), ("O=O", 2)],
    products=[("O=C=O", 1), ("O", 2)]
)

# Reactants and products
for species in rxn.reactants:
    print(f"Reactant: {species.formula}, coeff: {species.coefficient}")

for species in rxn.products:
    print(f"Product: {species.formula}, coeff: {species.coefficient}")

# All species
print(f"All species: {[s.formula for s in rxn.all_species]}")
```

### Coefficients

```python
# As dictionary
print(rxn.coefficients)
# {'CH4': 1.0, 'O2': 2.0, 'CO2': 1.0, 'H2O': 2.0}

# As stoichiometry vector (negative for reactants)
nu = rxn.stoichiometry_vector
print(nu)  # [-1, -2, 1, 2]
```

### Elements

```python
print(rxn.elements)  # {'C', 'H', 'O'}
```

## Thermodynamic Properties

### Enthalpy of Reaction

```python
rxn = Reaction.from_smiles(
    reactants=["CH4", "O=O"],
    products=["O=C=O", "O"]
)
rxn.balance()

# ΔH_rxn in kJ/mol
print(f"ΔH_rxn = {rxn.enthalpy:.1f} kJ/mol")
# or equivalently
print(f"ΔH_rxn = {rxn.delta_h:.1f} kJ/mol")
```

### Entropy of Reaction

```python
# ΔS_rxn in J/(mol·K)
print(f"ΔS_rxn = {rxn.entropy:.1f} J/(mol·K)")
print(f"ΔS_rxn = {rxn.delta_s:.1f} J/(mol·K)")
```

### Gibbs Free Energy

```python
# ΔG_rxn in kJ/mol at 298.15 K
print(f"ΔG_rxn = {rxn.gibbs_free_energy:.1f} kJ/mol")
print(f"ΔG_rxn = {rxn.delta_g:.1f} kJ/mol")
```

### Complete Example

```python
from phoenix import Reaction

# Combustion of methane
rxn = Reaction.from_smiles(
    reactants=["C", "O=O"],  # Methane + O2
    products=["O=C=O", "O"]  # CO2 + H2O
)
rxn.balance()

print(f"Reaction: {rxn}")
print(f"ΔH_rxn = {rxn.delta_h:.1f} kJ/mol")
print(f"ΔS_rxn = {rxn.delta_s:.1f} J/(mol·K)")
print(f"ΔG_rxn = {rxn.delta_g:.1f} kJ/mol")

if rxn.delta_g < 0:
    print("Reaction is spontaneous at 298 K")
```

**Expected output:**
```
Reaction: CH4 + 2 O2 -> CO2 + 2 H2O
ΔH_rxn = -802.6 kJ/mol
ΔS_rxn = 169.4 J/(mol·K)
ΔG_rxn = -853.1 kJ/mol
Reaction is spontaneous at 298 K
```

## ReactionSpecies

Each species in a reaction is wrapped in `ReactionSpecies`:

```python
from phoenix import Reaction

rxn = Reaction.from_smiles(
    reactants=[("CCO", 1)],
    products=[("C=C", 1), ("O", 1)]
)
rxn.balance()

for species in rxn.all_species:
    print(f"Formula: {species.formula}")
    print(f"Coefficient: {species.coefficient}")
    print(f"Is Auto: {species.is_auto}")
    print(f"Composition: {species.composition}")
```

## Error Handling

### OverconstrainedError

Raised when constraints are inconsistent:

```python
from phoenix import Reaction, OverconstrainedError

try:
    # Impossible reaction
    rxn = Reaction.from_smiles(
        reactants=[("CH4", 1), ("O=O", 1)],  # 1 O2 can't balance
        products=[("O=C=O", 1), ("O", 2)]    # Need 2 O2
    )
    rxn.balance()
except OverconstrainedError as e:
    print(f"Cannot balance: {e}")
    print(f"Imbalances: {e.imbalances}")
```

### UnderconstrainedError

Raised when multiple solutions exist:

```python
from phoenix import Reaction, UnderconstrainedError

try:
    # Multiple solutions possible
    rxn = Reaction.from_smiles(
        reactants=["C", "O=O"],
        products=["O=C=O", "C=O"]  # Both CO2 and CO
    )
    rxn.balance()
except UnderconstrainedError as e:
    print(f"Multiple solutions: {e}")
    print(f"Degrees of freedom: {e.degrees_of_freedom}")
    print(f"Suggestion: {e.suggestion}")
```

## Balancing Options

### Normalization

Control coefficient normalization:

```python
rxn = Reaction.from_smiles(
    reactants=["CH4", "O=O"],
    products=["O=C=O", "O"]
)

# Default: normalize so smallest = 1
rxn.balance(normalize=True)

# Without normalization
rxn.balance(normalize=False)
```

### Integer Preference

Control integer coefficient conversion:

```python
# Default: prefer integer coefficients
rxn.balance(prefer_integers=True)

# Allow fractional coefficients
rxn.balance(prefer_integers=False)
```

## String Representations

```python
rxn = Reaction.from_smiles(
    reactants=[("CH4", 1), ("O=O", 2)],
    products=[("O=C=O", 1), ("O", 2)]
)

# Human-readable
print(str(rxn))   # CH4 + 2 O2 -> CO2 + 2 H2O

# Repr
print(repr(rxn))  # Reaction(CH4 + O2 >> CO2 + H2O, balanced)

# Equation method
print(rxn.to_equation())  # CH4 + 2 O2 -> CO2 + 2 H2O
```

## Practical Examples

### Hydrogenation Reaction

```python
# Glycerol hydrogenation to 1,2-propanediol
rxn = Reaction.from_smiles(
    reactants=[("OCC(O)CO", 1), ("[H][H]", Auto)],  # Glycerol + H2
    products=[("CC(O)CO", 1), ("O", Auto)]          # Propanediol + H2O
)
rxn.balance()

print(f"Reaction: {rxn}")
print(f"ΔH_rxn = {rxn.delta_h:.1f} kJ/mol")
```

### Oxidation Reaction

```python
# Ethanol oxidation
rxn = Reaction.from_smiles(
    reactants=[("CCO", 1), ("O=O", Auto)],
    products=[("O=C=O", Auto), ("O", Auto)]
)
rxn.balance()

print(f"Reaction: {rxn}")
print(f"ΔH_rxn = {rxn.delta_h:.1f} kJ/mol")
```

### Decomposition Reaction

```python
# Hydrogen peroxide decomposition
rxn = Reaction.from_smiles(
    reactants=[("OO", 1)],  # H2O2
    products=[("O", 1), ("O=O", Auto)]
)
rxn.balance()

print(f"Reaction: {rxn}")
print(f"ΔH_rxn = {rxn.delta_h:.1f} kJ/mol")
```

## Troubleshooting

### Common Balancing Errors

#### OverconstrainedError: Impossible Constraints

This occurs when fixed coefficients make balancing mathematically impossible:

```python
from phoenix import Reaction, OverconstrainedError

try:
    # Fix O2=1 but need O2=2 for balance
    rxn = Reaction.from_smiles(
        reactants=[("CH4", 1), ("O=O", 1)],  # Fixed: 1 O2
        products=[("O=C=O", 1), ("O", 2)]    # Need: 2 O2 for 4 oxygen atoms
    )
    rxn.balance()
except OverconstrainedError as e:
    print(f"Error: {e}")
    # Error: Cannot balance reaction with given constraints
```

**Solution:** Remove or adjust fixed coefficients that cause the imbalance.

#### UnderconstrainedError: Multiple Solutions

This occurs when the system has more unknowns than constraints:

```python
from phoenix import Reaction, UnderconstrainedError

try:
    # Combustion producing both CO and CO2
    rxn = Reaction.from_smiles(
        reactants=["C", "O=O"],      # 1 carbon source
        products=["O=C=O", "C=O"]    # 2 carbon products
    )
    rxn.balance()
except UnderconstrainedError as e:
    print(f"Error: {e}")
    print(f"Degrees of freedom: {e.degrees_of_freedom}")
    # Multiple valid solutions exist
```

**Solution:** Fix one or more coefficients to constrain the system:

```python
rxn = Reaction.from_smiles(
    reactants=["C", "O=O"],
    products=[("O=C=O", 1), "C=O"]  # Fix CO2 coefficient
)
rxn.balance()  # Now has unique solution
```

### Common SMILES for Reactions

| Molecule | SMILES | Notes |
|----------|--------|-------|
| Methane | `C` | Single carbon = CH4 |
| Oxygen (O2) | `O=O` | Double bond |
| Water | `O` | Single oxygen = H2O |
| Carbon dioxide | `O=C=O` | |
| Carbon monoxide | `[C-]#[O+]` or `C=O` | |
| Hydrogen | `[H][H]` | Explicit H atoms |
| Nitrogen | `N#N` | Triple bond |
| Ammonia | `N` | Single nitrogen = NH3 |
| Hydrogen peroxide | `OO` | |

### Debugging Balance Issues

Check if a reaction can be balanced before calling `balance()`:

```python
from phoenix import Reaction

rxn = Reaction.from_smiles(
    reactants=["C", "O=O"],
    products=["O=C=O", "O"]
)

# Before balancing
print(f"All Auto: {rxn.all_auto}")  # True if all coefficients are Auto
print(f"Elements: {rxn.elements}")   # {'C', 'H', 'O'}

# After balancing
rxn.balance()
print(f"Balanced: {rxn.is_balanced}")  # True
print(f"Coefficients: {rxn.coefficients}")
```

### When to Use Each Option

| Scenario | Recommended Approach |
|----------|---------------------|
| Simple reactions (combustion, hydrolysis) | Full auto-balance with no fixed coefficients |
| Known stoichiometry | Explicit coefficients for all species |
| Partial information | Mix of fixed coefficients and `Auto` |
| Multiple valid products | Fix at least one product coefficient |

### Verifying Results

Always verify thermodynamic results are reasonable:

```python
from phoenix import Reaction

rxn = Reaction.from_smiles(
    reactants=["C", "O=O"],
    products=["O=C=O", "O"]
)
rxn.balance()

# Sanity checks
print(f"Reaction: {rxn}")
print(f"ΔH_rxn = {rxn.delta_h:.1f} kJ/mol")  # -802.6 kJ/mol (exothermic)
print(f"ΔS_rxn = {rxn.delta_s:.1f} J/(mol·K)")  # 169.4 J/(mol·K)
print(f"ΔG_rxn = {rxn.delta_g:.1f} kJ/mol")  # -853.1 kJ/mol (spontaneous)

# Combustion should be exothermic
assert rxn.delta_h < 0, "Combustion should release heat"
assert rxn.is_exothermic, "Expected exothermic reaction"
```

!!! warning "Validation Required"
    PHOENIX thermodynamic values are estimates based on Benson Group Additivity. For critical applications, validate against experimental data or literature values.

---

## Next Steps

- [Batch Processing](batch-processing.md): Screen multiple compounds
- [Thermodynamics](thermodynamics.md): Property details
- [API Reference](../api/reaction.md): Reaction API
