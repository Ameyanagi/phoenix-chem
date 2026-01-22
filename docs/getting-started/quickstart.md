# Quick Start

Get up and running with PHOENIX in 5 minutes.

!!! warning "Safety Notice"
    PHOENIX is for screening purposes only. Always validate results experimentally.

## Your First Compound

```python
from phoenix import Compound

# Create a compound from SMILES
ethanol = Compound.from_smiles("CCO")

# Basic properties
print(f"Formula: {ethanol.formula}")          # C2H6O
print(f"MW: {ethanol.molecular_weight:.2f}")  # 46.07 g/mol
print(f"Composition: {ethanol.composition}")  # {'C': 2, 'H': 6, 'O': 1}
```

## Thermodynamic Properties

```python
# Enthalpy of formation
hf = ethanol.enthalpy_of_formation
print(f"ΔHf° = {hf.value:.1f} {hf.unit}")  # ~-235 kJ/mol

# Entropy
s = ethanol.entropy
print(f"S° = {s.value:.1f} {s.unit}")  # ~290 J/(mol·K)

# Temperature-dependent properties
state = ethanol.thermo_at(T=400)  # at 400 K
print(f"H(400K) = {state.H.value:.1f} kJ/mol")
```

## Hazard Evaluation

```python
# Evaluate hazard for nitrobenzene
nitrobenzene = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")

# Quick hazard check
print(f"Oxygen Balance: {nitrobenzene.oxygen_balance:.1f}%")

# Full hazard evaluation
result = nitrobenzene.evaluate_hazard()
print(f"Hazard Class: {result.hazard_class}")           # HIGH, MEDIUM, or LOW
print(f"Max ΔHd: {result.max_decomposition_kJ_mol:.1f} kJ/mol")
```

## Decomposition Analysis

```python
# Analyze decomposition products
decomp = nitrobenzene.max_decomposition()

print(f"ΔHd = {decomp.delta_hd_kJ_mol:.1f} kJ/mol")
print(f"ΔHd = {decomp.delta_hd_cal_g:.1f} cal/g")

print("\nProducts:")
for product, moles in decomp.products.items():
    if moles > 0.01:
        print(f"  {product}: {moles:.2f} mol")
```

## Batch Screening

```python
from phoenix import screen

# Screen multiple compounds
smiles_list = [
    "CCO",                                    # Ethanol
    "CC(=O)C",                                # Acetone
    "c1ccccc1[N+](=O)[O-]",                   # Nitrobenzene
    "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",  # TNT
]

results = screen(smiles_list)

# Get DataFrame
df = results.dataframe
print(df[["formula", "hazard_class", "delta_hd_kJ_mol"]])

# Export to CSV
results.to_csv("screening_results.csv")
```

## Chemical Reactions

```python
from phoenix import Reaction, Auto

# Create a reaction with auto-balancing
rxn = Reaction.from_smiles(
    reactants=[("C", 1), ("O=O", Auto)],      # CH4 + O2
    products=[("O=C=O", Auto), ("O", Auto)],  # CO2 + H2O
)

print(f"Reaction: {rxn}")
print(f"ΔH_rxn = {rxn.delta_h:.1f} kJ/mol")
```

## Next Steps

- [Core Concepts](../user-guide/core-concepts.md): Understand oxygen balance, ΔHd, and CHETAH criteria
- [Hazard Evaluation](../user-guide/hazard-evaluation.md): Deep dive into hazard screening
- [Examples](../examples/index.md): More detailed examples
- [API Reference](../api/index.md): Complete API documentation
