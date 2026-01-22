# Working with Compounds

The `Compound` class is the primary entry point for PHOENIX.

## Creating Compounds

### From SMILES

The standard way to create a compound is from a SMILES string:

```python
from phoenix import Compound

# Simple molecules
ethanol = Compound.from_smiles("CCO")
acetone = Compound.from_smiles("CC(=O)C")

# Aromatic compounds
benzene = Compound.from_smiles("c1ccccc1")
nitrobenzene = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")

# Complex structures
tnt = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")
```

### With Phase Specification

Specify the phase for more accurate thermodynamic calculations:

```python
# Gas phase (default)
water_gas = Compound.from_smiles("O", phase="g")

# Liquid phase
water_liquid = Compound.from_smiles("O", phase="l")

# Solid phase
water_solid = Compound.from_smiles("O", phase="s")
```

## Basic Properties

### Molecular Information

```python
compound = Compound.from_smiles("CCO")

# Molecular formula (Hill notation)
print(compound.formula)          # C2H6O

# Molecular weight in g/mol
print(compound.molecular_weight)  # 46.07

# Elemental composition
print(compound.composition)       # {'C': 2, 'H': 6, 'O': 1}

# SMILES representations
print(compound.original_smiles)   # CCO (as provided)
print(compound.canonical_smiles)  # CCO (canonicalized)

# InChIKey for database lookup
print(compound.inchikey)          # LFQSCWFLJHTTHZ-UHFFFAOYSA-N
```

### Atom Count

```python
compound = Compound.from_smiles("Cc1ccccc1")  # Toluene

# Total atoms including hydrogens
print(compound.num_atoms)  # 15 (C7H8)
```

## Thermodynamic Properties

### Enthalpy of Formation

```python
compound = Compound.from_smiles("CCO")

# Access as property (at 298.15 K)
hf = compound.enthalpy_of_formation
print(f"ΔHf° = {hf.value:.1f} {hf.unit}")

# Access the raw float value
print(f"ΔHf° = {compound.delta_hf_kJ_mol:.1f} kJ/mol")

# At different temperatures
hf_400 = compound.enthalpy_of_formation(T=400)
hf_500 = compound.enthalpy_of_formation(T=500)
```

### Entropy

```python
compound = Compound.from_smiles("CCO")

# Standard entropy
s = compound.entropy
print(f"S° = {s.value:.1f} {s.unit}")  # J/(mol·K)

# At different temperature
s_500 = compound.entropy(T=500)
```

### Heat Capacity

```python
compound = Compound.from_smiles("CCO")

# Heat capacity at 298.15 K
cp = compound.heat_capacity()
print(f"Cp = {cp.value:.1f} {cp.unit}")

# At elevated temperature
cp_500 = compound.heat_capacity(temperature_K=500)
```

### Thermodynamic State

Get all properties at once using `thermo_at()`:

```python
compound = Compound.from_smiles("CCO")

# All properties at 500 K
state = compound.thermo_at(T=500)

print(f"H = {state.H.value:.1f} kJ/mol")
print(f"S = {state.S.value:.1f} J/(mol·K)")
print(f"Cp = {state.Cp.value:.1f} J/(mol·K)")
print(f"G = {state.G.value:.1f} kJ/mol")
```

## Hazard Properties

### Oxygen Balance

```python
compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")

# Oxygen balance as percentage
ob = compound.oxygen_balance
print(f"OB% = {ob:.1f}%")  # Negative = oxygen deficient
```

### Maximum Decomposition

```python
compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")

# Calculate decomposition
decomp = compound.max_decomposition()

print(f"ΔHd = {decomp.delta_hd_kJ_mol:.1f} kJ/mol")
print(f"ΔHd = {decomp.delta_hd_cal_g:.1f} cal/g")

# Decomposition products
for product, moles in decomp.products.items():
    if moles > 0.01:
        print(f"  {product}: {moles:.2f} mol")
```

### Full Hazard Evaluation

```python
compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")

# Complete hazard assessment
result = compound.evaluate_hazard()

print(f"Hazard Class: {result.hazard_class}")  # HIGH, MEDIUM, or LOW
print(f"Triggered Criteria: {result.triggered_criteria}")
print(f"Functional Group Alerts: {result.functional_group_alerts}")
```

## Warnings and Validation

### Automatic Validation

PHOENIX validates compounds during creation:

```python
from phoenix import (
    Compound,
    InvalidSmilesError,
    UnsupportedElementError,
    UnsupportedStructureError,
)

# Invalid SMILES raises InvalidSmilesError
try:
    Compound.from_smiles("not-a-smiles")
except InvalidSmilesError as e:
    print(f"Invalid SMILES: {e}")

# Unsupported elements raise UnsupportedElementError
try:
    Compound.from_smiles("[Fe]")  # Iron not supported
except UnsupportedElementError as e:
    print(f"Unsupported: {e.elements}")

# Charged species raise UnsupportedStructureError
try:
    Compound.from_smiles("[NH4+]")  # Charged species
except UnsupportedStructureError as e:
    print(f"Unsupported structure: {e.reason}")
```

### Large Molecule Warnings

PHOENIX warns for molecules with more than 100 atoms:

```python
import warnings

# Large molecules generate warnings
large_molecule = Compound.from_smiles("C" * 50)  # Very large
# UserWarning: Large molecule (150 atoms). Benson GA accuracy may be reduced...
```

### Accessing Warnings

```python
compound = Compound.from_smiles("C" * 50)

# Check for warnings generated during creation
for warning in compound.warnings:
    print(f"Warning: {warning}")
```

## Supported Elements

PHOENIX supports organic compounds containing:

| Element | Symbol |
|---------|--------|
| Carbon | C |
| Hydrogen | H |
| Nitrogen | N |
| Oxygen | O |
| Sulfur | S |
| Phosphorus | P |
| Fluorine | F |
| Chlorine | Cl |
| Bromine | Br |

Elements **not** supported: metals, iodine, silicon, boron, etc.

## RDKit Integration

Access the underlying RDKit molecule for advanced operations:

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

compound = Compound.from_smiles("CCO")

# Get RDKit Mol object
rdmol = compound.rdmol

# Use RDKit functions
logp = Descriptors.MolLogP(rdmol)
hbd = Descriptors.NumHDonors(rdmol)

print(f"LogP: {logp:.2f}")
print(f"H-bond donors: {hbd}")
```

## String Representations

```python
compound = Compound.from_smiles("c1ccccc1")

# repr shows canonical SMILES
print(repr(compound))  # Compound('c1ccccc1')

# str shows formula and MW
print(str(compound))   # C6H6 (MW=78.11)
```

## Equality and Hashing

Compounds are compared by canonical SMILES:

```python
# Same molecule, different SMILES
c1 = Compound.from_smiles("CCO")
c2 = Compound.from_smiles("OCC")

print(c1 == c2)  # True (same canonical SMILES)
print(hash(c1) == hash(c2))  # True

# Use in sets and dicts
compounds = {c1, c2}
print(len(compounds))  # 1 (deduplicated)
```

---

## Next Steps

- [Thermodynamics](thermodynamics.md): Deep dive into thermodynamic calculations
- [Hazard Evaluation](hazard-evaluation.md): Complete hazard assessment
- [Batch Processing](batch-processing.md): Screen multiple compounds
