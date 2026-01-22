# Examples

Practical examples demonstrating PHOENIX capabilities.

---

## Example Scripts

Example scripts are located in the `examples/` directory.

### Running Examples

```bash
# Navigate to project root
cd phoenix-chem

# Run an example
python examples/01_basic_compound.py
```

---

## Available Examples

### 01. Basic Compound Analysis

**File:** `examples/01_basic_compound.py`

Demonstrates:
- Creating compounds from SMILES
- Accessing molecular properties
- Getting thermodynamic data

```python
from phoenix import Compound

# Create compound
ethanol = Compound.from_smiles("CCO")

# Properties
print(f"Formula: {ethanol.formula}")
print(f"MW: {ethanol.molecular_weight:.2f} g/mol")
print(f"ΔHf°: {ethanol.delta_hf_kJ_mol:.1f} kJ/mol")
```

---

### 02. Hazard Evaluation

**File:** `examples/02_hazard_evaluation.py`

Demonstrates:
- Full hazard assessment
- Interpreting CHETAH criteria
- Accessing functional group alerts

```python
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
result = compound.evaluate_hazard()

print(f"Hazard Class: {result.hazard_class}")
print(f"Criteria: {result.triggered_criteria}")
print(f"ΔHd: {result.max_decomposition_cal_g:.1f} cal/g")
```

---

### 03. Decomposition Analysis

**File:** `examples/03_decomposition.py`

Demonstrates:
- Hierarchy vs LP methods
- Product distribution
- Gas generation calculations

```python
from phoenix import Compound

compound = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")

# Compare methods
comparison = compound.max_decomposition(method="both")
print(f"Hierarchy: {comparison.hierarchy_delta_hd:.1f} kJ/mol")
print(f"LP: {comparison.lp_delta_hd:.1f} kJ/mol")
print(f"Deviation: {comparison.deviation_percent:.2f}%")
```

---

### 04. Batch Screening

**File:** `examples/04_batch_screening.py`

Demonstrates:
- Processing multiple compounds
- Working with DataFrames
- Exporting results

```python
from phoenix import screen

smiles_list = ["CCO", "CC(=O)C", "c1ccccc1[N+](=O)[O-]"]
results = screen(smiles_list)

print(f"Screened: {results.successful}")
df = results.dataframe
print(df[["formula", "hazard_class", "delta_hd_cal_g"]])

results.to_csv("results.csv")
```

---

### 05. Reaction Thermodynamics

**File:** `examples/05_reactions.py`

Demonstrates:
- Creating balanced reactions
- Using Auto coefficients
- Calculating reaction thermodynamics

```python
from phoenix import Reaction, Auto

# Combustion
rxn = Reaction.from_smiles(
    reactants=["CH4", "O=O"],
    products=["O=C=O", "O"]
)
rxn.balance()

print(f"Reaction: {rxn}")
print(f"ΔH = {rxn.delta_h:.1f} kJ/mol")
print(f"ΔG = {rxn.delta_g:.1f} kJ/mol")
```

---

### 06. Temperature Dependence

**File:** `examples/06_temperature.py`

Demonstrates:
- Temperature-dependent properties
- Using ThermoState
- Vectorized calculations

```python
import numpy as np
from phoenix import Compound

compound = Compound.from_smiles("CCO")

# Single temperature
state = compound.thermo_at(T=500)
print(f"H(500K) = {state.H.value:.1f} kJ/mol")

# Temperature range
temps = np.linspace(300, 1000, 10)
hf_values = compound.enthalpy_of_formation(T=temps)
```

---

### 07. Error Handling

**File:** `examples/07_error_handling.py`

Demonstrates:
- Catching specific exceptions
- Graceful error recovery
- Batch processing with errors

```python
from phoenix import (
    Compound,
    InvalidSmilesError,
    UnsupportedElementError,
    PhoenixError,
)

smiles_list = ["CCO", "invalid", "[Fe]"]

for smiles in smiles_list:
    try:
        compound = Compound.from_smiles(smiles)
        print(f"{smiles}: {compound.formula}")
    except InvalidSmilesError:
        print(f"{smiles}: Invalid SMILES")
    except UnsupportedElementError as e:
        print(f"{smiles}: Unsupported elements {e.elements}")
```

---

### 08. NIST Reference Data

**File:** `examples/08_nist_data.py`

Demonstrates:
- Comparing estimates to reference data
- Accessing breakdown information
- Evaluating estimation accuracy

```python
from phoenix import Compound

compound = Compound.from_smiles("CCO")
hf = compound.enthalpy_of_formation

if hf.has_reference():
    print(f"Estimated: {hf.value:.1f} kJ/mol")
    print(f"Reference: {hf.reference_value.value:.1f} kJ/mol")
    print(f"Deviation: {hf.deviation:.1f} kJ/mol")
```

---

## Use Case Examples

### Screening a Chemical Library

```python
from phoenix import screen
import pandas as pd

# Read SMILES file
with open("compounds.smi") as f:
    smiles_list = [line.strip() for line in f]

# Screen all
results = screen(smiles_list)

# Filter high hazard
df = results.dataframe
high_hazard = df[df['hazard_class'] == 'HIGH']

print(f"High hazard: {len(high_hazard)}/{len(df)}")
high_hazard.to_csv("high_hazard_compounds.csv")
```

### Comparing Structural Isomers

```python
from phoenix import Compound

isomers = [
    ("CCO", "Ethanol"),
    ("COC", "Dimethyl ether"),
]

for smiles, name in isomers:
    c = Compound.from_smiles(smiles)
    print(f"\n{name} ({c.formula}):")
    print(f"  MW: {c.molecular_weight:.2f}")
    print(f"  ΔHf°: {c.delta_hf_kJ_mol:.1f} kJ/mol")
    print(f"  OB%: {c.oxygen_balance:.1f}%")
```

### Finding Optimal Oxygen Balance

```python
from phoenix import Compound

# Series of nitro compounds
compounds = [
    "c1ccccc1[N+](=O)[O-]",  # Nitrobenzene
    "c1cc([N+](=O)[O-])ccc1[N+](=O)[O-]",  # Dinitrobenzene
    "c1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",  # Trinitrobenzene
]

for smiles in compounds:
    c = Compound.from_smiles(smiles)
    print(f"{c.formula}: OB% = {c.oxygen_balance:.1f}%")
```

---

## Contributing Examples

To add new examples:

1. Create a new Python file in `examples/`
2. Follow the naming convention: `XX_description.py`
3. Include docstring explaining the example
4. Test that it runs successfully
5. Submit a pull request
