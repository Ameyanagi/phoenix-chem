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

## Intermediate Examples

These examples demonstrate combining multiple PHOENIX features for real-world workflows.

### Batch Screening with Hazard Summary

Screen multiple compounds and generate a summary report:

```python
from phoenix import screen, Compound
import pandas as pd

# Sample chemical library
smiles_list = [
    "CCO",                                    # Ethanol
    "CC(=O)C",                                # Acetone
    "CCCCCCCC",                               # Octane
    "c1ccccc1",                               # Benzene
    "c1ccccc1[N+](=O)[O-]",                   # Nitrobenzene
    "c1cc([N+](=O)[O-])ccc1[N+](=O)[O-]",     # 1,4-Dinitrobenzene
    "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",  # TNT
    "O=[N+]([O-])c1cc([N+](=O)[O-])c(N)c([N+](=O)[O-])c1", # Picramide
]

# Screen all compounds
results = screen(smiles_list)
df = results.dataframe

# Generate summary statistics
print("=== Hazard Screening Summary ===\n")
print(f"Total compounds screened: {len(df)}")
print(f"Successful: {results.successful}")
print(f"Failed: {results.failed}")

# Count by hazard class
print("\nHazard Class Distribution:")
hazard_counts = df['hazard_class'].value_counts()
for hazard_class, count in hazard_counts.items():
    pct = 100 * count / len(df)
    print(f"  {hazard_class}: {count} ({pct:.0f}%)")

# Filter high-hazard compounds
high_hazard = df[df['hazard_class'] == 'HIGH'].copy()
print(f"\n=== HIGH Hazard Compounds ({len(high_hazard)}) ===")
print(high_hazard[['formula', 'oxygen_balance', 'delta_hd_cal_g']].to_string(index=False))

# Find compounds in "dangerous" OB% range (-120 to +80)
df['ob_in_range'] = (df['oxygen_balance'] >= -120) & (df['oxygen_balance'] <= 80)
dangerous_ob = df[df['ob_in_range']]
print(f"\nCompounds with OB% in -120 to +80 range: {len(dangerous_ob)}")

# Export detailed report
df.to_csv("screening_report.csv", index=False)
print("\nFull report saved to screening_report.csv")
```

**Expected output:**
```
=== Hazard Screening Summary ===

Total compounds screened: 8
Successful: 8
Failed: 0

Hazard Class Distribution:
  HIGH: 4 (50%)
  LOW: 3 (38%)
  MEDIUM: 1 (12%)

=== HIGH Hazard Compounds (4) ===
  formula  oxygen_balance  delta_hd_cal_g
C6H5NO2          -162.5        -1072.0
C6H4N2O4          -96.5        -1055.9
C7H5N3O6          -74.0        -1362.4
C6H4N4O7          -36.9        -1498.3

Compounds with OB% in -120 to +80 range: 4

Full report saved to screening_report.csv
```

### Temperature Profile Analysis

Analyze thermodynamic properties across a temperature range:

```python
import numpy as np
from phoenix import Compound

# Create compound
ethanol = Compound.from_smiles("CCO")
print(f"Compound: {ethanol.formula} (Ethanol)")

# Define temperature range (in Kelvin)
temperatures = np.array([298.15, 400, 500, 600, 800, 1000])

# Calculate properties at each temperature
print("\n=== Temperature Profile ===")
print(f"{'T (K)':<10} {'ΔHf° (kJ/mol)':<15} {'S° (J/mol·K)':<15} {'Cp (J/mol·K)':<15} {'ΔGf° (kJ/mol)'}")
print("-" * 70)

for T in temperatures:
    state = ethanol.thermo_at(T=T)
    print(f"{T:<10.1f} {state.H.value:<15.1f} {state.S.value:<15.1f} {state.Cp.value:<15.1f} {state.G.value:.1f}")

# Vectorized calculation for plotting
T_range = np.linspace(300, 1200, 100)
hf_values = ethanol.enthalpy_of_formation(T=T_range)

print(f"\nEnthalpy range: {hf_values.min():.1f} to {hf_values.max():.1f} kJ/mol")
print(f"over temperature range: {T_range.min():.0f} to {T_range.max():.0f} K")

# Note: For plotting, you would use matplotlib:
# import matplotlib.pyplot as plt
# plt.plot(T_range, hf_values)
# plt.xlabel('Temperature (K)')
# plt.ylabel('ΔHf° (kJ/mol)')
# plt.savefig('ethanol_temperature_profile.png')
```

**Expected output:**
```
Compound: C2H6O (Ethanol)

=== Temperature Profile ===
T (K)      ΔHf° (kJ/mol)   S° (J/mol·K)    Cp (J/mol·K)    ΔGf° (kJ/mol)
----------------------------------------------------------------------
298.1      -235.1          289.9           81.0            -321.6
400.0      -227.8          319.8           92.3            -355.7
500.0      -218.7          346.2           109.5           -391.8
600.0      -207.9          370.3           126.0           -430.1
800.0      -181.8          413.3           151.9           -512.4
1000.0     -150.4          451.8           170.2           -602.2

Enthalpy range: -235.1 to -110.5 kJ/mol
over temperature range: 300 to 1200 K
```

### Comparing Decomposition Methods

Compare hierarchy and LP methods for a series of energetic compounds:

```python
from phoenix import Compound

# Energetic compounds
compounds = [
    ("c1ccccc1[N+](=O)[O-]", "Nitrobenzene"),
    ("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]", "TNT"),
    ("c1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]", "1,3,5-Trinitrobenzene"),
]

print("=== Decomposition Method Comparison ===\n")
print(f"{'Compound':<20} {'Hierarchy':<15} {'LP':<15} {'Deviation':<12} {'Preferred'}")
print("-" * 75)

for smiles, name in compounds:
    c = Compound.from_smiles(smiles)
    comparison = c.max_decomposition(method="both")

    hierarchy = comparison.hierarchy_delta_hd
    lp = comparison.lp_delta_hd
    deviation = comparison.deviation_percent

    # Choose method with more negative (more exothermic) result
    preferred = "Hierarchy" if hierarchy < lp else "LP"

    print(f"{name:<20} {hierarchy:<15.1f} {lp:<15.1f} {deviation:<12.2f}% {preferred}")

print("\n=== Product Comparison for TNT ===")
tnt = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")

# Hierarchy method products
hier_decomp = tnt.max_decomposition(method="hierarchy")
print("\nHierarchy method products:")
for product, moles in sorted(hier_decomp.products.items(), key=lambda x: -x[1]):
    if moles > 0.01:
        print(f"  {product}: {moles:.2f} mol")

# LP method products
lp_decomp = tnt.max_decomposition(method="lp")
print("\nLP method products:")
for product, moles in sorted(lp_decomp.products.items(), key=lambda x: -x[1]):
    if moles > 0.01:
        print(f"  {product}: {moles:.2f} mol")
```

### Reaction Thermodynamics Workflow

Calculate and compare thermodynamics for related reactions:

```python
from phoenix import Reaction, Compound

print("=== Combustion Series Analysis ===\n")

# Simple alkane combustion series
alkanes = [
    ("C", "Methane"),           # CH4
    ("CC", "Ethane"),           # C2H6
    ("CCC", "Propane"),         # C3H8
    ("CCCC", "Butane"),         # C4H10
]

print(f"{'Alkane':<12} {'ΔHcomb (kJ/mol)':<18} {'ΔHcomb/CH2 (kJ/mol)':<20} {'ΔGcomb (kJ/mol)'}")
print("-" * 70)

for smiles, name in alkanes:
    # Create combustion reaction
    rxn = Reaction.from_smiles(
        reactants=[smiles, "O=O"],
        products=["O=C=O", "O"]
    )
    rxn.balance()

    # Get compound for atom count
    compound = Compound.from_smiles(smiles)
    n_carbons = compound.composition.get('C', 0)

    # Per CH2 unit (approximate)
    dh_per_ch2 = rxn.delta_h / n_carbons

    print(f"{name:<12} {rxn.delta_h:<18.1f} {dh_per_ch2:<20.1f} {rxn.delta_g:.1f}")

# Compare exothermicity
print("\n=== Exothermicity Comparison ===")
print("All combustion reactions are highly exothermic (ΔH << 0)")
print("Heat release per carbon is approximately constant (~-655 kJ/mol C)")
```

**Expected output:**
```
=== Combustion Series Analysis ===

Alkane       ΔHcomb (kJ/mol)    ΔHcomb/CH2 (kJ/mol)  ΔGcomb (kJ/mol)
----------------------------------------------------------------------
Methane      -802.6             -802.6               -853.1
Ethane       -1427.9            -714.0               -1464.4
Propane      -2043.2            -681.1               -2073.5
Butane       -2658.3            -664.6               -2682.5

=== Exothermicity Comparison ===
All combustion reactions are highly exothermic (ΔH << 0)
Heat release per carbon is approximately constant (~-655 kJ/mol C)
```

---

## Contributing Examples

To add new examples:

1. Create a new Python file in `examples/`
2. Follow the naming convention: `XX_description.py`
3. Include docstring explaining the example
4. Test that it runs successfully
5. Submit a pull request
