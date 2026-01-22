# Functional Group Detection

PHOENIX detects reactive functional groups that may indicate chemical hazards.

## Overview

Functional group detection uses SMARTS patterns to identify structural features associated with:

- Explosive potential
- Thermal instability
- Shock sensitivity
- Reactivity concerns

!!! warning "Screening Only"
    Functional group alerts are for screening purposes. The presence of a functional group does not guarantee hazardous behavior.

## Quick Usage

```python
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")  # Nitrobenzene
result = compound.evaluate_hazard()

# Check alerts in hazard result
if result.functional_group_alerts:
    print("Alerts detected:")
    for alert in result.functional_group_alerts:
        print(f"  - {alert}")
```

## Detailed Detection

For more information about detected groups:

```python
from phoenix.hazard.functional_groups import detect_functional_groups
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
alerts = detect_functional_groups(compound)

for alert in alerts:
    print(f"Name: {alert.name}")
    print(f"  Count: {alert.count}")
    print(f"  SMARTS: {alert.smarts}")
    print(f"  Risk: {alert.description}")
```

## FunctionalGroupAlert

Each detection returns a `FunctionalGroupAlert`:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Name of the functional group |
| `smarts` | str | SMARTS pattern used |
| `count` | int | Number of occurrences |
| `description` | str | Hazard description |

## Detected Groups

### Nitro Compounds

| Name | SMARTS | Risk |
|------|--------|------|
| Nitro group | `[N+](=O)[O-]` | Shock-sensitive, explosive |
| Aromatic nitro | `c[N+](=O)[O-]` | Potential explosive |

```python
# Examples
nitrobenzene = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
nitromethane = Compound.from_smiles("C[N+](=O)[O-]")
```

### Peroxides

| Name | SMARTS | Risk |
|------|--------|------|
| Peroxide | `[OX2][OX2]` | Heat/shock/friction sensitive |
| Organic peroxide | `[CX4][OX2][OX2][CX4]` | Explosive decomposition |
| Peroxy acid | `[CX3](=O)[OX2][OX2]` | Oxidizer, explosive |
| Hydroperoxide | `[OX2][OX2H]` | Unstable |

```python
# Examples
h2o2 = Compound.from_smiles("OO")
mekp = Compound.from_smiles("CC(C)(OO)c1ccccc1")  # MEKP-like
```

### Azo and Diazo Compounds

| Name | SMARTS | Risk |
|------|--------|------|
| Azo group | `[NX2]=[NX2]` | Gas-releasing decomposition |
| Diazo group | `[CX3]=[N+]=[N-]` | Highly reactive |
| Diazonium | `[N+]#N` | Shock-sensitive explosive |

```python
# Examples
azobenzene = Compound.from_smiles("c1ccc(N=Nc2ccccc2)cc1")
```

### Azides

| Name | SMARTS | Risk |
|------|--------|------|
| Azide group | `[N-]=[N+]=[N-]` | Primary explosive |
| Organic azide | `[CX4][N-]=[N+]=[N-]` | Detonation risk |

```python
# Example
methyl_azide = Compound.from_smiles("CN=[N+]=[N-]")
```

### N-Oxides and Nitroso

| Name | SMARTS | Risk |
|------|--------|------|
| N-oxide | `[NX3+]([O-])` | Oxidizer potential |
| Nitroso group | `[NX2]=O` | Unstable |

### Nitrates and Nitrites

| Name | SMARTS | Risk |
|------|--------|------|
| Nitrate ester | `[OX2][N+](=O)[O-]` | Explosive (nitroglycerin) |
| Nitrite ester | `[OX2][NX2]=O` | Unstable |

```python
# Example (nitroglycerin-like)
nitrate = Compound.from_smiles("CON(=O)=O")  # Methyl nitrate
```

### Nitrogen-Nitrogen Bonds

| Name | SMARTS | Risk |
|------|--------|------|
| Hydrazine | `[NX3H2][NX3H2]` | Toxic, explosive |
| Hydrazide | `[CX3](=O)[NX3][NX3]` | Violent decomposition |

```python
# Examples
hydrazine = Compound.from_smiles("NN")
```

### Strained Rings

| Name | SMARTS | Risk |
|------|--------|------|
| Epoxide | `C1OC1` | Reactive, violent polymerization |
| Aziridine | `C1NC1` | Strained, reactive |

```python
# Examples
ethylene_oxide = Compound.from_smiles("C1CO1")
```

### Acyl Halides and Anhydrides

| Name | SMARTS | Risk |
|------|--------|------|
| Acyl halide | `[CX3](=O)[F,Cl,Br,I]` | Highly reactive |
| Acid anhydride | `[CX3](=O)[OX2][CX3](=O)` | Water-reactive |

```python
# Examples
acetyl_chloride = Compound.from_smiles("CC(=O)Cl")
acetic_anhydride = Compound.from_smiles("CC(=O)OC(=O)C")
```

### Other

| Name | SMARTS | Risk |
|------|--------|------|
| Acetylide | `[CX2]#[CX2-]` | Shock-sensitive |

## Multiple Groups

Compounds with multiple reactive groups are flagged individually:

```python
# TNT has multiple nitro groups
tnt = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")
alerts = detect_functional_groups(tnt)

for alert in alerts:
    print(f"{alert.name}: {alert.count} occurrence(s)")
```

## Using in Hazard Evaluation

Functional groups trigger CHETAH Criterion 4:

```python
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
result = compound.evaluate_hazard()

# Check if criterion 4 triggered
if 4 in result.triggered_criteria:
    print("Criterion 4 triggered by functional groups:")
    for alert in result.functional_group_alerts:
        print(f"  - {alert}")
```

## Custom SMARTS Queries

For advanced analysis, use RDKit directly:

```python
from rdkit import Chem
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
mol = compound.rdmol

# Custom SMARTS pattern
custom_pattern = "[c][N+](=O)[O-]"  # Aromatic nitro
pattern = Chem.MolFromSmarts(custom_pattern)

matches = mol.GetSubstructMatches(pattern)
print(f"Found {len(matches)} matches")
```

## Combining with Other Criteria

Cross-reference functional groups with other hazard indicators:

```python
from phoenix import Compound

compound = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")
result = compound.evaluate_hazard()

print(f"Hazard Class: {result.hazard_class}")
print(f"Criteria: {result.triggered_criteria}")
print(f"ΔHd: {result.max_decomposition_cal_g:.0f} cal/g")
print(f"OB%: {result.oxygen_balance_percent:.1f}%")
print(f"Alerts: {result.functional_group_alerts}")

# Correlation check
if 1 in result.triggered_criteria and 4 in result.triggered_criteria:
    print("WARNING: High ΔHd confirmed by functional group presence")
```

## Batch Analysis of Alerts

```python
from phoenix import screen

smiles_list = [
    "c1ccccc1[N+](=O)[O-]",  # Nitrobenzene
    "NN",                     # Hydrazine
    "C1CO1",                  # Ethylene oxide
    "CCO",                    # Ethanol (no alerts)
]

results = screen(smiles_list)
df = results.dataframe

# Count compounds with alerts
has_alerts = df['alerts'].apply(lambda x: len(x) > 0 if x else False)
print(f"Compounds with alerts: {has_alerts.sum()}/{len(df)}")

# Most common alerts
from collections import Counter
all_alerts = []
for alerts in df['alerts'].dropna():
    all_alerts.extend(alerts)

print("\nAlert frequency:")
for alert, count in Counter(all_alerts).most_common():
    print(f"  {alert}: {count}")
```

---

## Next Steps

- [Hazard Evaluation](hazard-evaluation.md): Complete hazard assessment
- [Batch Processing](batch-processing.md): Screen multiple compounds
- [Error Handling](error-handling.md): Handle exceptions
