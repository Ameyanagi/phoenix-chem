# Hazard Evaluation

PHOENIX implements CHETAH-style hazard screening based on ASTM E659.

!!! warning "Screening Only"
    Hazard evaluation is for screening purposes only. **Always validate results experimentally** before handling potentially hazardous materials.

## Quick Evaluation

```python
from phoenix import Compound

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")  # Nitrobenzene
result = compound.evaluate_hazard()

print(f"Hazard Class: {result.hazard_class}")
print(f"Triggered Criteria: {result.triggered_criteria}")
```

## HazardResult

The `evaluate_hazard()` method returns a `HazardResult` dataclass:

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
| `hazard_class` | str | HIGH, MEDIUM, or LOW |
| `triggered_criteria` | tuple[int] | List of triggered criteria (1-4) |
| `functional_group_alerts` | tuple[str] | Names of detected groups |
| `product_breakdown` | dict | Decomposition products |
| `gas_volume_L_g` | float | Gas generation in L/g |

### Example

```python
result = compound.evaluate_hazard()

# Access all properties
print(f"Formula: {result.formula}")
print(f"MW: {result.molecular_weight:.2f} g/mol")
print(f"ΔHf°: {result.delta_hf_kJ_mol:.1f} kJ/mol")
print(f"ΔHd: {result.max_decomposition_cal_g:.1f} cal/g")
print(f"OB%: {result.oxygen_balance_percent:.1f}%")
print(f"Hazard Class: {result.hazard_class}")
print(f"Criteria: {result.triggered_criteria}")
print(f"Alerts: {result.functional_group_alerts}")
print(f"Gas: {result.gas_volume_L_g:.3f} L/g")
```

## CHETAH Criteria

### Criterion 1: High Instability

**ΔHd < -300 cal/g**

Compounds with very high decomposition energy.

```python
# Example: Check for criterion 1
if result.max_decomposition_cal_g < -300:
    print("HIGH hazard: Criterion 1 triggered")
```

### Criterion 2: Medium Instability

**-300 ≤ ΔHd < -100 cal/g**

Compounds with moderate decomposition energy.

```python
# Example: Check for criterion 2
if -300 <= result.max_decomposition_cal_g < -100:
    print("MEDIUM hazard: Criterion 2 triggered")
```

### Criterion 3: Oxygen Balance

**-200% < OB% < +100%**

Compounds with reactive oxygen balance.

```python
# Example: Check for criterion 3
if -200 < result.oxygen_balance_percent < 100:
    print("OB% in reactive range: Criterion 3 triggered")
```

### Criterion 4: Functional Groups

Presence of known reactive functional groups.

```python
# Example: Check for criterion 4
if result.functional_group_alerts:
    print("Functional groups detected: Criterion 4 triggered")
    for alert in result.functional_group_alerts:
        print(f"  - {alert}")
```

## Hazard Classification Logic

```
IF Criterion 1 triggered:
    → HIGH

ELSE IF Criterion 2 triggered:
    → MEDIUM

ELSE IF Criteria 3 AND 4 triggered:
    → MEDIUM

ELSE IF ANY criterion triggered:
    → MEDIUM

ELSE:
    → LOW
```

## Examples by Hazard Class

### HIGH Hazard

```python
# TNT - Trinitrotoluene
tnt = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")
result = tnt.evaluate_hazard()
print(f"TNT: {result.hazard_class}")  # HIGH
print(f"ΔHd = {result.max_decomposition_cal_g:.0f} cal/g")
```

### MEDIUM Hazard

```python
# Nitrobenzene
nitrobenzene = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
result = nitrobenzene.evaluate_hazard()
print(f"Nitrobenzene: {result.hazard_class}")  # MEDIUM
```

### LOW Hazard

```python
# Ethanol
ethanol = Compound.from_smiles("CCO")
result = ethanol.evaluate_hazard()
print(f"Ethanol: {result.hazard_class}")  # LOW
```

## Interpreting Results

### High Hazard Compounds

- **Immediate concern**: High decomposition energy
- **Handling**: Extreme caution, specialized facilities
- **Validation**: Essential experimental characterization

### Medium Hazard Compounds

- **Potential concern**: One or more criteria triggered
- **Handling**: Standard chemical safety procedures
- **Validation**: Recommended for new compounds

### Low Hazard Compounds

- **General chemicals**: No CHETAH criteria triggered
- **Handling**: Normal laboratory safety
- **Note**: Still follow proper chemical safety

## Decomposition Products

Access the predicted decomposition products:

```python
result = compound.evaluate_hazard()

print("Decomposition products:")
for product, moles in result.product_breakdown.items():
    if moles > 0.01:
        print(f"  {product}: {moles:.2f} mol")
```

## Gas Generation

Predict gas generation from decomposition:

```python
result = compound.evaluate_hazard()

print(f"Gas volume: {result.gas_volume_L_g:.4f} L/g at STP")
```

This is useful for:

- Containment requirements
- Pressure buildup assessment
- Ventilation planning

## Working with Multiple Compounds

### Manual Loop

```python
smiles_list = [
    "CCO",                                    # Ethanol
    "c1ccccc1[N+](=O)[O-]",                   # Nitrobenzene
    "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",  # TNT
]

for smiles in smiles_list:
    compound = Compound.from_smiles(smiles)
    result = compound.evaluate_hazard()
    print(f"{result.formula}: {result.hazard_class}")
```

### Using Batch Screening

For many compounds, use batch processing:

```python
from phoenix import screen

results = screen(smiles_list)
df = results.dataframe

# Filter by hazard class
high_hazard = df[df['hazard_class'] == 'HIGH']
print(f"High hazard compounds: {len(high_hazard)}")
```

See [Batch Processing](batch-processing.md) for more details.

## Integration with Other Analysis

### Combining with Decomposition Analysis

```python
compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")

# Quick hazard check
result = compound.evaluate_hazard()

# Detailed decomposition analysis
if result.hazard_class != "LOW":
    decomp = compound.max_decomposition(method="both")
    print(f"Hierarchy: {decomp.hierarchy_result.delta_hd_cal_g:.1f} cal/g")
    print(f"LP: {decomp.lp_result.delta_hd_cal_g:.1f} cal/g")
    print(f"Deviation: {decomp.deviation_percent:.1f}%")
```

### Functional Group Details

```python
from phoenix.hazard.functional_groups import detect_functional_groups

compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
alerts = detect_functional_groups(compound)

for alert in alerts:
    print(f"{alert.name} (x{alert.count})")
    print(f"  SMARTS: {alert.smarts}")
    print(f"  Risk: {alert.description}")
```

---

## Best Practices

1. **Always validate experimentally** for new or unfamiliar compounds
2. **Use batch screening** for large datasets
3. **Check functional groups** for specific hazard mechanisms
4. **Consider context** - LOW doesn't mean "safe"
5. **Document results** for safety reviews

## Next Steps

- [Decomposition](decomposition.md): Detailed decomposition analysis
- [Functional Groups](functional-groups.md): Reactive group detection
- [Batch Processing](batch-processing.md): Screen multiple compounds
