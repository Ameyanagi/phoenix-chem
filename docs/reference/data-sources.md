# Data Sources

Thermodynamic data sources used by PHOENIX.

---

## Primary Sources

### pgradd (Benson Group Additivity)

**Library:** [pgradd](https://github.com/NREL/pgradd)

**Purpose:** Primary source for thermodynamic estimation

**Properties:**
- Enthalpy of formation (ΔHf°)
- Entropy (S°)
- Heat capacity (Cp)

**Data Coverage:**
- Organic compounds with common functional groups
- Temperature range: 298-1500 K (best accuracy)

**Method:** Benson Group Additivity with NASA polynomial coefficients

**Usage in PHOENIX:**
```python
from pgradd import ThermoChem
# PHOENIX wraps pgradd for Benson GA estimation
```

---

### chemicals (CalebBell/ChEDL)

**Library:** [chemicals](https://github.com/CalebBell/chemicals)

**Purpose:** Reference thermodynamic data

**Data Sources Included:**
- NIST WebBook
- Yaws' Handbook
- DIPPR database
- CRC Handbook
- Various literature sources

**Properties:**
- Formation enthalpies
- Standard entropies
- Heat capacities
- Vapor pressures
- Critical properties

**Usage in PHOENIX:**
```python
from chemicals import Hf, S0
# Used for reference value comparison
```

---

### NIST-JANAF Tables

**Source:** NIST Standard Reference Database 69

**Purpose:** High-accuracy reference data for common compounds

**Coverage:**
- ~2000 substances
- High-temperature data (0-6000 K for many)
- Experimental values with uncertainties

**Properties:**
- ΔHf°, S°, Cp, H°, G°
- Phase transition data

**Accuracy:** Generally ±1-5 kJ/mol for ΔHf°

---

## Decomposition Products

### Formation Enthalpies

Values used for decomposition calculations:

| Product | ΔHf° (kJ/mol) | Source |
|---------|---------------|--------|
| HF | -273.30 | NIST |
| N₂ | 0.00 | Reference |
| P₄O₁₀ | -2984.0 | CRC |
| H₂O (g) | -241.83 | NIST |
| CO₂ | -393.52 | NIST |
| SO₂ | -296.81 | NIST |
| CO | -110.53 | NIST |
| HCl | -92.31 | NIST |
| HBr | -36.29 | NIST |
| C (graphite) | 0.00 | Reference |
| H₂ | 0.00 | Reference |
| S (rhombic) | 0.00 | Reference |
| Br₂ (g) | 30.91 | NIST |

---

## Data Hierarchy

PHOENIX uses data in this priority order:

### For Thermodynamic Estimation

1. **pgradd (Benson GA)** - Primary estimation
2. **chemicals** - Fallback if Benson GA fails
3. **Hardcoded values** - For decomposition products

### For Reference Comparison

1. **NIST WebBook** (via chemicals)
2. **Yaws' Handbook** (via chemicals)
3. **Literature values**

---

## Accuracy Assessment

### Benson GA (pgradd)

| Compound Type | Typical Error (ΔHf°) |
|---------------|----------------------|
| Alkanes | ±2-5 kJ/mol |
| Alcohols, ethers | ±5-8 kJ/mol |
| Ketones, aldehydes | ±5-10 kJ/mol |
| Aromatic | ±5-10 kJ/mol |
| Nitro compounds | ±10-20 kJ/mol |
| Strained rings | ±10-30 kJ/mol |

### Reference Data (chemicals)

| Source | Typical Uncertainty |
|--------|---------------------|
| NIST experimental | ±1-5 kJ/mol |
| Yaws compilation | ±5-15 kJ/mol |
| Estimated values | ±10-30 kJ/mol |

---

## Data Limitations

### Benson GA Limitations

1. **Missing Groups**
   - Unusual bonding patterns
   - Rare functional groups
   - Some organometallic fragments

2. **Accuracy Degradation**
   - Strained ring systems
   - Multiple adjacent functional groups
   - Very large molecules (>100 atoms)

3. **Phase Corrections**
   - Primarily gas-phase data
   - Liquid/solid corrections less reliable

### Reference Data Limitations

1. **Coverage Gaps**
   - Many compounds lack experimental data
   - High-temperature data sparse

2. **Consistency Issues**
   - Different sources may disagree
   - Temperature/phase not always specified

---

## Validation

### Built-in Validation

PHOENIX can compare estimates to reference values:

```python
hf = compound.enthalpy_of_formation

if hf.has_reference():
    print(f"Estimated: {hf.value:.1f} kJ/mol")
    print(f"Reference: {hf.reference_value.value:.1f} kJ/mol")
    print(f"Deviation: {hf.deviation:.1f} kJ/mol")
```

### Recommended Validation

For safety-critical applications:

1. Compare Benson GA to experimental data
2. Check multiple reference sources
3. Validate with experimental measurements

---

## Adding Custom Data

PHOENIX does not currently support custom data injection.

For custom data needs:
1. Use RDKit for molecular processing
2. Implement custom estimation methods
3. Query external databases directly

---

## References

### pgradd

```bibtex
@software{pgradd,
  author = {NREL},
  title = {pgradd: Python Group Additivity},
  url = {https://github.com/NREL/pgradd}
}
```

### chemicals

```bibtex
@software{chemicals,
  author = {Bell, Caleb},
  title = {chemicals: Chemical properties component of ChEDL},
  url = {https://github.com/CalebBell/chemicals}
}
```

### NIST

```bibtex
@misc{nist_webbook,
  author = {NIST},
  title = {NIST Chemistry WebBook},
  url = {https://webbook.nist.gov/chemistry/}
}
```
