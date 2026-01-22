# Core Concepts

Understanding the fundamental concepts behind PHOENIX hazard evaluation.

!!! warning "Screening Tool"
    PHOENIX is for screening purposes only. All results require experimental validation.

## Oxygen Balance (OB%)

Oxygen balance indicates whether a compound has excess or deficient oxygen for complete combustion.

### Formula

$$
OB\% = \frac{-1600}{MW} \times \left( 2C + \frac{H}{2} + M - X - O \right)
$$

Where:
- $MW$ = Molecular weight (g/mol)
- $C$ = Number of carbon atoms
- $H$ = Number of hydrogen atoms
- $M$ = Number of metal atoms (assumed monovalent)
- $X$ = Number of halogen atoms
- $O$ = Number of oxygen atoms

### Interpretation

| OB% Range | Meaning | Risk |
|-----------|---------|------|
| OB% > 0 | Oxygen excess | Oxidizer potential |
| OB% ≈ 0 | Balanced | Maximum energy release |
| OB% < 0 | Oxygen deficient | Incomplete combustion |

### Example

```python
from phoenix import Compound

# TNT has significant oxygen deficiency
tnt = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")
print(f"TNT OB%: {tnt.oxygen_balance:.1f}%")  # ~-74%

# RDX is closer to balanced
rdx = Compound.from_smiles("O=[N+]([O-])N1CN([N+](=O)[O-])CN([N+](=O)[O-])C1")
print(f"RDX OB%: {rdx.oxygen_balance:.1f}%")  # ~-22%
```

Compounds with OB% between -200% and +100% trigger **CHETAH Criterion 3**.

---

## Enthalpy of Formation (ΔHf°)

The standard enthalpy of formation is the heat change when one mole of a compound is formed from its elements in their standard states.

### Sources

PHOENIX estimates ΔHf° using multiple methods:

1. **Benson Group Additivity (GA)** - Primary method
   - Decomposes molecules into functional groups
   - Sums pre-tabulated group contributions
   - Typical accuracy: ±5-10 kJ/mol for common organics

2. **Reference Data** - When available
   - NIST-JANAF thermochemical tables
   - Experimental literature values

### Usage

```python
from phoenix import Compound

ethanol = Compound.from_smiles("CCO")

# As property (at 298.15 K)
hf = ethanol.enthalpy_of_formation
print(f"ΔHf° = {hf.value:.1f} {hf.unit}")  # kJ/mol

# At different temperature
hf_500 = ethanol.enthalpy_of_formation(T=500)
print(f"ΔHf°(500K) = {hf_500.value:.1f} kJ/mol")
```

---

## Maximum Heat of Decomposition (ΔHd)

The maximum heat of decomposition represents the maximum possible energy release when a compound decomposes to its most stable products.

### Calculation

$$
\Delta H_d = \Delta H_f(\text{products}) - \Delta H_f(\text{reactant})
$$

Negative values indicate exothermic decomposition (energy release).

### Decomposition Hierarchy

PHOENIX uses thermodynamic priority rules to determine products:

| Priority | Reaction | Product |
|----------|----------|---------|
| 1 | F + H → HF | Hydrogen fluoride |
| 2 | N → ½N₂ | Nitrogen gas |
| 3 | P + O → ¼P₄O₁₀ | Phosphorus pentoxide |
| 4 | H + O → ½H₂O | Water |
| 5 | C + O → ½CO₂ | Carbon dioxide |
| 6 | S + O → SO₂ | Sulfur dioxide |
| 7 | C + ½O → CO | Carbon monoxide |
| 8 | Cl + H → HCl | Hydrogen chloride |
| 9 | Br + H → HBr | Hydrogen bromide |
| 10 | C → C(s) | Graphite |
| 11 | H → ½H₂ | Hydrogen gas |

### Example

```python
from phoenix import Compound

nitrobenzene = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
decomp = nitrobenzene.max_decomposition()

print(f"ΔHd = {decomp.delta_hd_kJ_mol:.1f} kJ/mol")
print(f"ΔHd = {decomp.delta_hd_cal_g:.1f} cal/g")
```

---

## CHETAH Criteria

PHOENIX implements the CHETAH hazard screening criteria from ASTM E659.

### Criterion 1: High Instability

**ΔHd < -300 cal/g** (-1255 kJ/kg)

Compounds meeting this criterion have high potential for violent decomposition.

### Criterion 2: Medium Instability

**-300 < ΔHd < -100 cal/g**

Compounds in this range have moderate decomposition energy.

### Criterion 3: Oxygen Balance Concern

**-200% < OB% < +100%**

Compounds in this OB range can participate in self-oxidation reactions.

### Criterion 4: Functional Group Alerts

Presence of known reactive functional groups:

- Nitro groups (-NO₂)
- Peroxides (-O-O-)
- Azides (-N₃)
- Azo compounds (-N=N-)
- Diazonium salts (-N₂⁺)
- Hydrazines (-NH-NH₂)
- Epoxides (strained rings)

### Hazard Classification

| Criteria Triggered | Classification |
|-------------------|----------------|
| Criterion 1 | **HIGH** |
| Criterion 2 | **MEDIUM** |
| Criteria 3 + 4 | **MEDIUM** |
| Any single criterion | **MEDIUM** |
| None | **LOW** |

---

## Phase States

PHOENIX supports calculations for different phases:

| Phase | Code | Description |
|-------|------|-------------|
| Gas | `g` | Default for Benson GA |
| Liquid | `l` | Liquid state corrections |
| Solid | `s` | Solid state corrections |

```python
from phoenix import Compound

# Gas phase (default)
water_g = Compound.from_smiles("O", phase="g")

# Liquid phase
water_l = Compound.from_smiles("O", phase="l")
```

!!! note "Phase Limitations"
    Phase corrections are available for some properties but not all.
    Most hazard calculations assume gas-phase products.

---

## Units

### Standard Units in PHOENIX

| Property | Unit |
|----------|------|
| Enthalpy | kJ/mol |
| Entropy | J/(mol·K) |
| Heat capacity | J/(mol·K) |
| Gibbs energy | kJ/mol |
| ΔHd (molar) | kJ/mol |
| ΔHd (specific) | cal/g |
| Oxygen balance | % |
| Molecular weight | g/mol |
| Gas volume | L/g |

### Unit Conversions

```python
# kJ/mol to cal/g
delta_hd_kJ_mol = -500
mw = 200  # g/mol
delta_hd_cal_g = (delta_hd_kJ_mol * 1000) / (4.184 * mw)
print(f"ΔHd = {delta_hd_cal_g:.1f} cal/g")  # -597.6 cal/g
```

---

## Next Steps

- [Compounds](compounds.md): Learn to create and analyze compounds
- [Hazard Evaluation](hazard-evaluation.md): Full hazard assessment workflow
- [Decomposition](decomposition.md): Detailed decomposition analysis
