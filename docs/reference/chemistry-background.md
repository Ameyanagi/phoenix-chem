# Chemistry Background

Scientific basis for PHOENIX calculations.

---

## Oxygen Balance

### Definition

Oxygen balance (OB%) indicates the excess or deficiency of oxygen in a compound for complete oxidation to CO₂, H₂O, SO₂, P₄O₁₀, etc.

### Formula

$$
OB\% = \frac{-1600}{MW} \times \left( 2a + \frac{b}{2} + M - c - d \right)
$$

Where for $C_a H_b O_c X_d$ (with halogen X):

| Symbol | Meaning |
|--------|---------|
| $a$ | Carbon atoms |
| $b$ | Hydrogen atoms |
| $c$ | Oxygen atoms |
| $d$ | Halogen atoms |
| $M$ | Metal atoms (assumed monovalent) |
| $MW$ | Molecular weight (g/mol) |

### Complete Formula (PHOENIX Implementation)

For compounds $C_a H_b N_c O_d S_e P_f F_g Cl_h Br_i$:

$$
OB\% = \frac{-1600}{MW} \times \left( 2a + \frac{b}{2} + \frac{5f}{2} + 2e - \frac{g}{2} - \frac{h}{2} - \frac{i}{2} - d \right)
$$

### Interpretation

| OB% | Condition | Example |
|-----|-----------|---------|
| OB% > 0 | Oxygen excess (oxidizer) | NH₄NO₃ (+20%) |
| OB% ≈ 0 | Balanced | RDX (-22%) |
| OB% < 0 | Oxygen deficient | TNT (-74%) |

Compounds with OB% close to zero release maximum energy on decomposition.

---

## Maximum Heat of Decomposition

### Definition

The maximum heat of decomposition (ΔHd) is the theoretical maximum energy release when a compound decomposes to its most thermodynamically stable products.

### Calculation

$$
\Delta H_d = \sum_i n_i \cdot \Delta H_f^{\circ}(\text{product}_i) - \Delta H_f^{\circ}(\text{reactant})
$$

Where:
- $n_i$ = moles of product $i$
- $\Delta H_f^{\circ}$ = standard enthalpy of formation

### Sign Convention

| Sign | Meaning |
|------|---------|
| ΔHd < 0 | Exothermic (energy release) |
| ΔHd > 0 | Endothermic (energy absorption) |

### Thermodynamic Hierarchy

PHOENIX uses a priority-based hierarchy to determine decomposition products:

| Priority | Reaction | ΔHf° Product | Reason |
|----------|----------|--------------|--------|
| 1 | F + H → HF | -273.3 kJ/mol | Most exothermic H-X bond |
| 2 | N → ½N₂ | 0.0 kJ/mol | Always forms N₂ |
| 3 | P + O → ¼P₄O₁₀ | -746.0 kJ/mol (per P) | Highly exothermic |
| 4 | H + O → ½H₂O | -241.8 kJ/mol | Water formation |
| 5 | C + O → ½CO₂ | -393.5 kJ/mol | Full oxidation |
| 6 | S + O → SO₂ | -296.8 kJ/mol | Sulfur oxidation |
| 7 | C + ½O → CO | -110.5 kJ/mol | Partial oxidation |
| 8 | Cl + H → HCl | -92.3 kJ/mol | After water |
| 9 | Br + H → HBr | -36.3 kJ/mol | After HCl |
| 10 | C → C(s) | 0.0 kJ/mol | Graphite |
| 11 | H → ½H₂ | 0.0 kJ/mol | Hydrogen gas |

### Unit Conversion

$$
\Delta H_d \text{ (cal/g)} = \frac{\Delta H_d \text{ (kJ/mol)} \times 1000}{4.184 \times MW}
$$

---

## CHETAH Criteria

ASTM E659 defines four hazard screening criteria.

### Criterion 1: High Instability

$$
\Delta H_d < -300 \text{ cal/g} \quad (-1255 \text{ kJ/kg})
$$

Indicates high potential for violent decomposition.

### Criterion 2: Medium Instability

$$
-300 \text{ cal/g} \leq \Delta H_d < -100 \text{ cal/g}
$$

Moderate decomposition energy.

### Criterion 3: Oxygen Balance Concern

$$
-200\% < OB\% < +100\%
$$

Compounds that can participate in self-oxidation reactions.

### Criterion 4: Functional Group Alerts

Presence of functional groups associated with:
- Explosive potential
- Thermal instability
- Shock sensitivity

### Classification Logic

```
IF Criterion 1: HIGH
ELSE IF Criterion 2: MEDIUM
ELSE IF Criteria 3 AND 4: MEDIUM
ELSE IF ANY criterion: MEDIUM
ELSE: LOW
```

---

## Benson Group Additivity

### Theory

Benson GA estimates thermodynamic properties by summing contributions from molecular fragments (groups).

$$
\Delta H_f^{\circ} = \sum_j n_j \cdot h_j
$$

Where:
- $n_j$ = occurrences of group $j$
- $h_j$ = enthalpy contribution of group $j$

### Group Notation

Groups are defined by central atom and neighbors:

| Notation | Meaning |
|----------|---------|
| C-(H)₃(C) | Carbon bonded to 3 H and 1 C |
| C-(H)₂(C)₂ | Carbon bonded to 2 H and 2 C |
| Cb-H | Benzene carbon bonded to H |
| O-(H)(C) | Oxygen in alcohol |

### Accuracy

| Property | Typical Uncertainty |
|----------|---------------------|
| ΔHf° | ±5-10 kJ/mol |
| S° | ±5-10 J/(mol·K) |
| Cp | ±5-10 J/(mol·K) |

### Limitations

- Strained ring corrections needed
- Unusual bonding patterns may lack data
- Large molecules (>100 atoms) less accurate

---

## Gas Generation

### Ideal Gas Law

Gas volume is calculated using:

$$
V = \frac{n \cdot R \cdot T}{P \cdot MW}
$$

Where:
- $V$ = volume per gram (L/g)
- $n$ = moles of gas per mole of compound
- $R$ = 0.08206 L·atm/(mol·K)
- $T$ = temperature (K)
- $P$ = 1 atm
- $MW$ = molecular weight (g/mol)

### At Standard Conditions (298.15 K)

$$
V_m = R \cdot T = 24.47 \text{ L/mol}
$$

### Gas Products

| Product | Phase |
|---------|-------|
| N₂ | Gas |
| H₂O | Gas |
| CO₂ | Gas |
| CO | Gas |
| HF, HCl, HBr | Gas |
| H₂ | Gas |
| SO₂ | Gas |
| O₂, F₂, Cl₂, Br₂ | Gas |

Non-gas products: C (graphite), S (rhombic), P₄O₁₀

---

## References

### Primary Sources

1. **ASTM E659** - Standard Test Method for Determining Thermodynamic Properties, Reactivities, and Initiation Hazards of Materials

2. **Benson, S.W.** - "Thermochemical Kinetics" 2nd Ed., Wiley, 1976

3. **NIST-JANAF** - Thermochemical Tables, 4th Ed., 1998

### Additional References

4. **Meyer et al.** - "Explosives" 6th Ed., Wiley-VCH, 2007

5. **Lide, D.R.** - CRC Handbook of Chemistry and Physics

6. **CalebBell/chemicals** - Python library for thermodynamic data
