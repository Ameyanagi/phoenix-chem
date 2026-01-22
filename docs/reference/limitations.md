# Limitations

Understanding PHOENIX's scope and constraints.

!!! warning "Critical Limitation"
    PHOENIX is a **screening tool only**. Results must not be the sole basis for safety decisions. Experimental validation is required.

---

## Supported Elements

### Included

| Element | Symbol | Notes |
|---------|--------|-------|
| Carbon | C | Full support |
| Hydrogen | H | Full support |
| Nitrogen | N | Full support |
| Oxygen | O | Full support |
| Sulfur | S | Full support |
| Phosphorus | P | Full support |
| Fluorine | F | Full support |
| Chlorine | Cl | Full support |
| Bromine | Br | Full support |

### Not Supported

| Category | Elements |
|----------|----------|
| Halogens | Iodine (I) |
| Metals | All (Fe, Cu, Zn, etc.) |
| Metalloids | Silicon (Si), Boron (B) |
| Other | Selenium (Se), Tellurium (Te) |

Attempting to use unsupported elements raises `UnsupportedElementError`.

---

## Structural Limitations

### Not Supported

| Structure | Example | Reason |
|-----------|---------|--------|
| Charged species | NH₄⁺, Cl⁻ | Benson GA not applicable |
| Radicals | CH₃• | Electronic effects |
| Coordination compounds | Fe(CO)₅ | Metal bonding |
| Polymers | Large repeating units | Size limitations |

### Partially Supported

| Structure | Limitation |
|-----------|------------|
| Strained rings | May lack correction data |
| Unusual bonding | Missing group contributions |
| Large molecules | Accuracy degrades >100 atoms |

---

## Accuracy Constraints

### Enthalpy of Formation (ΔHf°)

| Compound Type | Typical Uncertainty |
|---------------|---------------------|
| Simple hydrocarbons | ±5 kJ/mol |
| Common organics | ±5-10 kJ/mol |
| Nitro compounds | ±10-20 kJ/mol |
| Complex structures | ±15-30 kJ/mol |

### Maximum Decomposition (ΔHd)

| Source | Contribution to Error |
|--------|----------------------|
| Reactant ΔHf° uncertainty | ±5-20 kJ/mol |
| Product selection | Variable |
| Phase assumptions | ~5-10 kJ/mol |

**Combined uncertainty:** Typically ±10-30 kJ/mol for ΔHd

### Oxygen Balance

Calculation is exact for supported elements.

---

## Methodology Limitations

### Gas-Phase Assumption

PHOENIX calculations primarily assume:
- Gas-phase reactants (Benson GA)
- Gas-phase products for decomposition

**Impact:**
- Condensed-phase ΔHf° may differ by 10-50 kJ/mol
- Some products (C, S, P₄O₁₀) are solid

### Decomposition Model

**Assumptions:**
1. Complete decomposition to most stable products
2. Thermodynamic equilibrium (not kinetic)
3. Standard conditions (298.15 K, 1 atm)

**Not Modeled:**
- Kinetic barriers
- Detonation vs deflagration
- Partial decomposition pathways
- Confinement effects
- Impact/friction sensitivity

### Temperature Dependence

| Range | Reliability |
|-------|-------------|
| 200-1500 K | Good |
| <200 K | Extrapolation warning |
| >6000 K | Extrapolation warning |

---

## CHETAH Criteria Limitations

### What CHETAH Predicts

- Thermodynamic potential for energy release
- Structural features associated with hazards
- Relative hazard ranking

### What CHETAH Does Not Predict

- Initiation sensitivity (impact, friction, spark)
- Detonation velocity or pressure
- Critical diameter
- Storage stability
- Compatibility with other materials
- Environmental degradation

### False Negatives

Some hazardous compounds may be classified LOW:
- Peroxides after aging
- Shock-sensitive primary explosives
- Compounds with unusual decomposition pathways

### False Positives

Some compounds classified HIGH may be:
- Stable under normal conditions
- Only hazardous under specific circumstances

---

## Data Coverage Gaps

### Benson GA

Missing group data for:
- Some heterocyclic systems
- Unusual ring sizes
- Exotic functional groups
- Some organo-phosphorus compounds

When groups are missing, `MissingGroupError` is raised.

### Reference Data

Limited experimental data for:
- New compounds
- Energetic materials (safety restrictions)
- Unstable intermediates

---

## Computational Limitations

### Large Molecules

| Atoms | Performance | Accuracy |
|-------|-------------|----------|
| <50 | Fast, accurate | Good |
| 50-100 | Fast, accurate | Good |
| 100-200 | Moderate | Reduced |
| >200 | Slow | Uncertain |

Molecules >100 atoms generate `UserWarning`.

### LP Solver

Linear Programming decomposition may fail for:
- Unusual elemental compositions
- Numerical instabilities

Falls back to hierarchy method on failure.

---

## Validation Requirements

### When to Validate Experimentally

1. **Always** for new or unfamiliar compounds
2. HIGH hazard classifications
3. Safety-critical applications
4. Compounds near classification boundaries
5. Unusual structural features

### Recommended Validation Methods

| Property | Validation Method |
|----------|-------------------|
| ΔHf° | Bomb calorimetry |
| Decomposition | DSC/DTA |
| Sensitivity | Impact/friction testing |
| Stability | Accelerated aging |

---

## Use Case Boundaries

### Appropriate Uses

- Initial screening of compound libraries
- Prioritization for experimental testing
- Educational purposes
- Preliminary hazard assessment

### Inappropriate Uses

- Sole basis for safety decisions
- Regulatory compliance without validation
- Quantitative risk assessment
- Design of energetic materials

---

## Disclaimer

!!! danger "Safety Disclaimer"
    PHOENIX results are theoretical estimates based on thermodynamic calculations. They do not account for kinetic factors, sensitivity, or real-world conditions.

    **Always:**

    1. Validate results experimentally
    2. Consult qualified safety professionals
    3. Follow applicable regulations
    4. Use appropriate protective measures
    5. Have emergency procedures in place

    The authors and contributors of PHOENIX assume no liability for decisions made based on these calculations.
