## ADDED Requirements

### Requirement: Thermodynamic Estimation
The library SHALL provide accurate thermodynamic properties (ΔHf°, S°, Cp(T)) using Benson group additivity methods.

#### Scenario: Calculate enthalpy of formation
- **WHEN** a valid SMILES string is provided
- **THEN** return ΔHf° in kcal/mol with <5% error margin

### Requirement: Maximum Heat of Decomposition
The library SHALL calculate the maximum heat of decomposition using linear programming optimization.

#### Scenario: Optimize decomposition products
- **WHEN** molecular formula is analyzed
- **THEN** return maximum ΔHd and corresponding products

### Requirement: Hazard Classification
The library SHALL classify compounds based on energy density thresholds and safety criteria.

#### Scenario: Classify compound hazard
- **WHEN** thermodynamic data is available
- **THEN** return hazard class (low, medium, high risk)

### Requirement: Batch Screening
The library SHALL support processing multiple molecules efficiently.

#### Scenario: Screen SMILES list
- **WHEN** list of SMILES provided
- **THEN** return DataFrame with all calculated properties

### Requirement: Yoshida Correlation
The library SHALL implement Yoshida correlation for explosive properties from DSC data.

#### Scenario: Calculate EP from DSC
- **WHEN** Q_dsc and T_onset provided
- **THEN** return estimated EP value