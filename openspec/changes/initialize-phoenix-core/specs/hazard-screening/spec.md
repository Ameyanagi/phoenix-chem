## ADDED Requirements

### Requirement: Oxygen Balance
The system SHALL calculate the Oxygen Balance (OB%) of a molecule.

#### Scenario: OB of TNT
- **WHEN** `molecule.oxygen_balance()` is called for TNT
- **THEN** it returns -74.0%.

### Requirement: Yoshida Correlation
The system SHALL provide a function to calculate Yoshida correlations from DSC data.

#### Scenario: Yoshida Screen
- **WHEN** `phoenix.yoshida(Q_dsc=1000, T_onset=200)` is called
- **THEN** it returns hazard rankings (EP, SS).

### Requirement: Hazard Classification
The system SHALL classify molecules into hazard categories based on energy density and OB levels.

#### Scenario: High Energy Classification
- **WHEN** a molecule has Max $\Delta H_d$ > 3.0 kJ/g
- **THEN** it is classified as "High Hazard".
