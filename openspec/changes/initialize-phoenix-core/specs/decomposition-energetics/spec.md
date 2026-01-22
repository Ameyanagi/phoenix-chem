## ADDED Requirements

### Requirement: Maximum Heat of Decomposition
The system SHALL calculate the maximum heat of decomposition (Max $\Delta H_d$) by solving a Linear Programming problem for the most stable products.

#### Scenario: Decomposition of Nitrobenzene
- **WHEN** `molecule.max_decomposition()` is called for "C1=CC=C(C=C1)[N+](=O)[O-]"
- **THEN** it returns a value reflecting the energy released upon decomposition to the most stable products.

### Requirement: Gas Generation Volume
The system SHALL estimate the volume of gas generated per unit mass of substance at STP.

#### Scenario: Gas Volume Calculation
- **WHEN** `molecule.gas_generation()` is called
- **THEN** it returns the gas volume in L/g.
