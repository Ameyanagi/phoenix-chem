## ADDED Requirements

### Requirement: Benson GA Integration
The system SHALL estimate $\Delta H_f^\circ$, $S^\circ$, and $C_p(T)$ using Benson Group Additivity via `pgradd`.

#### Scenario: Gas-phase Enthalpy
- **WHEN** `molecule.delta_hf()` is called for a molecule with known groups (e.g., Methane)
- **THEN** it returns the estimated standard enthalpy of formation.

### Requirement: Missing Group Handling
The system SHALL warn the user if a molecule contains functional groups not supported by the Benson GA database.

#### Scenario: Unsupported Group
- **WHEN** a molecule with an exotic group is processed
- **THEN** the system issues a warning and indicates which atoms were not covered.
