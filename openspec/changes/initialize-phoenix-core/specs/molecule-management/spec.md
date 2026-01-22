## ADDED Requirements

### Requirement: SMILES Parsing
The system SHALL accept SMILES strings and convert them into internal `Molecule` objects using RDKit.

#### Scenario: Valid SMILES
- **WHEN** a valid SMILES string (e.g., "CC(=O)O") is provided
- **THEN** a `Molecule` object is successfully created.

### Requirement: RDKit Integration
The `Molecule` class SHALL allow access to the underlying RDKit `Mol` object.

#### Scenario: Accessing RDKit Mol
- **WHEN** a `Molecule` object is initialized
- **THEN** the `.rdmol` property returns a `rdkit.Chem.Mol` instance.

### Requirement: Basic Molecular Properties
The system SHALL provide molecular weight and formula for any initialized `Molecule`.

#### Scenario: Retrieving MW and Formula
- **WHEN** a `Molecule` (e.g., "C") is initialized
- **THEN** MW returns ~16.04 and formula returns "CH4".
