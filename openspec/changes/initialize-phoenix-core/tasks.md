## 1. Project Scaffolding
- [ ] 1.1 Update `pyproject.toml` with dependencies and project metadata.
- [ ] 1.2 Create the package structure: `phoenix/`, `phoenix/core/`, `phoenix/thermo/`, `phoenix/energetics/`, `phoenix/hazard/`, `phoenix/io/`.

## 2. Core Molecule Management
- [ ] 2.1 Implement `phoenix.core.molecule.Molecule` class with RDKit integration.
- [ ] 2.2 Add SMILES validation and normalization.

## 3. Thermodynamic Estimation
- [ ] 3.1 Integrate `pgradd` for Benson Group Additivity.
- [ ] 3.2 Implement `delta_hf`, `entropy`, and `heat_capacity` calculations.
- [ ] 3.3 Add phase correction estimations (Gas to Liquid/Solid).

## 4. Decomposition Energetics
- [ ] 4.1 Implement Max $\Delta H_d$ optimization using `scipy.optimize.linprog`.
- [ ] 4.2 Implement gas generation volume calculation.

## 5. Hazard Screening
- [ ] 5.1 Implement Oxygen Balance (OB%) calculation.
- [ ] 5.2 Implement Yoshida correlation for DSC data.
- [ ] 5.3 Implement hazard classification logic based on energy density.

## 6. Batch Processing and Reporting
- [ ] 6.1 Implement `phoenix.screen` for batch SMILES processing.
- [ ] 6.2 Implement reporting functionality (JSON/DataFrame).

## 7. Verification and Documentation
- [ ] 7.1 Write unit tests for all core components.
- [ ] 7.2 Validate against CHETAH reference cases.
- [ ] 7.3 Create initial usage documentation.
