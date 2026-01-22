## Context
The goal is to create a Pythonic, modular, and robust library for chemical hazard evaluation. The design must be extensible to support future features like mixture analysis and structural alerts.

## Goals
- Provide a clean `Molecule` centric API.
- Ensure high performance for batch screening.
- Maintain rigorous thermodynamic accuracy where Benson GA is applicable.

## Decisions

### 1. Molecule Object as Central Hub
- **Decision**: The `Molecule` class will hold the RDKit molecule object and cache computed properties.
- **Rationale**: This allows for lazy evaluation of expensive thermodynamic or energetics calculations and provides a single point of entry for the user.

### 2. Linear Programming for Max $\Delta H_d$
- **Decision**: Use `scipy.optimize.linprog` with a predefined set of stable decomposition products (e.g., $N_2, H_2O, CO_2, CO, C(s), HCl, HF$).
- **Rationale**: This is the standard approach for "worst-case" decomposition energy estimation.

### 3. pgradd for Benson GA
- **Decision**: Wrap `pgradd` to handle functional group decomposition and thermodynamic summation.
- **Rationale**: `pgradd` is an established library for Benson GA, avoiding the need to reinvent the group parsing logic.

### 4. Modular Internal Structure
- **Decision**: Separate concerns into `thermo`, `energetics`, and `hazard` modules.
- **Rationale**: Facilitates testing and future expansion.

## Risks / Trade-offs
- **Benson GA Coverage**: Many reactive functional groups might not be in standard Benson databases. We need a fallback or warning system.
- **LP Product Set**: The choice of products in the LP solver significantly affects the Max $\Delta H_d$. We should follow CHETAH's product hierarchy.

## Open Questions
- Should we support 3D conformation for any calculations? (Unlikely for MVP).
- How to handle charged species or radicals?
