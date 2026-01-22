# Change: Initialize PHOENIX Core Architecture and MVP Features

## Why
To establish the foundational architecture and feature set for PHOENIX as a modern successor to ASTM CHETAH, ensuring a robust, extensible, and user-friendly API for reactive chemical hazard evaluation.

## What Changes
- Define the core `Molecule` class API.
- Implement thermodynamic estimation using Benson Group Additivity.
- Implement maximum heat of decomposition ($\Delta H_d$) optimization via Linear Programming.
- Implement gas generation and oxygen balance calculations.
- Implement hazard classification and the Yoshida correlation.
- Establish a modular package structure.

## Impact
- Affected specs: `molecule-management`, `thermodynamic-estimation`, `decomposition-energetics`, `hazard-screening`.
- Affected code: New package structure under `phoenix/`.
