# Project Context

## Purpose

PHOENIX (Physicochemical Hazard Observation & Energetics Numerical Indexing eXpert) is a Python library for reactive chemical hazard evaluation. It serves as a modern, open-source successor to ASTM CHETAH, providing thermodynamic estimation, decomposition energetics, and hazard classification for process safety and chemical screening applications.

> **Warning**: This library is for screening purposes only. Results must not be used as the sole basis for safety decisions. Experimental validation is required before handling energetic materials. Consult qualified safety professionals and relevant regulations.

## License

MIT License - see [LICENSE](../LICENSE) file.

## Quick Start

### Installation

```bash
# Install from PyPI (once published)
pip install phoenix-chem

# Or using uv (recommended)
uv pip install phoenix-chem
```

### Development Setup

```bash
# Clone and setup development environment
git clone https://github.com/yourorg/phoenix-chem.git
cd phoenix-chem

# Create virtual environment and install with dev dependencies
uv venv --python 3.12
source .venv/bin/activate  # or `.venv\Scripts\activate` on Windows
uv pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Basic Usage

```python
from phoenix import Compound

# Analyze a compound from SMILES
compound = Compound.from_smiles("C1=CC=C(C=C1)[N+](=O)[O-]")  # nitrobenzene

# Get thermodynamic properties
print(f"ΔHf°: {compound.delta_hf_kJ_mol} kJ/mol")
print(f"Molecular Weight: {compound.molecular_weight} g/mol")

# Evaluate hazard
result = compound.evaluate_hazard()
print(f"Max ΔHd: {result.max_decomposition_kJ_mol} kJ/mol")
print(f"Max ΔHd: {result.max_decomposition_cal_g} cal/g")
print(f"Oxygen Balance: {result.oxygen_balance_percent}%")
print(f"Hazard Class: {result.hazard_class}")  # -> "Medium"

# Get full report with audit trail
report = compound.report()
print(report.to_json())
```

### Batch Screening

```python
from phoenix import screen

smiles_list = [
    "CC(=O)OOC(=O)C",      # peracetic acid
    "C1=CC=C(C=C1)[N+](=O)[O-]",  # nitrobenzene
    "CCO",                  # ethanol
]

results = screen(smiles_list)
print(results.to_dataframe())
```

### Running Tests

```bash
# Run all tests
pytest

# Run specific test tiers
pytest tests/unit/              # Unit tests only
pytest tests/integration/       # Integration tests
pytest tests/validation/ -v     # Validation suite (accuracy tests)

# Run tests for specific compound
pytest -k "test_tnt" -v

# Run with coverage
pytest --cov=phoenix --cov-report=html
```

## Goals

1. **Thermodynamic Property Estimation**: Calculate ΔHf°, S°, Cp(T) using Benson Group Additivity
2. **Maximum Heat of Decomposition**: Determine worst-case decomposition energy via analytical hierarchy or LP optimization
3. **Hazard Screening**: Classify compounds by energy release potential, oxygen balance, and functional group alerts
4. **Batch Processing**: Screen large SMILES datasets for reactive hazards
5. **Reproducibility**: Match or exceed CHETAH accuracy with transparent, documented algorithms

### Non-Goals (Out of Scope)

- Kinetic modeling or reaction rate prediction
- Detonation velocity or detonability prediction
- Chemical compatibility matrices
- Synthesis route suggestions
- Real-time process monitoring

## Tech Stack

### Runtime Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| Python | >= 3.12 | Runtime |
| rdkit | >= 2024.03 | Molecular parsing, SMILES canonicalization, atom counting |
| pgradd | >= 2.0 | Benson Group Additivity (VlachosGroup/PythonGroupAdditivity) |
| scipy | >= 1.11 | Linear programming for max ΔHd optimization |
| numpy | >= 1.26 | Numerical operations |
| pandas | >= 2.0 | Data handling, batch results |
| pydantic | >= 2.0 | Data validation, result models |
| chemicals | >= 1.1 | Thermodynamic data (70,000+ compounds) |

### Development Dependencies

| Package | Purpose |
|---------|---------|
| uv | Package and environment management |
| ruff | Linting and formatting (replaces black, isort, flake8) |
| ty | Type checking |
| pytest | Testing framework |
| pytest-cov | Coverage reporting |
| hypothesis | Property-based testing |
| pre-commit | Git hooks for code quality |

### Optional Dependencies (Future)

| Package | Purpose |
|---------|---------|
| janaf | Validation against NIST-JANAF tables |
| thermo (CalebBell) | Property lookup, validation against known compounds |
| Cantera | Advanced kinetics, equilibrium calculations |

## Project Conventions

### Code Style

- Follow PEP 8 with line length 100
- Use type hints for all public APIs (strict ty compliance)
- Use Google-style docstrings with explicit units in parameter descriptions
- Prefer `@dataclass(frozen=True)` for result objects (immutable)

**Tooling:**
```bash
# Format code
ruff format .

# Lint and auto-fix
ruff check . --fix

# Type check
ty check

# Run all checks (via pre-commit)
pre-commit run --all-files
```

### Naming Conventions

| Type | Convention | Example |
|------|------------|---------|
| Classes | `PascalCase` | `Compound`, `DecompositionResult` |
| Functions/methods | `snake_case` | `calculate_oxygen_balance` |
| Constants | `UPPER_SNAKE_CASE` | `DELTA_HF_PRODUCTS` |
| Private | `_leading_underscore` | `_validate_smiles` |
| Units in names | Suffix with unit | `delta_hf_kJ_mol`, `entropy_J_mol_K` |

### Architecture Patterns

#### Package Structure

```
phoenix/
├── __init__.py              # Public API exports
├── core/
│   ├── __init__.py
│   ├── compound.py          # Compound class, molecular parsing
│   ├── composition.py       # Elemental composition handling
│   ├── validation.py        # Input validation, SMILES sanitization
│   └── units.py             # Unit handling utilities
├── thermo/
│   ├── __init__.py          # Thermo module exports
│   ├── benson.py            # Benson GA via pgradd
│   ├── properties.py        # ΔHf°, S°, Cp(T) estimation
│   └── data.py              # Thermodynamic data via `chemicals` library
├── hazard/
│   ├── __init__.py
│   ├── decomposition.py     # Max ΔHd (analytical + LP)
│   ├── oxygen_balance.py    # OB% calculations
│   ├── gas_generation.py    # Gas volume estimation
│   ├── classification.py    # Hazard class assignment
│   └── functional_groups.py # SMARTS-based hazard group detection
├── correlations/
│   ├── __init__.py
│   ├── yoshida.py           # Yoshida EP/SS from DSC
│   └── other.py             # Future correlations
├── batch/
│   ├── __init__.py
│   └── screening.py         # Batch processing, parallelization
└── exceptions.py            # Custom exception classes
```

#### Design Principles

1. **Immutable Results**: All result objects are frozen dataclasses
2. **Lazy Evaluation**: Heavy imports (rdkit, scipy) deferred until needed
3. **Explicit Units**: All numeric values paired with unit documentation
4. **Fail-Fast Validation**: Input validation at API boundaries
5. **Audit Trail**: Results include method, assumptions, warnings, data versions

#### Data Flow

```
SMILES Input
    ↓
┌─────────────────┐
│ Compound.from_  │ ← Validation, canonicalization (rdkit)
│ smiles()        │
└────────┬────────┘
         ↓
┌─────────────────┐
│ Elemental       │ ← Atom counting, formula extraction
│ Composition     │
└────────┬────────┘
         ↓
┌─────────────────┐
│ Thermo Props    │ ← Benson GA (pgradd) or chemicals lookup
│ ΔHf°, S°, Cp(T) │
└────────┬────────┘
         ↓
┌─────────────────┐
│ Hazard Eval     │ ← Max ΔHd, OB%, functional groups
│                 │
└────────┬────────┘
         ↓
┌─────────────────┐
│ HazardResult    │ ← Classification + audit trail
└─────────────────┘
```

### Error Handling

PHOENIX uses explicit exceptions for error conditions:

| Exception | When Raised |
|-----------|-------------|
| `InvalidSmilesError` | SMILES string cannot be parsed by rdkit |
| `UnsupportedElementError` | Compound contains metals or unsupported elements |
| `UnsupportedStructureError` | Radicals, ions, or charged species |
| `MissingGroupError` | Benson GA lacks group contribution data |
| `DecompositionInfeasibleError` | LP solver finds no feasible product distribution |

```python
from phoenix import Compound
from phoenix.exceptions import InvalidSmilesError, UnsupportedElementError

try:
    compound = Compound.from_smiles("invalid_smiles")
except InvalidSmilesError as e:
    print(f"Parse error: {e}")
except UnsupportedElementError as e:
    print(f"Element not supported: {e}")
```

### Testing Strategy

- **Pytest** for all tests
- **Validation Set**: 50+ compounds with known CHETAH/NIST values
- **Tiers**:
  1. Unit tests: Individual functions, edge cases
  2. Integration tests: Full compound analysis pipelines
  3. Validation tests: Accuracy against reference data (TNT, RDX, PETN, etc.)
  4. Property-based tests: SMILES parsing robustness (hypothesis)
- **Tolerance**: ΔHf° within ±10 kJ/mol, OB% exact, hazard class matches CHETAH

### Git Workflow

- Feature branching from `main`
- Spec-driven development via OpenSpec
- Conventional commits: `feat:`, `fix:`, `docs:`, `refactor:`, `test:`
- PR requires:
  - All tests passing
  - `ruff check` clean
  - `ty check` clean
  - Pre-commit hooks passed

## Domain Context

### Key Concepts

| Term | Definition |
|------|------------|
| **SMILES** | Simplified Molecular-Input Line-Entry System. Text notation for molecular structures (e.g., `CCO` = ethanol) |
| **SMARTS** | SMILES arbitrary target specification. Pattern language for substructure matching |
| **Benson Group Additivity (GA)** | Method to estimate thermodynamic properties by summing contributions of functional groups. Accuracy: typically ±4 kJ/mol for ΔHf° |
| **Maximum Heat of Decomposition (ΔHd)** | Thermodynamic worst-case energy release assuming optimal product distribution. Not kinetically constrained |
| **Oxygen Balance (OB%)** | Measure of oxygen excess/deficiency for complete oxidation. OB≈0 maximizes energy density; negative = fuel-rich |
| **DSC** | Differential Scanning Calorimetry. Measures heat flow vs temperature to detect thermal events |
| **Yoshida Correlation** | Empirical relationship predicting explosive propagation (EP) and shock sensitivity (SS) from DSC data |
| **EP** | Explosive Propagation. Likelihood of self-sustained decomposition |
| **SS** | Shock Sensitivity. Susceptibility to initiation by mechanical impact |
| **CHETAH Criteria** | ASTM E659 hazard thresholds for chemical instability screening |

### CHETAH Hazard Criteria

| Criterion | Threshold | Interpretation |
|-----------|-----------|----------------|
| Criterion 1 | ΔHd < -300 cal/g | High instability potential |
| Criterion 2 | ΔHd < -100 cal/g | Medium instability potential |
| Criterion 3 | OB > -200% and < +100% | Oxidizer/fuel balance concern |
| Criterion 4 | Functional group alerts | Known reactive moieties present |

**Unit Conversions:**
- -300 cal/g = -1.255 kJ/g = -1255 J/g
- To convert kJ/mol to cal/g: `cal_g = (kJ_mol × 1000) / (4.184 × MW_g_mol)`
- PHOENIX uses **kJ/mol internally** and provides both units in results

### Thermodynamic Hierarchy (Max ΔHd)

Priority order for analytical decomposition calculation (maximizes energy release):

| Priority | Reaction | Rationale |
|----------|----------|-----------|
| 1 | F + H → HF | Strongest H-X bond (570 kJ/mol) |
| 2 | N → ½N₂ | Non-competitive, always forms N₂ |
| 3 | P + O → ¼P₄O₁₀ | Highly exothermic (-2984 kJ/mol) |
| 4 | H + O → H₂O | Primary hydrogen oxidation |
| 5 | C + O → CO₂ | Full carbon oxidation (if O available) |
| 6 | S + O → SO₂ | Sulfur oxidation |
| 7 | C + O → CO | Partial oxidation (O-limited) |
| 8 | Cl + H → HCl | After H₂O formation |
| 9 | Br + H → HBr | After HCl formation |
| 10 | C → C(s) | Graphite (very O-deficient) |
| 11 | H → ½H₂ | Excess hydrogen |
| 12 | S → S(s) | Excess sulfur |
| 13 | Excess halogens → X₂ | Elemental form |

**When is LP used vs. analytical hierarchy?**
- **Analytical hierarchy**: Default method, fast, deterministic, CHETAH-compatible
- **LP optimization**: Available for validation, research, or when analytical fails
- Results from both methods should agree for standard CHNO compounds

### Oxygen Balance Formula

For compounds with formula C_a H_b N_c O_d S_e P_f Cl_g Br_h:

```
OB% = -1600/MW × (2a + b/2 + 2e + 2.5f - g/2 - h/2 - d)
```

Where:
- `a, b, c, d, e, f, g, h` = atom counts for C, H, N, O, S, P, Cl, Br
- `MW` = molecular weight (g/mol)
- Assumes complete oxidation to CO₂, H₂O, SO₂, P₄O₁₀
- Halogens form HX (consume H, effectively "release" O equivalent)

### Data Sources

#### Python Thermochemical Libraries

| Library | PyPI | Content | Best For |
|---------|------|---------|----------|
| [chemicals](https://github.com/CalebBell/chemicals) | `pip install chemicals` | 70,000+ compounds, ΔHf°(g/l/s), DIPPR data | Formation enthalpies via `Hfg()`, `Hfl()`, `Hfs()` |
| [thermo](https://github.com/CalebBell/thermo) | `pip install thermo` | Full thermodynamic modeling, phase equilibria | Temperature-dependent Cp, phase transitions |
| [janaf](https://github.com/n-takumasa/py-janaf) | `pip install janaf` | NIST-JANAF tables, ~1800 species | High-accuracy reference data, Cp(T) polynomials |
| [CoolProp](https://github.com/CoolProp/CoolProp) | `pip install CoolProp` | 110+ pure fluids, transport properties | Refrigerants, industrial fluids |
| [Cantera](https://cantera.org/) | `pip install cantera` | Chemical kinetics, thermodynamics | Combustion, reaction equilibrium |

#### Primary Data Sources

| Source | Content | Access |
|--------|---------|--------|
| [NIST-JANAF Tables](https://janaf.nist.gov/) | ΔHf°, S°, Cp(T) for ~1800 species | `janaf` package |
| [NIST Chemistry WebBook](https://webbook.nist.gov/) | ΔHf° for 40,000+ compounds | `chemicals` library |
| [DIPPR 801](https://www.aiche.org/dippr) | Industrial property database | Via `chemicals` library |
| CHETAH Documentation | Hazard criteria, algorithms | ASTM E659 standard |

#### Data Access Strategy

1. **Runtime**: Use `chemicals` library directly via `phoenix.thermo.data` module
2. **Fallback**: Built-in NIST-JANAF values for core decomposition products (offline support)
3. **Caching**: `@lru_cache` for repeated lookups within a session
4. **Versioning**: Data source version recorded in result audit trail

## Important Constraints

### Accuracy Requirements

| Property | Target | Basis |
|----------|--------|-------|
| ΔHf° (gas) | ±10 kJ/mol | Benson GA typical error (±4 kJ/mol ideal, ±10 kJ/mol conservative) |
| ΔHd | ±15% | Acceptable for screening applications |
| OB% | Exact | Stoichiometric calculation (no estimation) |
| Hazard Class | Match CHETAH | Validated against reference compounds |

### Scope Limitations (MVP)

| Limitation | Details |
|------------|---------|
| **Supported Elements** | C, H, N, O, S, P, F, Cl, Br (no metals) |
| **Molecular Size** | Warning for >100 atoms (GA accuracy degrades) |
| **Phase State** | Gas-phase reference state only (see note below) |
| **No Radicals/Ions** | Reject charged species (formal charge ≠ 0) and open-shell radicals |
| **Single Components** | No mixture hazard analysis |

**Conservative Gas-Phase Assumption**: MVP uses gas-phase ΔHf° for reactants. Since ΔHf°(gas) > ΔHf°(solid/liquid) for most organic compounds, this yields larger (more conservative) decomposition energy estimates—appropriate for safety screening.

### Safety Considerations

1. **Disclaimer**: All outputs include screening-only caveat
2. **No Synthesis Routes**: Library does not provide preparation methods
3. **Audit Trail**: Results traceable to input, method, data version, and assumptions
4. **Conservative Estimates**: Err toward higher hazard classification when uncertain
5. **No Optimization for Energetics**: Library identifies hazards, not optimal explosive formulations

## External References

### Primary Standards

- **ASTM E659**: CHETAH hazard criteria and methodology
- **ASTM E1231**: Yoshida correlation for DSC data interpretation
- **ASTM E968**: DSC calibration standards

### Thermodynamic Data

- NIST-JANAF Thermochemical Tables, 4th Ed. (Chase, 1998)
- NIST Chemistry WebBook, SRD 69
- ATcT (Active Thermochemical Tables) for high-accuracy values

### Literature

- Meyer, Köhler, Homburg: "Explosives" 6th Ed. (Wiley-VCH, 2007) - ISBN 978-3-527-31656-4
- Cooper: "Explosives Engineering" (Wiley, 1996)
- Akhavan: "The Chemistry of Explosives" 3rd Ed. (RSC, 2011)
- pgradd paper: [doi.org/10.1016/j.cpc.2021.108221](https://doi.org/10.1016/j.cpc.2021.108221)

## Validation Compounds

### Primary Test Set

| Compound | Formula | ΔHf° (kJ/mol) | OB% | Hazard | Source |
|----------|---------|---------------|-----|--------|--------|
| TNT | C₇H₅N₃O₆ | -67 | -74% | High | NIST |
| RDX | C₃H₆N₆O₆ | +70 | -22% | High | NIST |
| PETN | C₅H₈N₄O₁₂ | -538 | -10% | High | NIST |
| Nitroglycerin | C₃H₅N₃O₉ | -370 | +3.5% | High | NIST |
| Ammonium nitrate | H₄N₂O₃ | -365 | +20% | High | NIST |
| Peracetic acid | C₂H₄O₃ | -390 | -105% | Medium | NIST |
| Benzoyl peroxide | C₁₄H₁₀O₄ | -369 | -155% | Medium | NIST |
| Ethylene oxide | C₂H₄O | -52.6 | -182% | Medium | NIST |
| Nitrobenzene | C₆H₅NO₂ | +12 | -163% | Medium | NIST |
| Acetone | C₃H₆O | -248 | -207% | Low | NIST |
| Water | H₂O | -242 | +89% | None | Reference |
| Urea | CH₄N₂O | -333 | -27% | Low | NIST |
| Carbon dioxide | CO₂ | -394 | +73% | None | Reference |

### Validation Criteria

| Property | Pass Criteria |
|----------|---------------|
| ΔHf° | Within ±10 kJ/mol of reference |
| OB% | Exact match (calculation is deterministic) |
| Hazard Class | Match CHETAH classification |
| Max ΔHd | Within ±15% of CHETAH value |

## Roadmap

### MVP (v0.1)

- [x] Core compound parsing and validation
- [x] Thermodynamic data access via `chemicals`
- [ ] Benson GA integration via `pgradd`
- [ ] Analytical max ΔHd calculation
- [ ] Oxygen balance calculation
- [ ] Hazard classification (Criteria 1-4)
- [ ] Functional group alerts
- [ ] Batch screening
- [ ] Validation suite (50+ compounds)

### v0.2 (Post-MVP)

- [ ] LP-based max ΔHd optimization
- [ ] Yoshida correlation (EP/SS from DSC)
- [ ] JSON/CSV export formats
- [ ] CLI interface

### v0.3 (Future)

- [ ] Phase corrections (liquid/solid ΔHf via Trouton's rule)
- [ ] Uncertainty quantification
- [ ] Extended element support (I, metals)

### v1.0 (Long-term)

- [ ] Mixture hazard compatibility screening
- [ ] Integration with process safety tools
- [ ] Regulatory report generation (UN GHS alignment)

## Contributing

### Getting Started

1. Fork the repository
2. Create a feature branch: `git checkout -b feat/my-feature`
3. Make changes following the code style guidelines
4. Run tests and checks: `pre-commit run --all-files && pytest`
5. Commit with conventional format: `git commit -m "feat: add new feature"`
6. Push and create a Pull Request

### Adding a Validation Compound

1. Find reference ΔHf° value (prefer NIST sources)
2. Add entry to `tests/validation/data/compounds.json`
3. Add test case to `tests/validation/test_accuracy.py`
4. Run validation suite: `pytest tests/validation/ -v`

### Adding a Functional Group Alert

1. Create SMARTS pattern for the reactive group
2. Add to `phoenix/hazard/functional_groups.py`
3. Add test case with known positive and negative examples
4. Document the hazard rationale

### Reporting Issues

- Use GitHub Issues for bug reports and feature requests
- Include SMILES string, expected vs actual output, and PHOENIX version
- For security concerns, email maintainers directly
