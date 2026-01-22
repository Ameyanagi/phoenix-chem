# Installation

## Requirements

- Python 3.12 or higher
- pip or uv package manager

## Install from PyPI

=== "pip"

    ```bash
    pip install phoenix-chem
    ```

=== "uv"

    ```bash
    uv pip install phoenix-chem
    ```

## Development Installation

For contributing or development:

```bash
# Clone the repository
git clone https://github.com/Ameyanagi/phoenix-chem.git
cd phoenix-chem

# Create virtual environment
uv venv --python 3.12
source .venv/bin/activate  # Linux/macOS
# or: .venv\Scripts\activate  # Windows

# Install with dev dependencies
uv pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

## Dependencies

PHOENIX automatically installs these dependencies:

| Package | Purpose |
|---------|---------|
| `rdkit` | Molecular parsing, SMILES handling |
| `pgradd` | Benson Group Additivity estimation |
| `chemicals` | Thermodynamic reference data |
| `scipy` | Linear programming for decomposition |
| `numpy` | Numerical operations |
| `pandas` | Data handling for batch results |
| `pydantic` | Data validation |

## Verify Installation

```python
import phoenix
print(phoenix.__version__)

# Quick test
from phoenix import Compound
ethanol = Compound.from_smiles("CCO")
print(f"Ethanol: {ethanol.formula}, MW = {ethanol.molecular_weight:.2f}")
```

## Troubleshooting

### RDKit Installation Issues

RDKit requires specific system libraries. If installation fails:

=== "Linux (Ubuntu/Debian)"

    ```bash
    sudo apt-get install libboost-all-dev
    pip install rdkit
    ```

=== "macOS"

    ```bash
    brew install boost
    pip install rdkit
    ```

=== "Windows"

    Use conda for easiest installation:
    ```bash
    conda install -c conda-forge rdkit
    pip install phoenix-chem
    ```

### pgradd Issues

The pgradd library (Benson GA) may have issues on some systems:

```bash
# Try installing with specific version
pip install pgradd>=2.9.13
```

If pgradd fails, PHOENIX falls back to the `chemicals` library for thermodynamic data.

## Optional Dependencies

For documentation development:

```bash
pip install phoenix-chem[docs]
```

For all optional features:

```bash
pip install phoenix-chem[all]
```
