# API Reference

Complete API documentation for PHOENIX.

## Module Structure

```
phoenix/
├── __init__.py          # Main exports
├── core/
│   ├── compound.py      # Compound class
│   └── reaction.py      # Reaction, ReactionSpecies, Auto
├── thermo/
│   ├── models.py        # ThermoProperty, ThermoState
│   ├── benson.py        # Benson GA estimation
│   ├── data.py          # NIST reference data
│   └── references.py    # Reference citations
├── hazard/
│   ├── classification.py # HazardResult, evaluate_hazard
│   ├── decomposition.py  # DecompositionResult, calculate_max_decomposition
│   ├── functional_groups.py # FunctionalGroupAlert
│   └── oxygen_balance.py # OB% calculation
├── batch/
│   ├── screening.py     # screen(), BatchResult
│   └── __init__.py
└── exceptions.py        # All exception classes
```

## Import Patterns

### Standard Imports

```python
# Core classes
from phoenix import Compound, Reaction, Auto, ReactionSpecies

# Batch processing
from phoenix import screen, BatchResult

# Results
from phoenix import HazardResult, ThermoProperty

# Exceptions
from phoenix import (
    PhoenixError,
    InvalidSmilesError,
    UnsupportedElementError,
    UnsupportedStructureError,
    MissingGroupError,
    DecompositionError,
    BalanceError,
    OverconstrainedError,
    UnderconstrainedError,
)
```

### Advanced Imports

```python
# Detailed functional group detection
from phoenix.hazard.functional_groups import (
    detect_functional_groups,
    FunctionalGroupAlert,
)

# Direct decomposition access
from phoenix.hazard.decomposition import (
    calculate_max_decomposition,
    DecompositionResult,
    DecompositionComparison,
)

# Thermodynamic data
from phoenix.thermo.data import get_formation_enthalpy

# Benson GA estimation
from phoenix.thermo.benson import (
    estimate_delta_hf,
    estimate_entropy,
    estimate_heat_capacity,
)
```

## Type Hints

PHOENIX uses Python type hints throughout:

```python
from typing import Literal

# Hazard class type
HazardClass = Literal["HIGH", "MEDIUM", "LOW"]

# Coefficient specification
CoeffSpec = float | int | _AutoType | None
```

## API Sections

<div class="grid cards" markdown>

-   **Compound**

    Core class for molecular representation

    [:octicons-arrow-right-24: Compound API](compound.md)

-   **Reaction**

    Reaction balancing and thermodynamics

    [:octicons-arrow-right-24: Reaction API](reaction.md)

-   **Hazard**

    Hazard evaluation and decomposition

    [:octicons-arrow-right-24: Hazard API](hazard.md)

-   **Thermo**

    Thermodynamic properties

    [:octicons-arrow-right-24: Thermo API](thermo.md)

-   **Batch**

    Batch screening

    [:octicons-arrow-right-24: Batch API](batch.md)

-   **Exceptions**

    Error handling

    [:octicons-arrow-right-24: Exceptions](exceptions.md)

</div>
