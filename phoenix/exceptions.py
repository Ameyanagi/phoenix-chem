"""Custom exception hierarchy for PHOENIX."""

from __future__ import annotations


class PhoenixError(Exception):
    """Base exception for all PHOENIX errors."""

    pass


class InvalidSmilesError(PhoenixError):
    """Raised when a SMILES string cannot be parsed."""

    def __init__(self, smiles: str, message: str | None = None):
        self.smiles = smiles
        if message is None:
            message = f"Invalid SMILES string: '{smiles}'"
        super().__init__(message)


class UnsupportedElementError(PhoenixError):
    """Raised when a compound contains unsupported elements."""

    def __init__(self, elements: list[str], message: str | None = None):
        self.elements = elements
        if message is None:
            message = f"Unsupported elements: {', '.join(elements)}"
        super().__init__(message)


class UnsupportedStructureError(PhoenixError):
    """Raised for radicals, ions, or other unsupported molecular structures."""

    def __init__(self, reason: str, smiles: str | None = None):
        self.reason = reason
        self.smiles = smiles
        message = f"Unsupported structure: {reason}"
        if smiles:
            message += f" (SMILES: {smiles})"
        super().__init__(message)


class MissingGroupError(PhoenixError):
    """Raised when Benson GA lacks group contribution data for a molecule."""

    def __init__(self, groups: list[str] | None = None, message: str | None = None):
        self.groups = groups or []
        if message is None:
            if groups:
                message = f"Missing Benson GA group data for: {', '.join(groups)}"
            else:
                message = "Missing Benson GA group data"
        super().__init__(message)


class DecompositionError(PhoenixError):
    """Raised when decomposition calculation fails."""

    def __init__(self, reason: str, formula: str | None = None):
        self.reason = reason
        self.formula = formula
        message = f"Decomposition calculation failed: {reason}"
        if formula:
            message += f" (formula: {formula})"
        super().__init__(message)


# =============================================================================
# Reaction Balancing Exceptions
# =============================================================================


class BalanceError(PhoenixError):
    """Base exception for reaction balancing errors."""

    pass


class OverconstrainedError(BalanceError):
    """
    Raised when reaction constraints are inconsistent.

    This occurs when user-specified coefficients violate atom conservation,
    making it impossible to balance the reaction.

    Attributes
    ----------
    imbalances : dict[str, float]
        Mapping of element to imbalance amount (positive = excess product)
    """

    def __init__(
        self,
        message: str = "Reaction is over-constrained: atom conservation cannot be satisfied",
        imbalances: dict[str, float] | None = None,
    ):
        self.imbalances = imbalances or {}
        if imbalances:
            details = ", ".join(f"{elem}: {val:+.4f}" for elem, val in imbalances.items())
            message = f"{message}. Imbalances: {details}"
        super().__init__(message)


class UnderconstrainedError(BalanceError):
    """
    Raised when multiple valid balanced solutions exist.

    This occurs when there are too few constraints (known coefficients)
    relative to the number of unknowns.

    Attributes
    ----------
    degrees_of_freedom : int
        Number of additional coefficients needed for unique solution
    suggestion : str | None
        Hint for how to resolve the issue
    """

    def __init__(
        self,
        message: str = "Reaction is under-constrained: multiple solutions exist",
        degrees_of_freedom: int = 0,
        suggestion: str | None = None,
    ):
        self.degrees_of_freedom = degrees_of_freedom
        self.suggestion = suggestion
        if degrees_of_freedom > 0:
            message = f"{message} ({degrees_of_freedom} degree(s) of freedom)"
        if suggestion:
            message = f"{message}. {suggestion}"
        super().__init__(message)
