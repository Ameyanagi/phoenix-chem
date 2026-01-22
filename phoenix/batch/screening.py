"""
Batch screening for multiple compounds.

Provides efficient processing of SMILES lists with DataFrame output,
partial failure handling, and progress reporting.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import pandas as pd

from phoenix.core.compound import Compound
from phoenix.exceptions import PhoenixError


@dataclass
class BatchResult:
    """
    Result of batch screening operation.

    Attributes
    ----------
    dataframe : pd.DataFrame
        DataFrame with all results
    successful : int
        Number of successfully processed compounds
    failed : int
        Number of failed compounds
    """

    dataframe: pd.DataFrame
    successful: int
    failed: int

    def to_csv(self, path: str, **kwargs: Any) -> None:
        """Export results to CSV file."""
        self.dataframe.to_csv(path, index=False, **kwargs)

    def to_json(self, **kwargs: Any) -> str:
        """Export results to JSON string."""
        return self.dataframe.to_json(orient="records", **kwargs)

    def summary(self) -> dict[str, Any]:
        """
        Get summary statistics for the batch.

        Returns
        -------
        dict
            Summary including counts, hazard class distribution, and ΔHd range
        """
        df = self.dataframe

        # Handle empty DataFrame
        if len(df) == 0:
            return {
                "total_compounds": 0,
                "successful": 0,
                "failed": 0,
                "hazard_class_counts": {},
                "delta_hd_cal_g_min": None,
                "delta_hd_cal_g_max": None,
            }

        # Hazard class distribution (only for successful rows)
        successful_df = df[df["error"].isna()]
        hazard_counts = (
            successful_df["hazard_class"].value_counts().to_dict() if len(successful_df) > 0 else {}
        )

        # ΔHd range
        delta_hd_values = successful_df["delta_hd_cal_g"].dropna()
        delta_hd_min = float(delta_hd_values.min()) if len(delta_hd_values) > 0 else None
        delta_hd_max = float(delta_hd_values.max()) if len(delta_hd_values) > 0 else None

        return {
            "total_compounds": len(df),
            "successful": self.successful,
            "failed": self.failed,
            "hazard_class_counts": hazard_counts,
            "delta_hd_cal_g_min": delta_hd_min,
            "delta_hd_cal_g_max": delta_hd_max,
        }


def screen(
    smiles_list: list[str],
    progress_callback: Callable[[int, int], None] | None = None,
) -> BatchResult:
    """
    Screen multiple compounds for reactive hazards.

    Parameters
    ----------
    smiles_list : list[str]
        List of SMILES strings to process
    progress_callback : Callable[[int, int], None], optional
        Callback function called with (current_index, total_count)
        for progress reporting

    Returns
    -------
    BatchResult
        Results with DataFrame and statistics

    Examples
    --------
    >>> results = screen(["CCO", "c1ccccc1[N+](=O)[O-]"])
    >>> print(results.dataframe.columns.tolist())
    ['smiles', 'canonical_smiles', 'formula', 'mw', ...]
    >>> print(results.summary())
    {'total_compounds': 2, 'successful': 2, ...}
    """
    results: list[dict[str, Any]] = []
    total = len(smiles_list)
    successful = 0
    failed = 0

    for i, smiles in enumerate(smiles_list):
        if progress_callback is not None:
            progress_callback(i, total)

        row = _process_single(smiles)
        results.append(row)

        if row.get("error") is None:
            successful += 1
        else:
            failed += 1

    # Final callback
    if progress_callback is not None:
        progress_callback(total, total)

    # Create DataFrame
    df = pd.DataFrame(results)

    # Ensure column order
    column_order = [
        "smiles",
        "canonical_smiles",
        "formula",
        "mw",
        "delta_hf_kJ_mol",
        "delta_hd_kJ_mol",
        "delta_hd_cal_g",
        "ob_percent",
        "hazard_class",
        "triggered_criteria",
        "alerts",
        "gas_volume_L_g",
        "error",
        "error_message",
    ]
    # Only include columns that exist
    df = df[[c for c in column_order if c in df.columns]]

    return BatchResult(dataframe=df, successful=successful, failed=failed)


def _process_single(smiles: str) -> dict[str, Any]:
    """Process a single SMILES string and return result dict."""
    row: dict[str, Any] = {
        "smiles": smiles,
        "canonical_smiles": None,
        "formula": None,
        "mw": None,
        "delta_hf_kJ_mol": None,
        "delta_hd_kJ_mol": None,
        "delta_hd_cal_g": None,
        "ob_percent": None,
        "hazard_class": None,
        "triggered_criteria": None,
        "alerts": None,
        "gas_volume_L_g": None,
        "error": None,
        "error_message": None,
    }

    try:
        # Parse compound
        compound = Compound.from_smiles(smiles)
        row["canonical_smiles"] = compound.canonical_smiles
        row["formula"] = compound.formula
        row["mw"] = compound.molecular_weight

        # Get thermodynamic properties
        row["delta_hf_kJ_mol"] = compound.delta_hf_kJ_mol

        # Evaluate hazard
        result = compound.evaluate_hazard()

        row["delta_hd_kJ_mol"] = result.max_decomposition_kJ_mol
        row["delta_hd_cal_g"] = result.max_decomposition_cal_g
        row["ob_percent"] = result.oxygen_balance_percent
        row["hazard_class"] = result.hazard_class
        row["triggered_criteria"] = list(result.triggered_criteria)
        row["alerts"] = list(result.functional_group_alerts)
        row["gas_volume_L_g"] = result.gas_volume_L_g

    except PhoenixError as e:
        row["error"] = type(e).__name__
        row["error_message"] = str(e)
    except Exception as e:
        row["error"] = "UnexpectedError"
        row["error_message"] = str(e)

    return row
