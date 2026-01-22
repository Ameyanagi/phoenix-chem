"""Integration tests for batch screening."""

from phoenix import screen
from phoenix.batch.screening import BatchResult


class TestBatchScreening:
    """Tests for batch screening functionality."""

    def test_screen_returns_batch_result(self):
        """Test screen() returns BatchResult."""
        smiles_list = ["CCO", "C"]
        result = screen(smiles_list)
        assert isinstance(result, BatchResult)

    def test_batch_result_has_dataframe(self):
        """Test BatchResult has DataFrame."""
        result = screen(["CCO"])
        assert hasattr(result, "dataframe")
        assert len(result.dataframe) == 1

    def test_dataframe_columns(self):
        """Test DataFrame has expected columns."""
        result = screen(["CCO"])
        expected_columns = [
            "smiles",
            "canonical_smiles",
            "formula",
            "mw",
            "hazard_class",
            "ob_percent",
        ]
        for col in expected_columns:
            assert col in result.dataframe.columns

    def test_successful_batch(self):
        """Test successful batch processing."""
        smiles_list = ["CCO", "C", "CC"]
        result = screen(smiles_list)

        assert result.successful == 3
        assert result.failed == 0
        assert len(result.dataframe) == 3

    def test_partial_failure(self):
        """Test partial failure handling."""
        smiles_list = ["CCO", "invalid_smiles", "C"]
        result = screen(smiles_list)

        assert result.successful == 2
        assert result.failed == 1
        assert len(result.dataframe) == 3

        # Check error row
        error_rows = result.dataframe[result.dataframe["error"].notna()]
        assert len(error_rows) == 1
        assert error_rows.iloc[0]["smiles"] == "invalid_smiles"

    def test_unsupported_element_error(self):
        """Test unsupported element error in batch."""
        smiles_list = ["CCO", "[Fe]", "C"]
        result = screen(smiles_list)

        assert result.failed == 1
        error_rows = result.dataframe[result.dataframe["error"].notna()]
        assert "UnsupportedElementError" in error_rows.iloc[0]["error"]

    def test_progress_callback(self):
        """Test progress callback is called."""
        smiles_list = ["CCO", "C", "CC"]
        progress_calls = []

        def callback(current, total):
            progress_calls.append((current, total))

        screen(smiles_list, progress_callback=callback)

        # Should be called for each compound plus final
        assert len(progress_calls) == 4
        assert progress_calls[-1] == (3, 3)

    def test_empty_list(self):
        """Test empty SMILES list."""
        result = screen([])
        assert result.successful == 0
        assert result.failed == 0
        assert len(result.dataframe) == 0


class TestBatchResultMethods:
    """Tests for BatchResult methods."""

    def test_to_json(self):
        """Test JSON export."""
        result = screen(["CCO"])
        json_str = result.to_json()
        assert isinstance(json_str, str)
        assert "CCO" in json_str or "OCC" in json_str  # canonical SMILES

    def test_summary(self):
        """Test summary statistics."""
        result = screen(["CCO", "invalid", "C"])
        summary = result.summary()

        assert summary["total_compounds"] == 3
        assert summary["successful"] == 2
        assert summary["failed"] == 1
        assert "hazard_class_counts" in summary

    def test_summary_empty_batch(self):
        """Test summary for empty batch."""
        result = screen([])
        summary = result.summary()

        assert summary["total_compounds"] == 0
        assert summary["delta_hd_cal_g_min"] is None


class TestBatchDataFrameOperations:
    """Tests for DataFrame operations on batch results."""

    def test_filter_by_hazard_class(self):
        """Test filtering DataFrame by hazard class."""
        # Use compounds that we know will have different classifications
        smiles_list = [
            "CCO",  # ethanol - likely LOW
            "c1ccccc1[N+](=O)[O-]",  # nitrobenzene - likely MEDIUM/HIGH
        ]
        result = screen(smiles_list)

        # Filter for specific class
        filtered = result.dataframe[result.dataframe["hazard_class"].notna()]
        assert len(filtered) == 2

    def test_sort_by_delta_hd(self):
        """Test sorting by decomposition energy."""
        smiles_list = ["CCO", "C", "CC"]
        result = screen(smiles_list)

        sorted_df = result.dataframe.sort_values("delta_hd_cal_g")
        assert len(sorted_df) == 3

    def test_successful_rows_have_values(self):
        """Test successful rows have all numeric values."""
        result = screen(["CCO"])
        row = result.dataframe.iloc[0]

        assert row["mw"] is not None
        assert row["delta_hf_kJ_mol"] is not None
        assert row["ob_percent"] is not None
        assert row["hazard_class"] is not None

    def test_error_rows_have_null_values(self):
        """Test error rows have null numeric values."""
        result = screen(["invalid"])
        row = result.dataframe.iloc[0]

        assert row["error"] is not None
        assert row["mw"] is None
        assert row["hazard_class"] is None


class TestBatchWithVariousCompounds:
    """Tests for batch processing with various compound types."""

    def test_batch_with_nitro_compounds(self):
        """Test batch processing of nitro compounds."""
        smiles_list = [
            "c1ccccc1[N+](=O)[O-]",  # nitrobenzene
            "CC(=O)O",  # acetic acid (for comparison)
        ]
        result = screen(smiles_list)

        assert result.successful == 2
        # Nitrobenzene should have functional group alerts
        nitro_row = result.dataframe[result.dataframe["formula"] == "C6H5NO2"].iloc[0]
        assert len(nitro_row["alerts"]) > 0

    def test_batch_with_halogens(self):
        """Test batch with halogenated compounds."""
        smiles_list = ["CF", "CCl", "CBr"]
        result = screen(smiles_list)

        assert result.successful == 3
        # All should have valid oxygen balance
        assert result.dataframe["ob_percent"].notna().all()

    def test_batch_reproducibility(self):
        """Test batch results are reproducible."""
        smiles_list = ["CCO", "C", "CC"]

        result1 = screen(smiles_list)
        result2 = screen(smiles_list)

        # Results should be identical
        assert result1.successful == result2.successful
        assert list(result1.dataframe["mw"]) == list(result2.dataframe["mw"])
