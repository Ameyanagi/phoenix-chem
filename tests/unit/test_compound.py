"""Unit tests for Compound class."""

import pytest

from phoenix import Compound, InvalidSmilesError, UnsupportedElementError, UnsupportedStructureError


class TestCompoundCreation:
    """Tests for Compound.from_smiles()."""

    def test_valid_smiles_ethanol(self):
        """Test parsing valid SMILES for ethanol."""
        compound = Compound.from_smiles("CCO")
        assert compound.formula == "C2H6O"
        assert abs(compound.molecular_weight - 46.07) < 0.1

    def test_valid_smiles_methane(self):
        """Test parsing valid SMILES for methane."""
        compound = Compound.from_smiles("C")
        assert compound.formula == "CH4"
        assert abs(compound.molecular_weight - 16.04) < 0.1

    def test_valid_smiles_nitrobenzene(self):
        """Test parsing aromatic nitrobenzene."""
        compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
        assert compound.formula == "C6H5NO2"
        assert compound.composition == {"C": 6, "H": 5, "N": 1, "O": 2}

    def test_valid_smiles_tnt(self):
        """Test parsing TNT."""
        # 2,4,6-trinitrotoluene
        compound = Compound.from_smiles("Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]")
        assert compound.composition["C"] == 7
        assert compound.composition["H"] == 5
        assert compound.composition["N"] == 3
        assert compound.composition["O"] == 6

    def test_invalid_smiles_raises_error(self):
        """Test that invalid SMILES raises InvalidSmilesError."""
        with pytest.raises(InvalidSmilesError) as exc_info:
            Compound.from_smiles("invalid_smiles")
        assert "invalid_smiles" in str(exc_info.value)

    def test_unsupported_element_metal(self):
        """Test that metals raise UnsupportedElementError."""
        with pytest.raises(UnsupportedElementError) as exc_info:
            Compound.from_smiles("[Fe]")
        assert "Fe" in exc_info.value.elements

    def test_unsupported_element_iodine(self):
        """Test that iodine raises UnsupportedElementError in MVP."""
        with pytest.raises(UnsupportedElementError) as exc_info:
            Compound.from_smiles("CI")  # iodomethane
        assert "I" in exc_info.value.elements
        assert "Iodine" in str(exc_info.value)

    def test_charged_species_rejected(self):
        """Test that charged species are rejected."""
        with pytest.raises(UnsupportedStructureError) as exc_info:
            Compound.from_smiles("[NH4+]")  # ammonium cation
        assert "Charged" in str(exc_info.value)

    def test_radical_rejected(self):
        """Test that radicals are rejected."""
        with pytest.raises(UnsupportedStructureError):
            Compound.from_smiles("[CH3]")  # methyl radical


class TestCompoundProperties:
    """Tests for Compound property access."""

    def test_canonical_smiles(self):
        """Test canonical SMILES generation."""
        # Different representations of same molecule
        c1 = Compound.from_smiles("CCO")
        c2 = Compound.from_smiles("OCC")
        assert c1.canonical_smiles == c2.canonical_smiles

    def test_composition_ethanol(self):
        """Test elemental composition extraction."""
        compound = Compound.from_smiles("CCO")
        assert compound.composition == {"C": 2, "H": 6, "O": 1}

    def test_num_atoms(self):
        """Test atom count."""
        compound = Compound.from_smiles("CCO")
        assert compound.num_atoms == 9  # 2C + 6H + 1O

    def test_rdmol_access(self):
        """Test RDKit Mol access."""
        compound = Compound.from_smiles("CCO")
        assert compound.rdmol is not None
        assert compound.rdmol.GetNumAtoms() == 9

    def test_original_smiles(self):
        """Test original SMILES preservation."""
        smiles = "CCO"
        compound = Compound.from_smiles(smiles)
        assert compound.original_smiles == smiles


class TestCompoundEquality:
    """Tests for Compound equality and hashing."""

    def test_equal_compounds(self):
        """Test that compounds with same canonical SMILES are equal."""
        c1 = Compound.from_smiles("CCO")
        c2 = Compound.from_smiles("OCC")
        assert c1 == c2

    def test_hash_equal_compounds(self):
        """Test that equal compounds have same hash."""
        c1 = Compound.from_smiles("CCO")
        c2 = Compound.from_smiles("OCC")
        assert hash(c1) == hash(c2)

    def test_compound_in_set(self):
        """Test that compounds work in sets."""
        c1 = Compound.from_smiles("CCO")
        c2 = Compound.from_smiles("OCC")
        c3 = Compound.from_smiles("C")
        s = {c1, c2, c3}
        assert len(s) == 2  # c1 and c2 are same


class TestLargeMoleculeWarning:
    """Tests for large molecule warning."""

    def test_large_molecule_warning(self):
        """Test warning for molecules >100 atoms."""
        # Create a large molecule (long alkane chain)
        smiles = "C" * 50  # 50 carbons = 152 atoms with H
        with pytest.warns(UserWarning, match=r"Large molecule.*atoms"):
            compound = Compound.from_smiles(smiles)
        assert len(compound.warnings) > 0

    def test_normal_molecule_no_warning(self):
        """Test no warning for normal-sized molecules."""
        import warnings as w_module

        with w_module.catch_warnings(record=True) as record:
            w_module.simplefilter("always")
            compound = Compound.from_smiles("CCO")
        # Filter for UserWarning about large molecules
        large_mol_warnings = [w for w in record if "Large molecule" in str(w.message)]
        assert len(large_mol_warnings) == 0
        assert len(compound.warnings) == 0


class TestSupportedElements:
    """Tests for supported element coverage."""

    def test_supported_elements_chno(self):
        """Test CHNO compounds are supported."""
        compound = Compound.from_smiles("CC(=O)O")  # acetic acid
        assert set(compound.composition.keys()).issubset(
            {"C", "H", "N", "O", "S", "P", "F", "Cl", "Br"}
        )

    def test_supported_element_sulfur(self):
        """Test sulfur is supported."""
        compound = Compound.from_smiles("CS")  # methanethiol
        assert "S" in compound.composition

    def test_supported_element_phosphorus(self):
        """Test phosphorus is supported."""
        compound = Compound.from_smiles("CP")  # methylphosphine
        assert "P" in compound.composition

    def test_supported_element_fluorine(self):
        """Test fluorine is supported."""
        compound = Compound.from_smiles("CF")  # fluoromethane
        assert "F" in compound.composition

    def test_supported_element_chlorine(self):
        """Test chlorine is supported."""
        compound = Compound.from_smiles("CCl")  # chloromethane
        assert "Cl" in compound.composition

    def test_supported_element_bromine(self):
        """Test bromine is supported."""
        compound = Compound.from_smiles("CBr")  # bromomethane
        assert "Br" in compound.composition

    def test_multi_element_compound(self):
        """Test compound with multiple supported elements."""
        # Hypothetical multi-element compound
        compound = Compound.from_smiles("CC(F)(Cl)Br")
        assert compound.composition["C"] == 2
        assert compound.composition["F"] == 1
        assert compound.composition["Cl"] == 1
        assert compound.composition["Br"] == 1
