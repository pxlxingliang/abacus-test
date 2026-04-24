import unittest
import numpy as np
import os
import tempfile
from pathlib import Path

from abacustest.lib_data.nao import AbacusNAO


class TestAbacusNAO(unittest.TestCase):
    """Test AbacusNAO class"""

    @classmethod
    def setUpClass(cls):
        cls.orb_file = os.path.join(
            os.path.dirname(__file__), "data", "orb", "H_gga_6au_100Ry_2s1p.orb"
        )

    def test_read_from_file(self):
        """Test reading NAO file"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertIsInstance(nao, AbacusNAO)

    def test_element(self):
        """Test element attribute"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(nao.element, "H")

    def test_energy_cutoff(self):
        """Test energy_cutoff attribute"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(nao.energy_cutoff, 100.0)

    def test_radius(self):
        """Test radius attribute"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(nao.radius, 6.0)

    def test_lmax(self):
        """Test lmax attribute"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(nao.lmax, 1)

    def test_l_orbs(self):
        """Test l_orbs attribute"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(nao.l_orbs, [2, 1])

    def test_mesh(self):
        """Test mesh attribute"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(nao.mesh, 601)

    def test_dr(self):
        """Test dr attribute"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(nao.dr, 0.01)

    def test_orbs_count(self):
        """Test number of orbitals"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(len(nao.orbs), 3)

    def test_orbs_data_length(self):
        """Test orbital data length"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        for orb in nao.orbs:
            self.assertEqual(len(orb["data"]), nao.mesh)

    def test_orbs_l_values(self):
        """Test orbital l values"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        expected_l = [0, 0, 1]
        actual_l = [orb["l"] for orb in nao.orbs]
        self.assertEqual(actual_l, expected_l)

    def test_orbs_n_values(self):
        """Test orbital n values"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        expected_n = [0, 1, 0]
        actual_n = [orb["n"] for orb in nao.orbs]
        self.assertEqual(actual_n, expected_n)

    def test_basis_num(self):
        """Test basis_num property"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(nao.basis_num, 5)

    def test_basis_num_calculation(self):
        """Test basis_num calculation matches expected formula"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        nbas = 0
        for i, norb in enumerate(nao.l_orbs):
            nbas += (2 * i + 1) * norb
        self.assertEqual(nao.basis_num, nbas)

    def test_get_orbital_label(self):
        """Test get_orbital_label method"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        self.assertEqual(nao.get_orbital_label(0, 0), "1s")
        self.assertEqual(nao.get_orbital_label(0, 1), "2s")
        self.assertEqual(nao.get_orbital_label(1, 0), "1p")
        self.assertEqual(nao.get_orbital_label(2, 0), "1d")

    def test_repr(self):
        """Test __repr__ method"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        repr_str = repr(nao)
        self.assertIn("AbacusNAO", repr_str)
        self.assertIn("element=H", repr_str)
        self.assertIn("energy_cutoff=100.0", repr_str)

    def test_direct_initialization(self):
        """Test direct initialization with all parameters"""
        orbs = [
            {"l": 0, "n": 0, "data": np.array([1.0, 2.0, 3.0])},
            {"l": 1, "n": 0, "data": np.array([0.5, 0.6, 0.7])},
        ]
        nao = AbacusNAO(
            element="Test",
            energy_cutoff=50.0,
            radius=5.0,
            lmax=1,
            l_orbs=[1, 1],
            mesh=3,
            dr=0.02,
            orbs=orbs,
        )
        self.assertEqual(nao.element, "Test")
        self.assertEqual(nao.energy_cutoff, 50.0)
        self.assertEqual(nao.radius, 5.0)
        self.assertEqual(nao.lmax, 1)
        self.assertEqual(nao.l_orbs, [1, 1])
        self.assertEqual(nao.mesh, 3)
        self.assertEqual(nao.dr, 0.02)
        self.assertEqual(len(nao.orbs), 2)

    def test_direct_initialization_with_default_orbs(self):
        """Test direct initialization with default orbs (None)"""
        nao = AbacusNAO(
            element="Test",
            energy_cutoff=50.0,
            radius=5.0,
            lmax=1,
            l_orbs=[1, 1],
            mesh=3,
            dr=0.02,
        )
        self.assertEqual(nao.orbs, [])

    def test_file_not_found(self):
        """Test reading non-existent file raises FileNotFoundError"""
        with self.assertRaises(FileNotFoundError):
            AbacusNAO.read_from_file("/non/existent/file.orb")

    def test_invalid_file_format(self):
        """Test reading invalid file format raises ValueError"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".orb", delete=False) as f:
            f.write("Invalid content\n")
            temp_file = f.name

        try:
            with self.assertRaises(ValueError):
                AbacusNAO.read_from_file(temp_file)
        finally:
            os.unlink(temp_file)

    def test_orbs_data_is_numpy_array(self):
        """Test that orbital data is stored as numpy array"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        for orb in nao.orbs:
            self.assertIsInstance(orb["data"], np.ndarray)

    def test_orbs_data_values(self):
        """Test orbital data values are reasonable"""
        nao = AbacusNAO.read_from_file(self.orb_file)
        for orb in nao.orbs:
            data = orb["data"]
            self.assertFalse(np.any(np.isnan(data)))
            self.assertFalse(np.any(np.isinf(data)))


if __name__ == "__main__":
    unittest.main()
