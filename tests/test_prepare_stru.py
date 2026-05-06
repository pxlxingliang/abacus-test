import unittest
import numpy as np
import os
import shutil
import tempfile
from pathlib import Path

from abacustest.lib_prepare.stru import AbacusAtomType, AbacusATOM, AbacusSTRU


class TestAbacusAtomType(unittest.TestCase):
    """Test AbacusAtomType class"""

    def setUp(self):
        self.atomtype1 = AbacusAtomType(label="Fe", element="Fe", mass=55.845, pp="Fe.upf")
        self.atomtype2 = AbacusAtomType(label="Fe", element="Fe", mass=55.845, pp="Fe.upf")
        self.atomtype3 = AbacusAtomType(label="Fe", element="Fe", mass=55.846, pp="Fe.upf")
        self.atomtype4 = AbacusAtomType(label="Fe", element="Fe", mass=55.845, pp="Fe.upf", type_mag=2.0)

    def test_atomtype_equality_same(self):
        self.assertEqual(self.atomtype1, self.atomtype2)

    def test_atomtype_equality_different_mass(self):
        self.assertNotEqual(self.atomtype1, self.atomtype3)

    def test_atomtype_equality_different_type_mag(self):
        self.assertNotEqual(self.atomtype1, self.atomtype4)

    def test_atomtype_default_list(self):
        default_list = self.atomtype1._default_list()
        self.assertEqual(default_list[0], "Fe")
        self.assertEqual(default_list[1], "Fe")
        self.assertEqual(default_list[2], 55.845)
        self.assertEqual(default_list[3], "Fe.upf")

    def test_atomtype_default_list_none_values(self):
        atomtype = AbacusAtomType(label="Fe")
        default_list = atomtype._default_list()
        self.assertEqual(default_list[2], 0.0)
        self.assertEqual(default_list[3], "")
        self.assertEqual(default_list[4], "")

    def test_atomtype_sort(self):
        atomtype_a = AbacusAtomType(label="A")
        atomtype_b = AbacusAtomType(label="B")
        self.assertTrue(atomtype_a < atomtype_b)


class TestAbacusATOM(unittest.TestCase):
    """Test AbacusATOM class"""

    def setUp(self):
        self.atom = AbacusATOM(label="Fe1", coord=(0.0, 0.0, 0.0))
        self.atom2 = AbacusATOM(label="Fe2", coord=(0.5, 0.5, 0.5), mag=2.0)
        self.atom3 = AbacusATOM(label="O1", coord=(1.0, 0.0, 0.0), angle1=90.0, angle2=45.0)

    def test_infer_element_from_label(self):
        self.assertEqual(self.atom.element, "Fe")
        self.assertEqual(self.atom2.element, "Fe")
        atom_h = AbacusATOM(label="H1", coord=(0, 0, 0))
        self.assertEqual(atom_h.element, "H")

    def test_mass_from_dict(self):
        self.assertGreater(self.atom.mass, 0)

    def test_noncolinear_false(self):
        self.assertFalse(self.atom.noncolinear)
        atom_mag = AbacusATOM(label="Fe1", coord=(0, 0, 0), mag=1.0)
        self.assertFalse(atom_mag.noncolinear)

    def test_noncolinear_true_angle(self):
        self.assertTrue(self.atom3.noncolinear)

    def test_noncolinear_true_vector(self):
        atom_vec = AbacusATOM(label="Fe1", coord=(0, 0, 0), mag=(1.0, 0.0, 0.0))
        self.assertTrue(atom_vec.noncolinear)

    def test_atommag_magnitude_no_mag(self):
        self.assertEqual(self.atom.atommag_magnitude, 0.0)

    def test_atommag_magnitude_with_mag(self):
        self.assertAlmostEqual(self.atom2.atommag_magnitude, 2.0)

    def test_atommag_magnitude_vector(self):
        atom_vec = AbacusATOM(label="Fe1", coord=(0, 0, 0), mag=(1.0, 0.0, 0.0))
        self.assertAlmostEqual(atom_vec.atommag_magnitude, 1.0)

    def test_set_atommag_scalar(self):
        self.atom.set_atommag(3.0)
        self.assertEqual(self.atom.mag, 3.0)

    def test_set_atommag_vector(self):
        self.atom.set_atommag((1.0, 1.0, 1.0))
        self.assertEqual(self.atom.mag, (1.0, 1.0, 1.0))

    def test_set_atommag_with_angles(self):
        self.atom.set_atommag(2.0, angle1=90.0, angle2=0.0)
        self.assertEqual(self.atom.mag, 2.0)
        self.assertEqual(self.atom.angle1, 90.0)
        self.assertEqual(self.atom.angle2, 0.0)

    def test_atomtype_property(self):
        atomtype = self.atom.atomtype
        self.assertIsInstance(atomtype, AbacusAtomType)
        self.assertEqual(atomtype.label, "Fe1")
        self.assertEqual(atomtype.element, "Fe")

    def test_sort_keep_order(self):
        atoms = [
            AbacusATOM(label="A1", coord=(0, 0, 0)),
            AbacusATOM(label="B1", coord=(1, 1, 1)),
            AbacusATOM(label="A2", coord=(2, 2, 2)),
            AbacusATOM(label="B2", coord=(3, 3, 3)),
        ]
        sorted_atoms, indices = AbacusATOM.sort(atoms, keep_first_order=True)
        self.assertEqual(sorted_atoms[0].label, "A1")
        # Test keeps original order (first appearance)

    def test_sort_new_order(self):
        atoms = [
            AbacusATOM(label="B1", coord=(0, 0, 0)),
            AbacusATOM(label="A1", coord=(1, 1, 1)),
        ]
        sorted_atoms, indices = AbacusATOM.sort(atoms, keep_first_order=False)
        self.assertEqual(sorted_atoms[0].label, "A1")

    def test_find_uniq_atomtypes(self):
        atoms = [
            AbacusATOM(label="Fe1", coord=(0, 0, 0)),
            AbacusATOM(label="Fe1", coord=(1, 1, 1)),
            AbacusATOM(label="O1", coord=(2, 2, 2)),
        ]
        unique_types = AbacusATOM.find_uniq_atomtypes(atoms)
        self.assertEqual(len(unique_types), 2)
        self.assertEqual(unique_types[0].natom, 2)
        self.assertEqual(unique_types[1].natom, 1)

    def test_same_type_true(self):
        atoms = [
            AbacusATOM(label="Fe1", coord=(0, 0, 0)),
            AbacusATOM(label="Fe1", coord=(1, 1, 1)),
        ]
        self.assertTrue(AbacusATOM.same_type(atoms))

    def test_same_type_false(self):
        atoms = [
            AbacusATOM(label="Fe1", coord=(0, 0, 0)),
            AbacusATOM(label="O1", coord=(1, 1, 1)),
        ]
        self.assertFalse(AbacusATOM.same_type(atoms))

    def test_same_type_empty(self):
        self.assertTrue(AbacusATOM.same_type([]))
        self.assertTrue(AbacusATOM.same_type([self.atom]))

    def test_rotate(self):
        rot_mat = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        atom = AbacusATOM(label="Fe1", coord=(1.0, 2.0, 3.0))
        atom.rotate(rot_mat)
        self.assertAlmostEqual(atom.coord[0], 2.0)
        self.assertAlmostEqual(atom.coord[1], 1.0)
        self.assertAlmostEqual(atom.coord[2], 3.0)


class TestAbacusSTRU(unittest.TestCase):
    """Test AbacusSTRU class"""

    def setUp(self):
        self.atom1 = AbacusATOM(label="Fe1", coord=(0.0, 0.0, 0.0))
        self.atom2 = AbacusATOM(label="Fe2", coord=(1.435, 1.435, 1.435))
        self.cell = [[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]]
        self.stru = AbacusSTRU(cell=self.cell, atoms=[self.atom1, self.atom2])
        self.work_dir = tempfile.mkdtemp()

    def tearDown(self):
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)

    def test_natoms(self):
        self.assertEqual(self.stru.natoms, 2)

    def test_labels(self):
        self.assertEqual(self.stru.labels, ["Fe1", "Fe2"])

    def test_elements(self):
        self.assertEqual(self.stru.elements, ["Fe", "Fe"])

    def test_coords(self):
        coords = self.stru.coords
        self.assertEqual(len(coords), 2)
        self.assertEqual(coords[0], (0.0, 0.0, 0.0))

    def test_coords_direct(self):
        coords_direct = self.stru.coords_direct
        self.assertEqual(len(coords_direct), 2)

    def test_append(self):
        new_atom = AbacusATOM(label="O1", coord=(0.5, 0.5, 0.0))
        self.stru.append(new_atom)
        self.assertEqual(self.stru.natoms, 3)

    def test_extend(self):
        new_atoms = [
            AbacusATOM(label="O1", coord=(0.5, 0.5, 0.0)),
            AbacusATOM(label="O2", coord=(0.6, 0.6, 0.0)),
        ]
        self.stru.extend(new_atoms)
        self.assertEqual(self.stru.natoms, 4)

    def test_insert(self):
        new_atom = AbacusATOM(label="O1", coord=(0.5, 0.5, 0.0))
        self.stru.insert(1, new_atom)
        self.assertEqual(self.stru.natoms, 3)
        self.assertEqual(self.stru.labels[1], "O1")

    def test_delitem(self):
        del self.stru[0]
        self.assertEqual(self.stru.natoms, 1)
        self.assertEqual(self.stru.labels[0], "Fe2")

    def test_getitem(self):
        atom = self.stru[0]
        self.assertEqual(atom.label, "Fe1")
        atoms = self.stru[0:2]
        self.assertEqual(len(atoms), 2)

    def test_setitem(self):
        new_atom = AbacusATOM(label="O1", coord=(0.5, 0.5, 0.0))
        self.stru[0] = new_atom
        self.assertEqual(self.stru.labels[0], "O1")

    def test_sort(self):
        atom_a = AbacusATOM(label="A1", coord=(0, 0, 0))
        atom_b = AbacusATOM(label="B1", coord=(1, 1, 1))
        stru = AbacusSTRU(cell=self.cell, atoms=[atom_b, atom_a])
        indices = stru.sort(keep_first_order=False)
        self.assertEqual(stru.labels[0], "A1")

    def test_set_pp(self):
        self.stru.set_pp({"Fe": "Fe_pbe.upf"})
        self.assertEqual(self.stru.pps[0], "Fe_pbe.upf")

    def test_set_orb(self):
        self.stru.set_orb({"Fe": "Fe.orb"})
        self.assertEqual(self.stru.orbs[0], "Fe.orb")

    def test_pp_dict(self):
        pp_dict = self.stru.pp_dict()
        self.assertIn("Fe1", pp_dict)

    def test_orb_dict(self):
        orb_dict = self.stru.orb_dict()
        self.assertIn("Fe1", orb_dict)

    def test_get_cell_param(self):
        a, b, c, alpha, beta, gamma = self.stru.get_cell_param()
        self.assertAlmostEqual(a, 2.87, delta=0.01)
        self.assertAlmostEqual(alpha, 90.0)
        self.assertAlmostEqual(beta, 90.0)
        self.assertAlmostEqual(gamma, 90.0)

    def test_write_read_stru(self):
        stru_file = os.path.join(self.work_dir, "test.STRU")
        self.stru.write(stru_file, fmt="stru")
        self.assertTrue(os.path.exists(stru_file))

        read_stru = AbacusSTRU.read(stru_file, fmt="stru")
        self.assertEqual(read_stru.natoms, self.stru.natoms)

    def test_write_read_poscar(self):
        poscar_file = os.path.join(self.work_dir, "POSCAR")
        self.stru.write(poscar_file, fmt="poscar")
        self.assertTrue(os.path.exists(poscar_file))

        read_poscar = AbacusSTRU.read(poscar_file, fmt="poscar")
        self.assertEqual(read_poscar.natoms, self.stru.natoms)

    def test_fix_atom_by_index(self):
        self.stru.fix_atom_by_index([0])
        self.assertEqual(self.stru.moves[0], (False, False, False))
        self.assertEqual(self.stru.moves[1], (True, True, True))

    def test_fix_atom_by_coord(self):
        self.stru.fix_atom_by_coord(min=0, max=0.5, direction=2)
        self.assertEqual(self.stru.moves[0], (False, False, False))

    def test_from_ase(self):
        from ase import Atoms
        ase_atoms = Atoms("Fe2", positions=[(0, 0, 0), (1.5, 1.5, 1.5)], cell=self.cell)
        stru = AbacusSTRU.from_ase(ase_atoms)
        self.assertEqual(stru.natoms, 2)

    def test_to_ase(self):
        ase_stru = self.stru.to("ase")
        self.assertEqual(len(ase_stru), 2)

    def test_create_subset(self):
        subset = self.stru.create_subset([0])
        self.assertEqual(subset.natoms, 1)
        self.assertEqual(subset.labels[0], "Fe1")

    def test_atomtypes(self):
        atomtypes = self.stru.atomtypes
        self.assertEqual(len(atomtypes), 2)

    def test_uniq_atomtypes(self):
        unique_types = self.stru.uniq_atomtypes
        self.assertIsInstance(unique_types[0], AbacusAtomType)

    def test_masses(self):
        masses = self.stru.masses
        self.assertEqual(len(masses), 2)
        self.assertGreater(masses[0], 0)

    def test_pps(self):
        self.assertEqual(len(self.stru.pps), 2)

    def test_orbs(self):
        self.assertEqual(len(self.stru.orbs), 2)

    def test_supercell_2x2x2(self):
        stru_super = self.stru.supercell([2, 2, 2])
        self.assertEqual(stru_super.natoms, 16)
        new_cell = stru_super.cell
        self.assertAlmostEqual(new_cell[0][0], 5.74)
        self.assertAlmostEqual(new_cell[1][1], 5.74)
        self.assertAlmostEqual(new_cell[2][2], 5.74)

    def test_supercell_1x2x1(self):
        stru_super = self.stru.supercell([1, 2, 1])
        self.assertEqual(stru_super.natoms, 4)

    def test_supercell_preserves_magnetic_moment(self):
        atom_mag = AbacusATOM(label="Fe1", coord=(0.0, 0.0, 0.0), mag=2.0)
        stru_mag = AbacusSTRU(cell=self.cell, atoms=[atom_mag])
        stru_super = stru_mag.supercell([2, 1, 1])
        self.assertEqual(stru_super.natoms, 2)
        self.assertEqual(stru_super[0].mag, 2.0)
        self.assertEqual(stru_super[1].mag, 2.0)


if __name__ == "__main__":
    unittest.main()
