import unittest
import numpy as np
import os
import shutil
import tempfile
from pathlib import Path

from abacustest.lib_prepare import abacus


class TestAbacusStru(unittest.TestCase):
    """Test AbacusStru class"""

    def setUp(self):
        self.work_dir = tempfile.mkdtemp()
        self.cell = [[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]]
        self.coord = [[0, 0, 0], [1.435, 1.435, 1.435]]
        self.stru = abacus.AbacusStru(
            label=["Fe"],
            cell=self.cell,
            coord=self.coord,
            atom_number=[2],
            pp=["Fe.upf"],
            element=["Fe"],
            magmom=[1.0],
            lattice_constant=1.0,
        )

    def tearDown(self):
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)

    def test_get_natoms(self):
        self.assertEqual(self.stru.get_natoms(), 2)

    def test_get_pp(self):
        pp = self.stru.get_pp()
        self.assertEqual(pp, ["Fe.upf"])

    def test_get_pp_total(self):
        pp_total = self.stru.get_pp(total=True)
        self.assertEqual(len(pp_total), 2)

    def test_get_orb(self):
        orb = self.stru.get_orb()
        self.assertIsNone(orb)

    def test_get_paw(self):
        paw = self.stru.get_paw()
        self.assertIsNone(paw)

    def test_get_label(self):
        label = self.stru.get_label(total=True)
        self.assertEqual(label, ["Fe", "Fe"])

    def test_get_label_not_total(self):
        label = self.stru.get_label(total=False)
        self.assertEqual(label, ["Fe"])

    def test_get_mass(self):
        mass = self.stru.get_mass()
        self.assertEqual(len(mass), 1)
        self.assertGreater(mass[0], 0)

    def test_get_element(self):
        element = self.stru.get_element(total=True, number=False)
        self.assertEqual(element, ["Fe", "Fe"])

    def test_get_element_with_number(self):
        element = self.stru.get_element(total=True, number=True)
        self.assertEqual(element[0], 26)

    def test_get_cell(self):
        cell = self.stru.get_cell(bohr=False)
        self.assertEqual(len(cell), 3)

    def test_get_cell_bohr(self):
        cell = self.stru.get_cell(bohr=True)
        self.assertEqual(len(cell), 3)

    def test_get_cell_param(self):
        a, b, c, alpha, beta, gamma = self.stru.get_cell_param()
        self.assertGreater(a, 0)
        self.assertAlmostEqual(alpha, 90.0)

    def test_get_volume(self):
        volume = self.stru.get_volume()
        self.assertGreater(volume, 0)

    def test_get_coord(self):
        coord = self.stru.get_coord(bohr=False, direct=False)
        self.assertEqual(len(coord), 2)

    def test_get_coord_direct(self):
        coord = self.stru.get_coord(bohr=False, direct=True)
        self.assertEqual(len(coord), 2)

    def test_get_mag(self):
        mag = self.stru.get_mag()
        self.assertEqual(mag, [1.0])

    def test_get_move(self):
        move = self.stru.get_move()
        self.assertEqual(len(move), 6)

    def test_globalidx2labelidx(self):
        label, idx = self.stru.globalidx2labelidx(0)
        self.assertEqual(label, "Fe")
        self.assertEqual(idx, 0)

    def test_globalidx2labelidx_second_atom(self):
        label, idx = self.stru.globalidx2labelidx(1)
        self.assertEqual(label, "Fe")
        self.assertEqual(idx, 1)

    def test_globalidx2labelidx_invalid_index(self):
        with self.assertRaises(ValueError):
            self.stru.globalidx2labelidx(10)

    def test_globalidx2labelidx_negative_index(self):
        with self.assertRaises(ValueError):
            self.stru.globalidx2labelidx(-1)

    def test_labelidx2globalidx(self):
        idx = self.stru.labelidx2globalidx("Fe", 1)
        self.assertEqual(idx, 1)

    def test_labelidx2globalidx_invalid_label(self):
        with self.assertRaises(AssertionError):
            self.stru.labelidx2globalidx("Invalid", 0)

    def test_set_pp(self):
        self.stru.set_pp(["NewFe.upf"])
        pp = self.stru.get_pp()
        self.assertEqual(pp, ["NewFe.upf"])

    def test_set_orb(self):
        self.stru.set_orb(["Fe.orb"])
        orb = self.stru.get_orb()
        self.assertEqual(orb, ["Fe.orb"])

    def test_set_mass(self):
        self.stru.set_mass([56.0])
        mass = self.stru.get_mass()
        self.assertEqual(mass[0], 56.0)

    def test_set_element(self):
        self.stru.set_element(["Co"])
        element = self.stru.get_element(total=False, number=True)
        self.assertEqual(element[0], 27)

    def test_set_atommag(self):
        self.stru.set_atommag([2.0, 2.0])
        mag = self.stru.get_atommag()
        self.assertEqual(mag, [2.0, 2.0])

    def test_set_atommag_none(self):
        self.stru.set_atommag(None)
        mag = self.stru.get_atommag()
        self.assertEqual(mag, [1.0, 1.0])

    def test_set_angle1(self):
        self.stru.set_angle1([90.0, 90.0])
        angle1 = self.stru.get_angle1()
        self.assertEqual(angle1, [90.0, 90.0])

    def test_set_angle2(self):
        self.stru.set_angle2([45.0, 45.0])
        angle2 = self.stru.get_angle2()
        self.assertEqual(angle2, [45.0, 45.0])

    def test_set_constrain(self):
        self.stru.set_constrain([[True, True, False], [False, False, False]])
        constrain = self.stru.get_constrain()
        self.assertEqual(constrain[0], [True, True, False])

    def test_get_constrain(self):
        constrain = self.stru.get_constrain()
        self.assertEqual(len(constrain), 2)

    def test_get_isconstrain(self):
        self.stru.set_constrain([True, False])
        is_constrain = self.stru.get_isconstrain()
        self.assertTrue(is_constrain[0])

    def test_get_lambda(self):
        self.stru.set_constrain([True, False])
        lambda_vals = self.stru.get_lambda()
        self.assertEqual(len(lambda_vals), 2)

    def test_get_atommag(self):
        mag = self.stru.get_atommag()
        self.assertEqual(len(mag), 2)

    def test_get_atommag_norm(self):
        mag = self.stru.get_atommag(norm=True)
        self.assertEqual(len(mag), 2)

    def test_set_dpks(self):
        self.stru.set_dpks("jle.orb")
        dpks = self.stru.get_dpks()
        self.assertEqual(dpks, "jle.orb")

    def test_supercell_2x1x1(self):
        new_stru = self.stru.supercell([2, 1, 1])
        self.assertEqual(new_stru.get_natoms(), 4)

    def test_supercell_2x2x2(self):
        new_stru = self.stru.supercell([2, 2, 2])
        self.assertEqual(new_stru.get_natoms(), 16)

    def test_delete_atom_all(self):
        new_stru = self.stru.delete_atom("Fe")
        self.assertEqual(new_stru.get_natoms(), 0)

    def test_delete_atom_by_index(self):
        new_stru = self.stru.delete_atom("Fe", idx=0)
        self.assertEqual(new_stru.get_natoms(), 1)

    def test_set_empty_atom(self):
        self.stru.set_empty_atom("Fe", idx=0)
        label = self.stru.get_label(total=False)
        self.assertIn("Fe_empty", label)

    def test_set_empty_atom_all(self):
        self.stru.set_empty_atom("Fe")
        label = self.stru.get_label(total=False)
        self.assertIn("Fe_empty", label)

    def test_write_stru(self):
        stru_file = os.path.join(self.work_dir, "test.STRU")
        self.stru.write(stru_file)
        self.assertTrue(os.path.exists(stru_file))

    def test_write_stru_direct(self):
        stru_file = os.path.join(self.work_dir, "test_direct.STRU")
        self.stru.write(stru_file, direct=True)
        self.assertTrue(os.path.exists(stru_file))

    def test_write_stru_cartesian(self):
        stru_file = os.path.join(self.work_dir, "test_cart.STRU")
        self.stru.write(stru_file, direct=False)
        self.assertTrue(os.path.exists(stru_file))

    def test_write2poscar(self):
        poscar_file = os.path.join(self.work_dir, "POSCAR")
        self.stru.write2poscar(poscar_file)
        self.assertTrue(os.path.exists(poscar_file))

    def test_to_ase(self):
        atoms = self.stru.to_ase()
        self.assertEqual(len(atoms), 2)

    def test_to_ase_empty2x(self):
        self.stru.set_empty_atom("Fe")
        atoms = self.stru.to_ase(empty2x=True)
        self.assertEqual(len(atoms), 2)

    def test_split_list(self):
        alist = [1, 2, 3, 4, 5]
        indices = [2, 3]
        result = self.stru.split_list(alist, indices)
        self.assertEqual(result, [[1, 2], [3, 4, 5]])

    def test_set_coord_cartesian(self):
        new_coord = [[0.5, 0.5, 0.5], [1.0, 1.0, 1.0]]
        self.stru.set_coord(new_coord, direct=False, bohr=False)
        coord = self.stru.get_coord(bohr=False, direct=False)
        self.assertEqual(len(coord), 2)

    def test_set_coord_direct(self):
        new_coord = [[0, 0, 0], [0.5, 0.5, 0.5]]
        self.stru.set_coord(new_coord, direct=True, bohr=False)
        coord = self.stru.get_coord(bohr=False, direct=True)
        self.assertEqual(len(coord), 2)

    def test_set_cell(self):
        new_cell = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]]
        self.stru.set_cell(new_cell, bohr=False, change_coord=True)
        cell = self.stru.get_cell(bohr=False)
        self.assertEqual(len(cell), 3)

    def test_get_stru(self):
        stru_dict = self.stru.get_stru()
        self.assertIn("cell", stru_dict)
        self.assertIn("coord", stru_dict)

    def test_perturb_stru_no_perturb(self):
        new_stru_list = self.stru.perturb_stru(0)
        self.assertEqual(len(new_stru_list), 1)

    def test_perturb_stru_cell(self):
        new_stru_list = self.stru.perturb_stru(2, cell_pert_frac=0.1)
        self.assertEqual(len(new_stru_list), 2)

    def test_perturb_stru_atom(self):
        new_stru_list = self.stru.perturb_stru(2, atom_pert_dist=0.1)
        self.assertEqual(len(new_stru_list), 2)

    def test_check_passes(self):
        self.assertTrue(self.stru._check())

    def test_get_dpks_none(self):
        dpks = self.stru.get_dpks()
        self.assertIsNone(dpks)


class TestAbacusStruMagToAngle(unittest.TestCase):
    """Test mag_to_angle and angle_to_mag methods"""

    def setUp(self):
        self.stru = abacus.AbacusStru(
            label=["Fe"],
            cell=[[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]],
            coord=[[0, 0, 0]],
            element=["Fe"],
        )

    def test_mag_to_angle_z_direction(self):
        angle1, angle2 = self.stru.mag_to_angle(0.0, 0.0, 1.0)
        self.assertAlmostEqual(angle1, 0.0, places=5)
        self.assertAlmostEqual(angle2, 0.0, places=5)

    def test_mag_to_angle_x_direction(self):
        angle1, angle2 = self.stru.mag_to_angle(1.0, 0.0, 0.0)
        self.assertAlmostEqual(angle1, 90.0, places=5)

    def test_angle_to_mag_z(self):
        mag = self.stru.angle_to_mag(1.0, 0.0, 0.0)
        self.assertAlmostEqual(mag[2], 1.0, places=5)

    def test_angle_to_mag_xy(self):
        mag = self.stru.angle_to_mag(1.0, 90.0, 0.0)
        self.assertAlmostEqual(mag[0], 1.0, places=5)


class TestAbacusStruMultipleAtomTypes(unittest.TestCase):
    """Test AbacusStru with multiple atom types"""

    def setUp(self):
        self.stru = abacus.AbacusStru(
            label=["Fe", "O", "Fe"],
            cell=[[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]],
            coord=[
                [0, 0, 0],
                [5, 5, 5],
                [2.5, 2.5, 2.5],
                [7.5, 7.5, 7.5],
                [3, 3, 3],
                [6, 6, 6],
            ],
            atom_number=[1, 2, 3],
            pp=["Fe.upf", "O.upf", "Fe.upf"],
            element=["Fe", "O", "Fe"],
        )

    def test_get_natoms(self):
        self.assertEqual(self.stru.get_natoms(), 6)

    def test_get_label_total(self):
        labels = self.stru.get_label(total=True)
        self.assertEqual(len(labels), 6)

    def test_get_element_total(self):
        elements = self.stru.get_element(total=True, number=False)
        self.assertEqual(elements.count("Fe"), 4)
        self.assertEqual(elements.count("O"), 2)

    def test_delete_one_type(self):
        new_stru = self.stru.delete_atom("O")
        self.assertEqual(new_stru.get_natoms(), 4)


class TestAbacusStruDirectCoordinates(unittest.TestCase):
    """Test AbacusStru with direct coordinates"""

    def setUp(self):
        self.stru = abacus.AbacusStru(
            label=["Fe"],
            cell=[[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]],
            coord=[[0, 0, 0], [0.5, 0.5, 0.5]],
            atom_number=[2],
            element=["Fe"],
            cartesian=False,
            lattice_constant=1.0,
        )

    def test_get_coord_direct(self):
        coord = self.stru.get_coord(bohr=False, direct=True)
        self.assertEqual(len(coord), 2)

    def test_get_coord_cartesian_from_direct(self):
        coord = self.stru.get_coord(bohr=False, direct=False)
        self.assertEqual(len(coord), 2)


if __name__ == "__main__":
    unittest.main()
