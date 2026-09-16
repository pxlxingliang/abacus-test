import unittest
import os
import shutil
import tempfile
import argparse

from abacustest.lib_prepare.stru import AbacusATOM, AbacusSTRU
from abacustest.lib_tools.tool_006_bsse import BsseTool, _parse_indices

BSSE_DIRS = ["BSSE_AB", "BSSE_A_ghostB", "BSSE_B_ghostA", "BSSE_A", "BSSE_B"]


class TestParseIndices(unittest.TestCase):
    def test_simple_list(self):
        self.assertEqual(_parse_indices("1,2,3"), [1, 2, 3])

    def test_ranges(self):
        self.assertEqual(_parse_indices("1,3-5"), [1, 3, 4, 5])

    def test_spaces(self):
        self.assertEqual(_parse_indices("1, 3 - 5"), [1, 3, 4, 5])

    def test_dedupe(self):
        self.assertEqual(_parse_indices("2,2,3-4"), [2, 3, 4])

    def test_invalid_non_positive(self):
        with self.assertRaises(ValueError):
            _parse_indices("0,1")
        with self.assertRaises(ValueError):
            _parse_indices("-1,2")
        with self.assertRaises(ValueError):
            _parse_indices("5-2")

    def test_invalid_syntax(self):
        with self.assertRaises(ValueError):
            _parse_indices("a,b")


class TestBsseTool(unittest.TestCase):
    def setUp(self):
        self.work_dir = tempfile.mkdtemp()
        cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        atoms = [
            AbacusATOM(label="H", coord=(0.0, 0.0, 0.0)),
            AbacusATOM(label="H", coord=(0.7, 0.0, 0.0)),
            AbacusATOM(label="O", coord=(1.2, 0.0, 0.0)),
            AbacusATOM(label="C", coord=(2.0, 0.0, 0.0)),
        ]
        self.stru = AbacusSTRU(cell=cell, atoms=atoms)
        self.input_file = os.path.join(self.work_dir, "STRU")
        self.stru.write(self.input_file, fmt="stru")
        self.input_file_src = os.path.join(self.work_dir, "INPUT")
        with open(self.input_file_src, "w") as f:
            f.write("suffix        test\n")
            f.write("ecutwfc       100\n")

    def tearDown(self):
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)

    def _make_params(self, frag1=None, frag2=None, output_dir=None,
                     stru_input=None, input_file_src=None):
        return argparse.Namespace(
            input=stru_input if stru_input is not None else self.input_file,
            output_dir=output_dir if output_dir is not None else self.work_dir,
            input_file=input_file_src,
            frag1=frag1,
            frag2=frag2,
        )

    def _read_labels(self, path):
        return AbacusSTRU.read(path, fmt="stru").labels

    def _assert_5_dirs(self, out_dir):
        for name in BSSE_DIRS:
            self.assertTrue(
                os.path.isfile(os.path.join(out_dir, name, "STRU")),
                f"Expected {name}/STRU in {out_dir}",
            )

    def test_frag1_only(self):
        tool = BsseTool()
        tool.run(self._make_params(frag1="1,2"))
        out = self.work_dir
        self._assert_5_dirs(out)

        # A = H(1,2), B = O(3), C(4)
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_AB", "STRU")),
                         ["H", "H", "O", "C"])
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_A_ghostB", "STRU")),
                         ["H", "H", "O_empty", "C_empty"])
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_B_ghostA", "STRU")),
                         ["H_empty", "H_empty", "O", "C"])
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_A", "STRU")),
                         ["H", "H"])
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_B", "STRU")),
                         ["O", "C"])

    def test_input_file_copied(self):
        tool = BsseTool()
        tool.run(self._make_params(frag1="1,2", input_file_src=self.input_file_src))
        out = self.work_dir
        self._assert_5_dirs(out)
        for name in BSSE_DIRS:
            self.assertTrue(
                os.path.isfile(os.path.join(out, name, "INPUT")),
                f"Expected INPUT in {name}",
            )
            with open(os.path.join(out, name, "INPUT")) as f:
                self.assertEqual(f.read(), "suffix        test\necutwfc       100\n")

    def test_input_file_detected_in_cwd(self):
        # Create a separate directory with an INPUT file and run from there
        cwd_dir = os.path.join(self.work_dir, "cwd")
        os.makedirs(cwd_dir)
        cwd_input = os.path.join(cwd_dir, "INPUT")
        with open(cwd_input, "w") as f:
            f.write("suffix        cwd\n")

        old_cwd = os.getcwd()
        try:
            os.chdir(cwd_dir)
            tool = BsseTool()
            tool.run(self._make_params(frag1="1,2"))
        finally:
            os.chdir(old_cwd)

        out = self.work_dir
        self._assert_5_dirs(out)
        for name in BSSE_DIRS:
            self.assertTrue(os.path.isfile(os.path.join(out, name, "INPUT")))
            with open(os.path.join(out, name, "INPUT")) as f:
                self.assertEqual(f.read(), "suffix        cwd\n")

    def test_no_input_copied_when_missing(self):
        tool = BsseTool()
        # input_file=None, and run from a directory without an INPUT file
        cwd_dir = os.path.join(self.work_dir, "noinput_cwd")
        os.makedirs(cwd_dir)
        old_cwd = os.getcwd()
        try:
            os.chdir(cwd_dir)
            tool.run(self._make_params(frag1="1,2"))
        finally:
            os.chdir(old_cwd)

        out = self.work_dir
        self._assert_5_dirs(out)
        for name in BSSE_DIRS:
            self.assertFalse(os.path.exists(os.path.join(out, name, "INPUT")))

    def test_both_fragments(self):
        tool = BsseTool()
        tool.run(self._make_params(frag1="1,2", frag2="3-4"))
        out = self.work_dir
        self._assert_5_dirs(out)
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_A_ghostB", "STRU")),
                         ["H", "H", "O_empty", "C_empty"])
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_B_ghostA", "STRU")),
                         ["H_empty", "H_empty", "O", "C"])

    def test_interleaved_fragments_grouped(self):
        # frag1 = atom2 (H) + atom3 (O), frag2 = atom1 (H) + atom4 (C)
        tool = BsseTool()
        tool.run(self._make_params(frag1="2,3"))
        out = self.work_dir
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_A_ghostB", "STRU")),
                         ["H_empty", "H", "O", "C_empty"])
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_B_ghostA", "STRU")),
                         ["H", "H_empty", "O_empty", "C"])
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_A", "STRU")),
                         ["H", "O"])
        self.assertEqual(self._read_labels(os.path.join(out, "BSSE_B", "STRU")),
                         ["H", "C"])

    def test_custom_output_dir(self):
        out_dir = os.path.join(self.work_dir, "bsse_out")
        tool = BsseTool()
        tool.run(self._make_params(frag1="1-2", output_dir=out_dir))
        self._assert_5_dirs(out_dir)
        self.assertFalse(os.path.exists(os.path.join(self.work_dir, "BSSE_A")))

    def test_pp_orb_preserved_for_empty_atoms(self):
        self.stru.set_pp({"H": "H.upf", "O": "O.upf", "C": "C.upf"})
        self.stru.set_orb({"H": "H.orb", "O": "O.orb", "C": "C.orb"})
        self.stru.write(self.input_file, fmt="stru")

        tool = BsseTool()
        tool.run(self._make_params(frag1="1-2"))
        stru_ghost = AbacusSTRU.read(os.path.join(self.work_dir, "BSSE_A_ghostB", "STRU"), fmt="stru")
        self.assertEqual(stru_ghost.pp_dict()["O_empty"], "O.upf")
        self.assertEqual(stru_ghost.pp_dict()["C_empty"], "C.upf")
        self.assertEqual(stru_ghost.orb_dict()["O_empty"], "O.orb")
        self.assertEqual(stru_ghost.orb_dict()["C_empty"], "C.orb")

    def test_no_fragment_error(self):
        out_dir = os.path.join(self.work_dir, "err_out")
        tool = BsseTool()
        tool.run(self._make_params(output_dir=out_dir))
        self.assertFalse(os.path.exists(os.path.join(out_dir, "BSSE_A")))

    def test_out_of_range_error(self):
        tool = BsseTool()
        tool.run(self._make_params(frag1="1,2,99"))
        self.assertFalse(os.path.exists(os.path.join(self.work_dir, "BSSE_A")))

    def test_overlap_error(self):
        tool = BsseTool()
        tool.run(self._make_params(frag1="1-3", frag2="3-4"))
        self.assertFalse(os.path.exists(os.path.join(self.work_dir, "BSSE_A")))

    def test_input_file_missing_error(self):
        tool = BsseTool()
        tool.run(self._make_params(frag1="1,2",
                                   input_file_src=os.path.join(self.work_dir, "no_such_INPUT")))
        self.assertFalse(os.path.exists(os.path.join(self.work_dir, "BSSE_A")))


if __name__ == "__main__":
    unittest.main()
