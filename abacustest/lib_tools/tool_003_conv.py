from abacustest.lib_tools.tool import Tool
from abacustest.lib_prepare.stru import AbacusSTRU, _guess_format
from abacustest.lib_prepare.comm import collect_pp
from abacustest.lib_prepare.abacus2vasp import gen_potcar
import os, shutil


class StructureConvertTool(Tool):
    @staticmethod
    def add_args(parser):
        parser.description = "Convert structure files between different formats (POSCAR, CIF, ABACUS STRU)."
        parser.add_argument("-i", "--input", required=True, help="Input structure file")
        parser.add_argument("-f", "--from", dest="from_fmt", default=None,
                            choices=["stru", "poscar", "vasp", "cif"],
                            help="Input format (auto-detect from filename if omitted)")
        parser.add_argument("-o", "--output", required=True, help="Output structure file")
        parser.add_argument("-t", "--to", dest="to_fmt", default=None,
                            choices=["stru", "poscar", "vasp", "cif"],
                            help="Output format (auto-detect from filename if omitted)")
        parser.add_argument("--pp", default=None, type=str,
                            help="Path to pseudopotential library (or use ABACUS_PP_PATH env). Only used when output is STRU.")
        parser.add_argument("--orb", default=None, type=str,
                            help="Path to orbital library (or use ABACUS_ORB_PATH env). Only used when output is STRU.")
        parser.add_argument("--potcar", default=None, type=str,
                            help="Path to VASP POTCAR library (or use VASP_POTCAR_PATH env). Only used when output is POSCAR.")
        parser.add_argument("--direct", action="store_true", default=None,
                            help="Write atomic positions in direct coordinates")
        parser.add_argument("--pporb-type", default=1, type=int, choices=[1, 2, 3],
                            help="How to handle pp/orb files: 1 = set absolute path in STRU file (default), "
                                 "2 = symlink pp/orb files to output dir, 3 = copy pp/orb files to output dir")
        return parser

    @staticmethod
    def _setup_pporb(stru, lib, out_dir, pporb_type, kind="pp"):
        """Handle pp/orb files according to pporb_type and return a dict mapping element to the value to write into STRU.

        Args:
            stru: AbacusSTRU object.
            lib (dict): Element -> source file path, returned by collect_pp.
            out_dir (str): Output directory where the STRU file is located.
            pporb_type (int): 1 = set absolute path in STRU file, 2 = symlink to out_dir, 3 = copy to out_dir.
            kind (str): "pp" or "orb".

        Returns:
            dict: Element -> pp/orb value to write into STRU, or {} if nothing found.
        """
        what = "pseudopotential" if kind == "pp" else "orbital"
        names = []
        elements = list(dict.fromkeys(stru.elements))
        for el in elements:
            if el in lib:
                src = lib[el]
                if pporb_type == 1:
                    value = os.path.abspath(src)
                else:
                    dst = os.path.join(out_dir, os.path.basename(src))
                    if pporb_type == 2:
                        if os.path.exists(dst):
                            os.remove(dst)
                        os.symlink(os.path.abspath(src), dst)
                    else:
                        shutil.copy2(src, dst)
                    value = os.path.basename(src)
                names.append(value)
            else:
                print(f"Warning: no {what} found for element '{el}'")
                names.append(None)
        d = dict(zip(elements, names))
        return {k: v for k, v in d.items() if v is not None}

    def run(self, params):
        in_file = params.input
        out_file = params.output

        if not os.path.isfile(in_file):
            print(f"Error: input file '{in_file}' does not exist.")
            return

        stru = AbacusSTRU.read(in_file, fmt=params.from_fmt)
        if stru is None:
            return

        in_fmt = params.from_fmt or _guess_format(in_file)
        print(f"Read structure from '{in_file}' ({in_fmt}): {len(stru)} atoms")

        out_dir = os.path.dirname(os.path.abspath(out_file))
        os.makedirs(out_dir, exist_ok=True)

        out_fmt = params.to_fmt or _guess_format(out_file)
        if out_fmt == "stru":
            pp_path = params.pp
            orb_path = params.orb
            if pp_path is None:
                pp_path = os.environ.get("ABACUS_PP_PATH", None)
            if orb_path is None:
                orb_path = os.environ.get("ABACUS_ORB_PATH", None)

            if pp_path is not None:
                pp_lib = collect_pp(pp_path)
                pp_dict = self._setup_pporb(stru, pp_lib, out_dir, params.pporb_type, kind="pp")
                if pp_dict:
                    stru.set_pp(pp_dict, key_type="element")

            if orb_path is not None:
                orb_lib = collect_pp(orb_path)
                orb_dict = self._setup_pporb(stru, orb_lib, out_dir, params.pporb_type, kind="orb")
                if orb_dict:
                    stru.set_orb(orb_dict, key_type="element")

        direct = params.direct if params.direct is not None else None
        stru.sort()
        if not stru.write(out_file, fmt=params.to_fmt, direct=direct):
            return

        print(f"Wrote structure to '{out_file}': {len(stru)} atoms")

        if out_fmt == "poscar":
            potcar_path = params.potcar
            if potcar_path is None:
                potcar_path = os.environ.get("VASP_POTCAR_PATH", None)
            if potcar_path is not None:
                elements = list(dict.fromkeys(stru.elements))
                gen_potcar(potcar_path, elements, os.path.join(out_dir, "POTCAR"))
