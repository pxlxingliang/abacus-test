from abacustest.lib_tools.tool import Tool
from abacustest.lib_prepare.stru import AbacusSTRU
from abacustest.lib_prepare.comm import collect_pp
import os, shutil


_SUPPORTED_FORMATS = {
    "stru": ("stru", "abacus/stru"),
    "poscar": ("poscar", "vasp"),
    "vasp": ("poscar", "vasp"),
    "cif": ("cif",),
}


def _guess_format(path: str) -> str:
    name = os.path.basename(path)
    if name == "POSCAR" or name.endswith(".vasp") or name.endswith(".poscar"):
        return "poscar"
    if name in ["STRU", "STRU_ION_D"] or name.endswith(".stru"):
        return "stru"
    if name.endswith(".cif"):
        return "cif"
    return None


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
        parser.add_argument("--direct", action="store_true", default=None,
                            help="Write atomic positions in direct coordinates")
        parser.add_argument("--copy-pp-orb", action="store_true",
                            help="Copy pp/orb files instead of symlinking them")
        return parser

    def run(self, params):
        in_file = params.input
        out_file = params.output

        if not os.path.isfile(in_file):
            print(f"Error: input file '{in_file}' does not exist.")
            return

        in_fmt = params.from_fmt or _guess_format(in_file)
        out_fmt = params.to_fmt or _guess_format(out_file)

        if in_fmt is None:
            print(f"Error: cannot detect input format from '{in_file}'. Please use --from.")
            return
        if out_fmt is None:
            print(f"Error: cannot detect output format from '{out_file}'. Please use --to.")
            return

        in_fmt = _SUPPORTED_FORMATS[in_fmt][0]
        out_fmt = _SUPPORTED_FORMATS[out_fmt][0]

        try:
            stru = AbacusSTRU.read(in_file, fmt=in_fmt)
        except Exception as e:
            print(f"Error reading '{in_file}' as {in_fmt}: {e}")
            return

        print(f"Read structure from '{in_file}' ({in_fmt}): {len(stru)} atoms")

        if out_fmt == "stru":
            pp_path = params.pp
            orb_path = params.orb
            if pp_path is None:
                pp_path = os.environ.get("ABACUS_PP_PATH", None)
            if orb_path is None:
                orb_path = os.environ.get("ABACUS_ORB_PATH", None)

            out_dir = os.path.dirname(os.path.abspath(out_file))
            os.makedirs(out_dir, exist_ok=True)

            if pp_path is not None:
                pp_lib = collect_pp(pp_path)
                pp_names = []
                elements = list(dict.fromkeys(stru.elements))
                for el in elements:
                    if el in pp_lib:
                        src = pp_lib[el]
                        dst = os.path.join(out_dir, os.path.basename(src))
                        if params.copy_pp_orb:
                            shutil.copy2(src, dst)
                        else:
                            if os.path.exists(dst):
                                os.remove(dst)
                            os.symlink(os.path.abspath(src), dst)
                        pp_names.append(os.path.basename(src))
                    else:
                        print(f"Warning: no pseudopotential found for element '{el}'")
                        pp_names.append(None)
                pp_dict = dict(zip(elements, pp_names))
                pp_dict = {k: v for k, v in pp_dict.items() if v is not None}
                if pp_dict:
                    stru.set_pp(pp_dict, key_type="element")

            if orb_path is not None:
                orb_lib = collect_pp(orb_path)
                orb_names = []
                elements = list(dict.fromkeys(stru.elements))
                for el in elements:
                    if el in orb_lib:
                        src = orb_lib[el]
                        dst = os.path.join(out_dir, os.path.basename(src))
                        if params.copy_pp_orb:
                            shutil.copy2(src, dst)
                        else:
                            if os.path.exists(dst):
                                os.remove(dst)
                            os.symlink(os.path.abspath(src), dst)
                        orb_names.append(os.path.basename(src))
                    else:
                        print(f"Warning: no orbital found for element '{el}'")
                        orb_names.append(None)
                orb_dict = dict(zip(elements, orb_names))
                orb_dict = {k: v for k, v in orb_dict.items() if v is not None}
                if orb_dict:
                    stru.set_orb(orb_dict, key_type="element")

        direct = params.direct if params.direct is not None else None
        try:
            stru.write(out_file, fmt=out_fmt, direct=direct)
        except Exception as e:
            print(f"Error writing to '{out_file}' as {out_fmt}: {e}")
            return

        print(f"Wrote structure to '{out_file}' ({out_fmt}): {len(stru)} atoms")
