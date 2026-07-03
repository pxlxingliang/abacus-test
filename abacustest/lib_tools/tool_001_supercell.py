from abacustest.lib_tools.tool import Tool
from abacustest.lib_prepare.stru import AbacusSTRU
import os


class SuperCellTool(Tool):
    @staticmethod
    def add_args(parser):
        parser.description = "Extend the unit cell to supercell. Supports ABACUS STRU, VASP POSCAR, and CIF formats."
        parser.add_argument("sc", help="the supercell size in a, b, c directions", nargs=3, type=int)
        parser.add_argument("-i", "--input", default=None, help="the input structure file.")
        parser.add_argument("-f", "--from-format", dest="from_fmt", default=None,
                            choices=["stru", "poscar", "vasp", "cif"],
                            help="Input format (auto-detect from filename if omitted)")
        parser.add_argument("-o", "--output", type=str, default=None,
                            help="the output structure file, default is {input}_{a}_{b}_{c}")
        parser.add_argument("-t", "--to-format", dest="to_fmt", default=None,
                            choices=["stru", "poscar", "vasp", "cif"],
                            help="Output format (auto-detect from filename if omitted)")
        parser.add_argument("--direct", action="store_true", default=None,
                            help="Write atomic positions in direct coordinates")
        return parser

    def run(self, params):
        if any(s <= 0 for s in params.sc):
            print("Error: supercell size must be positive integers.")
            return

        if params.input is None:
            print("Error: input file is required.")
            return

        if not os.path.isfile(params.input):
            print(f"Error: input file {params.input} does not exist.")
            return

        input_file = params.input
        stru = AbacusSTRU.read(input_file, fmt=params.from_fmt)
        if stru is None:
            return

        print("Original structure:")
        print("Atom numbers:", stru.natoms)
        print("Lattice parameter (Angstrom/Degree):\n", "%.2f %.2f %.2f %.1f %.1f %.1f" % tuple(stru.get_cell_param()))

        if params.output is None:
            output_file = f"{input_file}_{params.sc[0]}_{params.sc[1]}_{params.sc[2]}"
        else:
            output_file = params.output

        stru_super = stru.supercell(params.sc)
        direct = params.direct if params.direct is not None else None
        stru_super.sort()
        if not stru_super.write(output_file, fmt=params.to_fmt, direct=direct):
            return

        print(f"\nSupercell structure written to {output_file}")
        print("Supercell structure:")
        print("Atom numbers:", stru_super.natoms)
        print("Lattice parameter (Angstrom/Degree):\n", "%.2f %.2f %.2f %.1f %.1f %.1f" % tuple(stru_super.get_cell_param()))
