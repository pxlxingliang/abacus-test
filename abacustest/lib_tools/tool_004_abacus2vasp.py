from abacustest.lib_tools.tool import Tool
from abacustest.lib_prepare.abacus2vasp import Abacus2Vasp
import os


class Abacus2VaspTool(Tool):
    @staticmethod
    def add_args(parser):
        parser.description = "Transform ABACUS input files to VASP input files."
        parser.add_argument("-i", "--input", required=True,
                            help="Path to ABACUS working directory (containing INPUT, STRU, KPT)")
        parser.add_argument("-o", "--output", default=None,
                            help="Path to save VASP files (default: {input}.vasp)")
        parser.add_argument("--potcar", default=None, type=str,
                            help="Path to VASP POTCAR library, or read from VASP_POTCAR_PATH env variable")
        parser.add_argument("--set", default=[], action="append", nargs=2, metavar=("KEY", "VALUE"),
                            help="Additional VASP INCAR settings, e.g. --set ENCUT 500 --set EDIFF 1e-6")
        return parser

    def run(self, params):
        abacus_path = params.input
        if not os.path.isdir(abacus_path):
            print(f"Error: ABACUS path '{abacus_path}' does not exist or is not a directory.")
            return

        potcar = params.potcar
        if potcar is None:
            potcar = os.environ.get("VASP_POTCAR_PATH", None)

        vasp_setting = {}
        for k, v in params.set:
            vasp_setting[k] = v

        save_path = params.output or (abacus_path.rstrip("/") + ".vasp")

        print(f"Converting ABACUS files in '{abacus_path}' to VASP format...")
        os.makedirs(save_path, exist_ok=True)

        result = Abacus2Vasp(abacus_path, save_path=save_path, potcar=potcar, vasp_setting=vasp_setting)
        print(f"VASP input files written to '{save_path}'")
        return result
