from ..model import Model
from abacustest.lib_prepare.abacus import AbacusStru
import os

class SuperCellModel(Model):
    HAS_PREPARE_POST_COMMAND = False

    @staticmethod
    def model_name(): # type: ignore
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "supercell"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "extend the unit cell to supercell"
    
    @staticmethod
    def add_args(parser):
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand '''
        
        parser.description = "Prepare the ABACUS inputs file."
        parser.add_argument("sc", help="the supercell size in a, b, c directions", nargs=3, type=int)
        parser.add_argument('-i', '--input',default="STRU",help='the structure file.')
        parser.add_argument("-o", "--output", type=str, default=None, help="the output STRU file name, default is {i}_{a}_{b}_{c}, where {i} is input file, and {a}, {b}, {c} is the supercell size")
        return parser
    
    def run(self,params):
        '''
        Run the prepare subcommand with the given arguments
        The arguments are stored in params
        '''
        if any(s <= 0 for s in params.sc):
            print("Error: supercell size must be positive integers.")
            return
        
        if not os.path.isfile(params.input):
            print(f"Error: input file {params.input} does not exist.")
            return
        
        input_file = params.input
        stru = AbacusStru.ReadStru(input_file)
        if stru is None:
            print(f"Error: cannot read structure from {input_file}")
            return 
        print("Reading structure from", input_file)

        print("Original structure:")
        print("Atom numbers:", stru.get_natoms())
        print("Lattice parameter (Angstrom/Degree):\n", "%.2f %.2f %.2f %.1f %.1f %.1f" % tuple(stru.get_cell_param()))
        
        stru_super = stru.supercell(params.sc)
        if params.output is None:
            output_file = f"{input_file}_{params.sc[0]}_{params.sc[1]}_{params.sc[2]}"
        else:
            output_file = params.output
        stru_super.write(output_file)
        print(f"\nSupercell structure written to {output_file}")
        print("Supercell structure:")
        print("Atom numbers:", stru_super.get_natoms())
        print("Lattice parameter (Angstrom/Degree):\n", "%.2f %.2f %.2f %.1f %.1f %.1f" % tuple(stru_super.get_cell_param()))