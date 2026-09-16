"""
conv-out: Convert DFT software output files to other formats for visualization.

Currently supported:
  - Input:  abacus (ABACUS output)
  - Output: gview  (Gaussian View format)

Users should implement the actual conversion logic in the `run` method.
"""
from abacustest.lib_tools.tool import Tool
import os


_SUPPORTED_INPUT_FORMATS = ["abacus"]
_SUPPORTED_OUTPUT_FORMATS = ["gview", "xyz", "gro", "deepmd"]


class ConvOutTool(Tool):
    @staticmethod
    def add_args(parser):
        parser.description = (
            "Convert DFT software output to other formats for visualization. "
            "Currently supports ABACUS -> xyz conversion."
        )
        parser.add_argument(
            "-i", "--input", default=".",
            help="Input DFT output file or directory (e.g. ABACUS job path, default: current directory)"
        )
        parser.add_argument(
            "-f", "--from", dest="from_fmt", default="abacus",
            choices=_SUPPORTED_INPUT_FORMATS,
            help="Input format (default: abacus)"
        )
        parser.add_argument(
            "-o", "--output", required=True,
            help="Output file path"
        )
        parser.add_argument(
            "-t", "--to", dest="to_fmt", default="xyz",
            choices=_SUPPORTED_OUTPUT_FORMATS,
            help="Output format (default: xyz)"
        )
        return parser

    def run(self, params):
        in_file = params.input
        out_file = params.output
        in_fmt = params.from_fmt
        out_fmt = params.to_fmt

        import os
        if not os.path.exists(in_file):
            print(f"Error: input '{in_file}' does not exist.")
            return

        print(f"Input format : {in_fmt}")
        print(f"Output format: {out_fmt}")
        print(f"Input file   : {in_file}")
        print(f"Output file  : {out_file}")

        traj = None
        if in_fmt in ["abacus"]:
            traj = abacus_out2ase(in_file)
        else:
            raise NotImplementedError(f"Unsupported input format: {in_fmt}")

        if traj is None:
            print("Error: Failed to convert input file to ASE Atoms object.")
            return

        if out_fmt in ["xyz", "gro"]:
            from ase.io import write
            if out_fmt == "gro":
                write(out_file + ".gro", traj[0])
                #write(out_file + ".xtc", traj, format="xtc")
                write(out_file + ".dcd", traj, format="dcd")
            else:
                write(out_file, traj)

        else:
            raise NotImplementedError(f"Unsupported output format: {out_fmt}")


def abacus_out2ase(job_path):
    """
    Convert ABACUS output to ASE Trajectory format (ASE Atoms object).

    Parameters:
        job_path (str): Path to the ABACUS job directory containing output files.

    Returns:
        traj (list of ase.Atoms): List of ASE Atoms objects representing each step in the ABACUS calculation.
    """
    from abacustest import RESULT, ReadInput
    from ase import Atoms
    import numpy as np

    traj = []

    # Load the ABACUS job
    if os.path.isfile(job_path):
        if os.path.basename(job_path) in ["MD_dump"]:
            traj = parse_abacus_md_dump(job_path)
    elif os.path.isdir(job_path):
        input_param = ReadInput(os.path.join(job_path, "INPUT"))
        suffix = input_param.get("suffix", "ABACUS")
        if os.path.isfile(os.path.join(job_path, f"OUT.{suffix}", "MD_dump")):
            traj = parse_abacus_md_dump(os.path.join(job_path, f"OUT.{suffix}", "MD_dump"))
        else:
            result = RESULT(path=job_path, fmt="abacus")
            cells = result["cells"] # 
            coords = result["coordinates"]
            energies = result["energies"]
            forces = result["forces"]
            element = result["element"]

            # Create ASE Atoms object

            for i in range(len(cells)):
                cell = cells[i]
                pos = coords[i]
                energy = energies[i] if energies is not None else None
                force = forces[i] if forces is not None else None

                atoms = Atoms(
                    symbols=element,
                    positions=pos,
                    cell=cell,
                    pbc=True
                )
                if energy is not None:
                    atoms.info["energy"] = energy
                if force is not None:
                    atoms.arrays["forces"] = np.array(force).reshape(-1, 3)

                traj.append(atoms)
    
    return traj

def parse_abacus_md_dump(dump_file):
    """
    Parse the ABACUS MD_dump file and extract trajectory information.

    Parameters:
        dump_file (str): Path to the ABACUS MD_dump file.    
    
    Returns:
        traj (list of ase.Atoms): List of ASE Atoms objects representing each step in the ABACUS MD simulation.
    
    The MD_dump file is in format:
    MDSTEP:  0
LATTICE_CONSTANT: 0.529177000000 Angstrom
LATTICE_VECTORS
  26.512340425020  0.000000000000  0.000000000000
  0.000000000000  27.343803702790  0.000000000000
  0.000000000000  0.000000000000  27.419391273490
VIRIAL (kbar)
  212.157232346048  11.626082274312  -7.363588411683
  11.626082274312  190.382227993143  -7.806314800882
  -7.363588411683  -7.806314800882  194.941887052843
INDEX    LABEL    POSITION (Angstrom)    FORCE (eV/Angstrom)    VELOCITY (Angstrom/fs)
  0  Ti  6.386872470159  1.488970269005  9.731805680239  -1.224281570306  -1.016246093779  1.264566116164  0.002021379847  -0.000548275061  0.001284986167
  1  Ti  9.542809454027  7.129857634609  1.302973982869  3.422277570518  -0.022911406549  0.710183129739  0.000911249788  -0.001412225224  0.002583162709
  2  Ti  9.745805400701  1.195976119344  1.999960065801  0.520751851098  -1.133688688064  -1.062312771022  -0.004647938659  0.001248579237  -0.001089754698  
 ...
 

MDSTEP:  1
LATTICE_CONSTANT: 0.529177000000 Angstrom
...       
    """
    from ase import Atoms
    import numpy as np

    if not os.path.isfile(dump_file):
        raise FileNotFoundError(f"MD_dump file '{dump_file}' does not exist.")
    
    traj = []
    with open(dump_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("MDSTEP:"):
            # Start of a new MD step
            step_index = int(line.split()[1])
            lc = None
            cell = None
            virial = None
            positions = []
            forces = []
            velocities = []
            element = []
            virial = []

            # Read until the next MDSTEP or end of file
            j = i + 1
            while j < len(lines) and not lines[j].strip().startswith("MDSTEP:"):
                if lines[j].startswith("LATTICE_VECTORS"):
                    # Read lattice vectors
                    cell = []
                    for k in range(3):
                        j += 1
                        cell.append(list(map(float, lines[j].split())))
                    cell = np.array(cell)
                elif lines[j].startswith("LATTICE_CONSTANT"):
                    # Read lattice constant
                    lc = float(lines[j].split()[1])
                    j += 1
                elif lines[j].startswith("INDEX"):
                    # Read atom data
                    j += 1
                    while j < len(lines) and lines[j].strip():
                        parts = lines[j].split()
                        index = int(parts[0])
                        label = parts[1]
                        pos = list(map(float, parts[2:5]))
                        force = list(map(float, parts[5:8]))
                        vel = list(map(float, parts[8:11]))

                        element.append(label)
                        positions.append(pos)
                        forces.append(force)
                        velocities.append(vel)

                        j += 1
                elif lines[j].startswith("VIRIAL"):
                    # Read virial stress tensor
                    j += 1
                    virial = []
                    for k in range(3):
                        virial.append(list(map(float, lines[j].split())))
                        j += 1
                else:
                    j += 1
            
            # Create ASE Atoms object for this step
            if cell is not None:
                cell *= lc  # Scale cell by lattice constant
            atoms = Atoms(
                symbols=element,
                positions=positions,
                cell=cell,
                pbc=True
            )
            atoms.arrays["forces"] = np.array(forces)
            atoms.arrays["velocities"] = np.array(velocities)
            #atoms.arrays["virial"] = np.array(virial)
            traj.append(atoms)
            i = j  # Move to the next MDSTEP or end of file
        else:
            i += 1  # Skip lines until the next MDSTEP

    return traj

def unramp_traj(traj):
    """
    Unramp the trajectory to remove periodic boundary conditions.

    Parameters:
        traj (list of ase.Atoms): List of ASE Atoms objects representing the trajectory.

    Returns:
        unramped_traj (list of ase.Atoms): List of ASE Atoms objects with unramped positions.
    """
    from ase.geometry import wrap_positions
    unramped_traj = []
    for atoms in traj:
        positions = atoms.get_positions()
        cell = atoms.get_cell()
        unramped_positions = wrap_positions(positions, cell)
        atoms.set_positions(unramped_positions)
        unramped_traj.append(atoms)
    return unramped_traj    