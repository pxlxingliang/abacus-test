from ..model import Model
import json, os, shutil
from pathlib import Path
from typing import List, Dict, Any, Literal, Optional, Tuple

import numpy as np

from abacustest import AbacusSTRU, WriteInput, ReadInput, RESULT
from abacustest.constant import (
    RECOMMAND_IMAGE,
    RECOMMAND_COMMAND,
    RECOMMAND_MACHINE,
    RECOMMAND_VASP_IMAGE,
    RECOMMAND_VASP_COMMAND,
    RECOMMAND_VASP_MACHINE,
)
from abacustest.lib_model.comm import copy_abacusjob, get_largest_vacuum_dir
from abacustest.lib_data.grid import Potential
from abacustest.lib_prepare.abacus2vasp import Abacus2Vasp
from abacustest.lib_collectdata.collectdata import RESULT as CollectDataRESULT

EFIELD_DIRECTION_MAP = ["a", "b", "c"]


class WorkFuncModel(Model):
    @staticmethod
    def model_name():  # type: ignore
        """
        Name of the model, which will be used as the subcommand
        """
        return "workfunc"

    @staticmethod
    def description():  # type: ignore
        """
        Description of the model
        """
        return "Prepare and postprocess the work function calculation"

    @staticmethod
    def prepare_args(parser):
        """
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand
        """
        parser.description = "Prepare the inputs for the work function calculation."
        parser.add_argument(
            "-j",
            "--job",
            default=[],
            action="extend",
            nargs="*",
            help="the paths of ABACUS jobs, should contain INPUT, STRU, or KPT, and pseudopotential and orbital files",
        )
        parser.add_argument(
            "--vacuum",
            default="auto",
            choices=["a", "b", "c", "auto"],
            help="The direction of the vacuum for work function calculation, default is auto. If set to auto, the direction will try to be determined automatically.",
        )
        parser.add_argument(
            "--dipole-corr",
            action="store_true",
            help="Whether to use dipole correction for the work function calculation, default is not using.",
        )
        parser.add_argument(
            "--image",
            type=str,
            default=None,
            help="The image to use for the Bohrium job. Default is none, which means determine by used DFT software."
        )
        parser.add_argument(
            "--machine",
            type=str,
            default=None,
            help="The machine to use for the Bohrium job. Default is none, which means determine by used DFT software.",
        )
        parser.add_argument(
            "--dft-command",
            type=str,
            default=None,
            help=f"The command to run the dft job. Default is none, which means determine by used DFT software.",
        )
        parser.add_argument(
            "--dft",
            default="abacus",
            choices=["abacus", "vasp"],
            help="DFT software to use for work function calculation (ABACUS or VASP). Default is ABACUS.",
        )
        parser.add_argument(
            "--potcar_path",
            type=str,
            default=None,
            help="Path to VASP POTCAR files directory (required for VASP input generation)",
        )
        return parser

    def run_prepare(self, params):
        """
        Run the model with the given parameters
        """
        if not params.job:
            raise ValueError("No job specified, please use -j or --job to specify the job paths.")

        paths = prep_all_workfunc_jobs(params.job, params.dft, params.vacuum, params.dipole_corr, params.potcar_path)
        
        if params.image is None:
            image = RECOMMAND_IMAGE if params.dft == "abacus" else RECOMMAND_VASP_IMAGE if params.dft == "vasp" else None
        else:
            image = params.image
        if params.machine is None:
            machine = RECOMMAND_MACHINE if params.dft == "abacus" else RECOMMAND_VASP_MACHINE if params.dft == "vasp" else None
        else:
            machine = params.machine
        if params.dft_command is None:
            dft_command = RECOMMAND_COMMAND if params.dft == "abacus" else RECOMMAND_VASP_COMMAND if params.dft == "vasp" else None
        else:
            dft_command = params.dft_command

        # Create setting file
        if paths:
            setting = {"save_path": "results", "run_dft": []}
            setting["run_dft"].append(
                {
                    "ifrun": True,
                    "example": paths,
                    "command": dft_command,
                    "image": image,
                    "bohrium": {
                        "scass_type": machine,
                        "job_type": "container",
                        "platform": "ali",
                    },
                }
            )

        setting_file = "setting.json"
        json.dump(setting, open(setting_file, "w"), indent=4)
        print(f"\nInputs are generated in: {', '.join([str(p) for p in paths])}")
        print(f"You can modify '{setting_file}', and execute below command to run the abacustest to submit all jobs to bohrium:\n\tabacustest submit -p {setting_file}\n")

        # Print postprocess instructions
        print(f"\nAfter finishing the calculations, you can postprocess the work function results:")
        if params.dft == "abacus":
            print(f"  For ABACUS results: abacustest model workfunc post -j {' '.join(params.job)} --dft abacus")
        elif params.dft == "vasp":
            print(f"  For VASP results: abacustest model workfunc post -j {' '.join(params.job)} --dft vasp")

    @staticmethod
    def postprocess_args(parser):
        """
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand
        """
        parser.description = "Postprocess the work function calculation results."
        parser.add_argument(
            "-j",
            "--job",
            default=[],
            action="extend",
            nargs="*",
            help="the paths of the job directories, should contain the results of work function calculations generated by the prepare step.",
        )
        parser.add_argument(
            "--dft",
            default="abacus",
            choices=["abacus", "vasp"],
            help="DFT software used for work function calculation (ABACUS or VASP)",
        )
        return parser

    def run_postprocess(self, params):
        """
        Parse the parameters and run the postprocess process"""
        if not params.job:
            raise ValueError(
                "No job specified, please use -j or --job to specify the job paths."
            )

        print("Work function results (units in eV):")

        results_all = {}
        for job in params.job:
            if not os.path.isdir(job):
                raise ValueError("The job path is not a directory: %s" % job)

            # Determine job type based on parameters
            jobtype = "vasp" if params.dft == "vasp" else "abacus"
            results, plot_path, pot_file, plot_data_file = post_workfunc_calc(
                job, jobtype=jobtype
            )

            work_functions = ""
            for result in results:
                work_functions += f"{result['work_function']:4f} eV, "
            print(f"{job}: {work_functions}")

            results_all[job] = results

        json.dump(results_all, open("metrics.json", "w"), indent=4)
        print("\nThe postprocess is done. The metrics are saved in 'metrics.json'.")


def prep_all_workfunc_jobs(jobs: List[str], 
                           mode: Literal['abacus', 'vasp'], 
                           vacuum_dir: Literal['a', 'b', 'c', 'auto'] = 'auto', 
                           dipole_corr: bool = False, 
                           potcar_path: Optional[str] = None):
    paths = []
    for job in jobs:
        if not os.path.isdir(job):
            raise ValueError("The job path is not a directory: %s" % job)

        # Prepare ABACUS calculation if ABACUS is selected
        if mode == "abacus":
            abacus_path = prep_abacus_workfunc_calc(job, vacuum_dir, dipole_corr, os.path.join(job, "workfunc_job"))
            paths.append(abacus_path)

        # Prepare VASP calculation if VASP is selected
        elif mode == "vasp":
            vasp_path = prep_vasp_workfunc_calc(job, vacuum_dir, dipole_corr, potcar_path)
            paths.append(vasp_path)
    
    return paths
        

def prep_abacus_workfunc_calc(
    job: Path,
    vacuum_dir: Literal["a", "b", "c", "auto"] = "auto",
    dipole_corr: bool = False,
    workfunc_dir: Optional[str] = "workfunc_job"
) -> Path:
    workfunc_dir = Path(workfunc_dir).absolute()
    if workfunc_dir in os.listdir(job):
        print("Old workfunc_job directory found, remove and recreate it.")
        shutil.rmtree(Path(job) / "workfunc_job")

    copy_abacusjob(job, workfunc_dir, out_dir=False)

    input_params = ReadInput(os.path.join(workfunc_dir, "INPUT"))
    if input_params.get("nspin", 1) not in [1, 2]:
        raise ValueError("Only 1 or 2 spin is supported.")

    input_params["calculation"] = "scf"
    input_params["out_pot"] = 2

    # Find vacuum direction in automatic mode
    if vacuum_dir == "auto":
        stru_file = os.path.join(workfunc_dir, input_params.get("stru_file", "STRU"))
        stru = AbacusSTRU.read(stru_file)
        direction, max_vacuum_cart = get_largest_vacuum_dir(coords=stru.coords, cell=stru.cell)
        prep_info = {'vacuum_dir': direction}
        print(f"Automatically identified the vacuum along {direction} direction in {stru_file}")
    else:
        prep_info = {"vacuum_dir": vacuum_dir}

    # Dump vacuum direction for postprocess
    with open(os.path.join(workfunc_dir, "prep_info.json"), "w") as f:
        json.dump(prep_info, f, indent=2)

    if dipole_corr:
        input_params["efield_flag"] = 1
        input_params["dip_cor_flag"] = 1
        if vacuum_dir in ["a", "b", "c"]:
            input_params["efield_dir"] = EFIELD_DIRECTION_MAP.index(vacuum_dir)
        elif vacuum_dir == "auto":
            input_params["efield_dir"] = EFIELD_DIRECTION_MAP.index(direction)
        else:
            raise ValueError("Invalid vacuum direction: %s" % vacuum_dir)
        input_params['efield_amp'] = 0.0 # Automatically determine electric field strength
        input_params['efield_pos_max'] = None # Automatically determine electric field position
        input_params['efield_pos_dec'] = None # Automatically determine range of electric field
    
    WriteInput(input_params, os.path.join(workfunc_dir, "INPUT"))

    return str(workfunc_dir)


def prep_vasp_workfunc_calc(
    job: Path,
    vacuum_dir: Literal["a", "b", "c", "auto"] = "auto",
    dipole_corr: bool = False,
    potcar_path: Optional[str] = None,
) -> Path:
    """
    Prepare VASP work function calculation from ABACUS input.

    Args:
        job: Path to ABACUS job directory
        vacuum_dir: Vacuum direction for work function calculation
        dipole_corr: Whether to use dipole correction
        potcar_path: Path to VASP POTCAR files directory

    Returns:
        Path to prepared VASP work function calculation directory
    """
    import uuid

    # Check and remove existing VASP directory
    vasp_dir_name = "workfunc_vasp_job"
    if vasp_dir_name in os.listdir(job):
        print(f"Old {vasp_dir_name} directory found, remove and recreate it.")
        shutil.rmtree(Path(job) / vasp_dir_name)

    # Check and remove existing ABACUS directory if it exists
    abacus_dir_name = "workfunc_job" + str(uuid.uuid4())[:8] # avoid confilct with existed ABACUS job
    if abacus_dir_name in os.listdir(job):
        print(f"Existing {abacus_dir_name} directory found, removing it.")
        shutil.rmtree(Path(job) / abacus_dir_name)

    # First create ABACUS workfunc job to get proper INPUT parameters
    abacus_workfunc_dir = prep_abacus_workfunc_calc(job, vacuum_dir, dipole_corr, abacus_dir_name)

    # Read prep_info to get vacuum direction
    with open(os.path.join(abacus_workfunc_dir, "prep_info.json"), "r") as f:
        prep_info = json.load(f)

    # Create VASP directory
    vasp_workfunc_dir = Path(os.path.join(job, vasp_dir_name)).absolute()
    os.makedirs(vasp_workfunc_dir, exist_ok=True)

    # Copy prep_info to VASP directory
    with open(os.path.join(vasp_workfunc_dir, "prep_info.json"), "w") as f:
        json.dump(prep_info, f, indent=2)

    # Convert ABACUS input to VASP input
    try:
        if potcar_path is None:
            potcar_path = os.environ.get("VASP_PP_PATH", None)
        if potcar_path is None:
            print("Warning: No POTCAR provided")
        Abacus2Vasp(
            abacus_path=str(abacus_workfunc_dir),
            save_path=str(vasp_workfunc_dir),
            potcar=potcar_path,
            vasp_setting={"emax_coef": 1.5},
        )
    except Exception as e:
        raise RuntimeError(f"Failed to convert ABACUS to VASP input: {e}")

    # Remove the temporary ABACUS directory
    shutil.rmtree(abacus_workfunc_dir)

    return str(vasp_workfunc_dir)


def post_workfunc_calc(
    job: Path, jobtype: Literal["abacus", "vasp"] = "abacus"
) -> Tuple[List[Dict[str, Any]], Path, str, str]:
    # Postprocess calculation to obtain work function data
    if jobtype == "abacus":
        workfunc_job = os.path.join(job, "workfunc_job")
        results = RESULT(fmt="abacus", path=workfunc_job)
        efermi = results["efermi"]

        if results["normal_end"] is not True or results["converge"] is not True:
            raise RuntimeError(
                "ABACUS calculation didn't end normally or didn't reached SCF converge."
            )

        input_params = ReadInput(os.path.join(workfunc_job, "INPUT"))
        pot_file = os.path.join(
            workfunc_job,
            f"OUT.{input_params.get('suffix', 'ABACUS')}/ElecStaticPot.cube",
        )
        pot = Potential.from_cube(pot_file)

        prep_info = json.load(open(os.path.join(workfunc_job, "prep_info.json")))
        vacuum_dir = prep_info["vacuum_dir"]
        ave_elec_stat, coord_direct = pot.profile1d(axis=vacuum_dir, average=True)
        ave_elec_stat = -ave_elec_stat  # consider negative charge of electron
    elif jobtype == "vasp":
        workfunc_job = os.path.join(job, "workfunc_vasp_job")

        # Use lib_collectdata to get VASP calculation results
        try:
            result = CollectDataRESULT("vasp", path=workfunc_job)
            # Accessing keys will automatically call the registered methods
            efermi = result["efermi"]

            if efermi is None:
                # Check if OUTCAR exists
                outcar_file = os.path.join(workfunc_job, "OUTCAR")
                if not os.path.exists(outcar_file):
                    raise RuntimeError(
                        f"VASP OUTCAR file not found in {workfunc_job}. Please run the VASP calculation first."
                    )
                else:
                    raise RuntimeError(
                        "Could not find Fermi energy in VASP results. Calculation may have failed."
                    )

        except Exception as e:
            raise RuntimeError(f"Failed to collect VASP results: {str(e)}")

        # Read LOCPOT file using Potential.from_locpot
        locpot_file = os.path.join(workfunc_job, "LOCPOT")
        if not os.path.exists(locpot_file):
            raise RuntimeError(f"VASP LOCPOT file not found in {workfunc_job}")

        pot = Potential.from_locpot(locpot_file)

        # Read prep_info to get vacuum direction
        prep_info = json.load(open(os.path.join(workfunc_job, "prep_info.json")))
        vacuum_dir = prep_info["vacuum_dir"]

        # Get averaged electrostatic potential profile
        ave_elec_stat, coord_direct = pot.profile1d(axis=vacuum_dir, average=True)
        ave_elec_stat = -ave_elec_stat  # consider negative charge of electron

        pot_file = locpot_file
    else:
        raise ValueError(f"Unsupported job type: {jobtype}")

    profile_data_file = os.path.join(workfunc_job, "profiled.dat")
    with open(profile_data_file, "w") as f:
        f.write(f"{'fraction_coord':>16s}{'ave_elec_stat':>16s}\n")
        for i in range(len(coord_direct)):
            f.write(f"{coord_direct[i]:16.8f}{ave_elec_stat[i]:16.8f}\n")

    work_function_results = calculate_work_functions(ave_elec_stat, fermi_energy=efermi)

    plot_path = plot_averaged_elecstat_pot(
        coord_direct,
        ave_elec_stat,
        efermi,
        axis=vacuum_dir,
        work_function_results=work_function_results,
        plot_filename=os.path.join(workfunc_job, "elecstat_pot_profile.png"),
        title=os.path.basename(job),
    )

    return work_function_results, plot_path, pot_file, profile_data_file


def identify_potential_plateaus(averaged_potential: List, threshold: float = 0.01):
    """
    Identify plateaus in the electrostatic potential using derivative of the averaged electrostatic potential.
    """
    pot_derivatives = []
    n_points = len(averaged_potential)
    stepsize = 1 / (n_points - 1)
    for i in range(n_points):
        backward_idx = i - 1
        forward_idx = 0 if i == n_points - 1 else i + 1
        pot_derivative = (
            averaged_potential[forward_idx] - averaged_potential[backward_idx]
        ) / (stepsize * 2.0)
        pot_derivatives.append(pot_derivative)

    is_plateau = [abs(deriv) < threshold for deriv in pot_derivatives]

    plateau_ranges = []
    in_plateau, start_idx = False, None
    for i in range(n_points):
        if is_plateau[i] and not in_plateau:
            in_plateau = True
            start_idx = i
        elif not is_plateau[i] and in_plateau:
            in_plateau = False
            if not start_idx == i - 1:
                plateau_ranges.append((start_idx, i - 1))
                start_idx = None
        elif is_plateau[i] and i == n_points - 1:
            if not start_idx == i - 1:
                plateau_ranges.append((start_idx, i))
                in_plateau, start_idx = False, None

    if len(plateau_ranges) > 0:
        if plateau_ranges[-1][1] == n_points - 1:
            if len(plateau_ranges) > 1:
                combined_plateau = (
                    plateau_ranges[-1][0] - n_points,
                    plateau_ranges[0][1],
                )
                plateau_ranges[0] = combined_plateau
                plateau_ranges.pop()

    return plateau_ranges


def calculate_work_functions(averaged_potential: List, fermi_energy):
    """
    Calculate the work function from the averaged electrostatic potential.
    Dipole correction is suppoted and multiple plateau of electrostatic potential can be identified.
    """
    work_function_results = []
    plateau_ranges = identify_potential_plateaus(averaged_potential, threshold=0.01)
    for plateau_range in plateau_ranges:
        plateau_start, plateau_end = plateau_range
        plateau_averaged_potential = np.average(
            averaged_potential[plateau_start:plateau_end]
        )
        work_function_results.append(
            {
                "work_function": plateau_averaged_potential - fermi_energy,
                "plateau_start_fractional": plateau_start / len(averaged_potential),
                "plateau_end_fractional": plateau_end / len(averaged_potential),
            }
        )

    return work_function_results


def plot_averaged_elecstat_pot(
    coord: np.ndarray,
    averaged_elecstat_data: np.ndarray,
    efermi: float,
    axis: Literal["a", "b", "c"] = "c",
    work_function_results: Optional[List[Dict[str, float]]] = None,
    plot_filename: Optional[str] = "elecstat_pot_profile.png",
    title: Optional[str] = None,
) -> Path:
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 4))
    plt.plot(coord, averaged_elecstat_data, label="Average electrostatic potential")
    plt.axhline(
        y=efermi,
        linestyle="--",
        color="gray",
        alpha=0.5,
        label=f"Fermi Energy (={efermi:.2f} eV)",
    )
    plt.xlim(min(coord), max(coord))
    plt.xlabel("Fractional coordinate along " + axis)
    plt.ylabel("Electrostatic potential (eV)")

    # Add title if provided
    if title:
        plt.title(title)

    if work_function_results is not None:
        for result in work_function_results:
            # Calculate position of center of plateau
            start_idx = int(result["plateau_start_fractional"] * (len(coord) - 1))
            end_idx = int(result["plateau_end_fractional"] * (len(coord) - 1))
            mid_idx = (start_idx + end_idx) // 2
            x_mid = coord[mid_idx]

            plateau_potential = averaged_elecstat_data[mid_idx]
            work_func = result["work_function"]

            # Draw arrow between plateau and Fermi energy
            plt.annotate(
                "",
                xy=(x_mid, efermi),
                xytext=(x_mid, plateau_potential),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.5),
            )

            # Decide position of annotated text
            x_range = max(coord) - min(coord)
            offset = x_range * 0.02
            if x_mid + offset <= max(coord):
                text_x = x_mid + offset
                ha_align = "left"
            else:
                text_x = x_mid - offset
                ha_align = "right"

            text_y = (efermi + plateau_potential) / 2

            plt.text(
                text_x,
                text_y,
                f"{work_func:.2f} eV",
                fontsize=10,
                ha=ha_align,
                va="center",
            )

    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(plot_filename, dpi=300)
    plt.close()

    return Path(plot_filename).absolute()
