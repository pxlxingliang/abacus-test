from ..model import Model
import json, os
from typing import List, Dict, Any, Literal
import re
import shutil, glob, copy
from pathlib import Path
from itertools import groupby

from ase.build import bulk

from abacustest.lib_prepare.abacus import WriteKpt, WriteInput, ReadInput, AbacusStru, ReadKpt
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE, ELEMENT_CRYSTAL_STRUCTURES, A2BOHR
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.outresult import pandas_out

class VacancyModel(Model):
    @staticmethod
    def model_name(): # type: ignore
        """
        Name of the model, which will be used as the subcommand
        """
        return "vacancy"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "Calculate the vacancy formation energy for uncharged systems"
    
    @staticmethod
    def prepare_args(parser):
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand 
        '''
        parser.description = "Prepare the inputs for vacancy formation energy calculation."
        parser.add_argument('-j', '--job', default=[], action="extend", nargs="*", help='the paths of ABACUS jobs, should contain INPUT, STRU, or KPT, and pseudopotential and orbital files')
        parser.add_argument('-s', '--supercell', type=int, default=[1, 1, 1], nargs=3, help='the supercell size, default is [1, 1, 1]')
        parser.add_argument("--index", type=int, required=True, nargs="*", help="Index of the atom to be removed. Start from 0")
        parser.add_argument("--cal-element", action="store_true", help="Whether to do cell-relax calculation for elemental crystals. If not, a file containing element and its energy should be provided.")
        parser.add_argument("--force-thr-ev", type=float, default=0.01, help="The threshold of force convergence")
        parser.add_argument("--stress-thr-kbar", type=float, default=0.5, help="The threshold of stress convergence")
        parser.add_argument("--save-dir", type=str, default="vac_form_calcs", help="The folder to save the prepared directories of ABACUS jobs.")
        parser.add_argument("--image", type=str, default=RECOMMAND_IMAGE, help="The image to use for the Bohrium job, default is %s" % RECOMMAND_IMAGE)
        parser.add_argument("--machine", type=str, default=RECOMMAND_MACHINE, help="The machine to use for the Bohrium job, default is 'c32_m64_cpu'.")
        parser.add_argument("--abacus_command", type=str, default=RECOMMAND_COMMAND, help=f"The command to run the Abacus job, default is '{RECOMMAND_COMMAND}'.")

        return parser
    
    def run_prepare(self, params):
        '''
        Run the model with the given parameters
        '''
        if not params.job:
            raise ValueError("No job specified, please use -j or --job to specify the job paths.")
        
        folders = prepare_vacancy_jobs(params.job,
                                       params.supercell,
                                       params.index,
                                       params.cal_element,
                                       params.force_thr_ev,
                                       params.stress_thr_kbar,
                                       params.save_dir)
        
        # Write run scripts
        with open("run.sh", "w") as f1:
            f1.write(f"{params.abacus_command} | tee log\n")
        
        setting = {
            "save_path": "results",
            "run_dft": [{
                "ifrun": True,
                "example": folders,
                "extra_files":["run.sh"],
                "command": "bash run.sh",
                "image": params.image,
                "bohrium": {
                    "scass_type": params.machine,
                    "job_type": "container",
                    "platform": "ali"
                }
            }]
        } 
        
        setting_file = "setting.json"
        json.dump(setting, open(setting_file, "w"), indent=4)
        print("\nThe inputs are generated in", ", ".join(params.job))
        print(f"You can modify '{setting_file}' and 'run.sh', and execute below command to run the abacustest to submit all jobs to bohrium:\n\tabacustest submit -p setting.json\n")
        print(f"After finishing the calculations, you can enter the results directory, \nand run below command below to postprocess the vacancy formation energy results:\n\tabacustest model vacancy post -j {' '.join(params.job)}\n")
        

    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand 
        '''
        parser.description = "Postprocess the vacancy formation energy calculation results."
        parser.add_argument('-j', '--job', default=[], action="extend", nargs="*", help='the paths of the job directories, should contain the results of deformed structures generated by the prepare step.')
        parser.add_argument('--element-dir', default=None, help="The ABACUS job directory containing the elemental crystal structures. If not provided, reference atom energies will be read from provided reference file.")
        parser.add_argument("--ref-file", type=str, default='ref_energy.txt', help="The file to store the atom energy reference. Will be generated after a vacancy formation energy calculation." + \
                            "Has two columns: atom label, and atom energy. If not provided or selected atom does not contained in the atom energy reference file, the atom energy will be automatically calculated.")
        return parser   

    def run_postprocess(self, params):
        '''
        Parse the parameters and run the postprocess process
        '''
        metrics, ref_atom_energies = postprocess_vacancy(params.job, params.element_dir, params.ref_file)
        write_ref_atom_energies(params.ref_file, ref_atom_energies)
        json.dump(metrics, open("metrics_vacancy.json", "w"), indent=4)
        pandas_out(metrics)
        print(f"\nThe vacancy formation results are saved in 'metrics_vacancy.json'\nThe reference atom energies are saved in '{params.ref_file}'\n")

def read_ref_atom_energies(ref_atom_energy_file: str) -> Dict[str, float]:
    '''
    Read the reference atom energies from the file.
    '''
    ref_atom_energies = {}
    with open(ref_atom_energy_file, "r") as f:
        for line in f:
            label, energy = line.strip().split()
            ref_atom_energies[label] = float(energy)
    return ref_atom_energies

def write_ref_atom_energies(ref_atom_energy_file: str, ref_atom_energies: Dict[str, float]):
    '''
    Write the reference atom energies to the file.
    '''
    with open(ref_atom_energy_file, "w") as f:
        for label, energy in ref_atom_energies.items():
            f.write(f"{label} {energy}\n")

def prepare_vacancy_jobs(
    jobs: str,
    supercell: List[int],
    vacancy_indices: List[int],
    cal_element: bool = True,
    force_thr_ev: float = 0.01,
    stress_thr_kbar: float = 0.5,
    save_dir: str = "vac_form_calcs"
) -> None:
    '''
    Prepare the vacancy calculation jobs.
    '''
    #ref_atom_energies = read_ref_atom_energies(ref_file)
    prepared_elements = []
    vacancy_indices = remove_redundant_indices(vacancy_indices)

    old_folders = [os.path.join(save_dir, f) for f in os.listdir(save_dir) if os.path.isdir(os.path.join(save_dir, f))]
    if os.path.join(save_dir, 'element_crystals') in old_folders:
        shutil.rmtree(os.path.join(save_dir, 'element_crystals'))

    folders = []
    for job in jobs:
        job_dir = os.path.join(save_dir, job)
        print("Preparing vacancy formation energy calculation for job path:", job)
        # clean old folders
        if len(old_folders) > 0:
            if job_dir in old_folders:
                shutil.rmtree(job_dir)
        
        if not os.path.isfile(os.path.join(job, "STRU")):
            raise FileNotFoundError(f"STRU file not found in job path: {job}")
        if not os.path.isfile(os.path.join(job, "INPUT")):
            raise FileNotFoundError(f"INPUT file not found in job path: {job}")
        
        input_params = ReadInput(os.path.join(job, "INPUT"))
        original_stru = AbacusStru.ReadStru(os.path.join(job, input_params.get('stru_file', 'STRU')))

        input_params["suffix"] = "ABACUS"
        input_params['calculation'] = 'cell-relax'
        input_params['relax_method'] = 'cg'
        input_params['stru_file'] = None
        input_params['ntype'] = None
        input_params['force_thr_ev'] = force_thr_ev
        input_params['stress_thr'] = stress_thr_kbar

        pp_orb_files_fullname = glob.glob(os.path.join(job, "*.upf")) + glob.glob(os.path.join(job, "*.UPF")) \
                       + glob.glob(os.path.join(job, "*.orb"))
        pp_orb_files = [os.path.basename(i) for i in pp_orb_files_fullname]
        original_kpt_file = os.path.join(job, input_params.get('kpt_file', 'KPT'))

        # Prepare input files for the supercell
        supercell_stru = original_stru.supercell(supercell)
        supercell_jobpath = os.path.join(job_dir, f"vacancy_supercell_{supercell[0]}_{supercell[1]}_{supercell[2]}")
        os.makedirs(supercell_jobpath, exist_ok=True)
        copy_pp_orb_kpt_file(job, supercell_jobpath, pp_orb_files, original_kpt_file)
        write_inputs(supercell_jobpath, input_params, supercell_stru)
        folders.append(supercell_jobpath)

        # Prepare input files for the supercell with defect
        for idx in vacancy_indices:
            vacancy_element, labelidx = original_stru.globalidx2labelidx(idx)
            defect_supercell_stru = copy.deepcopy(supercell_stru)
            defect_supercell_stru.set_empty_atom(vacancy_element, labelidx)
            defect_supercell_jobpath = os.path.join(job_dir, f"vacancy_defect_{idx}_{vacancy_element}_{labelidx}_{supercell[0]}_{supercell[1]}_{supercell[2]}")
            os.makedirs(defect_supercell_jobpath, exist_ok=True)
            copy_pp_orb_kpt_file(job, defect_supercell_jobpath, pp_orb_files, original_kpt_file)
            write_inputs(defect_supercell_jobpath, input_params, defect_supercell_stru)
            folders.append(defect_supercell_jobpath)

            element_crys_path = os.path.join(save_dir, "element_crystals")
            os.makedirs(element_crys_path, exist_ok=True)
            # Prepare input files for the crystal structure of the vacancy element
            if cal_element and vacancy_element not in prepared_elements:
                element_type_index = [key for key, group in groupby(original_stru.get_element(number=False))].index(vacancy_element)
                vacancy_element_crys_jobpath = os.path.join(element_crys_path, f"{vacancy_element}")
                vacancy_element_pp = original_stru.get_pp()[element_type_index]
                vacancy_element_orb = original_stru.get_orb()[element_type_index]
                vacancy_element_crys_stru = build_most_stable_elementary_crys_stru(vacancy_element, vacancy_element_pp, vacancy_element_orb)
                os.makedirs(vacancy_element_crys_jobpath, exist_ok=True)
                copy_pp_orb_kpt_file(job, vacancy_element_crys_jobpath, pp_orb_files, original_kpt_file)
                write_inputs(vacancy_element_crys_jobpath, input_params, vacancy_element_crys_stru)
                folders.append(vacancy_element_crys_jobpath)
                prepared_elements.append(vacancy_element)
    
    return folders

def remove_redundant_indices(indices: List[int]):
    """
    Remove the same atom indices in the list.
    """
    new_indices = []
    for ind in indices:
        if ind not in new_indices:
            new_indices.append(ind)
        else:
            print(f"The atom index {ind} is already in the list, skip")
    
    return new_indices

def postprocess_vacancy(jobs: List[str],
                        element_dir: str,
                        ref_file: str) -> Dict[str, Any]:
    '''
    Postprocess the vacancy calculation results.
    
    Args:
        jobs: List of job paths.
    
    Returns:
        A dictionary containing the vacancy results.
    '''
    # Get reference atom energies from file
    if not os.path.isfile(ref_file):
        ref_atom_energies = {}
    else:
        ref_atom_energies = read_ref_atom_energies(ref_file)
    
    # Get element crystal results from calculation results. Will overwrite the reference atom energies from file
    if element_dir is not None and os.path.isdir(element_dir):
        for element_crys_job in glob.glob(os.path.join(element_dir, "*")):
            vac_ele_crys_stru = AbacusStru.ReadStru(os.path.join(element_crys_job, "STRU"))
            vacancy_element = vac_ele_crys_stru.get_element(number=False)[0]
            element_crys_job_result = read_relax_metrics(element_crys_job)
            e_vac_elem_crys = element_crys_job_result["energies"][-1]
            ref_atom_energies[vacancy_element] = e_vac_elem_crys / vac_ele_crys_stru.get_natoms()

    metrics = {}
    for job in jobs:
        print("Postprocessing vacancy calculation for job path:", job)

        sub_folders = [f for f in glob.glob(os.path.join(job, "vacancy_*")) if os.path.isdir(f)]

        defect_supercell_job_results = {}
        for sub_folder in sub_folders:
            if os.path.basename(sub_folder).startswith("vacancy_supercell"):
                supercell_job_results = read_relax_metrics(sub_folder)
            if os.path.basename(sub_folder).startswith("vacancy_defect"):
                words = sub_folder.split("_")
                vacancy_element, idx = words[3], int(words[2])
                defect_supercell_job_results[f'{vacancy_element}{idx}'] = read_relax_metrics(sub_folder)
        
        e_supercell = supercell_job_results["energies"][-1]

        for site in defect_supercell_job_results.keys():
            e_defect_supercell = defect_supercell_job_results[site]["energies"][-1]
            e_vac_form = (e_defect_supercell + ref_atom_energies[vacancy_element]) - e_supercell

            results = {
                'vac_formation_energy': e_vac_form,
                'supercell_job_relax_converge': supercell_job_results['relax_converge'],
                'supercell_job_normal_end': supercell_job_results['normal_end'],
                'supercell_relaxed_lattice_constant': supercell_job_results['lattice_constant'],
                'defect_supercell_job_relax_converge': defect_supercell_job_results[site]['relax_converge'],
                'defect_supercell_job_normal_end': defect_supercell_job_results[site]['normal_end'],
                'defect_supercell_relaxed_lattice_constant': defect_supercell_job_results[site]['lattice_constant'],
            }
            metrics[f"{job.rstrip('/')}-{site}"] = results

    return metrics, ref_atom_energies

def copy_pp_orb_kpt_file(src_dir: str, dst_dir: str, pp_orb_files: str, kpt_file: str):
    '''
    Copy the pp, orb, and kpt file from source directory to the destination directory.
    '''
    for file in pp_orb_files:
        os.symlink(os.path.abspath(os.path.join(src_dir, file)), os.path.join(dst_dir, file))
    
    if os.path.isfile(kpt_file):
        shutil.copy(kpt_file, os.path.join(dst_dir, os.path.basename(kpt_file)))

def write_inputs(folder, input_params: dict, stru: AbacusStru):
    # Write INPUT and STRU file
    WriteInput(input_params, os.path.join(folder, "INPUT"))
    stru.write(os.path.join(folder, "STRU"))

def build_most_stable_elementary_crys_stru(element: str, pp: str, orb: str) -> AbacusStru:
    """
    Build the most stable crystal structure of an element.
    Args:
        element (str): The chemical symbol of the element.
        pp (dict): Pseudopotential information from the original structure.
        orb (dict): Orbital information from the original structure.
    Returns:
        AbacusStru: The most stable crystal structure of the element.
    Raises:
        ValueError: If the element is not supported.
    """
    if element not in ELEMENT_CRYSTAL_STRUCTURES:
        raise ValueError(f"Element {element} not supported.")
    else:
        crystal_structure = ELEMENT_CRYSTAL_STRUCTURES[element]
        if crystal_structure['crystal'] == 'bcc':
            stru = bulk(element, 'bcc', a=crystal_structure['a'])
        elif crystal_structure['crystal'] == 'fcc':
            stru = bulk(element, 'fcc', a=crystal_structure['a'])
        elif crystal_structure['crystal'] == 'hcp':
            stru = bulk(element, 'hcp', a=crystal_structure['a'], c=crystal_structure['c'])
        elif crystal_structure['crystal'] == 'diamond':
            stru = bulk(element, 'diamond', a=crystal_structure['a'])
        else:
            raise ValueError(f"Crystal structure {crystal_structure['crystal']} not supported.")
    
    stru_abacus = AbacusStru(label=stru.get_chemical_symbols(),
                             cell=stru.get_cell(),
                             coord=stru.get_positions(),
                             lattice_constant=A2BOHR,
                             cartesian=False)
    stru_abacus.set_pp([pp])
    stru_abacus.set_orb([orb])

    return stru_abacus

def read_relax_metrics(job: str) -> List[str]:
    """
    Read the relaxation metrics from the log file.
    Args:
        job (str): The path to the job directory.
    Returns:
        List[str]: The relaxation metrics.
    """
    results = RESULT(path=job, fmt='abacus')
    relax_metrics = ["normal_end", "relax_steps", "largest_gradient",
                     "largest_gradient_stress", "relax_converge", "energies", "lattice_constant"]
    return {k: results[k] for k in relax_metrics}
