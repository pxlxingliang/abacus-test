from ..model import Model
import json, os
from typing import List, Dict, Any, Literal, Tuple, Optional, Union
import re
import shutil, glob, copy
from pathlib import Path
from itertools import groupby

import numpy as np
from ase.build import bulk, make_supercell
from ase import Atoms
from ase.io import read, write

from abacustest.lib_prepare.abacus import WriteKpt, WriteInput, ReadInput, AbacusStru, ReadKpt
from abacustest.lib_prepare.comm import kpt2kspacing
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE, ELEMENT_CRYSTAL_STRUCTURES, A2BOHR
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.lib_model.model_013_inputs import PrepInput, InputsModel
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
        parser.add_argument('-j', '--job', default=[], action="extend", nargs="*", help='the paths of structure files or ABACUS jobs. If structure files, should be in CIF, VASP POSCAR or ABACUS STRU format. If ABACUS job dirs, should contain INPUT, STRU, or KPT, and pseudopotential and orbital files')
        parser.add_argument('-s', '--supercell', type=int, default=[1, 1, 1], nargs=3, help='the supercell size, default is [1, 1, 1]')
        parser.add_argument('-i', "--index", type=int, default=None, nargs="*", help="Index of the atom to be removed. Start from 1")
        parser.add_argument("--index-file", type=str, default=None, help="The file containing the indices of atoms to be removed for each job. Each line corresponds to a job and the indices are separated by space. If provided, will override the --index argument.")
        parser.add_argument("--cal-reference", type=bool, default=True, help="Whether to do cell-relax calculation for elemental crystals. If not, a file containing element and its energy should be provided.")
        parser.add_argument("--ref-dir", type=str, default='ref_element', help="The directory containing the jobs of elemental crystal structures.")
        parser.add_argument("--max-step", type=int, default=100, help="The maximum number of steps to relax calculation. Default is 100.")
        parser.add_argument("--force-thr-ev", type=float, default=0.01, help="The threshold of force convergence")
        parser.add_argument("--stress-thr-kbar", type=float, default=0.5, help="The threshold of stress convergence")
        parser.add_argument("--image", type=str, default=RECOMMAND_IMAGE, help="The image to use for the Bohrium job, default is %s" % RECOMMAND_IMAGE)
        parser.add_argument("--machine", type=str, default=RECOMMAND_MACHINE, help="The machine to use for the Bohrium job, default is 'c32_m64_cpu'.")
        parser.add_argument("--abacus_command", type=str, default=RECOMMAND_COMMAND, help=f"The command to run the Abacus job, default is '{RECOMMAND_COMMAND}'.")
        parser.add_argument('--ftype', type=str, default='cif', help='The format of the structure files. Can be "cif", "poscar" or "stru". Default is cif.')
        parser.add_argument("--pp",default=None,type=str,help="the path of pseudopotential library, or read from enviroment variable ABACUS_PP_PATH")
        parser.add_argument("--orb",default=None,type=str,help="the path of orbital library, or read from enviroment variable ABACUS_ORB_PATH")
        parser.add_argument("--input",default=None,type=str,help="the template of input file, if not specified, the default input will be generated")
        parser.add_argument("--relax-kspacing", default=None, type=float, help="the kspacing for kpoint generation of relax calculation, default is None, which means kspacing in original input file will be used.")
        parser.add_argument("--scf-kspacing", default=None, type=float, help="the kspacing for scf calculation after relax, default is None, which means kspacing in original input file will be used.")
        parser.add_argument("--lcao", action="store_true", help="whether to use lcao basis, default is pw basis")
        parser.add_argument("--nspin", default=1, type=int, choices=[1, 2, 4], help="the number of spins, can be 1 (no spin), 2 (spin polarized), or 4 (non-collinear spin). Default is 1.")
        parser.add_argument("--soc", action="store_true", help="whether to use spin-orbit coupling, if True, nspin should be 4.")
        parser.add_argument("--dftu", action="store_true", help="whether to use DFT+U, default is False.")
        parser.add_argument("--dftu_param", default=None, nargs="+", help="the DFT+U parameters, should be element symbol and U value pairs like 'Fe 4 Ti 1'. If dftu is set, but dftu_param is not set, the default U values will be used: 4 eV for d orbital elements and 6 eV for f orbital elements.")
        parser.add_argument("--init_mag", default=None, nargs="+", help="the initial magnetic moment for magnetic elements, should be element symbol and magnetic moment pairs like 'Fe 4 Ti 1'.")
        parser.add_argument("--afm", action="store_true", help="whether to use antiferromagnetic calculation, default is False. Only valid when init_mag is set.")
        parser.add_argument("--copy_pp_orb", action="store_true", help="whether to copy the pseudopotential and orbital files to each job directory or link them. Default is False, which means linking the files.")

        return parser
    
    def run_prepare(self, params):
        '''
        Run the model with the given parameters
        '''
        if not params.job:
            raise ValueError("No job specified, please use -j or --job to specify the job paths.")
        if params.dftu_param is not None:
            dftu_param = InputsModel.parse_dftu_param(params.dftu_param)
        else:
            dftu_param = None
            
        if params.init_mag is not None:
            init_mag = InputsModel.parse_init_mag(params.init_mag)
        else:
            init_mag = None
        
        if params.index_file is None and params.index is None:
            raise ValueError("No vacancy index specified, please use --index or --index-file to specify the atom indices to be removed.")
        
        folders = prepare_vacancy_jobs(jobs=params.job,
                                       supercell=params.supercell,
                                       vacancy_indices=params.index_file if params.index_file is not None else params.index,
                                       cal_reference=params.cal_reference,
                                       ref_dir=params.ref_dir,
                                       max_step=params.max_step,
                                       force_thr_ev=params.force_thr_ev,
                                       stress_thr_kbar=params.stress_thr_kbar,
                                       ftype=params.ftype,
                                       pp=params.pp,
                                       orb=params.orb,
                                       input=params.input,
                                       relax_kspacing=params.relax_kspacing,
                                       scf_kspacing=params.scf_kspacing,
                                       lcao=params.lcao,
                                       nspin=params.nspin,
                                       soc=params.soc,
                                       dftu=params.dftu,
                                       dftu_param=dftu_param,
                                       init_mag=init_mag,
                                       afm=params.afm,
                                       copy_pp_orb=params.copy_pp_orb)
        
        run_sh = f"""{params.abacus_command} | tee log
if [ -d "final_scf" ]; then
    cp OUT.*/STRU_ION_D final_scf/STRU
    cd final_scf
    {params.abacus_command} | tee log
    cd ..
fi
"""
        # Write run scripts
        with open("run.sh", "w") as f1:
            f1.write(run_sh)
        
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
        print(f"After finishing the calculations, you can enter the results directory, \nand run below command below to postprocess the vacancy formation energy results:\n\tabacustest model vacancy post -j ...\n")
        

    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand 
        '''
        parser.description = "Postprocess the vacancy formation energy calculation results."
        parser.add_argument('-j', '--job', default=[], action="extend", nargs="*", help='the paths of the job directories, should contain the results of deformed structures generated by the prepare step.')
        parser.add_argument('--ref-dir', default='ref_element', help="The ABACUS job directory containing the elemental crystal structures. If not provided, reference atom energies will be read from provided reference file.")
        parser.add_argument("--ref-file", type=str, default='ref_energy.txt', help="The file to store the atom energy reference. Will be generated after a vacancy formation energy calculation." + \
                            "Has two columns: atom label, and atom energy. If not provided or selected atom does not contained in the atom energy reference file, the atom energy will be automatically calculated.")
        return parser   

    def run_postprocess(self, params):
        '''
        Parse the parameters and run the postprocess process
        '''
        metrics, ref_atom_energies = postprocess_vacancy(params.job, params.ref_dir, params.ref_file)
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

def parser_indices_file(index_file: str) -> Dict[str, List[int]]:
    """
    Parse the index file to get the vacancy indices for each job.
    The file should contain multiple lines, each line corresponds to a job and contains the indices of atoms to be removed, separated by space.
    Returns a dictionary with job names as keys and list of indices as values.

    Config file example:
    job1 1 2 3
    job2 4 5 6
    """
    job_indices = {}
    with open(index_file, "r") as f:
        for line in f:
            if line.strip() == "":
                continue
            parts = line.strip().split()
            job_name = parts[0]
            indices = [int(idx) for idx in parts[1:]]
            assert all(idx > 0 for idx in indices), f"All indices should be positive integers in line: {line.strip()}"
            job_indices[job_name] = indices
    return job_indices

def prepare_vacancy_jobs(
    jobs: str,
    supercell: Tuple[int, int, int],
    vacancy_indices: Union[str, List[int]],
    cal_reference: bool = True,
    ref_dir: str = "ref_element",
    max_step: int = 100,
    force_thr_ev: float = 0.01,
    stress_thr_kbar: float = 0.5,
    ftype: str = "cif",
    pp: Optional[str] = None,
    orb: Optional[str] = None,
    input: Optional[str] = None,
    relax_kspacing: Optional[float] = None,
    scf_kspacing: Optional[float] = None,
    lcao: bool = False,
    nspin: int = 1,
    soc: bool = False,
    dftu: bool = False,
    dftu_param: Optional[Dict[str, Union[float, Tuple[Literal["s", 'p', "d"], float]]]] = None,
    init_mag: Optional[Dict[str, float]] = None,
    afm: bool = False,
    copy_pp_orb: bool = False
) -> None:
    '''
    Prepare the vacancy calculation jobs.
    '''
    #ref_atom_energies = read_ref_atom_energies(ref_file)
    if pp is None:
        pp = os.environ.get("ABACUS_PP_PATH", None)
    if orb is None:
        orb = os.environ.get("ABACUS_ORB_PATH", None)
    prepared_elements = []
    if cal_reference:
        os.makedirs(ref_dir, exist_ok=True)

    if isinstance(vacancy_indices, list):
        job_indices_dict = {job: vacancy_indices for job in jobs}
    else:
        try:
            job_indices_dict = parser_indices_file(vacancy_indices)
        except Exception as e:
            raise ValueError(f"Error parsing vacancy indices file: {e}")

    folders = []
    for job in jobs:
        print("Job: ", job)
        if job not in job_indices_dict:
            raise ValueError(f"Vacancy indices not specified for job {job} in index file.")
        original_vacancy_indices = remove_redundant_indices(job_indices_dict[job])
        real_vacancy_indices = copy.deepcopy(original_vacancy_indices)
        if not os.path.isdir(job):
            # Create structure with vacancy
            #print((job, ftype, "scf", pp, orb, input, kpt, lcao, nspin, soc, dftu, dftu_param, init_mag, afm, copy_pp_orb))
            setting, job_path = PrepInput(files=job, 
                                          filetype=ftype, 
                                          jobtype="scf", 
                                          pp_path=pp, 
                                          orb_path=orb, 
                                          input_file=input, 
                                          kpt=None, 
                                          lcao=lcao, 
                                          nspin=nspin,
                                          soc=soc,
                                          dftu=dftu,
                                          dftu_param=dftu_param,
                                          init_mag=init_mag,
                                          afm=afm,
                                          copy_pp_orb=copy_pp_orb).run()
            
            if ftype != "stru":
                # Transform the vacancy indices to those in ABACUS STRU file
                categorized_idx = get_categorized_idx(job, ftype)
                for i in range(len(original_vacancy_indices)):
                    # Get global index in ABACUS STRU file
                    assert original_vacancy_indices[i]-1 in categorized_idx, f"Vacancy index {original_vacancy_indices[i]} is out of range for job {job}."
                    real_vacancy_indices[i] = categorized_idx.index(original_vacancy_indices[i]-1) + 1 # Use atom index starting from 1
            
            new_path = job.lstrip("./").replace("/", "_") + "_VACANCY"
            if os.path.exists(new_path):
                shutil.rmtree(new_path)
            shutil.move(job_path[0], new_path)

            job = new_path
            

        print("Preparing vacancy formation energy calculation for job path:", job)
        # clean old folders
        old_folders = [f for f in glob.glob(os.path.join(job, "vacancy_*")) if os.path.isdir(f)]
        if len(old_folders) > 0:
            print(f"Remove old vacancy_* folders in job path: {job}")
            for f in old_folders:
                shutil.rmtree(f)
        
        if not os.path.isfile(os.path.join(job, "STRU")):
            raise FileNotFoundError(f"STRU file not found in job path: {job}")
        if not os.path.isfile(os.path.join(job, "INPUT")):
            raise FileNotFoundError(f"INPUT file not found in job path: {job}")
        
        input_params = ReadInput(os.path.join(job, "INPUT"))
        original_stru = AbacusStru.ReadStru(os.path.join(job, input_params.get('stru_file', 'STRU')))

        # all indices should be valid
        for idx in real_vacancy_indices:
            if idx < 1 or idx > original_stru.get_natoms():
                raise ValueError(f"Vacancy index {idx} is out of range for job {job} with {original_stru.get_natoms()} atoms.")

        input_params["suffix"] = "ABACUS"
        input_params['calculation'] = 'cell-relax'
        input_params['relax_method'] = 'cg'
        input_params['symmetry'] = 0
        input_params['stru_file'] = None
        input_params['ntype'] = None
        input_params['relax_nmax'] = max_step
        input_params['force_thr_ev'] = force_thr_ev
        input_params['stress_thr'] = stress_thr_kbar
        if relax_kspacing is not None:
            input_params['kspacing'] = relax_kspacing

        pp_orb_files_fullname = glob.glob(os.path.join(job, "*.upf")) + glob.glob(os.path.join(job, "*.UPF")) \
                       + glob.glob(os.path.join(job, "*.orb"))
        
        pp_orb_files = [os.path.basename(i) for i in pp_orb_files_fullname]
        original_kpt_file = os.path.join(job, input_params.get('kpt_file', 'KPT'))

        # Prepare input files for the supercell
        # cell-relax job for the original structure
        original_stru_jobpath = os.path.join(job, "vacancy_original_stru")
        os.makedirs(original_stru_jobpath, exist_ok=True)
        create_new_abacus_job(job,
                              original_stru_jobpath,
                              pp_orb_files,
                              extra_input=input_params,
                              stru=original_stru, 
                              kpt_file=original_kpt_file)
        # scf job after cell-relax job for the original structure
        if scf_kspacing is None:
            if 'kspacing' not in input_params.keys():
                original_kpt = ReadKpt(os.path.join(original_stru_jobpath, "KPT"))
                original_kspacing = kpt2kspacing(original_kpt[0][:3], original_stru.get_cell())
            else:
                original_kspacing = input_params['kspacing']
            scf_kspacing = min(original_kspacing*0.9, 0.10)
        
        create_new_abacus_job(job,
                              os.path.join(original_stru_jobpath, "final_scf"),
                              pp_orb_files,
                              extra_input={'calculation': 'scf', 'kspacing': scf_kspacing},
                              stru=None, # Don't prepare STRU file
                              kpt_file=None)
        
        folders.append(original_stru_jobpath)

        # Prepare input files for the supercell with defect
        for i, idx in enumerate(real_vacancy_indices):
            # If ABACUS inputs directory is provided
            vacancy_element, labelidx = original_stru.globalidx2labelidx(idx-1) # Use atom index staring from 0 in globalidx2labelidx
            supercell_stru = original_stru.supercell(supercell)
            defect_supercell_stru = copy.deepcopy(supercell_stru)
            defect_supercell_stru.set_empty_atom(vacancy_element, labelidx)
            defect_supercell_jobpath = os.path.join(job, f"vacancy_defect_{original_vacancy_indices[i]}_{vacancy_element}_{supercell[0]}_{supercell[1]}_{supercell[2]}")
            os.makedirs(defect_supercell_jobpath, exist_ok=True)
            # cell-relax job for the defect structure
            create_new_abacus_job(job,
                                  defect_supercell_jobpath,
                                  pp_orb_files,
                                  extra_input=input_params,
                                  stru=defect_supercell_stru,
                                  kpt_file=original_kpt_file)
            # scf job after cell-relax job for the defect structure
            create_new_abacus_job(job,
                                  os.path.join(defect_supercell_jobpath, "final_scf"),
                                  pp_orb_files,
                                  extra_input={'calculation': 'scf', 'kspacing': scf_kspacing},
                                  stru=None, # Don't prepare STRU file
                                  kpt_file=None) 
            defect_supercell_stru.write2cif(os.path.join(defect_supercell_jobpath, "STRU.cif"), empty2x=True)
            folders.append(defect_supercell_jobpath)

            # Prepare input files for the crystal structure of the vacancy element
            if cal_reference and vacancy_element not in prepared_elements:
                element_type_index = [key for key, group in groupby(original_stru.get_element(number=False))].index(vacancy_element)
                vacancy_element_crys_jobpath = os.path.join(ref_dir, f"{vacancy_element}")
                # Remove old folders when cal_reference is True
                if os.path.isdir(vacancy_element_crys_jobpath):
                    shutil.rmtree(vacancy_element_crys_jobpath)
                vacancy_element_pp = original_stru.get_pp()[element_type_index]
                vacancy_element_orb = original_stru.get_orb()[element_type_index]
                vacancy_element_crys_stru = build_most_stable_elementary_crys_stru(vacancy_element, vacancy_element_pp, vacancy_element_orb)

                create_new_abacus_job(job,
                                      vacancy_element_crys_jobpath,
                                      pp_orb_files,
                                      extra_input=input_params,
                                      stru=vacancy_element_crys_stru,
                                      kpt_file=None)
                
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
                        ref_dir: str,
                        ref_file: str = None) -> Dict[str, Any]:
    '''
    Postprocess the vacancy calculation results.
    
    Args:
        jobs: List of job paths.
    
    Returns:
        A dictionary containing the vacancy results.
    '''
    # Get reference atom energies from file
    if ref_file is None or not os.path.isfile(ref_file):
        ref_atom_energies = {}
    else:
        ref_atom_energies = read_ref_atom_energies(ref_file)
    
    # Get element crystal results from calculation results. Will overwrite the reference atom energies from file
    if ref_dir is not None and os.path.isdir(ref_dir):
        for element_crys_job in glob.glob(os.path.join(ref_dir, "*")):
            vac_ele_crys_stru = AbacusStru.ReadStru(os.path.join(element_crys_job, "STRU"))
            vacancy_element = vac_ele_crys_stru.get_element(number=False)[0]
            element_crys_job_result = read_metrics(element_crys_job, jobtype="relax")
            e_vac_elem_crys = element_crys_job_result["energies"][-1]
            ref_atom_energies[vacancy_element] = e_vac_elem_crys / vac_ele_crys_stru.get_natoms()

    metrics = {}
    for job in jobs:
        print("Postprocessing vacancy calculation for job path:", job)

        sub_folders = [f for f in glob.glob(os.path.join(job, "vacancy_*")) if os.path.isdir(f)]

        defect_supercell_job_results, defect_supercell_scf_results = {}, {}
        for sub_folder in sub_folders:
            if os.path.basename(sub_folder).startswith("vacancy_original_stru"):
                original_stru_job_results = read_metrics(sub_folder, jobtype="relax")
                original_stru_scf_results = read_metrics(os.path.join(sub_folder, 'final_scf'), jobtype="scf")
            if os.path.basename(sub_folder).startswith("vacancy_defect"):
                words = sub_folder.split('/')[-1].split("_")
                vacancy_element, idx = words[3], int(words[2])
                supercell = (int(words[4]), int(words[5]), int(words[6]))
                defect_supercell_job_results[f'{vacancy_element}{idx}'] = read_metrics(sub_folder, jobtype="relax")
                defect_supercell_scf_results[f'{vacancy_element}{idx}'] = read_metrics(os.path.join(sub_folder, 'final_scf'), jobtype="scf")

        for site in defect_supercell_job_results.keys():
            if original_stru_job_results["energies"] is None or defect_supercell_job_results[site]["energies"] is None:
                e_vac_form = None
            elif original_stru_scf_results["converge"] is False or defect_supercell_scf_results[site]["converge"] is False:
                e_vac_form = None
            else:
                e_vac_form = (defect_supercell_scf_results[site]["energy"] + ref_atom_energies[vacancy_element]) - original_stru_scf_results["energy"] * supercell[0] * supercell[1] * supercell[2]

            results = {
                'vac_formation_energy': e_vac_form,
                'original_stru_job_relax_converge': original_stru_job_results['relax_converge'],
                'original_stru_job_normal_end': original_stru_job_results['normal_end'],
                'original_stru_job_max_force': original_stru_job_results['largest_gradient'][-1],
                'original_stru_job_max_stress': original_stru_job_results['largest_gradient_stress'][-1],
                'original_stru_relaxed_lattice_constant': original_stru_job_results['lattice_constant'],
                'defect_supercell_job_relax_converge': defect_supercell_job_results[site]['relax_converge'],
                'defect_supercell_job_normal_end': defect_supercell_job_results[site]['normal_end'],
                'defect_supercell_job_max_force': defect_supercell_job_results[site]['largest_gradient'][-1],
                'defect_supercell_job_max_stress': defect_supercell_job_results[site]['largest_gradient_stress'][-1],
                'defect_supercell_relaxed_lattice_constant': defect_supercell_job_results[site]['lattice_constant'],
            }
            metrics[f"{job.rstrip('/')}-{site}"] = results

    return metrics, ref_atom_energies

def create_new_abacus_job(src_dir: str,
                          dst_dir: str,
                          pp_orb_files: str,
                          extra_input: dict=None,
                          stru: Tuple[AbacusStru, str]=None,
                          kpt_file: str=None):
    """
    Create new abacus job from a source directory and new settings.
    """
    if not os.path.isdir(dst_dir):
        os.makedirs(dst_dir, exist_ok=True)
    
    for file in pp_orb_files:
        os.symlink(os.path.abspath(os.path.join(src_dir, file)), os.path.join(dst_dir, file))
    
    if extra_input is None:
        shutil.copy(os.path.join(src_dir, "INPUT"), os.path.join(dst_dir, "INPUT"))
    else:
        input_params = ReadInput(os.path.join(src_dir, "INPUT"))
        for key, value in extra_input.items():
            input_params[key] = value
        WriteInput(input_params, os.path.join(dst_dir, "INPUT"))
    
    
    if stru is not None:
        if isinstance(stru, AbacusStru):
            stru.write(os.path.join(dst_dir, "STRU"))
        elif isinstance(stru, str):
            shutil.copy(stru, os.path.join(dst_dir, "STRU"))

    if kpt_file is not None and os.path.isfile(kpt_file):
        shutil.copy(kpt_file, os.path.join(dst_dir, os.path.basename(kpt_file)))

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
                             cartesian=True)
    stru_abacus.set_pp([pp])
    stru_abacus.set_orb([orb])

    return stru_abacus

def read_metrics(job: str, jobtype: Literal['relax', 'scf']) -> List[str]:
    """
    Read the relaxation metrics from the log file.
    Args:
        job (str): The path to the job directory.
    Returns:
        List[str]: The relaxation metrics.
    """
    results = RESULT(path=job, fmt='abacus')
    if jobtype == 'relax':
        relax_metrics = ["normal_end", "relax_steps", "largest_gradient",
                         "largest_gradient_stress", "relax_converge", "energies", "lattice_constant"]
        return {k: results[k] for k in relax_metrics}
    elif jobtype == 'scf':
        scf_metrics = ["normal_end", "energy", "converge"]
        return {k: results[k] for k in scf_metrics}
    else:
        raise ValueError(f"Job type {jobtype} not supported.")

def get_categorized_idx(original_stru_file: str, original_stru_format: str):
    """
    Get index correspondace between original structure file and ABACUS STRU file.
    Return a list, where its index is the index in the transformed ABACUS STRU file, and its value is the index in the original structure file.
    """
    original_stru = read(original_stru_file, format=original_stru_format)

    unique_labels = []
    for label in original_stru.get_chemical_symbols():
        if label not in unique_labels:
            unique_labels.append(label)
    
    original_indices = []
    for label in unique_labels:
        for i, atom_label in enumerate(original_stru.get_chemical_symbols()):
            if atom_label == label:
                original_indices.append(i)
    
    return original_indices
