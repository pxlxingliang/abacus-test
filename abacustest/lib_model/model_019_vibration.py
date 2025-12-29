from ..model import Model
import json, os
from abacustest.lib_prepare.abacus import WriteKpt, WriteInput, ReadInput, AbacusStru, ReadKpt
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

from typing import List, Dict, Any

import shutil, glob, copy
import numpy as np
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo
from abacustest.lib_prepare.abacus import AbacusStru, ReadInput, WriteInput
from abacustest.lib_prepare.stru import AbacusSTRU
from pathlib import Path
import math
import os
import copy
import json


class VibrationModel(Model):
    @staticmethod
    def model_name(): # type: ignore
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "vibration"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "Calculate the vibration frequency of selected atoms"
    
    @staticmethod
    def prepare_args(parser):
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand 
        '''
        parser.description = "Prepare the inputs for vibration frequency calculation."
        parser.add_argument('-j', '--job', default=[], action="extend", nargs="*", help='the paths of ABACUS jobs, should contain INPUT, STRU, or KPT, and pseudopotential and orbital files')
        parser.add_argument('-s', "--stepsize", default=0.01, type=float, help="The step size for finite difference calculation, default is 0.01 Angstrom, which should be a small positive value.")
        parser.add_argument('-i', "--index", default=None, nargs="*", type=int, help="The indices of atoms to calculate the vibration frequency, starts from 1. You can specify multiple indices separated by space.")
        parser.add_argument("--image", type=str, default=RECOMMAND_IMAGE, help="The image to use for the Bohrium job, default is %s" % RECOMMAND_IMAGE)
        parser.add_argument("--machine", type=str, default=RECOMMAND_MACHINE, help="The machine to use for the Bohrium job, default is 'c32_m64_cpu'.")
        parser.add_argument("--abacus_command", type=str, default=RECOMMAND_COMMAND, help=f"The command to run the Abacus job, default is '{RECOMMAND_COMMAND}'.")
        return parser
    
    
    def run_prepare(self,params):
        '''
        Run the model with the given parameters
        '''
        if not params.job:
            raise ValueError("No job specified, please use -j or --job to specify the job paths.")
        
        if params.stepsize <= 0:
            raise ValueError("Step size should be a positive value.")

        folders = []
        for i, job_path in enumerate(params.job):
            input_params = ReadInput(os.path.join(job_path, "INPUT"))
            stru_file = input_params.get('stru_file', "STRU")
            stru = AbacusStru.ReadStru(os.path.join(job_path, stru_file))
            input_params['calculation'] = 'scf'
            input_params['cal_force'] = 1

            vib_work_path = os.path.join(job_path, "vib")

            if os.path.exists(vib_work_path) and os.path.isdir(vib_work_path):
                shutil.rmtree(vib_work_path) # Clear the existed job folder

            disped_stru_job_paths = prepare_abacus_vibration_analysis(input_params,
                                                                      stru,
                                                                      work_path=os.path.join(job_path, "vib"),
                                                                      abacus_inputs_dir=Path(job_path),
                                                                      selected_atoms=params.index,
                                                                      stepsize=params.stepsize)

            folders.append(disped_stru_job_paths)
        
        run_sh = f"""#!/bin/bash
{params.abacus_command}
"""        
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
        print(f"After finishing the calculations, you can enter the results directory, \nand run below command below to postprocess the vibration results:\n\tabacustest model vibration post -j {' '.join(params.job)}\n")
        
    
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand 
        '''
        parser.description = "Postprocess the vibration frequency calculation results."
        parser.add_argument('-j', '--job', default=[], action="extend", nargs="*", help='the paths of the job directories, should contain the results of deformed structures generated by the prepare step.')
        parser.add_argument('-t', '--temperature', default=298.15, type=float, help="The temperature of the calculation, default is 298.15 K.")
        return parser   

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        if not params.job:
            raise ValueError("No job specified, please use -j or --job to specify the job paths.")
        
        results = {}
        for job_path in params.job:
            result = post_abacus_vibration_analysis_onejob(job_path, params.temperature)
            results[job_path] = result
        
        # Save the vibration analysis results
        json.dump(results, open("metrics_vibration.json", "w"), indent=4)
        print("\nThe vibration analysis results are saved in 'metrics_vibration.json'")


def copy_abacusjob(src_dir: str, dst_dir: str, input_file=True, stru=True, kpt=True, pp=True, orb=True, out_dir=False, link_pp_orb=True):
    """
    Copy ABACUS job files from src_dir to dst_dir.
    """
    os.makedirs(dst_dir, exist_ok=True)
    input_params = ReadInput(os.path.join(src_dir, "INPUT"))
    if input_file:
        shutil.copy2(os.path.join(src_dir, "INPUT"), os.path.join(dst_dir, "INPUT"))
    if stru:
        stru_file = input_params.get('stru_file', "STRU")
        shutil.copy2(os.path.join(src_dir, stru_file), os.path.join(dst_dir, stru_file))
    if kpt:
        kpt_file = input_params.get('kpt_file', "KPT")
        try:
            shutil.copy2(os.path.join(src_dir, kpt_file), os.path.join(dst_dir, kpt_file))
        except FileNotFoundError:
            pass
    if pp:
        if link_pp_orb:
            pp_files = glob.glob(os.path.join(src_dir, "*.upf")) + glob.glob(os.path.join(src_dir, "*.UPF"))
            for pp_file in pp_files:
                dst_path = os.path.join(dst_dir, os.path.basename(pp_file))
                if os.path.islink(dst_path) or os.path.exists(dst_path):
                    os.remove(dst_path)
                os.symlink(os.path.abspath(pp_file), dst_path)
        else:
            pp_files = glob.glob(os.path.join(src_dir, "*.upf")) + glob.glob(os.path.join(src_dir, "*.UPF"))
            for pp_file in pp_files:
                shutil.copy2(pp_file, os.path.join(dst_dir, os.path.basename(pp_file)))
    if orb:
        if link_pp_orb:
            orb_files = glob.glob(os.path.join(src_dir, "*.orb"))
            for orb_file in orb_files:
                dst_path = os.path.join(dst_dir, os.path.basename(orb_file))
                if os.path.islink(dst_path) or os.path.exists(dst_path):
                    os.remove(dst_path)
                os.symlink(os.path.abspath(orb_file), dst_path)
    if out_dir:
        out_dir_path = os.path.join(src_dir, f"OUT.{input_params.get('suffix', 'ABACUS')}")
        if os.path.exists(out_dir_path) and os.path.isdir(out_dir_path):
            target_out_dir = os.path.join(dst_dir, os.path.basename(out_dir_path))
            shutil.copytree(out_dir_path, target_out_dir, dirs_exist_ok=True)

def prepare_abacus_vibration_analysis(input_params: Dict[str, Any],
                                      stru: AbacusStru,
                                      work_path: Path,
                                      abacus_inputs_dir: Path,
                                      selected_atoms: List[int] = None,
                                      stepsize: float = 0.01,
):
    """
    Prepare ABACUS input files for vibrational analysis.
    """
    stru_file = input_params.get('stru_file', "STRU")
    displaced_stru = copy.deepcopy(stru)
    original_stru_coord = np.array(stru.get_coord(bohr=False, direct=False))

    if selected_atoms is None:
        selected_atoms = [i for i in range(stru.get_natoms())]
    else:
        selected_atoms = [i-1 for i in selected_atoms] # Use 0-based index
        selected_atoms.sort()
    
    DIRECTION_MAP = ['x', 'y', 'z']
    STEP_MAP = {'+': 1, '-': -1}
    job_paths = []

    # Prepare ABACUS input files for the given structure
    abacus_scf_work_path = os.path.join(work_path, "SCF")
    original_stru_job_path = os.path.join(abacus_scf_work_path, "eq")
    os.makedirs(original_stru_job_path)

    copy_abacusjob(abacus_inputs_dir, original_stru_job_path, input_file=False)
    WriteInput(input_params, os.path.join(original_stru_job_path, "INPUT"))
    job_paths.append(original_stru_job_path)

    # Dump selected atoms and stepsize info to a json file used in postprocess
    prepare_params = {
        "selected_atoms": selected_atoms,
        "stepsize": stepsize
    }
    with open(os.path.join(original_stru_job_path, "prepare_params.json"), "w") as f:
        json.dump(prepare_params, f, indent=2)

    # Prepare ABACUS input files for each displaced structure. nfree is assumed to be 2.
    for selected_atom in selected_atoms:
        for direction in range(3): # x, y and z directions
            displaced_stru_coord = copy.deepcopy(original_stru_coord)
            for step in STEP_MAP.keys(): # Two steps along one direction
                disped_stru_job_path = os.path.join(abacus_scf_work_path, f"disp_{selected_atom+1}_{DIRECTION_MAP[direction]}{step}")
                copy_abacusjob(abacus_inputs_dir, disped_stru_job_path, input_file=False, stru=False)

                displaced_stru_coord[selected_atom][direction] = original_stru_coord[selected_atom][direction] + stepsize * STEP_MAP[step]
                displaced_stru.set_coord(displaced_stru_coord, bohr=False, direct=False)
                WriteInput(input_params, os.path.join(disped_stru_job_path, "INPUT"))
                displaced_stru.write(os.path.join(disped_stru_job_path, stru_file))

                job_paths.append(disped_stru_job_path)

    return job_paths

def identify_complex_types(complex_array):
    real_part = np.real(complex_array)
    imag_part = np.imag(complex_array)

    is_real = np.isclose(imag_part, 0)
    is_pure_imag = np.isclose(real_part, 0) & ~np.isclose(imag_part, 0)
    is_general = ~is_real & ~is_pure_imag

    return is_real, is_pure_imag, is_general

def collect_force(abacusjob_dir):
    from abacustest.lib_collectdata.collectdata import RESULT
    
    abacusresult = RESULT(fmt="abacus", path=abacusjob_dir)
    metrics = {i: abacusresult[i] for i in ['force', 'normal_end', 'converge']}
    if metrics['normal_end'] is not True:
        print(f"ABACUS calculation in {abacusjob_dir} didn't end normally")
    elif metrics['converge'] is not True:
        print(f"ABACUS calculation in {abacusjob_dir} didn't reached SCF convergence")
    else:
        pass

    return metrics['force']

def dump_cache_forces_json(stru: AbacusSTRU,
                           work_dir: Path,
                           disped_stru_job_paths: List[str]
):
    """
    Dump cache forces json files for ase.vibration module needed.
    """
    vib_cache_dir = os.path.join(work_dir, "vib/cache")
    os.makedirs(vib_cache_dir, exist_ok=True)
    cache_forces_json = {"forces": {
        "__ndarray__": [[stru.get_natoms(), 3],
                        "float64",
                        []]
    }}
    for disped_stru_job_path in disped_stru_job_paths:
        cache_forces_json["forces"]["__ndarray__"][2] = collect_force(disped_stru_job_path)
        disped_stru_job_basename = os.path.basename(disped_stru_job_path)
        if os.path.basename(disped_stru_job_path) == "eq":
            cache_file = "eq"
        else:
            # get filename from job path
            infos = disped_stru_job_basename.split("_")
            atom_index, step = int(infos[1]) - 1, infos[2]
            cache_file = f"cache.{atom_index}{step}.json"
        with open(os.path.join(vib_cache_dir, cache_file), "w") as fin:
            json.dump(cache_forces_json, fin, indent=2)
    
    return vib_cache_dir

def post_abacus_vibration_analysis_onejob(work_dir: Path,
                                          temperature: float = 298.15):
    """
    Post-process ABACUS vibration analysis results for one job.
    """
    # Read parameters from prepare process
    with open(os.path.join(work_dir, "vib/SCF/eq/prepare_params.json"), "r") as fin:
        prepare_params = json.load(fin)
    
    # Save cached json file as ase.vibration needed
    stru = AbacusStru.ReadStru(os.path.join(work_dir, "vib/SCF/eq/STRU"))
    disped_stru_job_paths = glob.glob(os.path.join(work_dir, "vib/SCF/*"))
    vib_cache_dir = dump_cache_forces_json(stru, work_dir, disped_stru_job_paths)
    
    vib = Vibrations(stru.to_ase(),
                     name=vib_cache_dir,
                     indices=prepare_params['selected_atoms'],
                     delta=prepare_params['stepsize'],
                     nfree=2)
    
    # Do the vibration analysis
    vib.summary()
    # Generate list of frequencies in the return value
    frequencies = vib.get_frequencies()
    real_freq_mask, imag_freq_mask, complex_freq_mask = identify_complex_types(frequencies)
    real_freq, imag_freq = np.real(frequencies[real_freq_mask]).tolist(), frequencies[imag_freq_mask].tolist()
    for key, value in enumerate(imag_freq):
        imag_freq[key] = -math.fabs(value.imag) # Represent imaginary frequency with negative number
    freqs = imag_freq + real_freq

    # Write animations of normal modes in ASE traj format
    vib.write_mode()

    # Thermochemistry calculations
    # Vibrations.get_energies() gets `h \nu` for each mode, which is from the eigenvalues of force constant
    # matrix. The force constant matrix should be a real symmetric matrix mathematically, but due to numerical
    # errors during calculating its matrix element, it will deviate from symmetric matrix slightly, and its eigenvalue
    # will have quite small imaginary parts. Magnitude of imaginary parts will decrease as accuracy of force
    # increases, and it's safe to use norm of the complex eigenvalue as vibration energy if force is 
    # accurate enough.
    vib_energies = vib.get_energies()
    vib_energies_float = [float(np.linalg.norm(i)) for i in vib_energies]
    zero_point_energy = sum(vib_energies_float) / 2
    thermo = HarmonicThermo(vib_energies, ignore_imag_modes=True)
    entropy = thermo.get_entropy(temperature)
    free_energy = thermo.get_helmholtz_energy(temperature)

    return {'frequencies': freqs,
            'zero_point_energy': float(zero_point_energy),
            'vib_entropy': float(entropy),
            'vib_free_energy': float(free_energy)}
