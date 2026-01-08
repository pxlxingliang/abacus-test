from ..model import Model
import json, os
from abacustest.lib_prepare.abacus import WriteKpt, WriteInput, ReadInput, AbacusStru, ReadKpt
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

from typing import List, Dict, Any

import shutil, glob, copy
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.lib_model.comm import clean_none_list


class BECModel(Model):
    @staticmethod
    def model_name(): # type: ignore
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "bec"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "Calculate the Born effective charge by finite difference method."
    
    @staticmethod
    def prepare_args(parser):
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand 
        '''
        parser.description = "Prepare the inputs for Born effective charge calculation."
        parser.add_argument('-j', '--job', default=[], action="extend", nargs="*", help='the paths of ABACUS jobs, should contain INPUT, STRU, or KPT, and pseudopotential and orbital files')
        parser.add_argument("--stepsize", default=0.01, type=float, help="The step size for finite difference calculation, default is 0.01 Angstrom, which should be a small positive value.")
        parser.add_argument("--index", default=[0], nargs="*", type=int, help="The indices of atoms to calculate the Born effective charge, default is 0 (the first atom). You can specify multiple indices separated by space.")
        parser.add_argument("--dir", default=["x", "y", "z"], nargs="*", choices=["x", "y", "z"],type=str, help="The directions to displace the atoms, can be x, y, or z. Default is x y z. You can specify multiple directions separated by space.")
        parser.add_argument("--type", type=str, default="f", choices=["f", "b", "c"], help="The type of displacement, f: forward displacement, b: backward displacement, c: central difference (both forward and backward displacements). Default is f.")
        #parser.add_argument("--no-k-continuity", action="store_true", help="Whether to not set use_k_continuity in nscf calculations. Default is set. This parameter is used for ABACUS >= 3.10.1")
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

        folders = prepare_bec_jobs(params.job,
                                   params.stepsize,
                                   params.index,
                                   params.dir,
                                   params.type,
                                   k_continuity=False)
        # write run scripts
        with open("run.sh", "w") as f1:
            f1.write(f"cp INPUT.scf INPUT\ncp KPT.scf KPT\n{params.abacus_command} | tee scf.log\n\n")
            for i in range(1,4):
                f1.write(f"cp INPUT.nscf{i} INPUT\ncp KPT.nscf{i} KPT\n{params.abacus_command} | tee nscf{i}.log\n")
                f1.write(f"mv OUT.ABACUS/running_nscf.log OUT.ABACUS/running_nscf{i}.log\nrm OUT.ABACUS/*.restart\n")
        
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
        print(f"After finishing the calculations, you can enter the results directory, \nand run below command below to postprocess the BEC results:\n\tabacustest model bec post -j {' '.join(params.job)}\n")
        
    
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand 
        '''
        parser.description = "Postprocess the elastic calculation results."
        parser.add_argument('-j', '--job', default=[], action="extend", nargs="*", help='the paths of the job directories, should contain the results of deformed structures generated by the prepare step.')
        return parser   

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        if not params.job:
            raise ValueError("No job specified, please use -j or --job to specify the job paths.")
        
        metrics, becs = postprocess_bec(params.job)
        json.dump(becs, open("metrics_bec.json", "w"), indent=4)
        json.dump(metrics, open("metrics.json", "w"), indent=4)
        print("\nThe BEC results are saved in 'metrics_bec.json'")

        import pandas as pd
        c = "Born effective charge tensor for:"
        for k, v in becs.items():
            c += f"\n{k}:\n"
            bect, index = clean_none_list(v["bec_tensor (e)"], ["X", "Y", "Z"])
            pd_df = pd.DataFrame(bect, columns=['X', 'Y', 'Z'], index=index)
            # save 4 decimal places
            pd_df = pd_df.round(4)
            c += f"{pd_df}\n"

        c += """
NOTE: The BEC tensor is in unit of e, calculated by delta P / delta R, 
where P is polarization vector in cartesian coordinates (e/Volume * A), 
and R is atomic position in cartesian coordinates (A). 
Displacement is the actual atomic displacement in cartesian coordinates (A).
The p_vec_disp and p_vec_org are polarization vectors along the three cell vectors.
The p1(C/m^2) is the polarization in C/m^2 along the three cell vectors, and p2(C/m^2) is the polarization in C/m^2 in the cartesian coordinates.
        """
        with open("bec_summary.txt", "w") as f1:
            f1.write(c)
        print("\nSummary of Born effective charge tensors:")
        print(c)
            
        
        
def prepare_bec_jobs(jobs: str,
                     stepsize: float,
                     index: List[int],
                     directions: List[str],
                     disp_type: str = "f",
                     k_continuity: bool = True
                    ) -> None:
    '''
    Prepare the BEC calculation jobs.
    
    Args:
        jobs: List of job paths.
        stepsize: Step size for finite difference.
        index: List of atom indices to displace.
        directions: List of directions to displace ('x', 'y', 'z').
        disp_type: Type of displacement, 'f' for forward, 'b' for backward, 'c' for central difference.
    '''
    if disp_type not in ["f", "b", "c"]:
        raise ValueError(f"Invalid displacement type: {disp_type}. Should be 'f', 'b', or 'c'.")

    if stepsize == 0:
        raise ValueError("Step size should be a non-zero value.")
    if stepsize < 0:
        stepsize = -stepsize
    
    directions = list(set(directions))
    directions_map = {"x": [stepsize, 0.0, 0.0],
                      "y": [0.0, stepsize, 0.0],
                      "z": [0.0, 0.0, stepsize]}
    if set(directions) - set(directions_map.keys()):
        raise ValueError(f"Invalid directions specified: {directions}. Valid options are 'x', 'y', 'z'.")

    folders = []
    for job in jobs:
        # clean org folders
        print("Preparing BEC calculation for job path:", job)
        old_folders = [f for f in glob.glob(os.path.join(job, "bec_*")) if os.path.isdir(f)]
        if len(old_folders) > 0:
            print(f"Remove old bec_* folders in job path: {job}")
        for f in old_folders:
            shutil.rmtree(f)
        
        if not os.path.isfile(os.path.join(job, "STRU")):
            raise FileNotFoundError(f"STRU file not found in job path: {job}")
        if not os.path.isfile(os.path.join(job, "INPUT")):
            raise FileNotFoundError(f"INPUT file not found in job path: {job}")
        
        org_files = [f for f in os.listdir(job) if f not in ["INPUT", "STRU", "KPT"]]
        input_param = ReadInput(os.path.join(job, "INPUT"))
        stru = AbacusStru.ReadStru(os.path.join(job, "STRU"))
        kpt, kpt_model = ReadKpt(job)
        input_param.pop("kspacing", None) # we should remove kspacing, but use KPT instead
        input_param["suffix"] = "ABACUS"
        # if SPIN1_CHG.cube exists, use it as initial charge density
        if os.path.isfile(os.path.join(job, "OUT.ABACUS/SPIN1_CHG.cube")):
            input_param["init_chg"] = "file"
        
        if disp_type != "c":
            org_path = os.path.join(job, f"bec_org")
            os.makedirs(org_path, exist_ok=True)
            for f in org_files:
                os.symlink(os.path.abspath(os.path.join(job, f)), os.path.join(org_path, f))
            write_inputs(org_path, stru, input_param, kpt, kpt_model, k_continuity=k_continuity)
            folders.append(org_path)
        

        for direction in directions:           
            # write displaced structure folders
            for idx in index:
                disp_vec = directions_map[direction]

                if disp_type in ["f", "c"]:
                    disp_stru = copy.deepcopy(stru)
                    atom_coord = disp_stru.get_coord(bohr=False, direct=False)
                    atom_coord[idx] = [atom_coord[idx][i] + disp_vec[i] for i in range(3)]
                    disp_stru.set_coord(atom_coord, bohr=False, direct=False)
                    disp_path = os.path.join(job, f"bec_disp_atom{idx}_{direction}")
                    os.makedirs(disp_path, exist_ok=True)
                    for f in org_files:
                        os.symlink(os.path.abspath(os.path.join(job, f)), os.path.join(disp_path, f))
                    write_inputs(disp_path, disp_stru, input_param, kpt, kpt_model, k_continuity=k_continuity)
                    folders.append(disp_path)
                
                if disp_type in ["b", "c"]:
                    disp_stru = copy.deepcopy(stru)
                    atom_coord = disp_stru.get_coord(bohr=False, direct=False)
                    atom_coord[idx] = [atom_coord[idx][i] - disp_vec[i] for i in range(3)]
                    disp_stru.set_coord(atom_coord, bohr=False, direct=False)
                    disp_path = os.path.join(job, f"bec_disp_atom{idx}_{direction}_back")
                    os.makedirs(disp_path, exist_ok=True)
                    for f in org_files:
                        os.symlink(os.path.abspath(os.path.join(job, f)), os.path.join(disp_path, f))
                    write_inputs(disp_path, disp_stru, input_param, kpt, kpt_model, k_continuity=k_continuity)
                    folders.append(disp_path)

    return folders
            
            
def write_inputs(folder, stru, input_param, kpt, kpt_model, k_continuity: bool = True):
    # 1. write INPUT
    input_param_copy = copy.deepcopy(input_param)
    input_param_copy["calculation"] = "scf"
    input_param_copy["out_chg"] = 1
    input_param_copy["out_bandgap"] = 1
    WriteInput(input_param_copy, os.path.join(folder, "INPUT.scf"))
    
    # 2. write KPT
    WriteKpt(kpt, os.path.join(folder, "KPT.scf"),kpt_model)
    
    # 3. write INPUT nscf
    input_param_copy["calculation"] = "nscf"
    input_param_copy["out_chg"] = None
    input_param_copy["init_chg"] = "file"
    input_param_copy["berry_phase"] = 1
    if k_continuity:
        input_param_copy["use_k_continuity"] = True
    for i in range(1,4):
        input_param_copy["gdir"] = i
        WriteInput(input_param_copy, os.path.join(folder, f"INPUT.nscf{i}"))
    
    # 4. write KPT nscf
    for i in range(3):
        kpt_nscf = copy.deepcopy(kpt)
        kpt_nscf[i] = 2 * kpt[i] # double the k-point grid in the displacement direction
        WriteKpt(kpt_nscf, os.path.join(folder, f"KPT.nscf{i+1}"),kpt_model)
    
    # 5. write STRU
    stru.write(os.path.join(folder, "STRU"))

    
def postprocess_bec(jobs: List[str]) -> Dict[str, Any]:
    '''
    Postprocess the BEC calculation results.
    
    Args:
        jobs: List of job paths.
    
    Returns:
        A dictionary containing the BEC results.
    '''
    xyz_idx = {"x":0, "y":1, "z":2}
    metrics = {}
    becs = {}
    for job in jobs:
        print("Postprocessing BEC calculation for job path:", job)
        if not os.path.isdir(job):
            print(f"Warning: job path {job} is not a directory, skip.")
            continue

        sub_folders = [f for f in glob.glob(os.path.join(job, "bec_*")) if os.path.isdir(f)]
        sub_folders.sort()
        label = AbacusStru.ReadStru(os.path.join(sub_folders[0], "STRU")).get_label(total=True)
        
        for folder in sub_folders:
            metrics[folder] = read_metrics(folder)
            stru = AbacusStru.ReadStru(os.path.join(folder, "STRU"))
            atom_coord = stru.get_coord(bohr=False, direct=False)
            cell = stru.get_cell(bohr=False)
            p_vec1, mod1, volume1, p_cm2_1 = read_p(os.path.join(folder, "OUT.ABACUS/running_nscf1.log"))
            p_vec2, mod2, volume2, p_cm2_2 = read_p(os.path.join(folder, "OUT.ABACUS/running_nscf2.log"))
            p_vec3, mod3, volume3, p_cm2_3 = read_p(os.path.join(folder, "OUT.ABACUS/running_nscf3.log"))
            # check volume consistency
                                           
            metrics[folder]["atom_coord"] = atom_coord
            metrics[folder]["cell"] = cell
            metrics[folder]["p_vec"] =  [p_vec1, p_vec2, p_vec3]  # p_vec1/2/3 is along cell vectors
            metrics[folder]["p1(C/m^2)"] = [p_cm2_1, p_cm2_2, p_cm2_3] # along a/b/c directions
            metrics[folder]["p2(C/m^2)"] = cal_polarization_cartesian([p_cm2_1, p_cm2_2, p_cm2_3], cell)  # in cartesian coordinates
            metrics[folder]["volume"] = volume1
            metrics[folder]["mod"] = [mod1, mod2, mod3]
            
        
        for folder in sub_folders:
            # bec_disp_atom{idx}_{direction}
            if "bec_disp" not in os.path.basename(folder):
                continue
            parts = os.path.basename(folder).split("_")
            idx = int(parts[2][4:])
            direction = parts[3]

            org_folder = os.path.join(job, f"bec_org")
            # find the org_folder
            if folder.endswith("_back") and folder[:-5] in sub_folders:
                continue  # will be handled in the positive displacement
            elif (not folder.endswith("_back")) and (folder + "_back" in sub_folders):
                org_folder = folder + "_back"

            if org_folder not in metrics:
                print(f"Warning: folder {org_folder} not found for displaced folder {folder}. Skipping.")
                continue
            
            p_vec_disp = metrics[folder]["p_vec"]
            p_vec_org = metrics[org_folder]["p_vec"]
            if p_vec_disp is None or p_vec_org is None:
                print(f"Warning: Polarization vector not found for folder {folder} or {org_folder}. Skipping.")
                continue
            
            coord_disp = metrics[folder]["atom_coord"][idx]
            coord_org = metrics[org_folder]["atom_coord"][idx]
            delta_r_atom = [(coord_disp[i] - coord_org[i]) for i in range(3)]
            disp_idx = xyz_idx[direction]
            delta_r_norm = sum([delta_r_atom[i]**2 for i in range(3)])**0.5 * (1 if delta_r_atom[disp_idx] >=0 else -1)
            delta_p = cal_delta_p(p_vec_org, p_vec_disp, metrics[org_folder]["mod"], metrics[folder]["mod"])
            delta_p = cal_polarization_cartesian(delta_p, metrics[org_folder]["cell"])  # transform to cartesian
            bec_tensor = [delta_p[i] / delta_r_norm for i in range(3)]
            
            ilabel = f"{job}:atom{idx}_{label[idx]}"

            if ilabel not in becs:
                becs[ilabel] = {
                    "bec_tensor (e)": [None, None, None],
                    "displacement (A)": [None, None, None],
                    "p_vec_org ((e/Volume)*A)": p_vec_org,
                    "mod_org ((e/Volume)*A)": metrics[org_folder]["mod"],
                    "p_vec_disp ((e/Volume)*A)": [None, None, None]
                }

            becs[ilabel]["bec_tensor (e)"][xyz_idx[direction]] = bec_tensor
            becs[ilabel]["displacement (A)"][xyz_idx[direction]] = delta_r_atom
            becs[ilabel]["p_vec_disp ((e/Volume)*A)"][xyz_idx[direction]] = p_vec_disp
    
    return metrics, becs

def cal_polarization_cartesian(p_vec, cell):
    '''
    Convert polarization vector along cell vectors to cartesian coordinates.
    
    Args:
        p_vec: Polarization vector in cell vector basis.
        cell: Cell vectors.
    '''
    import numpy as np
    cell_matrix = np.array(cell)
    p_cartesian = np.sum([p_vec[i] * cell_matrix[i] / np.linalg.norm(cell_matrix[i]) for i in range(3)], axis=0)

    return p_cartesian.tolist()

def cal_delta_p(p1, p2, mod1, mod2):
    '''
    Calculate the delta p considering the polarization quantum.
    
    Args:
        p1: Polarization vector 1.
        p2: Polarization vector 2.
        mod1: Modulus for polarization vector 1.
        mod2: Modulus for polarization vector 2.
    '''
    delta_p = [p2[i] - p1[i] for i in range(3)]
    for i in range(3):
        if delta_p[i] > 0:
            if abs(delta_p[i] - mod2[i]) < abs(delta_p[i]):
                delta_p[i] = delta_p[i] - mod2[i]
        else:
            if abs(delta_p[i] + mod2[i]) < abs(delta_p[i]):
                delta_p[i] = delta_p[i] + mod2[i]
    return delta_p

def read_metrics(job: str) -> List[str]:
    '''
    Read the available metrics from the job results.
    
    Args:
        job: Path to the job directory.
        
    Returns:
        A list of available metric names.
    '''
    shutil.copy(os.path.join(job, "INPUT.scf"), os.path.join(job, "INPUT"))
    r = RESULT(path=job, fmt="abacus",output=os.path.join(job, "scf.log"))
    keys = ["normal_end", "converge","energy", "scf_steps", "denergy_last", "drho_last", "band_gap"]
    return {k: r[k] for k in keys}

def read_p(logf):
    """Read polarization from log file.
    
    Args:
        logf: Path to the log file.
        
    Returns:
        p_vec: List of polarization vector components (P * V, with unit in e*A)
    
    """
    """
        The calculated polarization direction is in R1 direction

 P =    4.6793707  (mod    4.7318739)  (   4.6793707,   0.0000000,   0.0000000) (e/Omega).bohr

 P =    0.0212834  (mod    0.0215222)  (   0.0212834,   0.0000000,   0.0000000) e/bohr^2

 P =    1.2168069  (mod    1.2304597)  (   1.2168069,   0.0000000,   0.0000000) C/m^2
    """
    from abacustest.constant import BOHR2A
    with open(logf) as f1: lines = f1.readlines()  
    
    p_vec = None
    volume = None
    mod = None
    p_cm2 = None
    for i, line in enumerate(lines):
        if "The calculated polarization direction is " in line:
            p_vec = float(lines[i+2].split()[2]) * BOHR2A
            mod = float(lines[i+2].split("mod")[1].split(")")[0]) * BOHR2A
            p_cm2 = float(lines[i+6].split()[2])  # in C/m^2
        elif "Volume (A^3) =" in line:
            volume = float(line.split("=")[1])
    return p_vec, mod, volume, p_cm2        
                    
