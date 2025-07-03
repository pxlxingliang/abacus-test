from ..model import Model
from . import comm
import argparse,json, os
from abacustest.lib_prepare.abacus import WriteKpt, WriteInput, gen_stru, ReadInput
from pathlib import Path
from abacustest.lib_model.model_012_band import PrepBand
from abacustest.constant import RECOMMAND_IMAGE

from typing import List, Dict, Any
from pathlib import Path


JOB_TYPES = {"scf": {"calculation": "scf", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                    "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                    "mixing_beta": 0.7,  "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                    "#cal_force": 1, "#cal_stress": 1,
                    "kspacing": "0.1 # unit in 1/bohr"}, 
             "relax": {"calculation": "relax", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                       "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                       "mixing_beta": 0.7, "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                       "cal_force": 1, "#cal_stress": 1,"kspacing": "0.1 # unit in 1/bohr",
                       "relax_method": "cg # or bfgs bfgs_trad cg_bfgs sd fire",
                       "relax_nmax": 60, "force_thr_ev": "0.01  # unit in eV/A", "#stress_thr": "0.5 # unit in kbar",
                       }, 
             "cell-relax":{"calculation": "cell-relax", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                       "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                       "mixing_beta": 0.7, "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                       "cal_force": 1, "cal_stress": 1, "kspacing": "0.1 # unit in 1/bohr",
                       "relax_method": "cg # or bfgs, bfgs_trad, cg_bfgs, sd, fire",
                       "relax_nmax": 60, "force_thr_ev": "0.01  # unit in eV/A", "stress_thr": "0.5 # unit in kbar",
                       "fixed_axes": "None # or volume, shape, a, b, c, ab, ac, bc; only valid for cell-relax calculation to fix some axes",
                       }, 
             "md":{"calculation": "md", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                       "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                       "mixing_beta": 0.7, "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                       "#cal_force": 1, "#cal_stress": 1,
                       "kspacing": "0.1 # unit in 1/bohr", 
                       "md_type": "nvt  # or npt, nve, langevin, fire, msst",
                       "md_nstep": "10  # number of steps",
                       "md_dt": "1.0  # unit in fs", 
                       "md_tfirst": "100  # unit in K",
                       "md_tlast": "100  # unit in K",}, 
             "band":{"calculation": "scf", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                    "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                    "mixing_beta": 0.7,  "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                    "#cal_force": 1, "#cal_stress": 1,
                    "kspacing": "0.1 # unit in 1/bohr"}
        }

LCAO_PARAM = {
    "basis_type": "lcao",
    "ks_solver": "genelpa",
    "#gamma_only": 0,
    "ecutwfc": 100,
    "scf_thr": 1e-7,
}


class InputsModel(Model):
    @staticmethod
    def model_name(): # type: ignore
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "inputs"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "Prepare the ABACUS inputs of specified model"
    
    @staticmethod
    def add_args(parser):
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand '''
        
        parser.description = "Prepare the ABACUS inputs file."
        parser.add_argument('-f', '--file',default=[], action="extend",nargs="*" ,help='the structure files')
        parser.add_argument("--ftype",default="cif",type=str,help="the structure type, should be cif or dpdata supportted type like: poscar, abacus/stru, ..", )
        parser.add_argument("--jtype", default="scf", type=str, help=f"the job type, should be one of {list(JOB_TYPES.keys())}")
        parser.add_argument("--pp",default=None,type=str,help="the path of pseudopotential library, or read from enviroment variable ABACUS_PP_PATH")
        parser.add_argument("--orb",default=None,type=str,help="the path of orbital library, or read from enviroment variable ABACUS_ORB_PATH")
        parser.add_argument("--input",default=None,type=str,help="the template of input file, if not specified, the default input will be generated")
        parser.add_argument("--kpt", default=None, type=int, nargs="*", help="the kpoint setting, should be one or three integers")
        parser.add_argument("--lcao", action="store_true", help="whether to use lcao basis, default is pw basis")
        return parser
    
    def run(self,params):
        if params.kpt is not None and len(params.kpt) not in [1, 3]:
            raise ValueError("The kpoint setting should be one or three integers.")
        pinput = PrepInput(
            files=params.file,
            filetype=params.ftype,
            jobtype=params.jtype,
            pp_path=params.pp,
            orb_path=params.orb,
            input_file=params.input,
            kpt=params.kpt,
            lcao=params.lcao,
        )
        pinput.run()
        return 0

class PrepInput:
    def __init__(self, files, filetype, jobtype, pp_path=None, orb_path=None, input_file=None, kpt=None,
                 abacus_command="OMP_NUM_THREADS=1 mpirun -np 16 abacus", machine="c32_m64_cpu", 
                 image=RECOMMAND_IMAGE,
                 lcao=False):
        if jobtype not in JOB_TYPES:
            raise ValueError(f"Unsupported job type: {jobtype}.\nSupported job types are {list(JOB_TYPES.keys())}.")
        
        if not isinstance(files, list):
            self.files = [files]
        else:
            self.files = files
            
        self.filetype = filetype
        self.jobtype = jobtype
        self.pp_path = pp_path
        self.orb_path = orb_path
        self.input_file = input_file
        self.kpt = kpt
        self.abacus_command = abacus_command
        self.machine = machine
        self.image = image
        self.lcao = lcao
        
        if self.pp_path is None and os.environ.get("ABACUS_PP_PATH") is not None:
            self.pp_path = os.environ["ABACUS_PP_PATH"]
        if self.orb_path is None and os.environ.get("ABACUS_ORB_PATH") is not None:
            self.orb_path = os.environ["ABACUS_ORB_PATH"]
        
        print("Structure files:", self.files)
        print("Structure type:", self.filetype)
        print("Job type:", self.jobtype)
        print("Pseudopotential path:", self.pp_path)
        print("Orbital path:", self.orb_path)
        print("Input file:", self.input_file)
        print("Kpoint:", self.kpt)
        print("")
    
    def run(self):
        recommanded_param = JOB_TYPES[self.jobtype]
        auto_ecutwfc = True
        if self.lcao:
            recommanded_param.update(LCAO_PARAM)
            auto_ecutwfc = False
        
        if self.input_file is None:
            print("Automatically generate the input file.")
            input_param = recommanded_param
        else:
            input_param = ReadInput(self.input_file)
            if "calculation" in input_param and input_param.get("calculation") != recommanded_param["calculation"]:
                print(f"Warning: the calculation type in the input file is {input_param['calculation']}, but the job type is {self.jobtype}.")
                print(f"         Automatically set the calculation type to {recommanded_param['calculation']}.")
                input_param["calculation"] = recommanded_param["calculation"]
            if "ecutwfc" in input_param:
                auto_ecutwfc = False
            
            # update other required parameters
            for key in recommanded_param:
                if key in input_param or (key.startswith("#") and key[1:] in input_param):
                    continue
                input_param[key] = recommanded_param[key]
            
        jobs = self.gen_abacus_inputs(self.files, self.filetype, 
                                      self.pp_path, self.orb_path, 
                                      input_param, self.kpt,
                                      auto_ecutwfc=auto_ecutwfc)
        job_path = list(jobs.keys())
        if self.jobtype == "band":
            prep_band = PrepBand(job_path, run_command=self.abacus_command)
            _, run_script, extra_files = prep_band.run()
        else:
            run_script_file = "run.sh"
            Path(run_script_file).write_text(self.abacus_command)
            run_script = f"bash {run_script_file}"
            extra_files = [run_script_file]
            

        # generate the parameter setting for abacustest submit
        setting = {
            "save_path": "results",
            "run_dft": {
                "example": job_path,
                "extra_files": extra_files,
                "command": run_script,
                "image": self.image,
                "bohrium": {
                    "scass_type": self.machine,
                    "job_type": "container",
                    "platform": "ali"
                }
            }
        }
        
        setting_file = "setting.json"
        json.dump(setting, open(setting_file, "w"), indent=4)
        print("\nThe inputs are generated in", ", ".join(job_path))
        print(f"Please modify the ABACUS running command in script {extra_files[0]} if needed.")
        print(f"You can modify '{setting_file}', and execute the command 'abacustest submit -p setting.json' to run the abacustest to submit all jobs to bohrium.")
        print(f"Or you can 'cd' to each job path and execute '{extra_files[0]}' to run the job.")
        
        return setting, job_path
    
    @staticmethod
    def gen_abacus_inputs(stru_files, stru_type, pp_path, orb_path, input_param=None, kpoint=None, auto_ecutwfc=True):
        """
        Generate the abacus input files.

        Parameters
        ----------
        stru_files : list
            The structure files, should be a list of cif or dpdata supported data
        stru_type : str
            The structure type, cif or dpdata supported format
        pp_path : str
            The pseudopotential path
        orb_path : str
            The orbital path
        input_param : dict
            The specified input parameter, which will used in INPUT
        kpoint : list
            The kpoint setting, should be a list of three int, will generate a KPT file
        auto_ecutwfc : bool
            Whether to automatically find the ecutwfc from ecutwfc.json file, default is True

        Returns
        -------
        job_path : dict
            The job path, which is a dict, key is the job path, value is {"element": element, "pp": pp, "orb": orb}
        """
        job_path = gen_stru(stru_files, stru_type, pp_path, orb_path, tpath=".")
        # job_path is dict: key is the job path, value is {"element": element, "pp": pp, "orb": orb}
        # the pp and orb file has been linked to the job path

        if job_path is None or len(job_path) == 0:
            raise ValueError("No valid structure file found.")

        # write INPUT file and KPT file

        default_input = {
            "ecutwfc": 100,
            "calculation": "scf",
            "basis_type": "pw"
        }
        if input_param is not None:
            default_input.update(input_param)
        # will find the recommand ecutwfc from ecutwfc.json file
        if pp_path is not None and os.path.isfile(os.path.join(pp_path, "ecutwfc.json")):
            recommand_ecutwfc = json.load(open(os.path.join(pp_path, "ecutwfc.json"), "r"))
        else:
            recommand_ecutwfc = None

        if auto_ecutwfc and recommand_ecutwfc is not None:
            print("Based on recommended ecutwfc, automatically set ecutwfc for each job.")
            
        for path, job in job_path.items():
            element = job["element"]
            if auto_ecutwfc and recommand_ecutwfc is not None:
                ecutwfcs = []
                for e in element:
                    if e in recommand_ecutwfc:
                        ecutwfcs.append(recommand_ecutwfc[e])
                    else:
                        print(f"Warning: No recommended ecutwfc for element {e}.")
                if len(ecutwfcs) == 0:
                    print(f"Warning: No recommended ecutwfc for elements {element}. Using default ecutwfc {default_input['ecutwfc']}.")
                else:      
                    default_input["ecutwfc"] = max(ecutwfcs)
                    recutwfc = {i: recommand_ecutwfc.get(i,None) for i in element}
                    print(f"{path}: set ecutwfc to {default_input['ecutwfc']}, recommended ecutwfc: {recutwfc}")

            if kpoint is not None:
                default_input.pop("kspacing",None)
                WriteKpt(kpoint, os.path.join(path, "KPT"))
            
            WriteInput(default_input, os.path.join(path, "INPUT"))

        return job_path