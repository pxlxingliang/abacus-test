from ..model import Model
from . import comm
import argparse,json, os
from abacustest.lib_prepare.abacus import WriteKpt, WriteInput, gen_stru, ReadInput, AbacusStru
from pathlib import Path
from abacustest.lib_model.model_012_band import PrepBand
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

import warnings

from typing import List, Dict, Union, Optional, Tuple, Literal


JOB_TYPES = {"scf": {"calculation": "scf", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                    "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                    "mixing_beta": 0.8,  "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                    "pw_diag_ndim":    2, "pw_diag_nmax":    20,
                    "precision": "double  # or single",
                    "#cal_force": 1, "#cal_stress": 1,
                    "kspacing": "0.14 # unit in 1/bohr"}, 
             "relax": {"calculation": "relax", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                       "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                       "mixing_beta": 0.8, "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                       "pw_diag_ndim":    2, "pw_diag_nmax":    20,
                       "precision": "double  # or single",
                       "cal_force": 1, "#cal_stress": 1,"kspacing": "0.14 # unit in 1/bohr",
                       "relax_method": "cg # or bfgs bfgs_trad cg_bfgs sd fire",
                       "relax_nmax": 60, "force_thr_ev": "0.01  # unit in eV/A", "#stress_thr": "0.5 # unit in kbar",
                       }, 
             "cell-relax":{"calculation": "cell-relax", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                       "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                       "mixing_beta": 0.8, "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                       "pw_diag_ndim":    2, "pw_diag_nmax":    20,
                       "precision": "double  # or single",
                       "cal_force": 1, "cal_stress": 1, "kspacing": "0.14 # unit in 1/bohr",
                       "relax_method": "cg # or bfgs, bfgs_trad, cg_bfgs, sd, fire",
                       "relax_nmax": 60, "force_thr_ev": "0.01  # unit in eV/A", "stress_thr": "0.5 # unit in kbar",
                       "fixed_axes": "None # or volume, shape, a, b, c, ab, ac, bc; only valid for cell-relax calculation to fix some axes",
                       }, 
             "md":{"calculation": "md", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                       "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                       "mixing_beta": 0.8, "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                       "pw_diag_ndim":    2, "pw_diag_nmax":    20,
                       "precision": "double  # or single",
                       "#cal_force": 1, "#cal_stress": 1,
                       "kspacing": "0.14 # unit in 1/bohr", 
                       "md_type": "nvt  # or npt, nve, langevin, fire, msst",
                       "md_nstep": "10  # number of steps",
                       "md_dt": "1.0  # unit in fs", 
                       "md_tfirst": "100  # unit in K",
                       "md_tlast": "100  # unit in K",}, 
             "band":{"calculation": "scf", "symmetry": 1, "ecutwfc": 80, "scf_thr": 1e-8, "scf_nmax": 100,
                    "smearing_method": "gauss", "smearing_sigma": 0.015, "mixing_type": "broyden",
                    "mixing_beta": 0.8,  "basis_type": "pw  # or lcao", "ks_solver": "dav_subspace  # or genelpa for lcao basis",
                    "pw_diag_ndim":    2, "pw_diag_nmax":    20,
                    "precision": "double  # or single",
                    "#cal_force": 1, "#cal_stress": 1,
                    "kspacing": "0.14 # unit in 1/bohr"}
        }

LCAO_PARAM = {
    "basis_type": "lcao",
    "ks_solver": "genelpa",
    "#gamma_only": 0,
    "ecutwfc": 100,
    "scf_thr": 1e-7,
}

ELEMENT_DFTU_D = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                  "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
                  "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "La", "Ac", "Th"]
ELEMENT_DFTU_F = ["Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
                  "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]

MAG_ELEMENT = ELEMENT_DFTU_D + ELEMENT_DFTU_F


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
        parser.add_argument("--nspin", default=1, type=int, choices=[1, 2, 4], help="the number of spins, can be 1 (no spin), 2 (spin polarized), or 4 (non-collinear spin). Default is 1.")
        parser.add_argument("--soc", action="store_true", help="whether to use spin-orbit coupling, if True, nspin should be 4.")
        parser.add_argument("--dftu", action="store_true", help="whether to use DFT+U, default is False.")
        parser.add_argument("--dftu_param", default=None, nargs="+", help="the DFT+U parameters, should be element symbol and U value pairs like 'Fe 4 Ti 1'. If dftu is set, but dftu_param is not set, the default U values will be used: 4 eV for d orbital elements and 6 eV for f orbital elements.")
        parser.add_argument("--init_mag", default=None, nargs="+", help="the initial magnetic moment for magnetic elements, should be element symbol and magnetic moment pairs like 'Fe 4 Ti 1'.")
        parser.add_argument("--afm", action="store_true", help="whether to use antiferromagnetic calculation, default is False. Only valid when init_mag is set.")
        parser.add_argument("--copy_pp_orb", action="store_true", help="whether to copy the pseudopotential and orbital files to each job directory or link them. Default is False, which means linking the files.")
        parser.add_argument("--download-pporb", nargs="?", type=str,default=None,const="apns-v1",choices=["apns-v1"], help="Download the recommended pseudopotential and orbital files from AISquare.")
        return parser
    
    @staticmethod
    def parse_dftu_param(values):
        """
        Parse the values from the command line arguments.
        """
        if len(values) % 2 != 0:
            raise ValueError("The DFT+U parameters should be in pairs of element symbol and U value.")
        dftu_param = {}
        for i in range(0, len(values), 2):
            element = values[i]
            u_value = float(values[i + 1])
            if element not in dftu_param:
                dftu_param[element] = u_value
            else:
                warnings.warn(f"Element {element} already has a DFT+U value, overwriting it with {u_value}.")
        return dftu_param
     
    @staticmethod
    def parse_init_mag(values):
        """
        Parse the initial magnetic moment values from the command line arguments.
        """
        if len(values) % 2 != 0:
            raise ValueError("The initial magnetic moment parameters should be in pairs of element symbol and magnetic moment.")
        init_mag = {}
        for i in range(0, len(values), 2):
            element = values[i]
            mag_value = float(values[i + 1])
            if element not in init_mag:
                init_mag[element] = mag_value
            else:
                warnings.warn(f"Element {element} already has an initial magnetic moment, overwriting it with {mag_value}.")
        return init_mag   
    
    def download_pporb(self, version):
        """
        Download the recommended pseudopotential and orbital files from AISquare.
        """
        if version is None:
            return
        
        apns_v1 = "https://store.aissquare.com/datasets/af21b5d9-19e6-462f-ada1-532f47f165f2/ABACUS-APNS-PPORBs-v1.zip"
        
        if version == "apns-v1":
            # download the apns-v1 pseudopotential and orbital files
            import urllib.request
            import zipfile
            
            zip_file = "ABACUS-APNS-PPORBs-v1.zip"
            print(f"Downloading the recommended pseudopotential and orbital files from {apns_v1} ...")
            urllib.request.urlretrieve(apns_v1, zip_file)
            print(f"Unzipping the file {zip_file} ...")
            with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                zip_ref.extractall(".")
            os.remove(zip_file)
            print("Download and unzip completed.")
            # set the pp and orb path
        
            print("You can set the pseudopotential and orbital path by --pp and --orb, or set the environment variable ABACUS_PP_PATH and ABACUS_ORB_PATH.")
            
        
    
    def run(self,params):
        if params.download_pporb is not None:
            self.download_pporb(params.download_pporb)
            return 0
        
        if params.kpt is not None and len(params.kpt) not in [1, 3]:
            raise ValueError("The kpoint setting should be one or three integers.")
        
        if params.dftu_param is not None:
            dftu_param = self.parse_dftu_param(params.dftu_param)
        else:
            dftu_param = None
            
        if params.init_mag is not None:
            init_mag = self.parse_init_mag(params.init_mag)
        else:
            init_mag = None
        
        pinput = PrepInput(
            files=params.file,
            filetype=params.ftype,
            jobtype=params.jtype,
            pp_path=params.pp,
            orb_path=params.orb,
            input_file=params.input,
            kpt=params.kpt,
            lcao=params.lcao,
            nspin=params.nspin,
            soc=params.soc,
            dftu=params.dftu,
            dftu_param=dftu_param,
            init_mag=init_mag,
            afm=params.afm,
            abacus_command=RECOMMAND_COMMAND,
            machine=RECOMMAND_MACHINE,  # default Bohrium machine type for CPU jobs
            image=RECOMMAND_IMAGE,  # default recommended image for ABACUS jobs 
            copy_pp_orb=params.copy_pp_orb 
        )
        pinput.run()
        return 0

class PrepInput:
    """Prepare the ABACUS inputs for the specified model.
    This class will generate the ABACUS inputs files based on the specified structure files, job type, and other parameters.
    
    Parameters
    ----------
    files : list
        The structure files, should be a list of cif or dpdata supported data
    filetype : str
        The structure type, cif or dpdata supported format
    pp_path : str
        The pseudopotential path, if not specified, will read from the environment variable ABACUS_PP_PATH
    orb_path : str
        The orbital path, if not specified, will read from the environment variable ABACUS_ORB_PATH
    input_file : str
        The template of input file, if not specified, the default input will be generated
    kpt : list
        The kpoint setting, should be a list of three int, will generate a KPT file
    abacus_command : str
        The command to run ABACUS, default is RECOMMAND_COMMAND
    machine : str
        The machine type, default is RECOMMAND_MACHINE, which is the Bohrium machine type for CPU jobs.
    image : str
        The Docker image to use, default is RECOMMAND_IMAGE, which is the recommended image for ABACUS jobs.
    lcao : bool
        Whether to use LCAO basis, default is False, which means using PW basis.
    nspin : int
        The number of spins, can be 1 (no spin), 2 (spin polarized), or 4 (non-collinear spin). Default is 1.
    soc : bool
        Whether to use spin-orbit coupling, if True, nspin should be 4.
    dftu : bool
        Whether to use DFT+U, default is False.
    dftu_param : dict
        The DFT+U parameters, should be a dict like {"Fe": 4, "Ti": 1}, where the key is the element symbol and the value is the U value.
        Value can also be a list of two values, and the first value is the orbital (p, d, f) to apply DFT+U, and the second value is the U value.
        For example, {"Fe": ["d", 4], "O": ["p", 1]} means applying DFT+U to Fe 3d orbital with U=4 eV and O 2p orbital with U=1 eV.
        If dftu is True, but dftu_param is None, the default U values will be used: 4 eV for d orbital elements and 6 eV for f orbital elements.
    init_mag : dict or None
        the initial magnetic moment for magnetic elements, should be a dict like {"Fe": 4, "Ti": 1}, where the key is the element symbol and the value is the initial magnetic moment.
    afm : bool
        Whether to use antiferromagnetic calculation, default is False. If True, half of the magnetic elements will be set to negative initial magnetic moment.
    copy_pp_orb : bool
        Whether to copy the pseudopotential and orbital files to each job directory or link them. Default is False, which means linking the files.
    """
    
    def __init__(self, files: Union[str, List[str], Path, List[Path]],
                 filetype: str, 
                 jobtype: str, 
                 pp_path: Optional[Union[str, Path]]=None, 
                 orb_path: Optional[Union[str, Path]]=None, 
                 input_file: Optional[Union[str, Path]]=None, 
                 kpt: Optional[Tuple[int, int, int]]=None,
                 abacus_command: str=RECOMMAND_COMMAND, 
                 machine: str=RECOMMAND_MACHINE, 
                 image: str=RECOMMAND_IMAGE,
                 lcao: bool=False,
                 nspin: Literal[1,2,4]=1,  # can be 1 or 2 or 4, 1: no spin, 2: spin polarized, 4: non-collinear spin
                 soc: bool =False, # spin-orbit coupling, if True, nspin should be 4 
                 dftu: bool=False,
                 dftu_param: Optional[Dict[str, Union[float, Tuple[Literal["s", 'p', "d"], float]]]]=None,
                 init_mag: Optional[Dict[str, float]] =None,
                 afm: bool = False, 
                 copy_pp_orb: bool = False,
                 ):
        if jobtype not in JOB_TYPES:
            raise ValueError(f"Unsupported job type: {jobtype}.\nSupported job types are {list(JOB_TYPES.keys())}.")
        
        if not isinstance(files, list):
            self.files = [files]
        else:
            self.files = files
        
        # check the parameters
        if nspin not in [1, 2, 4]:
            raise ValueError(f"Invalid nspin value: {nspin}. It should be 1 (no spin), 2 (spin polarized), or 4 (non-collinear spin).")
        if soc and nspin != 4:
            warnings.warn("Spin-orbit coupling is set, but nspin is not 4. Setting nspin to 4 for non-collinear spin calculation.")
            nspin = 4

        self.filetype = filetype
        self.jobtype = jobtype
        self.input_file = input_file
        self.kpt = kpt
        self.abacus_command = abacus_command
        self.machine = machine
        self.image = image
        self.lcao = lcao
        self.nspin = nspin
        self.soc = soc
        self.dftu = dftu
        self.dftu_param = dftu_param
        self.init_mag = init_mag
        self.afm = afm
        self.copy_pp_orb = copy_pp_orb
        
        self.pp_path = pp_path
        self.orb_path = orb_path
        if pp_path is None and os.environ.get("ABACUS_PP_PATH") is not None:
            self.pp_path = os.environ["ABACUS_PP_PATH"]
        if orb_path is None and os.environ.get("ABACUS_ORB_PATH") is not None:
            self.orb_path = os.environ["ABACUS_ORB_PATH"]

        print("Structure files:", self.files)
        print("Structure type:", self.filetype)
        print("Job type:", self.jobtype)
        print("Pseudopotential path:", self.pp_path)
        print("Orbital path:", self.orb_path)
        print("Input file:", self.input_file)
        print("Kpoint:", self.kpt)
        print("LCAO basis:", self.lcao)
        print("Nspin:", self.nspin)
        print("SOC:", self.soc)
        print("DFTU:", self.dftu)
        print("DFTU parameters:", self.dftu_param)
        print("Initial magnetic moment:", self.init_mag)
        print("AFM:", self.afm)
        print("ABACUS command:", self.abacus_command)
        print("Machine:", self.machine)
        print("Image:", self.image)
        print("")
    
    def run(self):
        jobs = gen_stru(self.files, self.filetype, self.pp_path, self.orb_path, tpath=".", copy_pp_orb=self.copy_pp_orb)

        if self.input_file is not None:
            if not os.path.isfile(self.input_file):
                raise ValueError(f"Input file {self.input_file} not found.")
            else:
                input_template = ReadInput(self.input_file)
        else:
            input_template = None
        
        for job in jobs:
            recommand_ecutwfc = jobs[job].get("recommand_ecutwfc", None)
            element = jobs[job]["element"]
            
            input_param = self.generate_input(element)
            input_param = self.update_input(input_param, input_template, recommand_ecutwfc, element)
            
            WriteInput(input_param, os.path.join(job, "INPUT"))
            self.set_init_mag(job, input_param.get("nspin", 1))
            self.write_kpt(job)
        
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

    
    def write_kpt(self, job):
        if self.kpt is not None:
            WriteKpt(self.kpt, os.path.join(job, "KPT"))
    
    def set_init_mag(self, job, nspin):
        """
        Set the initial magnetic moment for the job.
        If afm is True, half of the magnetic elements will be set to negative initial magnetic moment.
        """
        if not self.init_mag and not self.afm:
            return
        
        noncolin = nspin == 4
        
        stru_file = os.path.join(job, "STRU")
        if not os.path.isfile(stru_file):
            raise ValueError(f"Structure file {stru_file} not found.")
        stru = AbacusStru.ReadStru(stru_file)
        element = stru.get_element(number=False, total=True)
        
        if self.afm and not self.init_mag:
            mag_element = list(set([e for e in element if e in MAG_ELEMENT]))
            if len(mag_element) == 0:
                raise ValueError(f"No magnetic elements found in the structure {stru_file}, but AFM is set to True. Please set init_mag or remove afm option.")
            else:
                print(f"Setting initial magnetic moment (3 uB) for antiferromagnetic calculation for elements: {mag_element}.")
            init_mag = {e: 3 for e in mag_element}  # default initial magnetic moment is 1 for all magnetic elements
        else:
            init_mag = self.init_mag
        
        atom_mag = [0 for _ in element]
        atommag_idx = {}
        for i, e in enumerate(element):
            if e in init_mag:
                if e not in atommag_idx:
                    atommag_idx[e] = []
                atommag_idx[e].append(i)
                atom_mag[i] = init_mag[e]
        
        if len(atommag_idx) == 0:
            print(f"No magnetic elements found in the structure {stru_file}, skipping initial magnetic moment setting.")
            return
        
        if self.afm:
            # set half of the magnetic elements to negative initial magnetic moment
            print("Setting initial magnetic moment for antiferromagnetic calculation.")
            for e in atommag_idx:
                idx = atommag_idx[e]
                if len(idx) > 1:
                    for i in idx[::2]:
                        atom_mag[i] = -atom_mag[i]
                    print(f"Odd indexed atoms of element {e} are set to negative magnetic moment.")

        if noncolin:
            atom_mag = [[0,0,i] for i in atom_mag]
            print(f"Non-collinear spin calculation, set initial magnetic moment with three components [0, 0, m] for each atom.")
        
        stru.set_atommag(atom_mag)
        stru.write(stru_file)
        
    
    def update_input(self, input_param, input_template, recommand_ecutwfc, element):
        """
        Update the input parameters based on the input template and recommended ecutwfc.
        """
        if input_template is not None: 
            if "calculation" in input_template:
                calculation = input_template["calculation"]
                if calculation != self.jobtype:
                    warnings.warn(f"The calculation type in the input template is {calculation}, but the job type is {self.jobtype}. The input template will be ignored.")
                input_template.pop("calculation")
            input_param.update(input_template)
        else:
            input_template = {}
        
        # only when input_template has not set ecutwfc, and the basis_type is pw, we will reset the ecutwfc based on the recommand_ecutwfc
        if "ecutwfc" not in input_template and input_param.get("basis_type").startswith("pw"):
            if recommand_ecutwfc is not None and list(set(recommand_ecutwfc)) != [None]:
                input_param["ecutwfc"] = max([i for i in recommand_ecutwfc if i is not None])
                print(f"Set ecutwfc to {input_param['ecutwfc']}, element: {element}, recommended ecutwfc: {recommand_ecutwfc}")
        
        if not input_param.get("basis_type", "pw").startswith("pw"):
            # remove pw specific parameters
            input_param.pop("pw_diag_ndim", None)
            input_param.pop("pw_diag_nmax", None)
            
        return input_param
    
    def generate_input(self, element):
        """
        Set the special input parameters for the job.
        """
        input_param = JOB_TYPES[self.jobtype].copy()
        pw_basis = input_param.get("basis_type", "pw").startswith("pw")
        
        # LCAO
        if self.lcao:
            input_param.update(LCAO_PARAM)
            input_param.pop("pw_diag_ndim", None)
            input_param.pop("pw_diag_nmax", None)
        
        # NSPIN
        if self.nspin == 2:
            input_param["nspin"] = 2
            input_param["mixing_beta"] = 0.4
            input_param["symmetry"] = 0
            input_param["onsite_radius"] = 3  # to generate the atomic magnetic moment
            input_param["out_mul"] = 1  # to output the mulliken charge
            
        if self.soc or self.nspin == 4:
            if self.nspin != 4:
                warnings.warn("nspin is set to 4 for non-collinear spin calculation.")
            input_param["nspin"] = 4
            input_param["noncolin"] = 1
            input_param["mixing_beta"] = 0.4
            input_param["symmetry"] = -1
            input_param["onsite_radius"] = 3
            input_param["out_mul"] = 1  # to output the mulliken charge
            if self.soc:
                input_param["lspinorb"] = 1
        
        # DFTU
        if self.dftu:
            if self.nspin == 1:
                warnings.warn("DFTU is set, but nspin is 1. We suggest to set nspin to 2 or 4 for DFTU calculation.")
            input_param["dft_plus_u"] = 1
            orbital_corr = []
            hubbard_u = []
            for e in element:
                # set u for element defined in dftu_param
                if self.dftu_param is not None and e in self.dftu_param:
                    if isinstance(self.dftu_param[e], (list,tuple)):
                        orbital = {"p":1, "d":2, "f":3}.get(self.dftu_param[e][0], None)  
                        if orbital is None:
                            raise ValueError(f"Unsupported orbital type {self.dftu_param[e][0]} for element {e}. Supported types are p, d, f.")
                        orbital_corr.append(orbital)
                        hubbard_u.append(self.dftu_param[e][1])
                    else:
                        if e in ELEMENT_DFTU_D:
                            orbital_corr.append(2)
                        elif e in ELEMENT_DFTU_F:
                            orbital_corr.append(3)
                        else:
                            orbital_corr.append(1)
                        hubbard_u.append(self.dftu_param[e])
                # if not define dftu_param, then set u for d and f orbital elements
                elif self.dftu_param is None:
                    if e in ELEMENT_DFTU_D:
                        orbital_corr.append(2)
                        hubbard_u.append(4)  # default U value for d orbital
                    elif e in ELEMENT_DFTU_F:
                        orbital_corr.append(3)
                        hubbard_u.append(6)  # default U value for f orbital
                    else:
                        orbital_corr.append(-1)
                        hubbard_u.append(0)  # no DFTU for other elements
                else:
                    orbital_corr.append(-1)
                    hubbard_u.append(0)  # no DFTU for other elements
            input_param["orbital_corr"] = orbital_corr
            input_param["hubbard_u"] = hubbard_u
            input_param["onsite_radius"] = 3
            if len(set(orbital_corr)) > 1 or orbital_corr[0] != -1:
                if pw_basis:
                    input_param["uramping"] = max(hubbard_u)
                else:
                    input_param["mixing_dmr"] = 1
            
            input_param["mixing_restart"] = 0.001
        
        if input_param.get("basis_type", "pw").startswith("pw") and input_param.get("out_mul", ""):
            # out_mul is only valid for lcao basis
            input_param.pop("out_mul")
            
        if input_param.get("basis_type", "pw").startswith("pw") and input_param.get("nspin",1) == 2 and "onsite_radius" in input_param:
            # onsite_radius is invalid for pw basis with nspin=2
            input_param.pop("onsite_radius")

        return input_param
    