from ..model import Model
import os, glob, json, sys
from . import comm


class AseRelax(Model):
    '''
    Base class for all models
    
    Some specific workflow will be implemented in the subclass.
    Such as EOS, BAND, etc.
    
    The model should have the following methods:
    '''
    @staticmethod
    def model_name():
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "aserelax"
    
    @staticmethod
    def description():
        '''
        Description of the model
        '''
        return "Prepare the inputs for relax by ASE + ABACUS"
    
    @staticmethod
    def add_args(parser):
        '''
        Add arguments for the model
        The arguments can not be command, model, modelcommand 
        '''
        parser.description = "This script is used to run relax job by ase-abacus!\n Need all inputs file of abaucs in one folder."
        parser.add_argument('-o', '--optimize', type=str, default="BFGS",help="the relax optimize method of ASE, can be: 'CG', 'BFGS', 'BFGSLineSearch', 'Berny', 'FIRE', 'GPMin', 'GoodOldQuasiNewton', 'LBFGS', 'LBFGSLineSearch', 'MDMin', 'ODE12r', 'QuasiNewton', 'RestartError'. Default is BFGS")
        parser.add_argument('-a', '--abacus', type=str,  default="abacus",help='the path of abacus executable, default is abacus')
        parser.add_argument('--fmax', type=float,  default=None,help='the fmax for the relax (eV/A). Default is None. and will read from INPUT file or 0.0257112 eV/A.') 
        parser.add_argument('--omp', type=int,  default=1,help='number of OMP parrallel, default 1 ')
        parser.add_argument('--mpi', type=int,  default=0,help='number of MPI parrallel, default 0, which means all cores/number of omp.' )
        parser.add_argument('-j','--job', type=str,  default=".",help='the path of abacus inputs, default is current folder.' )
    
    def run(self,params):
        '''
        Parse the parameters and run the model
        '''
        if params.mpi == 0:
            tcore = os.system("lscpu | grep 'Core(s) per socket' | awk '{print $4}'")
            mpi = int(int(tcore)/params.omp)
        else:
            mpi = params.mpi
        aserelax = ExeAseRelax(params.job, params.abacus, params.omp, mpi, params.optimize, params.fmax, "aserelax")
        aserelax.run()
    
    @staticmethod
    def prepare_args(parser):
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand '''
        
        parser.add_argument("-j","--jobs",default=[],type=str,help="the path of jobs to be tested",action="extend",nargs="*",)
        
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        pass
        
    
    @staticmethod
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand'''
        pass
    
    
    
    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        pass
    

class ExeAseRelax:
    def __init__(self, job, abacus, omp, mpi, optimize, fmax, wrok_path):
        self.job = job
        self.abacus = abacus
        self.omp = omp
        self.mpi = mpi
        self.optimize = optimize
        self.fmax = fmax
        self.work_path = wrok_path
    
    def run(self):
        os.makedirs(self.work_path,exist_ok=True)
        input_param, atoms = self.prepare_param()
        self.set_fmax(input_param)
        opt, opt_module = self.set_ase_opt()
        self.print_info(opt_module)
        qn_init = self.exe_ase(atoms, input_param, opt)
        
        metrics = self.read_metrics(self.work_path)
        metrics.update({"relax_steps": int(qn_init.nsteps) + 1,
                   "relax_converge": bool(qn_init.converged()) })
        json.dump(metrics,open("metrics.json","w"),indent=4)
            
    def prepare_param(self):
        # read inputs and return a dict of INPUT parameters
        # Should also set extra parameters: pp/orb/kpts
        from abacustest.lib_prepare.abacus import ReadInput, AbacusStru, ReadKpt
        
        init_path = self.job
        input_param = ReadInput(os.path.join(init_path,"INPUT"))
        input_param["calculation"] = "scf"  # only run scf
        stru = AbacusStru.ReadStru(os.path.join(init_path,"STRU"))
        if not stru:
            raise ValueError(f"STRU file not found in {init_path}")
        
        labels = stru.get_label(total=False)
        cell = stru.get_cell(bohr=False,)
        coord = stru.get_coord(bohr=False,direct=False)
        pp = stru.get_pp()
        orb = stru.get_orb()
        input_param["pp"] = {labels[i]:pp[i] for i in range(len(labels))}
        if orb:
            input_param["basis"]  = {labels[i]:orb[i] for i in range(len(labels))}

        kpt = ReadKpt(init_path)
        if kpt:
            input_param["kpts"] = kpt[0][:3]
            
        # construct an ASE Atoms object
        from ase import Atoms
        atoms = Atoms(symbols=stru.get_label(total=True), positions=coord, cell=cell)
        return input_param, atoms
    
    def set_fmax(self, input_param):
        if not self.fmax:
            if input_param.get("force_thr"):
                self.fmax = float(input_param.get("force_thr")) * 25.7112
            elif input_param.get("force_thr_ev"):
                self.fmax = float(input_param.get("force_thr_ev"))
            else:
                self.fmax = 0.0257112
    
    def print_info(self, opt_module):
        print("\n")
        print("SETTING OMP_NUM_THREADS:",self.omp)
        print("SETTING     mpi process:",self.mpi)
        print("SETTING          abacus:",self.abacus)
        print("SETTING      INPUT path:",os.path.abspath(self.job))
        print("SETTING    running path:",os.path.abspath(self.work_path))
        print("SETTING      fmax(eV/A):",self.fmax)
        print("SETTING    ASE OPTIMIZE:",opt_module)
        print("\n")
                
    def exe_ase(self, atoms, input_param, optimizer):
        from ase.calculators.abacus import Abacus, AbacusProfile
        
        os.environ['OMP_NUM_THREADS'] = f'{self.omp}'
        profile = AbacusProfile(
            argv=['mpirun', '-np', f'{self.mpi}', self.abacus])
        atoms.calc = Abacus(profile=profile, directory=self.work_path,
                         **input_param)
        qn_init = optimizer(atoms, trajectory='init_opt.traj')
        qn_init.run(fmax=self.fmax,steps = input_param.get("relax_nmax",100))
        print("RELAX STEPS:",qn_init.get_number_of_steps())

        return qn_init

    def set_ase_opt(self):
        ase_optimize = ['BFGS', 'BFGSLineSearch', 'Berny', 'FIRE', 'GPMin', 'GoodOldQuasiNewton', 'LBFGS', 'LBFGSLineSearch', 'MDMin', 'ODE12r', 'QuasiNewton', 'RestartError']
        ase_opt_lower = [i.lower() for i in ase_optimize]
        optimize = self.optimize.lower()
        if optimize == "cg":
            from ase.optimize.sciopt import  SciPyFminCG  # ,SciPyFminBFGS
            
            return SciPyFminCG, "ase.optimize.sciopt.SciPyFminCG"
        elif optimize in ase_opt_lower:
            idx = ase_opt_lower.index(optimize)
            opt = ase_optimize[idx]
            import importlib
            optimize_module = importlib.import_module("ase.optimize")
            return getattr(optimize_module, opt), f"ase.optimize.{opt}"
        else:
            print("Do not support optimize:",self.optimize)
            sys.exit(1)    
    
    def read_metrics(self,jobpath):
        from abacustest.lib_collectdata.collectdata import RESULT
        iresult = RESULT(fmt="abacus",path=jobpath)
        return {"force": iresult["force"],
                "stress": iresult["stress"],
                "energy": iresult["energy"],
                "energy_per_atom": iresult["energy_per_atom"],
                "version": iresult["version"]}