import ase
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
        parser.add_argument('--cellrelax', type=int,  default=None,help='if relax the box. 0: no, 1: yes. Default will read the INPUT, and set to 1 only when calculation is cell-relax' )
        parser.add_argument('-j','--job', type=str,  default=".",help='the path of abacus inputs, default is current folder.' )
    
    def run(self,params):
        '''
        Parse the parameters and run the model
        '''
        if params.mpi == 0:
            tcore = comm.get_physical_cores()
            mpi = int(int(tcore)/params.omp)
        else:
            mpi = params.mpi
        aserelax = ExeAseRelax(params.job, params.abacus, params.omp, mpi, params.optimize, params.fmax, params.cellrelax,"aserelax")
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
    def __init__(self, job, abacus, omp, mpi, optimize, fmax, relax_cell, wrok_path):
        self.job = job
        self.abacus = abacus
        self.omp = omp
        self.mpi = mpi
        self.optimize = optimize
        self.fmax = fmax
        self.work_path = wrok_path
        self.relax_cell = relax_cell
        self.logfile = "aserelax.log"

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
        if self.relax_cell == None:
            self.relax_cell = input_param.get("calculation","scf") == "cell-relax"
        input_param["calculation"] = "scf"  # only run scf
        stru = AbacusStru.ReadStru(os.path.join(init_path,"STRU"))
        if not stru:
            raise ValueError(f"STRU file not found in {init_path}")

        labels = stru.get_label(total=False)
        pp = stru.get_pp()
        orb = stru.get_orb()
        input_param["pp"] = {labels[i]:os.path.abspath(os.path.join(self.job,input_param.get("pseudo_dir",""),pp[i])) for i in range(len(labels))}
        if input_param.get("basis_type") in ["lcao"] and orb:
            input_param["basis"]  = {labels[i]:os.path.abspath(os.path.join(self.job,input_param.get("pseudo_dir",""),orb[i])) for i in range(len(labels))}

        kpt = ReadKpt(init_path)
        if kpt:
            input_param["kpts"] = kpt[0][:3]

        # construct an ASE Atoms object
        from ase import Atoms
        from ase.constraints import FixAtoms
        cell = stru.get_cell(bohr=False,)
        coord = stru.get_coord(bohr=False,direct=False)
        mag = stru.get_atommag()
        move = stru.get_move()
        c = FixAtoms(mask=[i==0 for i in move])

        atoms = Atoms(symbols=stru.get_label(total=True), positions=coord, cell=cell, pbc=[True,True,True], magmoms=mag, constraint=c)
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
        logs = ""
        logs += "SETTING OMP_NUM_THREADS: {}\n".format(self.omp)
        logs += "SETTING     mpi process: {}\n".format(self.mpi)
        logs += "SETTING     ABACUS path: {}\n".format(self.abacus)
        logs += "SETTING        job path: {}\n".format(os.path.abspath(self.job))
        logs += "SETTING       work path: {}\n".format(os.path.abspath(self.work_path))
        logs += "SETTING      fmax(eV/A): {}\n".format(self.fmax)
        logs += "SETTING      cell relax: {}\n".format(bool(self.relax_cell))
        logs += "SETTING    ASE OPTIMIZE: {}\n".format(opt_module)
        logs += "\n"
        print("\n" + logs)
        with open(self.logfile,"w") as f:
            f.write(logs)

    def exe_ase(self, atoms, input_param, optimizer):
        from ase.calculators.abacus import Abacus, AbacusProfile
        from ase.constraints import UnitCellFilter, ExpCellFilter

        os.environ['OMP_NUM_THREADS'] = f'{self.omp}'
        profile = AbacusProfile(
            argv=['mpirun', '-np', f'{self.mpi}', self.abacus])
        atoms.calc = Abacus(profile=profile, directory=self.work_path,
                         **input_param)

        if self.relax_cell:
            ucf = ExpCellFilter(atoms,)
            opt = optimizer(ucf, trajectory='init_opt.traj')
        else:
            opt = optimizer(atoms, trajectory='init_opt.traj', logfile="aserelax.log")

        opt.run(fmax=self.fmax,steps = input_param.get("relax_nmax",100))
        opt.log()
        print("RELAX STEPS:",opt.get_number_of_steps())

        return opt

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

    def read_metrics(self,jobpath,opt):
        from abacustest.lib_collectdata.collectdata import RESULT
        iresult = RESULT(fmt="abacus",path=jobpath)
        force = iresult["force"]
        stress = iresult["stress"]
        if force == None:
            max_force = None
            max_force_comp = None
        else:
            max_force = max([(sum([force[3*i+j]**2 for j in range(3)]))**0.5 for i in range(len(force)//3)])
            max_force_comp = max([abs(i) for i in force])

        if stress == None:
            pressure = None
        else:
            pressure = (stress[0] + stress[4] + stress[8])/3.0
            
        ase_force = opt.atoms.get_forces()
        ase_stress = opt.atoms.get_stress(voigt=False)
        
        return {
            "version": iresult["version"],
            "energy": iresult["energy"],
            "energy_per_atom": iresult["energy_per_atom"],
            "max_force": max_force,
            "max_force_comp": max_force_comp,
            "pressure": pressure,
            "force": force,
            "stress": stress,
            "ase_force": ase_force,
            "ase_stress": ase_stress,
        }
