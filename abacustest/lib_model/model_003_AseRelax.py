from ..model import Model
import os, json, sys, inspect
from . import comm, comm_plot


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
        parser.add_argument('--cellrelax', nargs='?',type=int, const=1, default=None,help='if relax the box. 0: no, 1: yes. Default will read the INPUT, and set to 1 only when calculation is cell-relax' )
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
        parser.add_argument("-j","--jobs",type=str,help="the path of jobs to be tested",action="extend",nargs="*",)
        parser.add_argument("-c", "--rundftcommand", type=str, default="abacustest model aserelax -o BFGS --mpi 32 --omp 1",help="the command to execute aserelax, default is 'abacustest model aserelax -o BFGS' ")
        parser.add_argument("-i","--image",default="registry.dp.tech/dptech/prod-471/abacus-ase:20240522",type=str,help="the used image. Should has ABACUS/ASE-ABACUS/abacustest in image", )
        parser.add_argument("--machine", default="c32_m128_cpu", help="the machine to run the abacus. Default is c32_m128_cpu")
        parser.add_argument("-r", "--run", default=0, help="if run the test. Default is 0.", type=int)
        
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        real_jobs = comm.get_job_list(params.jobs)
        setting = {
            "save_path": "results",
            "bohrium_group_name": "ase-abacus-relax",
            "run_dft": {
                "example": real_jobs,
                "command": params.rundftcommand,
                "image": params.image,
                "bohrium": {
                    "scass_type": params.machine,
                    "job_type": "container",
                    "platform": "paratera",
                },
            },
        }
        
        comm.dump_setting(setting)
        comm.doc_after_prepare("ASE-ABACUS relax", real_jobs, ["setting.json"],has_prepare=False)
        print("After finish the calculation, you can run below command to do the postprocess:")
        print(f"    abacustest model {self.model_name()} post -j {' '.join(real_jobs)} -r result.json\n")
        
        if params.run:
            bash_script = "aserelax.sh"
            with open(bash_script,"w") as f:
                f.write("abacustest submit -p setting.json\n")
                f.write("cd results\n")
                f.write(f"abacustest model {self.model_name()} post -j {' '.join(real_jobs)} -r result.json\n")
            os.system(f"bash {bash_script} &")

    @staticmethod
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand'''
        parser.add_argument("-j","--jobs",type=str,default=["."],help="The path of aserelax job. Should has metrics.json file generated after aserelax.",action="extend",nargs="*",)
        parser.add_argument("-o","--output",default="aserelax.png",help="the output picture name, default is aserelax.png",)
        parser.add_argument("-r","--result",default="result.json",help="Save the data of the plot to a json file. Default is result.json.",)
        parser.add_argument('--abacus', nargs='?',type=int, const=1, default=None,help='postprocess the job as an abacus relax job' )
        parser.add_argument('--metric', nargs='?',type=str, const="metrics.json", default=None,help='postprocess a metrics.json file. If set metric, then will only plot data in this file' )
        parser.add_argument('--noplot', nargs='?',type=int, const=1, default=None,help='If not plot the image' )

    def _post_gen_allmetrics(self,params):
        if params.metric:
            allmetrics = json.load(open(params.metric))
        else:
            allmetrics = {}
            jobs = params.jobs if len(params.jobs) == 1 else params.jobs[1:]
            jobs = comm.get_job_list(jobs)
            if len(jobs) == 0:
                print("No jobs are specified.")
                return {}

            if params.abacus:
                from abacustest.lib_collectdata.collectdata import RESULT
                for ijob in jobs:
                    result = RESULT(path=ijob,fmt="abacus")
                    forces = result["forces"]
                    stresses = result["stresses"]
                    fmax = []
                    smax = []
                    if forces:
                        for force in forces:
                            fmax.append(max([ (sum([j**2 for j in force[3*i:3*i+3]]))**0.5 for i in range(len(force)//3)]))
                    if stresses:
                        for stress in stresses:
                            smax.append(max([ abs(s) for s in stress ]))
                            
                    allmetrics[ijob.rstrip("/")] = {
                        "version": result["version"],
                        "energy": result["energy"],
                        "energy_traj": result["energies"],
                        "fmax_traj": fmax,
                        "smax_traj": smax,
                        "relax_steps": result["relax_steps"],
                        "relax_converge": result["relax_converge"]
                    }
            else:
                for ijob in jobs:
                    if os.path.isfile(os.path.join(ijob,"metrics.json")):
                        metrics = json.load(open(os.path.join(ijob,"metrics.json")))
                        allmetrics[ijob.rstrip("/")] = metrics
                    else:
                        print(f"metrics.json not found in {ijob}")
                if len(allmetrics) == 0:
                    print("No metrics.json found in all jobs.")
                    return {}
        
        return allmetrics
    
    def _post_gen_plot_data(self,allmetrics):
        plot_data = {}
        for ik,iv in allmetrics.items():
            enes = iv.get("energy_traj")
            fmax = iv.get("fmax_traj")
            if enes and fmax:
                plot_data[ik] = {
                    "enes": enes,
                    "fmax": fmax,
                    "smax": iv.get("smax_traj", None),
                }
        return plot_data 
    
    def _post_plot(self,plot_data,output):
        import matplotlib.pyplot as plt
        import numpy as np
        comm_plot.set_font("Times New Roman", font_size=12)
        examples = list(plot_data.keys())
        examples.sort()
        
        ncol = int((len(examples))**0.5)
        nrow = int(len(examples)//ncol)
        while nrow * ncol < len(examples):
            nrow += 1
        
        fig, axes = plt.subplots(nrow,ncol,figsize=(6*ncol,4*nrow))
        for i,example in enumerate(examples):
            if nrow * ncol == 1:
                ax0 = axes
            elif nrow == 1 or ncol == 1:
                ax0 = axes[i]
            else:
                ax0 = axes[i//ncol][i%ncol]
            ax0.plot(plot_data[example]["enes"],label="Energy (eV)",color="blue",linestyle="-",marker="o")
            ax1 = ax0.twinx()
            ax1.plot(plot_data[example]["fmax"],label="Fmax (eV/$\mathrm{\AA}$)",color="red",linestyle="--",marker="o")
            if plot_data[example].get("smax"):
                ax1.plot(np.array(plot_data[example]["smax"])/100,label="Stress_max (kbar/100)",color="green",linestyle="--",marker="x")
            ax0.set_title(example)
            ax0.set_xlabel("Step")
            ax0.set_ylabel("Energy (eV)")
            if plot_data[example].get("smax"):
                ax1.set_ylabel("Force (eV/$\mathrm{\AA}$)/Stress (kbar/100)")
            else:
                ax1.set_ylabel("Fmax (eV/$\mathrm{\AA}$)")
            ax0.legend(loc="upper left")
            ax1.legend(loc="upper right")
            ene_max = max(plot_data[example]["enes"])
            ene_min = min(plot_data[example]["enes"])
            if plot_data[example].get("smax"):
                fmax_max = max(max(plot_data[example]["fmax"]), max(plot_data[example]["smax"])/100)
                fmax_min = min(min(plot_data[example]["fmax"]), min(plot_data[example]["smax"])/100)
            else:
                fmax_max = max(plot_data[example]["fmax"])
                fmax_min = min(plot_data[example]["fmax"])
            ax0.set_ylim(ene_min,ene_max + (ene_max - ene_min)*0.2)
            ax1.set_ylim(fmax_min,fmax_max + (fmax_max - fmax_min)*0.2)
            # plot line y = 0.01
            ax1.axhline(y=0.01, color='gray', linestyle='--', linewidth=0.5)
            # set y log scale
            if fmax_max / fmax_min > 100:
                ax1.set_yscale("log")
        plt.tight_layout()
        plt.savefig(output,dpi=300)
    
    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        all_metrics = self._post_gen_allmetrics(params)
        plot_data = self._post_gen_plot_data(all_metrics)
        if not params.noplot:
            self._post_plot(plot_data,params.output)
            json.dump(comm.gen_supermetrics(params.output),open("supermetrics.json","w"),indent=4)
        if not params.metric: 
            json.dump(all_metrics,open(comm.bak_file("metrics.json",bak_org=False),"w"),indent=4)
        json.dump(plot_data,open(params.result,"w"),indent=4)
        

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
        self.cal_stress = False
        self.logfile = "aserelax.log"

    def run(self):
        os.makedirs(self.work_path,exist_ok=True)
        input_param, atoms = self.prepare_param()
        self.set_fmax(input_param)
        optimizer, opt_module = self.set_ase_opt()
        self.print_info(opt_module)
        qn_init = self.exe_ase(atoms, input_param, optimizer)

        metrics = self.read_metrics(self.work_path,qn_init)
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
        if input_param.get("cal_stress") or self.relax_cell:
            self.cal_stress = True
        input_param["calculation"] = "scf"  # only run scf
        stru = AbacusStru.ReadStru(os.path.join(init_path,"STRU"))
        if not stru:
            raise ValueError(f"STRU file not found in {init_path}")

        labels = stru.get_label(total=False)
        pp = stru.get_pp()
        orb = stru.get_orb()
        input_param["pp"] = {labels[i]:os.path.abspath(os.path.join(init_path,input_param.get("pseudo_dir",""),pp[i])) for i in range(len(labels))}
        if input_param.get("basis_type") in ["lcao"] and orb:
            input_param["basis"]  = {labels[i]:os.path.abspath(os.path.join(init_path,input_param.get("orbital_dir",""),orb[i])) for i in range(len(labels))}

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
        c = FixAtoms(mask=[True if set(i) == {0} else False for i in move])

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
        init_signature = inspect.signature(AbacusProfile.__init__)
        if "command" in init_signature.parameters:
            command = f'mpirun -np {self.mpi} {self.abacus}'
            profile = AbacusProfile(command=command)
        else:
            profile = AbacusProfile(
                argv=['mpirun', '-np', f'{self.mpi}', self.abacus])
        atoms.calc = Abacus(profile=profile, directory=self.work_path,
                         **input_param)

        if self.relax_cell:
            ucf = ExpCellFilter(atoms,)
            opt = optimizer(ucf, trajectory='init_opt.traj',logfile = self.logfile)
        else:
            opt = optimizer(atoms, trajectory='init_opt.traj', logfile=self.logfile)

        opt.run(fmax=self.fmax,steps = input_param.get("relax_nmax",100))
        #opt.log()
        print("RELAX STEPS:",opt.get_number_of_steps()+1)

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
        
        ase_force = opt.atoms.get_forces().tolist()
        
        if ase_force:
            max_force = max((sum([j**2 for j in i]))**0.5 for i in ase_force)
            max_force_comp = max(max([abs(j) for j in i]) for i in ase_force)
        else:
            max_force = None
            max_force_comp = None
        
        if os.path.isfile(self.logfile):
            with open(self.logfile) as f:
                logs = f.readlines()
            nsteps = opt.get_number_of_steps() + 1
            enes = []
            fmaxs = []
            for idx,line in enumerate(logs):
                if "Step     Time          Energy          fmax" in line:
                    i = idx + 1
                    while nsteps > 0:
                        if logs[i].strip() == "": continue
                        ene, fmax = logs[i].split()[3:5]
                        enes.append(float(ene))
                        fmaxs.append(float(fmax))
                        i += 1
                        nsteps -= 1
        else:
            enes = None
            fmaxs = None

        return {
            "version": iresult["version"],
            "energy": iresult["energy"],
            "energy_per_atom": iresult["energy_per_atom"],
            "max_force": max_force,
            "max_force_comp": max_force_comp,
            "energy_traj": enes,
            "fmax_traj": fmaxs,
            "force": ase_force,
            "abacus_force": iresult["force"],
            "abacus_stress": iresult["stress"],
        }
