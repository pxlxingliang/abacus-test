import os
import copy
import numpy as np
from ..model import Model
from . import comm
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,WriteInput,WriteKpt
from abacustest.lib_prepare.comm import kspacing2kpt

class Eos(Model):
    '''
    Prepare the input files for the EOS calculation
    '''
    @staticmethod
    def model_name():
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "eos"
    
    @staticmethod
    def description():
        '''
        Description of the model
        '''
        return "Prepare and postprocess the EOS calculation"
    
    @staticmethod
    def add_args(parser):
        '''
        Add arguments for the model
        The arguments can not be command, model, modelcommand 
        '''
        pass
    
    def run(self,params):
        '''
        Parse the parameters and run the model
        '''
        pass
    
    @staticmethod
    def prepare_args(parser):
        parser.description = "Prepare the input files for the EOS calculation. Will generate some eos folders in each example path."
        parser.add_argument('-s', '--start', type=float, default=0.9,help="the maximum scaling down of volume, should be less than 1, default is 0.9")
        parser.add_argument('-e', '--end', type=float,  default=1.1,help='the maximum scaling up of volume, should be larger than 1, default is 1.1')
        parser.add_argument('-d', '--step', type=float,  default=0.025,help='the step size of volume scaling, default is 0.025')
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is the current path.')
        parser.add_argument('--kspacing', default=0, type=int, help="if kspacing is defined in INPUT, then use the this kspacing in all EOS calculation. If set to 0, then will generate the KPT file by kspacing and use this KPT for all EOS calculation. 0: not use; 1: use. default 0")
        parser.add_argument('--relax',default=0, type=int,help='the calculation type. -1:use the setting in INPUT; 0: scf; 1: relax; 2:cell-relax with fixed volume; default -1')
        parser.add_argument('--clean',default=1, type=int,help='if clean the eos folder in each example path, 0: not clean; 1: clean. default 1')
        #parser.add_argument("-r", "--run", default=0, help="if run the test. Default is 0.", type=int)
    
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        jobs = params.job if len(params.job) == 1 else params.job[1:]
        real_jobs = comm.get_job_list(jobs)
        peos = PrepareEos(real_jobs,params.start,params.end,params.step,params.relax,params.clean,params.kspacing)
        subfolders = peos.run()
        
        setting = {
            "save_path": "results",
            "bohrium_group_name": "EOS",
            "run_dft": {
                "example": subfolders,
                "command": "OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log",
                "image": "registry.dp.tech/deepmodeling/abacus-intel:latest",
                "bohrium": {
                    "scass_type": "c32_m64_cpu",
                    "job_type": "container",
                    "platform": "ali",
                },
            },
        }
        comm.dump_setting(setting)
        comm.doc_after_prepare("EOS", subfolders, ["setting.json"],has_prepare=False)
    
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
    
    
class PrepareEos:
    '''
    A class to prepare the input files for the EOS calculation.
    '''
    
    def __init__(self,jobs,start,end,step,relax_type,clean,use_kspacing):
        self.jobs = jobs
        self.start = start
        self.end = end
        self.step = step # the step size of volume scaling, default is 0.025
        self.relax_type = relax_type  # -1:use the setting in INPUT; 0: scf; 1: relax; 2:cell-relax with fixed volume; default -1
        self.clean = clean # if clean the eos folder in each example path, 0: not clean; 1: clean. default 1
        self.use_kspacing = use_kspacing # if kspacing is defined in INPUT, then use the this kspacing in all EOS calculation. If set to 0, then will generate the KPT file by kspacing and use this KPT for all EOS calculation. 0: not use; 1: use. default 0
        self.eos_path = "eos"
        pass
    
    def run(self):
        print("prepare EOS inputs")
        subfolders = []
        for job in self.jobs:
            if self.clean:
                comm.clean_files(job,folder_list=[self.eos_path+"*"])
            subfolders += self.prepare_onejob(job,self.start,self.end,self.step,self.relax_type,job)
        return subfolders

    @staticmethod
    def prepare_onejob(ijob,t_path,start,end,step,relax_type,use_kspacing=False,eos_path="eos"):
        # read the input/STRU/kpt/pp/orb files
        abacus_input = comm.get_abacus_inputfiles(ijob)
        input_param = abacus_input["input"]
        stru = abacus_input["stru"]
        otherfiles = abacus_input["extra_files"]
        kpt = abacus_input["kpt"]
        if abacus_input["pp"]:
            otherfiles += abacus_input["pp"]
        if abacus_input["orb"]:
            otherfiles += abacus_input["orb"]
        
        if not input_param or not stru:
            print(f"ERROR: read input/STRU failed in {ijob}")
            return []

        # modify calculation type
        if relax_type == 0:
            input_param["calculation"] = "scf"
        elif relax_type == 1:
            input_param["calculation"] = "relax"
        elif relax_type == 2:
            input_param["calculation"] = "cell-relax"
            input_param["fixed_axes"] = "volume" 
        
        # set the KPT
        # if kspacing is defined in INPUT, then 
        if "kspacing" in input_param and not use_kspacing:
            kpt = kspacing2kpt(input_param["kspacing"],stru.get_cell(bohr=True))  # three integers
            input_param.pop("kpt_file","")  # remove the kpt_file defined in INPUT

        # write the input file
        subfolders = []
        # write eos0
        new_stru = copy.deepcopy(stru)
        cell = stru.get_cell(bohr=True)
        new_stru.set_cell(cell,bohr=True,change_coord=True)
        spath = os.path.join(t_path,eos_path+"0")
        PrepareEos.write_onejob(new_stru, input_param, kpt, otherfiles, spath)
        subfolders.append(spath)

        # write scaled up/down stru
        number = 1
        while True:
            lcd = 1 - number*step
            lcu = 1 + number*step
            if lcd < start and lcu > end:
                break
            if lcd >= start:
                new_stru = copy.deepcopy(stru)
                new_cell = np.array(cell) * lcd**(1/3)
                new_stru.set_cell(new_cell,bohr=True,change_coord=True)
                spath = os.path.join(t_path,f"{eos_path}{-number}")
                PrepareEos.write_onejob(new_stru, input_param, kpt, otherfiles, spath)
                subfolders.append(spath)
            if lcu <= end:
                new_stru = copy.deepcopy(stru)
                new_cell = np.array(cell) * lcu**(1/3)
                new_stru.set_cell(new_cell,bohr=True,change_coord=True)
                spath = os.path.join(t_path,f"{eos_path}{number}")
                PrepareEos.write_onejob(new_stru, input_param, kpt, otherfiles, spath)
                subfolders.append(spath)
            number += 1
        return subfolders
    
    @staticmethod
    def write_onejob(new_stru, input_param, kpt, extra_files, t_path):
        # wrtie the input file to the t_path
        new_stru.write_stru(os.path.join(t_path,"STRU"))
        WriteInput(input_param,os.path.join(t_path,"INPUT"))
        if isinstance(kpt,list):
            WriteKpt(kpt,os.path.join(t_path,"KPT"))
        elif isinstance(kpt,str) and os.path.isfile(kpt):
            os.symlink(kpt,os.path.join(t_path,os.path.basename(kpt)))
        for ifile in extra_files:
            if os.path.isfile(ifile):
                os.symlink(ifile,os.path.join(t_path,os.path.basename(ifile)))
    