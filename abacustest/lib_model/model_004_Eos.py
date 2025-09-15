import os,glob,json
import copy
import numpy as np
from ..model import Model
from . import comm
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,WriteInput,WriteKpt
from abacustest.lib_prepare.comm import kspacing2kpt
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

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
        parser.add_argument('--relax',default=-1, type=int,help='the calculation type. -1:use the setting in INPUT; 0: scf; 1: relax; 2:cell-relax with fixed volume; default -1')
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
                "command": RECOMMAND_COMMAND,
                "image": RECOMMAND_IMAGE,
                "bohrium": {
                    "scass_type": RECOMMAND_MACHINE,
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
        Catch some key values for several EOS jobs, and then do below things:
        1. plot the EOS curve: volume vs energy_per_atom
        2. fit the EOS curve: .
        3. calculate the Delta value between two EOS curves: .
        '''
        parser.description = "Postprocess the EOS calculation. Will plot the EOS curve and fit the EOS curve."
        parser.add_argument('-j', '--job',default=[], action="extend",nargs="*" ,help='the path of abacus inputs, default is the current path. There should be some eos folders in each example path.')
        parser.add_argument("-s", '--save', default="result.json", type=str, help="the file to save the EOS results. Default is result.json")
        parser.add_argument('-r', '--result',default=None, type=str,help='the result file that contains the volume and energy_per_atom of each example. If has defined job, then will deal with both job and result file.')
        parser.add_argument('--ref', default=[], action="extend",nargs="*" , help='the reference result file. Will calculate the Delta value between results and reference. The format should be like result.json')
        parser.add_argument('--sumfile', default=None,const="alleos.png" , type=str, nargs='?',help='if plot all the EOS curves in one figure, then save the figure to this file. Default file name is alleos.png')

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        jobs = params.job 
        result = params.result
        ref = params.ref
        savef = params.save
        sumfile = params.sumfile
        
        if not jobs and not result:
            print("No job or metric is specified.")
            return
        
        PostEos(jobs,result,savef,ref,sumfile).run()

    
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
            subfolders += PrepareEos.prepare_onejob(job,job,self.start,self.end,self.step,self.relax_type,self.use_kspacing,self.eos_path)
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
        os.makedirs(t_path,exist_ok=True)
        new_stru.write(os.path.join(t_path,"STRU"))
        WriteInput(input_param,os.path.join(t_path,"INPUT"))
        if isinstance(kpt,list):
            WriteKpt(kpt,os.path.join(t_path,"KPT"))
        elif isinstance(kpt,str) and os.path.isfile(kpt):
            os.symlink(kpt,os.path.join(t_path,os.path.basename(kpt)))
        for ifile in extra_files:
            if os.path.isfile(ifile):
                os.symlink(ifile,os.path.join(t_path,os.path.basename(ifile)))


class PostEos:
    '''
    A class to postprocess the EOS calculation
    '''
    
    def __init__(self,jobs,result,savef,refs,plot_summary_file=None):
        self.jobs = jobs
        self.result = result
        self.savef = savef # the file to save the EOS results which only volume/energy_per_atom
        self.refs = refs
        self.eos_path = "eos"
        self.metrics_file = "metrics.json" # the file to save the metrics which contains more details
        self.plot_summary_file = plot_summary_file # if plot the summary file
    
    def run(self):
        allresults,allmetrics = self.collectdata_result()
        json.dump(allresults,open(self.savef,"w"),indent=4)
        json.dump(allmetrics,open(self.metrics_file,"w"),indent=4)
        result_ref = self.collect_refs()
        PostEos.plot_eos(allresults,result_ref,self.plot_summary_file)
        
    def collectdata_result(self):
        allresults = {}  # key is {jobname}, value is volume and energy_per_atom
        allmetrics = {}  # key is {jobname/eos*}, value is volume, energy_per_atom, converge, normal_end, total_time, scf_steps, version...
        if self.result:
            if os.path.isfile(self.result):
                allresults = json.load(open(self.result))
        
        if self.jobs:
            for ijob in self.jobs:
                vs = []
                es = []
                alleos = glob.glob(f"{ijob}/{self.eos_path}*")
                for ieos in alleos:
                    if os.path.isdir(ieos):
                        print("collecting",ieos)
                        r = RESULT(path=ieos,fmt="abacus")
                        if None not in [r["volume"], r["energy_per_atom"]]:
                            vs.append(r["volume"])
                            es.append(r["energy_per_atom"])
                        allmetrics[ieos] = {
                            "volume": r["volume"],
                            "energy_per_atom": r["energy_per_atom"],
                            "natom": r["natom"],
                            "converge": r["converge"],
                            "normal_end": r["normal_end"],
                            "relax_converge": r["relax_converge"],
                            "relax_steps": r["relax_steps"],
                            "total_time": r["total_time"],
                            "scf_steps": r["scf_steps"],
                            "version": r["version"]
                        }
                if vs:
                    allresults[ijob.rstrip("/")] = {"volume":vs,"energy_per_atom":es}
        return allresults,allmetrics
    
    def collect_refs(self):
        result_ref = {}
        if self.refs:
            for iref in self.refs:
                if os.path.isfile(iref):
                    result_ref.update(json.load(open(iref)))
        return result_ref

    @staticmethod
    def plot_eos(results, ref_results={}, one_plot_file=None):
        '''
        Plot the EOS curve,
        
        Args:
            results (dict): The results of the EOS calculation. The key is the jobname, the value is a dict of volume and energy_per_atom.
            ref_results (dict): The reference results of the EOS calculation. The key is the jobname, the value is the volume and energy_per_atom.
            one_plot_file (str): The file to save the EOS curve. If not defined, then will show the EOS curve.
        
        If ref_results is defined, then will also plot the reference EOS curve, and calculate the Delta value between the results and reference.
        If one_plot_file is defined, then will also plot all the EOS curves in one figure.
        '''
        from . import comm_eos
        import matplotlib.pyplot as plt
        
        examples = []
        for k,v in results.items():
            if v and "volume" in v and "energy_per_atom" in v and len(v["volume"]) == len(v["energy_per_atom"]):
                examples.append(k)
            else:
                print(f"Example {k} has no valid results. Need keys: volume and energy_per_atom.")
                
        n = len(examples)
        if n == 0:
            print("No example to plot.")
            return
        
        # save plot data to a csv file
        plot_data = ""
        for ik, k in enumerate(examples):
            fig, ax = plt.subplots(1,1,figsize=(8,5))
            plot_r = comm_eos.plot_eos_one(ax,results[k], ref_results.get(k,None), label_size=16, legend_size=14,x_label="Volume ($A^3$)",y_label="Energy ($eV/atom$)")
            title = "EOS of " + k
            if plot_r and plot_r.get("delta",None):
                delta = plot_r["delta"]
                title += f" (Delta={delta:.2e})"
            ax.set_title(title,fontsize=16)
            
            filename = k.strip("/").replace("/","_") + ".png"
            plt.tight_layout()
            plt.savefig(filename)
            plt.close()
                      
            plot_data += f"Example {k}:\n"
            plot_data += "x," + ",".join([f"{x:.2f}" for x in results[k]["volume"]]) + "\n"
            plot_data += "y," + ",".join([f"{x:.6f}" for x in results[k]["energy_per_atom"]]) + "\n"
            if ref_results.get(k,None):
                plot_data += "x_ref," + ",".join([f"{x:.2f}" for x in ref_results[k]["volume"]]) + "\n"
                plot_data += "y_ref," + ",".join([f"{x:.6f}" for x in ref_results[k]["energy_per_atom"]]) + "\n"
            if plot_r:
                if plot_r["fit"]:
                    v0, e0, fitx, fity, b0, bp, res = plot_r["fit"]
                    plot_data += "fitx," + ",".join([f"{x:.2f}" for x in fitx]) + "\n"
                    plot_data += "fity," + ",".join([f"{x:.6f}" for x in fity]) + "\n"
                    plot_data += f"v0,{v0:.2f}\n"
                    plot_data += f"e0,{e0:.6f}\n"
                    plot_data += f"b0,{b0:.6f}\n"
                    plot_data += f"bp,{bp:.6f}\n"
                    plot_data += f"res,{res:.6e}\n"
                if plot_r["fit_ref"]:
                    v0_ref, e0_ref, fitx_ref, fity_ref, b0_ref, bp_ref, res_ref = plot_r["fit_ref"]
                    plot_data += "fitx_ref," + ",".join([f"{x:.2f}" for x in fitx_ref]) + "\n"
                    plot_data += "fity_ref," + ",".join([f"{x:.6f}" for x in fity_ref]) + "\n"
                    plot_data += f"v0_ref,{v0_ref:.2f}\n"
                    plot_data += f"e0_ref,{e0_ref:.6f}\n"
                    plot_data += f"b0_ref,{b0_ref:.6f}\n"
                    plot_data += f"bp_ref,{bp_ref:.6f}\n"
                    plot_data += f"res_ref,{res_ref:.6e}\n"
                if plot_r.get("delta",None):
                    plot_data += f"delta,{plot_r['delta']:.6e}\n"
            plot_data += "\n"
        
        with open("plot_data.csv","w") as f:
            f.write(plot_data)  
        
        def flat_list(l):
            a = []
            if not isinstance(l,list):
                return [l]
            for i in l:
                if isinstance(i,list):
                    a += flat_list(i)
                else:
                    a.append(i)
            return a
        
        if one_plot_file:
            ncol = int(n ** 0.5)   
            nrow = int(n / ncol) 
            while ncol * nrow < n:
                ncol += 1
            figs, axs = plt.subplots(nrow, ncol, figsize=(ncol*8, nrow*5))
            axs_flat = flat_list(axs)
            
            for ik, k in enumerate(examples):
                ax = axs_flat[ik]
                plot_r = comm_eos.plot_eos_one(ax,results[k], ref_results.get(k,None), label_size=16, legend_size=14,x_label="Volume ($A^3$)",y_label="Energy ($eV/atom$)")
                title = "EOS of " + k
                if plot_r:
                    delta = plot_r.get("delta",None)
                    if delta != None:
                        title += f" (Delta={delta:.2e})"
                ax.set_title(title)
            plt.tight_layout()
            plt.savefig(one_plot_file)
            plt.close()

        
            
            
            
    
    
        

    
        
        
        
        
        
        