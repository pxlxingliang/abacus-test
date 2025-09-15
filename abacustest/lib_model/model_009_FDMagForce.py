import os,sys,glob,json,shutil,copy,traceback
from abacustest import constant
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,WriteInput,WriteKpt
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.prepare import PrepareAbacus
import numpy as np
from ..model import Model
from . import comm
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

PREPARE_DOC="""
    This model is used to do the finite difference test for magnetic force. 
    Should prepare a file named 'mag_force_info.txt' in the folder of inputs, which contains the label, number of atom, and the xyz component of the magnetic force.
    '
    C 2 x y z
    H 1 x y 
    O 3 y
    '
    which means do the finite difference test of magnetic force 2-nd C on x y z component, and 1-st H on x y component, and 3-rd O on y component.
"""
POST_DOC="""
    the subfolder of the job should be like:
    fdf: the original input
    fdf_<LABEL_NAME>_<LABEL_IDX>_<x/y/z>_<+/->_<inum>: the input with cell vector changed by +/- diff
    there should be a file named "step.txt" in the fdf folder, which contains the step of finite difference test.
"""

class fdmagforce(Model):
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
        return "fdmagforce"
    
    @staticmethod
    def description():
        '''
        Description of the model
        '''
        return "finite difference of magenetic force"
    
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
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand '''
        
        parser.description = PREPARE_DOC
        parser.add_argument('-d', '--step', type=float, default=0.1,help="the step of finite difference (unit: uB), default 0.1")
        parser.add_argument('-n', '--number', type=int, default=5,help="the number of finite difference of one direction, default 5")
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is current folder')
        parser.add_argument('-t', '--type', type=int, default=0,help="the type of finite difference. 0: cartesian, 1: angle1, 2: angle2, 3: norm. Default 0")
        
        return parser
    
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        jobs = ["."] if len(params.job) == 1 else params.job[1:]
        infof = "mag_force_info.txt"
        subfolders = PrepareFDMagForce(infof,jobs,params.step,params.number).run()

        if subfolders:
            setting = {
                "save_path": "results",
                "bohrium_group_name": "FDMagForce",
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
            comm.doc_after_prepare("FDMagForce", subfolders, ["setting.json"],has_prepare=False)
    
    @staticmethod
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand'''
        
        parser.description = POST_DOC
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is current folder.')
        parser.add_argument('-t', '--type',default=0,type=int ,help='the job type. 0:abaus. Default 0')
        return parser

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        jobs = ["."] if len(params.job) == 1 else params.job[1:]
        PostProcessFDMagForce(jobs,params.type).run()
    
class PrepareFDMagForce:
    def __init__(self,infof,jobs,step,number):
        self.infof = infof
        self.jobs = jobs
        self.step = step
        self.number = number
    
    def run(self):
        step = self.step
        num = self.number
        print("step(uB):",self.step)

        subfolders = []
        for job in self.jobs:
            for ijob in glob.glob(job):
                if not os.path.isdir(ijob):
                    continue
                # remove fdf folders
                for ifile in glob.glob(os.path.join(ijob,"fdf*")):
                    if os.path.isdir(ifile):
                        shutil.rmtree(ifile)
                print("prepare inputs in ",ijob)
                infos = self.read_info(ijob)
                if infos:
                    isub = self.preapre_inputs(infos,ijob,step,num=num)
                    if isub:
                        subfolders.extend(isub)
                else:
                    print("ERROR: preapre failed in",ijob)
        # write relative path of subfolders to example.txt
        return subfolders
    
    def read_info(self,path):
        """
        the txt.info should be like:
        C 2 x y z
        H 1 x y
        
        or 
        
        C 2
        
        Return:
        [["C",2,["x","y","z"]],["H",1,["x","y"]]]
        or 
        [["C", 2]]
        The label, number of atom, and the xyz component of the magnetic force.
        """
        def print_error(path,line):
            print("ERROR: mag_force_info.txt format error in ",path)
            print("line:",line)
            print(PREPARE_DOC)
        if os.path.isfile(os.path.join(path,self.infof)):
            with open(os.path.join(path,self.infof)) as f:
                lines = f.readlines()
                info = []
                # each atom has a line append to info with a lit like: ["Fe", 1, ["x", "y", "z"]]
                # if not define the component, then only the label and number
                for line in lines:
                    if line.strip() == "":
                        continue
                    sline = line.strip().split()
                    if len(sline) < 2 or not sline[1].isdigit():
                        print_error(path,line)
                        continue
                    if len(sline) == 2:
                        info.append([sline[0],int(sline[1]),["x","y","z"]])
                    else:
                        xyz = [i.lower() for i in sline[2:]]
                        # if xyz has other character, then error
                        has_error = False
                        for i in xyz:
                            if i not in ["x","y","z"]:
                                print_error(path,line)
                                has_error = True
                                break
                        if not has_error:
                            info.append([sline[0],int(sline[1]),xyz])
                return info
        else:
            print("ERROR: mag_force_info.txt not found in ",path)
            return None
        
    def prepare_one_stru_abacus(self,istru,final_path,otherfiles,input_param):
        os.makedirs(final_path,exist_ok=True)
        istru.write(os.path.join(final_path,"STRU"))
        # copy other files to final_path
        for file in otherfiles:
            os.symlink(file,os.path.join(final_path,os.path.basename(file)))

        # write input
        PrepareAbacus.WriteInput(input_param,os.path.join(final_path,"INPUT"))

    def preapre_inputs(self,infos, init_path, step = 0.001, num=3) :
        subfolders = []
        input_param = {"calculation": "scf"}
        otherfiles = []
        for ifile in os.listdir(init_path):
            if ifile == "INPUT":
                input_param = PrepareAbacus.ReadInput(os.path.join(init_path,"INPUT"))
                input_param["calculation"] = "scf"
                continue
            elif ifile != "STRU" and not ifile.startswith(".") and os.path.isfile(os.path.join(init_path,ifile)):
                otherfiles.append(os.path.abspath(os.path.join(init_path,ifile)))

        stru = AbacusStru.ReadStru(os.path.join(init_path,"STRU"))
        if not stru:
            print(f"ERROR: read STRU failed in {init_path}")
            return None

        label = stru.get_label()
        atom_num = [label.count(i) for i in label] # the number of atoms for each element
        atom_mag = stru.get_atommag()
        new_stru = copy.deepcopy(stru)

        # write original input
        self.prepare_one_stru_abacus(stru,os.path.join(init_path,"fdf"),otherfiles,input_param)
        with open(os.path.join(init_path,"fdf","step.txt"),"w") as f: f.write(str(step))
        subfolders.append(os.path.join(init_path,"fdf"))

        xyz_idx = {"x":0,"y":1,"z":2}
        # write other inputs
        for info in infos:
            if len(info) == 2:
                label_name, label_idx = info
                xyz = [None]
            else:
                label_name, label_idx, xyz = info
            if label_name not in label:
                print(f"ERROR: label {label_name} not in STRU")
                continue
            if label_idx > atom_num[label.index(label_name)]: 
                print(f"ERROR: label {label_name} has only {atom_num[label.index(label_name)]} atoms")
                continue
            
            atom_idx = 0
            for i in range(len(label)):
                if label[i] == label_name:
                    atom_idx += 1
                    if atom_idx == label_idx:
                        break
            atom_idx = i
            atom_mag_norm = norm = np.linalg.norm(atom_mag[atom_idx])

            for ixyz in xyz:
                for inum in range(num):
                    if ixyz == None:
                        subfolder_name1 = f"fdf_{label_name}_{label_idx}_norm_+_{inum+1}"
                        subfolder_name2 = f"fdf_{label_name}_{label_idx}_norm_-_{inum+1}"
                    else:
                        subfolder_name1 = f"fdf_{label_name}_{label_idx}_{ixyz}_+_{inum+1}"
                        subfolder_name2 = f"fdf_{label_name}_{label_idx}_{ixyz}_-_{inum+1}"
                        
                    subfolders.append(os.path.join(init_path,subfolder_name1))
                    subfolders.append(os.path.join(init_path,subfolder_name2))

                    # write stru
                    new_atommag = copy.deepcopy(atom_mag)
                    if ixyz == None:
                        new_atommag[atom_idx] = (np.array(atom_mag[atom_idx]) / atom_mag_norm * (atom_mag_norm + step * (inum+1))).tolist()
                    else:
                        new_atommag[atom_idx][xyz_idx[ixyz]] += step * (inum+1)
                    new_stru = copy.deepcopy(stru)
                    new_stru.set_atommag(new_atommag)
                    self.prepare_one_stru_abacus(new_stru,
                            os.path.join(init_path,subfolder_name1),
                            otherfiles,
                            input_param)
                    
                    new_atommag = copy.deepcopy(atom_mag)
                    if ixyz == None:
                        new_atommag[atom_idx] = (np.array(atom_mag[atom_idx]) / atom_mag_norm * (atom_mag_norm - step * (inum+1))).tolist()
                    else:
                        new_atommag[atom_idx][xyz_idx[ixyz]] -= step * (inum+1)
                    new_stru = copy.deepcopy(stru)
                    new_stru.set_atommag(new_atommag)
                    self.prepare_one_stru_abacus(new_stru,
                            os.path.join(init_path,subfolder_name2),
                            otherfiles,
                            input_param)
        return subfolders
    
    
class PostProcessFDMagForce:
    def __init__(self,jobs,job_type=0):
        self.jobs = jobs
        
        job_type_dict = {0:"abacus",1:"qe",2:"vasp"}
        #if job_type not in job_type_dict:
        if job_type != 0:
            print(f"ERROR: job type {job_type} not supported")
            job_type = 0
        self.job_type_str = job_type_dict.get(job_type,"abacus")
        
        pass
    
    def run(self):
        allresults,cases = self.get_results()
        json.dump(allresults,open("metrics.json","w"),indent=4)
        
        plot_results = self.gen_plot_results(allresults,cases)
        allpngs = self.plot_force(plot_results)
        
        if allpngs:
            json.dump({i[:-4]: {"type": "image", "file": i} for i in allpngs},open("supermetrics.json","w"),indent=4) 
    
    def get_results(self):
        """
        return allresults, cases
        allresults is a dict with key is all fdf jobs and value is a dict of key metrics
        cases is a dict with key is the case name and value is a dict of key metrics
        """
        allresults = {}
        cases = {}
        for ijob in self.jobs:
            for iijob in glob.glob(ijob):
                if not os.path.isdir(iijob):
                    continue
                if not os.path.isdir(os.path.join(iijob,"fdf")):
                    print(f"ERROR: {iijob} does not have fdf folder")
                    continue
                if not os.path.isfile(os.path.join(iijob,"fdf","step.txt")):
                    print(f"ERROR: {iijob} does not have step.txt")
                    continue
                
                step = float(open(os.path.join(iijob,"fdf","step.txt")).readline()) 
                print(f"get results in {iijob}")
                for job in glob.glob(os.path.join(iijob,"fdf*")):
                    job_name = os.path.basename(job)
                    job_name_list = job_name.split("_")
                    # job_name should be like fdf or fdf_Fe_x_1_+_1
                    if os.path.isdir(job) and (job_name == "fdf" or len(job_name_list) == 6):
                        r = RESULT(path=job,fmt=self.job_type_str)
                        allresults[job] = {
                            "converge": r["converge"],
                            "energy": r["energy"],
                            "drho_last": r["drho_last"],
                            "denergy_last": r["denergy_last"],
                            "denergy_womix_last": r["denergy_womix_last"],
                            "atomlabel_list": r["atomlabel_list"],
                            "total_time": r["total_time"] ,
                            "mag": r["ds_mag"],
                            "mag_force": r["ds_mag_force"]     
                        }
                        if len(job_name_list) == 6:
                            icasename = os.path.join(iijob,"_".join(job_name_list[1:4]))
                            if icasename not in cases:
                                cases[icasename] = {"step": step, "number": int(job_name_list[-1]),"label": job_name_list[1], "idx": int(job_name_list[2]), "xyz": job_name_list[3],
                                                    "label_idx":None}
                            if r["atomlabel_list"] != None and cases[icasename]["label_idx"] == None:
                                label_idx = None
                                idx = 0
                                for i in range(len(r["atomlabel_list"])):
                                    if r["atomlabel_list"][i] == job_name_list[1]:
                                        idx += 1
                                    if idx == int(job_name_list[2]):
                                        label_idx = i
                                        break
                                cases[icasename]["label_idx"] = label_idx
                                
                            if int(job_name_list[-1]) > cases[icasename]["number"]:
                                cases[icasename]["number"] = int(job_name_list[-1])
        return allresults,cases        

    def gen_fdf_name(self,example_name,case_name,step_num):
        # case name should be like: Fe_1_x
        if step_num == 0:
            return os.path.join(example_name,"fdf")
        elif step_num > 0:
            return os.path.join(example_name,f"fdf_{case_name}_+_{step_num}")
        else:
            return os.path.join(example_name,f"fdf_{case_name}_-_{-step_num}")

    def gen_plot_results(self,allresults,cases):
        if not allresults:
            return {}
        plot_results = {}  # results for FD at different STRU with same step
        
        for case in cases:
            example_name = os.path.dirname(case)
            case_name = os.path.basename(case)
            input_r = allresults.get(os.path.join(example_name,"fdf"),{})  # the results of input structure
            step = cases[case]["step"]
            maxnumber = cases[case]["number"]
            label_idx = cases[case]["label_idx"]
            xyz_idx = {"x":0, "y":1, "z":2}[cases[case]["xyz"]]
            if label_idx == None:
                print(f"ERROR: {case} is not successful")
                continue
            
            plot_results[case] = {"r1":{"pos": [], "energy": [], "mag_force": [], "fd_mag_force": []},
                                  "r2":{"step": [], "fd_mag_force": [], "mag_force": None if input_r.get("mag_force") == None else input_r.get("mag_force")[label_idx][xyz_idx]}}
            
            # plot_results1
            for i in range(-maxnumber+1,maxnumber):
                fdf_name = self.gen_fdf_name(example_name,case_name,i)
                fdf_name_p = self.gen_fdf_name(example_name,case_name,i+1)
                fdf_name_m = self.gen_fdf_name(example_name,case_name,i-1)
                mag_force = allresults.get(fdf_name,{}).get("mag_force")
                energy_p = allresults.get(fdf_name_p,{}).get("energy")
                energy_m = allresults.get(fdf_name_m,{}).get("energy")
                
                plot_results[case]["r1"]["pos"].append(i*step)
                plot_results[case]["r1"]["energy"].append(allresults.get(fdf_name,{}).get("energy"))
                plot_results[case]["r1"]["mag_force"].append(None if mag_force == None else mag_force[label_idx][xyz_idx])
                if energy_p != None and energy_m != None:
                    plot_results[case]["r1"]["fd_mag_force"].append((energy_m - energy_p) / (2 * step))
                else:
                    plot_results[case]["r1"]["fd_mag_force"].append(None)
            
            # plot_results2
            for i in range(1,maxnumber+1):
                fdf_name_p = self.gen_fdf_name(example_name,case_name,i)
                fdf_name_m = self.gen_fdf_name(example_name,case_name,-i)
                energy_p = allresults.get(fdf_name_p,{}).get("energy")
                energy_m = allresults.get(fdf_name_m,{}).get("energy")
                plot_results[case]["r2"]["step"].append(i*step*2)
                if energy_p != None and energy_m != None:
                    plot_results[case]["r2"]["fd_mag_force"].append((energy_m - energy_p) / (2 * step * i))
                else:
                    plot_results[case]["r2"]["fd_mag_force"].append(None)
            
            json.dump(plot_results[case],open(os.path.join(example_name,f"{case_name}.json"),"w"),indent=4)
                
        return plot_results

    def cal_rmsd(self,list1,list2):
        l1,l2 = comm.clean_none_list(list1,list2)
        if len(l1) == 0:
            return None
        return (sum([(i-j)**2 for i,j in zip(l1,l2)]) / len(l1))**0.5
    
    def plot_force(self, results):
        import matplotlib.pyplot as plt
        
        ub_unit = "$\mu$B"
        
        allpngs = []
        for case in results:
            distance = results[case]["r1"]["pos"]
            energy = results[case]["r1"]["energy"]
            mag_force = results[case]["r1"]["mag_force"]
            fd_mag_force = results[case]["r1"]["fd_mag_force"]
            
            d,e,f,fd = comm.clean_none_list(distance,energy,mag_force,fd_mag_force)
            x, y = comm.clean_none_list(results[case]["r2"]["step"],results[case]["r2"]["fd_mag_force"])
            if len(d) == 0 or len(x) == 0:
                continue
            
            fig,ax = plt.subplots(1,2,figsize=(12,5))
            ax1 = ax[0]
            #fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()
            ax1.plot(d,e,"b+-",label="Energy (eV)")
            ax2.plot(d,f,"ro-",label=f"Analytic(eV/{ub_unit})")
            ax2.plot(d,fd,"m*-",label=f"Finite difference (eV/{ub_unit})")
            
            ymax1 = max(e)
            ymin1 = min(e)
            ymax2 = max(max(f),max(fd))
            ymin2 = min(min(f),min(fd))
            ax1.set_xlabel(f"Magnetic moment variation ({ub_unit})")
            ax1.set_ylabel("Energy (eV)",color="b")
            ax2.set_ylabel(f"Magnetic Force (eV/{ub_unit})",color="r")
            ax1.spines['left'].set_color('blue')
            ax1.tick_params(axis='y', colors='blue')
            ax2.spines['right'].set_color('red')
            ax2.tick_params(axis='y', colors='red')
            ax1.set_ylim(ymin1-(ymax1-ymin1)*0.1,ymax1+(ymax1-ymin1)*0.2)
            ax2.set_ylim(ymin2-(ymax2-ymin2)*0.1,ymax2+(ymax2-ymin2)*0.2)
            rmsd = self.cal_rmsd(f,fd)
            ax1.set_title(case + f"(FD-Ana RMSD={rmsd:.2e})")
            print(case + f"(FD-Ana RMSD={rmsd:.2e} eV/uB)")
            ax1.legend(loc="upper left")
            ax2.legend(loc="upper right")
            ax2.grid(True)

            ax3 = ax[1]
            # plot force0-FD vs distance0
            y_ref = results[case]["r2"]["mag_force"]
            ymin = min(y_ref, min(y))
            ymax = max(y_ref, max(y))
            # sort x,y by x
            ax3.plot(x,y,"m*-",label=f"Finite Difference(eV/{ub_unit})")
            if y_ref != None:
                ax3.axhline(y_ref,ls="--",color="red",label=f"Analytic(eV/{ub_unit})")
            ax3.set_xlabel(f"Step size of finite difference ({ub_unit})")
            ax3.set_ylabel(f"Magnetic force (eV/{ub_unit})")
            ax3.set_title(case + f"(FD VS step-size)")
            ax3.legend(loc="upper right")
            ax3.grid(True)
            ax3.set_ylim(ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.2)

            plt.subplots_adjust(right=0.85) 
            png = case + ".png"
            plt.tight_layout()
            plt.savefig(png)
            plt.close()
            allpngs.append(png)
        return allpngs