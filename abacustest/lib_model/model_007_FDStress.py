import os,sys,glob,json,shutil,copy
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,WriteInput,WriteKpt
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.prepare import PrepareAbacus
from ..model import Model
from . import comm
from abacustest.lib_prepare.comm import kspacing2kpt

import matplotlib.pyplot as plt
import numpy as np
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

class fdstress(Model):
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
        return "fdstress"
    
    @staticmethod
    def description():
        '''
        Description of the model
        '''
        return "finite difference of stress"
    
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
        
        parser.description = "Do the finite difference of stress test."
        parser.add_argument('-d', '--step', type=float, default=0.0001,help="the step to change the cell, default 0.0001. This is a coefficient, the real step is step*cell_vector")
        parser.add_argument('-n', '--number', type=int,  default=5,help='the number of points to calculate the stress tensor, default 5')
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is current folder.')
        parser.add_argument('--comp',type=str,default=["11","12","13","22","23","33"], nargs="*" ,help='the stress components to test, default is 11 12 13 22 23 33.')
        return parser
    
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        jobs = ["."] if len(params.job) == 1 else params.job[1:]

        subfolders = PrepareFDStress(jobs,params.step,params.number,params.comp).run()

        if subfolders:
            setting = {
                "save_path": "results",
                "bohrium_group_name": "FDStress",
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
            comm.doc_after_prepare("FDStress", subfolders, ["setting.json"],has_prepare=False)
    
    @staticmethod
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand'''
        parser.description = "Postprocess the finite difference of stress test. There should have several cell_* folders in the job folder."
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus results, default is current folder.')
        parser.add_argument('--read',type=int, const=1, default=0,nargs="?",help='if read the results.json in each job folder, default is 0, no read.')
        parser.add_argument('-t', '--type',default=0,type=int ,help='the job type. 0:abaus, 1:qe. Default 0')
        return parser

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        jobs = ["."] if len(params.job) == 1 else params.job[1:]
        PostProcessFDStress(jobs,params.read,params.type).run()
    
class PrepareFDStress:
    def __init__(self,jobs,step,number,comp=["11","12","13","22","23","33"]):
        self.jobs = jobs
        self.step = step
        self.number = number
        self.comp = comp
    
    def run(self):
        print("prepare inputs")
        print("step:",self.step)
        print("number:",self.number)
        subfolders = []
        for job in self.jobs:
            for ijob in glob.glob(job):
                # remove fdf folders
                for ifile in glob.glob(os.path.join(ijob,"cell_*")):
                    if os.path.isdir(ifile):
                        shutil.rmtree(ifile)
                print("prepare inputs in ",ijob)
                ifolders = self.preapre_stru(ijob)
                if ifolders:
                    subfolders.extend(ifolders)
        return subfolders
    
    def gen_kpt(self, init_path, kspacing, cell):
        kpt = kspacing2kpt(kspacing, cell)
                
    def preapre_stru(self,init_path):
        input_param = {"calculation": "scf","cal_stress":1}
        otherfiles = []
        for ifile in os.listdir(init_path):
            if ifile == "INPUT":
                input_param = PrepareAbacus.ReadInput(os.path.join(init_path,"INPUT"))
                input_param["calculation"] = "scf"
                input_param["cal_stress"] = 1
                continue
            elif ifile != "STRU" and not ifile.startswith(".") and os.path.isfile(os.path.join(init_path,ifile)):
                otherfiles.append(os.path.abspath(os.path.join(init_path,ifile)))
                
        stru = AbacusStru.ReadStru(os.path.join(init_path,"STRU"))
        if not stru:
            print(f"ERROR: read STRU failed in {init_path}")
            return None
                
        cell = stru.get_cell(bohr=True)  # need use cell as the input of AbacusStru
        #print(label,atom_num,pp,orb,mass,dpks,coord,cell)
        
        # we need to fix the kpt.
        kspacing = input_param.pop("kspacing",None)
        if kspacing is not None:
            kpt = kspacing2kpt(kspacing, cell)
            kptf = (os.path.abspath(os.path.join(init_path,"KPT")))
            WriteKpt(kpt,kptf)
            if kptf not in otherfiles:
                otherfiles.append(kptf)

        write_0 = False
        subfolders = []
        for i in range(3):
            for j in range(3):
                # i is the index of cell vector, j is the index of x/y/z component of cell vector
                for k in range(self.number*-1,self.number+1):
                    # k is the index of diff coef
                    if k == 0 and write_0:
                       continue
                    if f"{i+1}{j+1}" not in self.comp:
                        continue
                    sub_folder = os.path.join(init_path,f"cell_{i+1}_{j+1}_{k}")
                    if k == 0:
                        sub_folder = os.path.join(init_path,f"cell_0")
                    os.makedirs(sub_folder,exist_ok=True)
                    
                    if k == 0:
                        json.dump({"step":self.step,"number":self.number},open(os.path.join(sub_folder,"param.json"),"w"))
                        write_0 = True

                    new_cell = copy.deepcopy(cell)
                    for m in range(3):
                        new_cell[m][i] += k * self.step * cell[m][j]
                    new_stru = copy.deepcopy(stru)
                    new_stru.set_cell(new_cell,bohr=True,change_coord=True)
                    new_stru.write(os.path.join(sub_folder,"STRU"))
                    WriteInput(input_param,os.path.join(sub_folder,"INPUT"))
                    
                    # copy other files to sub_folder
                    for file in otherfiles:
                        if not os.path.isfile(os.path.join(sub_folder,os.path.basename(file))):
                            os.symlink(file,os.path.join(sub_folder,os.path.basename(file)))
                    subfolders.append(sub_folder)
        return subfolders
    
class PostProcessFDStress:
    def __init__(self,jobs,readr,jobtype=0):
        self.jobs = jobs
        self.readr = readr
        
        jobtype_dict = {0:"abacus",1:"qe"}
        if jobtype not in jobtype_dict:
            print(f"jobtype {jobtype} is not supported.")
            jobtype = 0
        self.jobtype = jobtype_dict[jobtype]
        
        unit_trans = {
            0: 1602.1757722389546,  # in abacus, the stress transfer to kbar is multiply a coef = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
                 # = 4.35974394e-18 / 2 / pow(0.529177e-10, 3) * 1.0e-8
                 # Ry to eV = 13.605698
                 # bohr3 to angstrom3 = pow(ModuleBase::BOHR_TO_A, 3) = 0.5291770**3
                 # so, when I get energy in eV and volume in angstrom3, I can get the stress = energy / 13.605698 / (volume/0.5291770**3)  * 4.35974394e-18 / 2 / pow(0.529177e-10, 3) * 1.0e-8
                 # = energy / volume * 1602.1757722389546
            1: 1602.1777023087203  # in QE, the energy and stress are transfered by abacustest
                                    # KBAR2HARTREEPERBOHR3 = 3.398927420868445E-6
                                    # BOHR2A = 0.52917721092, HARTREE2EV = 27.211396132
                                    # so, stress = energy / 27.211396132 / (volume / 0.52917721092**3) / 3.398927420868445E-6 = energy / volume * 1602.1777023087203
        }
        self.unitcoef = unit_trans[jobtype]
        
    def run(self):
        alljobs = []
        for job in self.jobs:
            for ijob in glob.glob(job):
                if os.path.isdir(ijob) and ijob not in alljobs:
                    alljobs.append(ijob) 
        sm = {}   
        metrics = {}       
        for ijob in alljobs:
            if os.path.isdir(ijob):
                print("postprocess in ",ijob)
                jsonf = os.path.join(ijob,"results.json")
                if self.readr and os.path.isfile(jsonf):
                    values = json.load(open(jsonf,"r"))
                else:
                    if not os.path.isfile(os.path.join(ijob,"cell_0","param.json")):
                        print(f"param.json is not exist in {ijob}/cell_0")
                        continue
                    param = json.load(open(os.path.join(ijob,"cell_0","param.json"),"r"))
                    values = {"param": param}
                    for cell in glob.glob(f"{ijob}/cell_*"):
                        values[os.path.basename(cell)] = self.get_value(cell)
                    json.dump(values,open(jsonf,"w"),indent=4)
                
                for k,v in values.items():
                    if k != "param":
                        metrics[os.path.join(ijob,k)] = v
                    
                if "step" not in values.get("param",{}) or "number" not in values.get("param",{}):
                    print(f"param.json is not correct in {ijob}. Need step and number.")
                    continue
                diff = values.get("param",{}).get("step",0.0001)
                number = values.get("param",{}).get("number",6)

                results1 = self.gen_result1(values,diff,number)
                results2 = self.gen_result2(values,diff,number)
                json.dump(results1,open(f"{ijob}/results1.json","w"),indent=4)
                json.dump(results2,open(f"{ijob}/results2.json","w"),indent=4)
                self.plot1(results1,f"{ijob}/stress_tensor1.png")
                self.plot2(results2,f"{ijob}/stress_tensor2.png")
                sm[f"{ijob}/stress_tensor1"] = {
                    "type": "image",
                    "file": f"{ijob}/stress_tensor1.png"
                }
                sm[f"{ijob}/stress_tensor2"] = {
                    "type": "image",
                    "file": f"{ijob}/stress_tensor2.png"
                }
                self.write2csv(ijob,results1,results2)

        # write the supermetrics.json
        json.dump(metrics,open("metrics.json","w"),indent=4)
        json.dump(sm,open("supermetrics.json","w"),indent=4)
        
    def get_value(self,ipath):
        if not os.path.isdir(ipath):
            return {}

        iresult = RESULT(fmt=self.jobtype,path=ipath)
        if iresult["stress"] != None:
            stress= [[iresult["stress"][0],iresult["stress"][1],iresult["stress"][2]],
                      [iresult["stress"][3],iresult["stress"][4],iresult["stress"][5]],
                      [iresult["stress"][6],iresult["stress"][7],iresult["stress"][8]]]
        else:
            stress = None
        return {"stress": stress,
                "energy": iresult["energy"],
                "version": iresult["version"],
                "volume": iresult["volume"],
                "converge": iresult["converge"],
                "drho_last": iresult["drho_last"],
                "denergy_last": iresult["denergy_last"],
                "total_time": iresult["total_time"]}

    def gen_result1(self,value_dict,diff,number):
        # calculate the stress tensor of -number+1, ..., number-1,
        # and compare the difference between the calculated stress tensor and the stress tensor from abacus.
        # will return a dict
        results = {
            "x":[i*diff for i in range(-number+1,number)],
            "11":{"analytical":[],"numerical":[]},
            "12":{"analytical":[],"numerical":[]},
            "13":{"analytical":[],"numerical":[]},
            "21":{"analytical":[],"numerical":[]},
            "22":{"analytical":[],"numerical":[]},
            "23":{"analytical":[],"numerical":[]},
            "31":{"analytical":[],"numerical":[]},
            "32":{"analytical":[],"numerical":[]},
            "33":{"analytical":[],"numerical":[]}
        }

        for i in range(3):
            for j in range(3):
                ij = f"{i+1}{j+1}"
                for k in range(number*-1 + 1,number):
                    cell_name = f"cell_{i+1}_{j+1}_{k}"
                    if k == 0:
                        cell_name = "cell_0"
                    volume = value_dict.get(cell_name,{}).get("volume",None)
                    stress_abacus = value_dict.get(cell_name,{}).get("stress",None)
                    if stress_abacus != None:
                        stress_abacus = stress_abacus[i][j]
                    results[ij]["analytical"].append(stress_abacus)

                    cell1 = f"cell_{i+1}_{j+1}_{k-1}"
                    cell2 = f"cell_{i+1}_{j+1}_{k+1}"

                    if k -1 == 0:
                        cell1 = "cell_0"
                    if k +1 == 0:
                        cell2 = "cell_0"
                    v1 = value_dict.get(cell1,{}).get("volume",None)
                    v2 = value_dict.get(cell2,{}).get("volume",None)
                    #print(i,j,v1,v2)
                    energy_k1 = value_dict.get(cell1,{}).get("energy",None)
                    energy_k2 = value_dict.get(cell2,{}).get("energy",None)
                    
                    stress_finite_diff = None
                    if energy_k1 == None:
                        print(f"Get energy from cell_{i+1}_{j+1}_{k-1} failed.")
                    if energy_k2 == None:
                        print(f"Get energy from cell_{i+1}_{j+1}_{k+1} failed.")
                    if volume == None:
                        print(f"Get volume from cell_{i+1}_{j+1}_{k} failed.")
                        
                    if not (energy_k1 == None or energy_k2 == None or volume == None):
                        if i == j:
                            rdiff = diff / (1 + k*diff)
                        else:
                            rdiff = diff
                        #stress_finite_diff = (energy_k1 - energy_k2) / (2*diff) / volume * 1602.17662
                        stress_finite_diff = (energy_k1 - energy_k2) / (2*rdiff) / volume * self.unitcoef

                    results[ij]["numerical"].append(stress_finite_diff)
        return results              

    def gen_result2(self,value_dict,diff,number):
        # calculate the stress tensor with diff coef 0, and based on different diff coef k and -k.
        # will return the deviation of stress tensor (numerical - analytical).
        results = {
            "x": [i*diff*2 for i in range(1,number+1)],
            "11":[],"12":[],"13":[],
            "21":[],"22":[],"23":[],
            "31":[],"32":[],"33":[],
        }
        volume = value_dict.get("cell_0",{}).get("volume",None)
        stress_analytical = value_dict.get("cell_0",{}).get("stress",None)

        for i in range(3):
            for j in range(3):
                ij = f"{i+1}{j+1}"
                for k in range(1,number+1):

                    e1 = value_dict.get(f"cell_{i+1}_{j+1}_{-k}",{}).get("energy",None)
                    e2 = value_dict.get(f"cell_{i+1}_{j+1}_{k}",{}).get("energy",None)
                    v1 = value_dict.get(f"cell_{i+1}_{j+1}_{-k}",{}).get("volume",None)
                    v2 = value_dict.get(f"cell_{i+1}_{j+1}_{k}",{}).get("volume",None)

                    if e1 == None or e2 == None or volume == None or stress_analytical == None:
                        deviation = None
                    else:
                        stress_finite_diff = (e1 - e2) / (2*k*diff) / volume * self.unitcoef
                        #stress_finite_diff = (e1 - e2) /(v2 - v1) * 1602.17662
                        deviation = stress_finite_diff - stress_analytical[i][j]

                    results[ij].append(deviation)
        return results

    def plot1(self,results,fname):
        # plot two figures, one is the both analytical and numerical stress tensor, the other is the deviation.
        # x axis is the diff coef, y axis is the stress tensor, each component of stress tensor is one legend.
        # Plot two figures in one figure, and the first figure has 9*2 legends, the second figure has 9 legends.
        fig, axs = plt.subplots(2,1,figsize=(12,10),sharex=True)
        axs[0].set_title("stress tensor",fontsize=14)
        axs[1].set_title("deviation (F.D. - Ana.))",fontsize=14)
        #axs[0].set_xlabel("diff coef")
        axs[1].set_xlabel("configurations (cell deformation degree)(F.D.=Finite Difference, Ana.=Analytical)",fontsize=14)
        axs[0].set_ylabel("stress tensor (kbar)",fontsize=14)
        axs[1].set_ylabel("deviation (kbar)",fontsize=14)
        colors = ["red","orange","purple","brown","green","pink","gray","olive","blue"]
        for i in range(3):
            for j in range(3):
                ij = f"{i+1}{j+1}"
                x,y1,y2 = comm.clean_none_list(results["x"],results[ij]["analytical"],results[ij]["numerical"])
                if len(x) > 0:
                    axs[0].plot(x,y1,label=f"Ana. {ij}",linestyle="--",color=colors[i*3+j],marker="+")
                    axs[0].plot(x,y2,label=f"F.D. {ij}",linestyle="-",color=colors[i*3+j],marker="o")
                    deviation = np.array(y1) - np.array(y2)
                    axs[1].plot(x,deviation,label=f"Dev. {ij}",linestyle="-",color=colors[i*3+j],marker="o")
        axs[0].legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.)
        axs[1].legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.)
        plt.savefig(fname)
        plt.close()
        return


    def plot2(self,results,fname):
        # plot one figure, x axis is the diff coef, y axis is the deviation of stress tensor, each component of stress tensor is one legend.
        # Plot one figure, and the figure has 9 legends.

        fig, axs = plt.subplots(1,1,figsize=(12,5))
        axs.set_title("deviation (F.D. - Ana.))",fontsize=14)
        axs.set_xlabel("finite difference step (cell deformation degree) (F.D.=Finite Difference, Ana.=Analytical)",fontsize=14)
        axs.set_ylabel("deviation (kbar)",fontsize=14)
        colors = ["red","orange","purple","brown","green","pink","gray","olive","blue"]
        for i in range(3):
            for j in range(3):
                ij = f"{i+1}{j+1}"
                x,y = comm.clean_none_list(results["x"],results[ij])
                if len(x) > 0:
                    axs.plot(x,y,label=f"Dev. {ij}",linestyle="-",color=colors[i*3+j],marker="o")
        axs.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.)
        plt.savefig(fname)
        plt.close()
        return    

    def write2csv(self,ijob,results1,results2):
        csvf = os.path.join(ijob,"stress_tensor1.csv")
        with open(csvf,"w") as f:
            f.write("The stress tensor calculated by analytical and finite difference method. Each x related to one configuration, and x is offset of cell vector.\n")
            f.write("x,")
            for i in range(3):
                for j in range(3):
                    f.write(f"analytical_{i+1}{j+1},")
                    f.write(f"numerical_{i+1}{j+1},")
            f.write("\n")
            for i in range(len(results1["x"])):
                f.write(f"{results1['x'][i]},")
                for j in range(3):
                    for k in range(3):
                        f.write(f"{results1[f'{j+1}{k+1}']['analytical'][i]},")
                        f.write(f"{results1[f'{j+1}{k+1}']['numerical'][i]},")
                f.write("\n")

        # save results2 to csv
        csvf = os.path.join(ijob,"stress_tensor2.csv")
        with open(csvf,"w") as f:
            f.write("The stress tensor calculated by finite difference method. x is the offset of cell vector.\n")
            f.write("x,")
            for i in range(3):
                for j in range(3):
                    f.write(f"{i+1}{j+1},")
            f.write("\n")
            for i in range(len(results2["x"])):
                f.write(f"{results2['x'][i]},")
                for j in range(3):
                    for k in range(3):
                        f.write(f"{results2[f'{j+1}{k+1}'][i]},")
                f.write("\n")        