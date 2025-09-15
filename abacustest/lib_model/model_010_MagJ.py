import os,sys,glob,json,shutil,copy,traceback
from abacustest import constant
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,WriteInput,WriteKpt
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.prepare import PrepareAbacus
import numpy as np
from ..model import Model
from . import comm,comm_magj
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

PREPARE_DOC="""
    This model is used to prepare the calculation of the magnetic exchange interactions (J) of two atom.
    You should prepare a file named magj.txt in the folder, which define the atoms to calculate the J.
    Each line of the file defines a pair of atoms with format: <Atomlabel1> <AtomIndex1> <Atomlabel2> <AtomIndex2>
    The AtomIndex is start from 1 in each type of atom.  
"""

POST_DOC="""
    This model is used to postprocess the calculation of the magnetic exchange interactions (J) of two atom.
    You can get the J from the output files.
"""

class magj(Model):
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
        return "magj"
    
    @staticmethod
    def description():
        '''
        Description of the model
        '''
        return "Calculate the magnetic exchange interactions of two atom"
    
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
        parser.add_argument('-d', '--step', type=float, default=1,help="the tilt step. Default 1 degree.")
        parser.add_argument('-n', '--number', type=int, default=5,help="the tilt times. Default is 5 times. The maximum tilt degree is step*number")
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is current folder')
        parser.add_argument('-f', '--file', type=str, default="magj.txt",help="the file name of atom info, which define which two atoms to calculate the J. Default is magj.txt")
        
        return parser
    
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        jobs = ["."] if len(params.job) == 1 else params.job[1:]
        infof = "magj.txt"
        subfolders = PrepareMagJ(infof,jobs,params.step,params.number).run()

        if subfolders:
            setting = {
                "save_path": "results",
                "bohrium_group_name": "MagJ",
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
            comm.doc_after_prepare("MagJ", subfolders, ["setting.json"],has_prepare=False)
    
    @staticmethod
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand'''
        
        parser.description = POST_DOC
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is current folder.')
        return parser

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        jobs = ["."] if len(params.job) == 1 else params.job[1:]
        PostProcessMagJ(jobs,"magj").run()


def Condition2folder(subfolder, label1, label1idx, label2, label2idx, angle, case):   
    # CASE0: magj_Fe1_Fe2_angle1_1
    # CASE1: magj_Fe1_Fe2_angle1_2
    # CASE2: magj_Fe1_Fe2_angle1_3
    # 1: tilt atom1, 2: tilt atom2, 3: tilt atom1 and atom2 together
    # check: abacustest/lib_model/comm_magj.py
    angle = str(angle).rstrip("0").rstrip(".")
    return f"{subfolder}_{label1}{label1idx+1}_{label2}{label2idx+1}_angle{angle}_{case+1}"  

def Folder2Condition(folder):
    def split_label(label):
        if label[1].isdigit():
            return label[0],int(label[1])-1
        else:
            return label[:2],int(label[2:])-1
    try:
        ks = folder.split("_")
        subfolder = ks[0]
        label1, idx1 = split_label(ks[1])
        label2, idx2 = split_label(ks[2])
        angle = float(ks[3][5:])
        case = int(ks[4])-1
        return (subfolder, label1, idx1, label2, idx2, angle, case)
    except:
        traceback.print_exc()
        print("ERROR: folder name format error in",folder)
        return None

    
class PrepareMagJ:
    def __init__(self,infof,jobs,step,number):
        self.infof = infof
        self.jobs = jobs
        self.step = step
        self.number = number
        self.subfolder = "magj"
    
    def run(self):
        step = self.step
        num = self.number
        print("step (degree):",self.step)
        print("number of steps:",self.number)

        subfolders = []
        for job in self.jobs:
            for ijob in glob.glob(job):
                if not os.path.isdir(ijob):
                    continue
                # remove magj folders
                for ifile in glob.glob(os.path.join(ijob,self.subfolder+"*")):
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
        Fe 1 Fe 2
        
        Return:
        [["Fe", 0, "Fe", 1]]
        """
        def print_error(path,line):
            print(f"ERROR: {self.infof} format error in ",path)
            print("line:",line)
            print(PREPARE_DOC)
        if os.path.isfile(os.path.join(path,self.infof)):
            with open(os.path.join(path,self.infof)) as f:
                lines = f.readlines()
                info = []
                for line in lines:
                    if line.strip() == "":
                        continue
                    sline = line.strip().split()
                    if len(sline) == 4:
                        info.append([sline[0],int(sline[1])-1,sline[2],int(sline[3])-1])
                        if info[-1][1] < 0 or info[-1][3] < 0:
                            print_error(path,line)
                            return None
                    else:
                        print_error(path,line)
                        return None
                return info
        else:
            print(f"ERROR: {self.infof} not found in ",path)
            return None
        
    def prepare_one_stru_abacus(self,istru,final_path,otherfiles,input_param):
        os.makedirs(final_path,exist_ok=True)
        istru.write(os.path.join(final_path,"STRU"))
        # copy other files to final_path
        for file in otherfiles:
            os.symlink(file,os.path.join(final_path,os.path.basename(file)))

        # write input
        PrepareAbacus.WriteInput(input_param,os.path.join(final_path,"INPUT"))

    def preapre_inputs(self,infos, init_path, step = 1.0, num=3) :
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

        label = stru.get_label(total=True)
        atom_mag = stru.get_atommag()
        new_stru = copy.deepcopy(stru)

        # write original input
        self.prepare_one_stru_abacus(stru,os.path.join(init_path,self.subfolder),otherfiles,input_param)
        #with open(os.path.join(init_path,subfolder,"step.txt"),"w") as f: f.write(str(step))
        subfolders.append(os.path.join(init_path,self.subfolder))
        
        for info in infos:
            label1, index1, label2, index2 = info
            if label1 not in label or label2 not in label:
                print(f"ERROR: label {label1} or {label2} not found in STRU in {init_path}")
                continue
            global_idx1 = 0 
            global_idx2 = 0
            for i in range(len(label)):
                if label[i] == label1 and label[:i+1].count(label1) == index1+1:
                    global_idx1 = i
                if label[i] == label2 and label[:i+1].count(label2) == index2+1:
                    global_idx2 = i
            mag1 = atom_mag[global_idx1]
            mag2 = atom_mag[global_idx2]
            mags = comm_magj.prepare_atom_mag(mag1, mag2, step, num)
            for i, mag in enumerate(mags):
                angle = (i+1)*step
                
                for j in range(3):
                    atom_mag[global_idx1] = mag[j][0]
                    atom_mag[global_idx2] = mag[j][1]
                    new_stru.set_atommag(atom_mag)
                    subpath = os.path.join(init_path, Condition2folder(self.subfolder, label1, index1, label2, index2, angle,j ))
                    self.prepare_one_stru_abacus(new_stru,
                                                 subpath,
                                                 otherfiles,input_param)
                    subfolders.append(subpath)

        return subfolders
    
class PostProcessMagJ:
    def __init__(self,job,subfolder="magj"):
        self.job = job
        self.subfolder = "magj"
    
    def run(self):
        all_metrics = {}
        all_results = {}
        for ijob in self.job:
            if not os.path.isdir(ijob):
                continue
            print("postprocess in ",ijob)
            metrics = self.get_metrics(ijob)
            json.dump(metrics,open(os.path.join(ijob,"metrics.json"),"w"),indent=4)
            all_metrics.update(metrics)
            
            results = self.metrics2results(metrics)
            json.dump(results,open(os.path.join(ijob,"results.json"),"w"),indent=4)
            all_results[ijob] = results
        
        json.dump(all_metrics,open("metrics.json","w"),indent=4)
        json.dump(all_results,open("results.json","w"),indent=4)
        # print the J
        for k,v in all_results.items():
            for kk,vv in v.items():
                if kk == "original_energy (eV)":
                    continue
                print(f"{k} {kk} J:",vv["j (meV)"])
                
    
    def get_metrics(self, path):
        # find the subfolders in path
        # the subfolder name is like "magj_He1_He2_angle1_1"
        metrics = {}
        for ipath in glob.glob(os.path.join(path, self.subfolder+"*")):
            if not os.path.isdir(ipath):
                continue
            basename = os.path.basename(ipath)
            if not (basename == self.subfolder or len(basename.split("_")) == 5):
                continue
            
            r = RESULT(path=ipath, fmt="abacus")
            metrics[ipath] = {
                "converge": r["converge"],
                "energy": r["energy"],
                "drho_last": r["drho_last"],
                "denergy_last": r["denergy_last"],
                "denergy_womix_last": r["denergy_womix_last"],
                "atomlabel_list": r["atomlabel_list"],
                "total_time": r["total_time"] ,
                "mag": r["ds_mag"]
            }
        return metrics
    
    def get_atom_mag(self, atomlabel_list, mag, atomlabel, atomindex):
        """
        Get the magnetic moment of the atoms in atomlabel_list
        
        Parameters
        ----------
        atomlabel_list : list of str
            The list of atom label of all atoms
        mag : list of lists
            The magnetic moment of all atoms, which should respect the order of atomlabel_list
        atomlabel : str
            The atom label to get the magnetic moment
        atomindex : int
            The atom index of this atom type
        """
        idx = 0
        for i, label in enumerate(atomlabel_list):
            if label == atomlabel:
                if idx == atomindex:
                    return mag[i]
                idx += 1
        return None
    
    def metrics2results(self, metrics):
        """
        Convert the metrics to results
        
        Parameters  
        ----------
        metrics : dict
            The metrics of the calculation. The key is the subfolder name, and the value is the metrics of the calculation.
            
        Returns
        -------
        results : dict
            The results that is convient to calculate the J
        
        results = {
            "original_energy": -100,  # the energy of the original structure
            "Fe1_Fe2":{
                "angle": [1,2,3,4,5], # the target tilt angle
                "angle_cos": [[0.1,0.2,0.3],...] # the cos of the tilt angle, which is the average of the angle_cos_detail
                "angle_detail": [[1,1,1],...] # the detail of the tilt angle of atom1, atom2, and atom1 and atom2 together
                "angle_cos_detail": [[0.1,0.2,0.3],...] # the cos of the tilt angle of atom1, atom2, and atom1 and atom2 together
                "energy": [[-100,-101,-102],...] # the energy of tilt angle of atom1, atom2, and atom1 and atom2 together
                "delta_energy": [1,2,3,4,5], # (E_togather - E_original) - (E_atom1 - E_original) - (E_atom2 - E_original)
            }
        """
        metrics_new = {os.path.basename(k):v for k,v in metrics.items()}
        
        results = {
            "original_energy (eV)": metrics_new.get(self.subfolder,{}).get("energy",None)
        }
        for k,v in metrics_new.items():
            if k == self.subfolder:
                continue
            
            conditions = Folder2Condition(k)
            if conditions is None:
                continue
            
            _, label1, idx1, label2, idx2, tilt_angle, job_idx = conditions
            # job_idx: 0: tilt atom1, 1: tilt atom2, 2: tilt atom1 and atom2 together
            
            case_name = f"{label1}{idx1+1}_{label2}{idx2+1}"
            if case_name not in results:
                results[case_name] = {
                    "angle": [],
                    "angle_detail": [],
                    "angle_cos_detail": [],
                    "energy (eV)": [],
                }
            if tilt_angle not in results[case_name]["angle"]:
                results[case_name]["angle"].append(tilt_angle)
                results[case_name]["angle_detail"].append([None,None,None])
                results[case_name]["angle_cos_detail"].append([None,None,None])
                results[case_name]["energy (eV)"].append([None,None,None])
                
            angle_idx = results[case_name]["angle"].index(tilt_angle)
            mag1 = self.get_atom_mag(v["atomlabel_list"], v["mag"], label1, idx1)
            mag2 = self.get_atom_mag(v["atomlabel_list"], v["mag"], label2, idx2)
            if None not in [mag1,mag2]:
                results[case_name]["angle_detail"][angle_idx][job_idx] = np.arccos(np.dot(mag1,mag2)/np.linalg.norm(mag1)/np.linalg.norm(mag2))*180/np.pi
                results[case_name]["angle_cos_detail"][angle_idx][job_idx] = np.dot(mag1,mag2)/np.linalg.norm(mag1)/np.linalg.norm(mag2)
                
            results[case_name]["energy (eV)"][angle_idx][job_idx] = v["energy"]
        
        e_original = results["original_energy (eV)"]
        for k,v in results.items():
            if k == "original_energy (eV)":
                continue
            results[k]["angle_cos"] = []
            results[k]["delta_energy (eV)"] = []
            for i in range(len(v["angle"])):
                # calcualte the averge of the angle_cos_detail
                angle_cos_detail = v["angle_cos_detail"][i]
                if None in angle_cos_detail:
                    results[k]["angle_cos"].append(None)
                else:
                    results[k]["angle_cos"].append(np.mean(angle_cos_detail))
                
                # calculate the delta_energy    
                e1 = v["energy (eV)"][i][0]
                e2 = v["energy (eV)"][i][1]
                etot = v["energy (eV)"][i][2]
                if None in [e1,e2,etot,e_original]:
                    results[k]["delta_energy (eV)"].append(None)
                else:
                    results[k]["delta_energy (eV)"].append((etot - e_original) - (e1 - e_original) - (e2 - e_original))
                    
            results[k]["j (meV)"] = self.cal_j(results[k]["delta_energy (eV)"], results[k]["angle_cos"])
        return results
    
    def cal_j(self, deltae, cos_theta):
        """
        Calculate the J from the results.
        Will fit: \delta E = J * (1 - cos(\theta))
        
        Parameters
        ----------
        results : dict
            The results that is convient to calculate the J
        
        Returns
        -------
        j : dict
            The J of the atoms
        """
        
        # clean the None value
        de = []
        ct = []
        for i in range(len(deltae)):
            if deltae[i] is not None and cos_theta[i] is not None:
                de.append(deltae[i])
                ct.append(1 - cos_theta[i])
        
        if len(de) < 1:
            return None
        else:
            # fit the J by polyfit
            j = np.polyfit(ct,de,1)[0]
            return j * 1000
            
            
            

        
        