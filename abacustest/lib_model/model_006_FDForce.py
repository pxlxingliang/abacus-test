import os,sys,glob,json,shutil,copy
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,WriteInput,WriteKpt
from abacustest.prepare import PrepareAbacus
from ..model import Model
from . import comm

PREPARE_DOC="""
    This model is used to prepare the inputs to do the finite difference test of force. 
    Should prepare a file named 'info.txt' in the folder of inputs, which contains the label, number of atom, x/y/z component of cell vector, such as:
    '
    C 2 x y z
    H 1 x y 
    O 3 y
    '
    which means do the finite difference test of force 2-nd C on x y z component, and 1-st H on x y component, and 3-rd O on y component.
"""

class fdforce(Model):
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
        return "fdforce"
    
    @staticmethod
    def description():
        '''
        Description of the model
        '''
        return "finite difference of force"
    
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
        parser.add_argument('-d', '--step', type=float, default=0.01,help="the step of finite difference (unit: bohr), default 0.01")
        parser.add_argument('-n', '--number', type=int, default=2,help="the number of finite difference of one direction, default 2")
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is current folder')
        
        return parser
    
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        jobs = ["."] if len(params.job) == 1 else params.job[1:]
        infof = "info.txt"
        subfolders = PrepareFDForce(infof,jobs,params.step,params.number).run()

        if subfolders:
            setting = {
                "save_path": "results",
                "bohrium_group_name": "FDForce",
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
            comm.doc_after_prepare("FDForce", subfolders, ["setting.json"],has_prepare=False)
    
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
    
class PrepareFDForce:
    def __init__(self,infof,jobs,step,number):
        self.infof = infof
        self.jobs = jobs
        self.step = step
        self.number = number
    
    def run(self):
        step = self.step
        num = self.number
        print("step(bohr):",self.step)

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
        """
        def print_error(path,line):
            print("ERROR: info.txt format error in ",path)
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
                    if len(sline) < 3 or not sline[1].isdigit():
                        print_error(path,line)
                        continue
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
            print("ERROR: info.txt not found in ",path)
            return None
        
    def prepare_one_stru_abacus(self,istru,final_path,otherfiles,input_param):
        os.makedirs(final_path,exist_ok=True)
        istru.write(os.path.join(final_path,"STRU"))
        # copy other files to final_path
        for file in otherfiles:
            shutil.copy(file,final_path)

        # write input
        PrepareAbacus.WriteInput(input_param,os.path.join(final_path,"INPUT"))

    def preapre_inputs(self,infos, init_path, step = 0.001, num=3) :
        subfolders = []
        input_param = {"calculation": "scf","cal_force":1}
        otherfiles = []
        for ifile in os.listdir(init_path):
            if ifile == "INPUT":
                input_param = PrepareAbacus.ReadInput(os.path.join(init_path,"INPUT"))
                input_param["calculation"] = "scf"
                input_param["cal_force"] = 1
                continue
            elif ifile != "STRU" and not ifile.startswith(".") and os.path.isfile(os.path.join(init_path,ifile)):
                otherfiles.append(os.path.join(init_path,ifile))

        stru = AbacusStru.ReadStru(os.path.join(init_path,"STRU"))
        if not stru:
            print(f"ERROR: read STRU failed in {init_path}")
            return None

        label = stru.get_label()
        atom_num = [label.count(i) for i in label] # the number of atoms for each element
        coord = stru.get_coord(direct=False,bohr=True)
        cell = stru.get_cell(bohr=True)  # need use cell as the input of AbacusStru
        new_stru = copy.deepcopy(stru)
        new_stru.set_coord(coord,bohr=True,direct=False)

        # write original input
        self.prepare_one_stru_abacus(new_stru,os.path.join(init_path,"fdf"),otherfiles,input_param)
        with open(os.path.join(init_path,"fdf","step.txt"),"w") as f: f.write(str(step))
        subfolders.append(os.path.join(init_path,"fdf"))

        xyz_idx = {"x":0,"y":1,"z":2}
        # write other inputs
        for info in infos:
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

            for ixyz in xyz:
                for inum in range(num):
                    subfolder_name1 = f"fdf_{label_name}_{label_idx}_{ixyz}_+_{inum+1}"
                    subfolder_name2 = f"fdf_{label_name}_{label_idx}_{ixyz}_-_{inum+1}"
                    subfolders.append(os.path.join(init_path,subfolder_name1))
                    subfolders.append(os.path.join(init_path,subfolder_name2))

                    # write stru
                    new_coord = copy.deepcopy(coord)

                    new_coord[atom_idx][xyz_idx[ixyz]] += step * (inum+1)
                    new_stru = copy.deepcopy(stru)
                    new_stru.set_coord(new_coord,bohr=True,direct=False)
                    self.prepare_one_stru_abacus(new_stru,
                            os.path.join(init_path,subfolder_name1),
                            otherfiles,
                            input_param)

                    new_coord = copy.deepcopy(coord)
                    new_coord[atom_idx][xyz_idx[ixyz]] -= step * (inum+1)
                    new_stru = copy.deepcopy(stru)
                    new_stru.set_coord(new_coord,bohr=True,direct=False)
                    self.prepare_one_stru_abacus(new_stru,
                            os.path.join(init_path,subfolder_name2),
                            otherfiles,
                            input_param)
        return subfolders