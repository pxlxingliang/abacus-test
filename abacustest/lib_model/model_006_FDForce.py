import os,sys,glob,json,shutil,copy,traceback
from abacustest import constant
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,WriteInput,WriteKpt
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.prepare import PrepareAbacus
from ..model import Model
from . import comm

PREPARE_DOC="""
    This model is used to do the finite difference test of force. 
    Should prepare a file named 'info.txt' in the folder of inputs, which contains the label, number of atom, x/y/z component of cell vector, such as:
    '
    C 2 x y z
    H 1 x y 
    O 3 y
    '
    which means do the finite difference test of force 2-nd C on x y z component, and 1-st H on x y component, and 3-rd O on y component.
"""
POST_DOC="""
    the subfolder of the job should be like:
    fdf: the original input
    fdf_<LABEL_NAME>_<LABEL_IDX>_<x/y/z>_<+/->_<inum>: the input with cell vector changed by +/- diff
    there should be a file named "step.txt" in the fdf folder, which contains the step of finite difference test.
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
        
        parser.description = POST_DOC
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is current folder.')
        return parser

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        jobs = ["."] if len(params.job) == 1 else params.job[1:]
        PostProcessFDForce(jobs).run()
    
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
    
    
class PostProcessFDForce:
    def __init__(self,jobs):
        self.jobs = jobs
        pass
    
    def run(self):
        allresults = {}

        for ijob in self.jobs:
            for iijob in glob.glob(ijob):
                if not os.path.isdir(iijob):
                    continue
                if not os.path.isdir(os.path.join(iijob,"fdf")):
                    print(f"ERROR: {iijob} does not have fdf folder")
                    continue
                print(f"get results in {iijob}")
                allresults[iijob.rstrip("/")] = self.get_onejob(iijob)

        allpngs = self.plot_force(allresults)
        metrics = {}
        for ik,iv in allresults.items():
            if iv == None:
                continue
            for jk,jv in iv.items():
                if jk in ["force","converge","labels","step"]:
                    continue
                
                metrics[ik + "/" + jk] = {
                    "Ana-force(eV/A)": jv["abacus"],
                    "FD-force(eV/A)": jv["fd"],
                    "Ana-FD(eV/A)": None if jv["abacus"] == None or jv["fd"] == None else jv["abacus"] - jv["fd"],
                    "distance(A)": jv["distance(A)"],
                    "energy(eV)": jv["energy"],
                    "energy+(eV)": jv["energy+"],
                    "energy-(eV)": jv["energy-"],
                    "force+(eV/A)": jv["force+"],
                    "force-(eV/A)": jv["force-"],
                    "converge(input)": iv["converge"],
                    "converge-FD+": jv["converge+"],
                    "converge-FD-": jv["converge-"],
                    "step(bohr)": iv["step"]  
                }
        json.dump(metrics,open("metrics.json","w"),indent=4)   
        if len(allpngs) > 0:
            json.dump({i[:-4]: {"type": "image", "file": i} for i in allpngs},open("supermetrics.json","w"),indent=4)      
    
    def get_onejob(self,job):
        # read results from fdf
        result = RESULT(path=os.path.join(job,"fdf"),fmt="abacus")
        force = result["force"]
        converge = result["converge"]
        labels = result["atomlabel_list"]
        energy = result["energy"]
        allresults = {"force": force, "converge": converge, "labels": labels}

        if not os.path.isfile(os.path.join(job,"fdf","step.txt")):
            print(f"ERROR: {job} does not have step.txt")
            return None
        step = float(open(os.path.join(job,"fdf","step.txt")).readline())
        allresults["step"] = step

        for i in glob.glob(os.path.join(job,"fdf_*_+_*")):
            ii = i.split("/")[-1].split("_")

            if len(ii) != 6:
                continue
            try:
                ilabel = ii[1]
                iidx = int(ii[2])
                xyz = ii[3]
                xyz_idx = {"x":0,"y":1,"z":2}[xyz]
                atom_idx = 0
                for j in range(len(labels)):
                    if labels[j] == ilabel:
                        atom_idx += 1
                    if atom_idx == iidx:
                        break
                atom_idx = j
                step_num = int(ii[5])

                result1 = RESULT(path=i,fmt="abacus")
                energy1 = result1["energy"]
                converge1 = result1["converge"]
                force1 = result1["force"]

                if not os.path.isdir(os.path.join(job,f"fdf_{ilabel}_{iidx}_{xyz}_-_{step_num}")):
                    print(f"ERROR: {job} does not have fdf_{ilabel}_{iidx}_{xyz}_-_{step_num}")
                    continue
                result2 = RESULT(path=os.path.join(job,f"fdf_{ilabel}_{iidx}_{xyz}_-_{step_num}"),fmt="abacus")
                energy2 = result2["energy"]
                converge2 = result2["converge"]
                force2 = result2["force"]

                if energy1 != None and energy2 != None:
                    fd = (energy2 - energy1) / (2 * step * constant.BOHR2A * step_num)
                else:
                    fd = None
                abacus = None if force ==None else force[atom_idx*3 + xyz_idx]
                allresults[f"{ilabel}_{iidx}_{xyz}_{step_num}"] = {"abacus": abacus,
                                                        "fd": fd,
                                                        "distance(A)": step * constant.BOHR2A * step_num,
                                                        "energy": energy,
                                                        "energy+": energy1,
                                                        "energy-": energy2,
                                                        "force+": None if force1 ==None else force1[atom_idx*3 + xyz_idx],
                                                        "force-": None if force2 == None else force2[atom_idx*3 + xyz_idx],
                                                        "converge+": converge1,
                                                        "converge-": converge2,}

            except:
                traceback.print_exc()
                continue
        return allresults

    def plot_force(self,allresults):
        # plot the force calculated by abacus and finite difference
        # plot two figures for each atom direction, one is the energy of each step, the other is the analytic and finite difference force of each step
        results = {}
        labels = []
        for ik,iv in allresults.items():
            if iv == None:
                continue
            for jk,jv in iv.items():
                if jk in ["force","converge","labels","step"]:
                    continue
                ilabel = ik + "/" + "_".join(jk.split("_")[:3])
                if ilabel not in results:
                    labels.append(ilabel)
                    results[ilabel] = {
                        "distance": [0],
                        "energy":[jv["energy"]],
                        "force":[jv["abacus"]],
                        "force0-Ana": jv["abacus"],
                        "distance0": [],
                        "force0-FD": []
                    }
                results[ilabel]["distance"].append(jv["distance(A)"])
                results[ilabel]["energy"].append(jv["energy+"])
                results[ilabel]["force"].append(jv["force+"])
                results[ilabel]["distance"].append(-1*jv["distance(A)"])
                results[ilabel]["energy"].append(jv["energy-"])
                results[ilabel]["force"].append(jv["force-"])
                results[ilabel]["distance0"].append(2 * jv["distance(A)"])
                results[ilabel]["force0-FD"].append((jv["energy-"] - jv["energy+"]) / (2 * jv["distance(A)"]))

        labels.sort()
        allpngs = []
        for ilabel in labels:
            energy = results[ilabel]["energy"]
            distance = results[ilabel]["distance"]
            force = results[ilabel]["force"]
            if None in energy or None in force:
                print(f"ERROR: {ilabel} has None in energy or force")
                continue
            # sort the data by distance
            results[ilabel]["energy"] = energy = [x for _,x in sorted(zip(distance,energy))]
            results[ilabel]["force"] = force = [x for _,x in sorted(zip(distance,force))]
            distance.sort()
            fd_force = []
            for i in range(1,len(distance)-1):
                fd_force.append((energy[i+1] - energy[i-1]) / (distance[i-1] - distance[i+1]))
            results[ilabel]["fd_force"] = fd_force

            # calculate the RMSD of the force
            rmsd = 0
            for i in range(len(force)-2):
                rmsd += (force[i+1] - fd_force[i])**2
            rmsd = (rmsd / len(force))**0.5

            # plot the two figures, one in left, the other in right
            import matplotlib.pyplot as plt
            fig,ax = plt.subplots(1,2,figsize=(12,5))
            ax1 = ax[0]
            #fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()
            ax1.plot(distance,energy,"b+-",label="Energy(eV)")
            ax2.plot(distance,force,"ro-",label="Analytic(eV/A)")
            ax2.plot(distance[1:-1],fd_force,"m*-",label="Finite Difference(eV/A)")
            ymax1 = max(energy)
            ymin1 = min(energy)
            ymax2 = max(max(force),max(fd_force))
            ymin2 = min(min(force),min(fd_force))
            ax1.set_xlabel("distance(A)")
            ax1.set_ylabel("Energy(eV)",color="b")
            ax2.set_ylabel("Force(eV/A)",color="r")
            ax1.spines['left'].set_color('blue')
            ax1.tick_params(axis='y', colors='blue')
            ax2.spines['right'].set_color('red')
            ax2.tick_params(axis='y', colors='red')
            ax1.set_ylim(ymin1-(ymax1-ymin1)*0.1,ymax1+(ymax1-ymin1)*0.2)
            ax2.set_ylim(ymin2-(ymax2-ymin2)*0.1,ymax2+(ymax2-ymin2)*0.2)
            ax1.set_title(ilabel + f"(FD-Ana RMSD={rmsd:.2e})")
            ax1.legend(loc="upper left")
            ax2.legend(loc="upper right")
            ax2.grid(True)

            ax3 = ax[1]
            # plot force0-FD vs distance0
            x = results[ilabel]["distance0"]
            y = results[ilabel]["force0-FD"]
            y_ref = results[ilabel]["force0-Ana"]
            ymin = min(y_ref, min(y))
            ymax = max(y_ref, max(y))
            # sort x,y by x
            x,y = zip(*sorted(zip(x,y)))
            ax3.plot(x,y,"m*-",label="Finite Difference(eV/A)")
            ax3.axhline(y_ref,ls="--",color="red",label="Analytic(eV/A)")
            ax3.set_xlabel("Step size of finite difference(A)")
            ax3.set_ylabel("Force(eV/A)")
            ax3.set_title(ilabel + f"(FD VS step-size)")
            ax3.legend(loc="upper right")
            ax3.grid(True)
            ax3.set_ylim(ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.2)


            plt.subplots_adjust(right=0.85) 
            png = ilabel + ".png"
            plt.tight_layout()
            plt.savefig(png)
            plt.close()
            allpngs.append(png)

        json.dump(results,open("results.json","w"),indent=4)    
        return allpngs