from ..model import Model
import os, glob, json, traceback, copy
from . import comm
import numpy as np
from abacustest.lib_prepare.abacus import ReadInput, WriteInput, ReadKpt, WriteKpt, AbacusStru
from abacustest.constant import RY2EV
from abacustest.lib_collectdata.comm import cal_band_gap

class BandModel(Model):
    @staticmethod
    def model_name(): # type: ignore
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "band"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "Prepare and postprocess the calcualtion of band structure"
    
    @staticmethod
    def prepare_args(parser):
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand '''
        
        parser.description = "Prepare the band structure calculation. Will generate the scf/nscf/kpath input files."
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is current folder')
        parser.add_argument("-c", "--rundftcommand", type=str, default="OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log",help="the command to execute aabcus, default is 'OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log' ")
        parser.add_argument("-i","--image",default="registry.dp.tech/deepmodeling/abacus-intel:latest",type=str,help="the used image. Default is: registry.dp.tech/deepmodeling/abacus-intel:latest", )
        parser.add_argument("--machine", default="c32_m64_cpu", help="the machine to run the abacus. Default is c32_m64_cpu")
        parser.add_argument("-r", "--run", default=0, help="if run the test. Default is 0.", type=int)
        return parser
    
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        real_jobs = comm.get_job_list(params.job)
        real_jobs, run_script = PrepBand(real_jobs,params.rundftcommand).run()
        if len(real_jobs) == 0:
            print("No valid job is found")
            return
        
        setting = {
            "save_path": "results",
            "bohrium_group_name": "band-strcutre",
            "run_dft": {
                "example": real_jobs,
                "command": run_script,
                "image": params.image,
                "bohrium": {
                    "scass_type": params.machine,
                    "job_type": "container",
                    "platform": "ali",
                },
            },
        }
        
        comm.dump_setting(setting)
        comm.doc_after_prepare("band", real_jobs, ["setting.json"],has_prepare=False)
        print("After finish the calculation, you can run below command to do the postprocess:")
        print(f"    abacustest model {self.model_name()} post -j {' '.join(real_jobs)} -r result.json\n")
        
        if params.run:
            bash_script = "cal_band_structure.sh"
            with open(bash_script,"w") as f:
                f.write("abacustest submit -p setting.json\n")
            os.system(f"bash {bash_script} &")

    @staticmethod
    def postprocess_args(parser):
        '''
        PostProcess the band structure calculation. Will plot the band structure and calculate the band gap.
        '''
        parser.description = "PostProcess the band structure calculation. Will plot the band structure and calculate the band gap."
        parser.add_argument('-j', '--job',default=[], action="extend",nargs="*" ,help='the path of abacus inputs, default is the current path. There should be some eos folders in each example path.')
        parser.add_argument("--input", default="INPUT.nscf", type=str, help="the input file name, default is INPUT.nscf")
        parser.add_argument("--kpt", default="KPT.nscf", type=str, help="the kpoint file name, default is KPT.nscf")
        parser.add_argument("--range", default=[-5,5], type=float, help="The range of the band structure, default is efermi-5, efermi+5 eV", nargs=2)
        return parser
        
        
    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        jobs = comm.get_job_list(params.job)
        PostBand(jobs, input_file=params.input, kpt_file=params.kpt, band_range=params.range).run()
        print("The band structure is plotted in the band.png file")
                
class PrepBand:
    def __init__(self, jobs,run_command=None,job_type="abacus"):
        self.jobs = jobs
        self.job_type = job_type
        self.run_command = self.modify_command(run_command)
        self.scf_input = "INPUT.scf"
        self.scf_kpt = "KPT.scf"
        self.nscf_input = "INPUT.nscf"
        self.nscf_kpt = "KPT.nscf"
        self.run_script = "run.sh"
    
    def run(self):
        if self.job_type == "abacus":
            jobs = self.run_abacus()
        else:
            print(f"Error: the job type {self.job_type} is not supported")
            jobs = []
        
        return jobs, f"bash {self.run_script}"
                
    def modify_command(self,command):
        # if the command is end with | tee out.log, remove it
        if command is None:
            return "abacus"
        commands = command.split()
        if "|" in commands and "tee" in commands and commands.index("|") == commands.index("tee") - 1:
            return " ".join(commands[:commands.index("|")])  
        else:
            return command        
                   
    def run_abacus(self):
        valid_jobs = []
        for job in self.jobs:
            if not os.path.isdir(job):
                print(f"Warning: {job} is not a valid folder")
                continue
            inputfile = os.path.join(job,"INPUT")
            strufile = os.path.join(job,"STRU")
            if not os.path.isfile(inputfile):
                print(f"Warning: {job} is not a valid folder, can not find {inputfile}")
                continue
            if not os.path.isfile(strufile):
                print(f"Warning: {job} is not a valid folder, can not find {strufile}")
                continue

            try:
                kpt = ReadKpt(job)
            except:
                traceback.print_exc()
                print(f"Warning: {job} is not a valid folder, can not get the K points information from KPT file or INPUT file")
                continue

            input_param = ReadInput(inputfile)
            scf_input = copy.deepcopy(input_param)
            nscf_input = copy.deepcopy(input_param)
            
            stru = AbacusStru.ReadStru(strufile)
            if stru is None:
                print(f"Warning: {job} is not a valid folder, can not get the structure information from STRU file")
                continue

            # generate the scf input
            scf_input["calculation"] = "scf"
            scf_input["out_chg"] = 1
            scf_input["kpoint_file"] = self.scf_kpt
            scf_input.pop("kspacing",None)
            WriteInput(scf_input,os.path.join(job,self.scf_input))
            WriteKpt(kpt[0], os.path.join(job,self.scf_kpt),model=kpt[1])
            
            # generate the nscf input
            nscf_input["calculation"] = "nscf"
            nscf_input["init_chg"] = "file"
            nscf_input["out_band"] = 1
            nscf_input["kpoint_file"] = self.nscf_kpt
            nscf_input.pop("kspacing",None)
            WriteInput(nscf_input,os.path.join(job,self.nscf_input))
            stru.get_kline(kpt_file=os.path.join(job,self.nscf_kpt),point_number=20)
            
            # gen run script
            with open(os.path.join(job,self.run_script),"w") as f:
                f.write(f"cp {self.scf_input} INPUT\n")
                f.write(f"{self.run_command} | tee scf.log\n")
                f.write(f"cp {self.nscf_input} INPUT\n")
                f.write(f"{self.run_command} | tee nscf.log\n")
            valid_jobs.append(job)
        return valid_jobs
        

class PostBand:
    def __init__ (self, jobs, input_file=None, kpt_file=None, band_range=None):
        self.jobs = jobs
        self.input_file = input_file 
        self.kpt_file = kpt_file  
        self.band_range = [-5, 5 ] if band_range is None else band_range
                    
    def run(self):
        for job in self.jobs:
            print(job)
            inputfile = self.find_input(job)
            kpt = self.get_kpoints(inputfile, job)
            suffix = "ABACUS"
            if inputfile is not None:
                input_param = ReadInput(inputfile)
                suffix = input_param.get("suffix","ABACUS")
            band_file = os.path.join(job, "OUT."+input_param.get("suffix","ABACUS"), "BANDS_1.dat")
            if not os.path.isfile(band_file):
                print(f"Warning: can not find the band file {band_file}")
                continue
            bands = self.get_band(band_file)
            self.plot_band(bands, kpt, os.path.join(job,"band.png"), efermi=self.get_efermi(os.path.join(job,"OUT."+suffix,"running_nscf.log")))
    
    def get_efermi(self, logfile):
        with open(logfile) as f:
            lines = f.readlines()
            
        for line in lines[::-1]:
            if "E_Fermi" in line:
                return float(line.split()[-1])
            elif "read in fermi energy" in line:
                return float(line.split()[-1]) * RY2EV
        return None
                
    @staticmethod
    def get_band(band_file):
        bands = np.loadtxt(band_file)
        return bands[:,2:].T # the first two columns are index
                 
    def find_input(self, job):
        if self.input_file is not None and os.path.isfile(os.path.join(job,self.input_file)):
            return os.path.join(job,self.input_file)
        else:
            for ifile in ["INPUT.nscf", "INPUT"]:
                if os.path.isfile(os.path.join(job,ifile)):
                    return os.path.join(job,ifile)
        return None
    
    def get_kpoints(self, inputfile, job):
        if self.kpt_file is not None and os.path.isfile(os.path.join(job,self.kpt_file)):
            return ReadKpt(os.path.join(job,self.kpt_file))

        if inputfile is not None:
            input_param = ReadInput(inputfile)
            kpt_file = input_param.get("kpoint_file","KPT")
            if os.path.isfile(os.path.join(job,kpt_file)):
                return ReadKpt(os.path.join(job,kpt_file))
        
        for ifile in ["KPT.nscf", "KPT", "KLINES","LINES"]:
            if os.path.isfile(os.path.join(job,ifile)):
                return ReadKpt(os.path.join(job,ifile))
            
        return None

    @staticmethod
    def rearrange_plotdata(band, line_points):
        # line_points is a list of the high symmetry k points and contarns the k points and symbols, 
        # which is generated by the ReadKpt.get_kline, line [0.0, 0.0, 0.0, 20, '#G']
        npoint = []
        symbols = []
        for line in line_points:
            npoint.append(line[3])
            symbols.append(line[4].lstrip("#"))
        
        assert sum(npoint) == len(band[0]), "The number of k points is not equal to the band data"
        
        # if one point has only 1 k point, then we should merge it with the next point
        x = sum(npoint) - npoint.count(1)
        if npoint[-1] == 1: x += 1
        
        symbol_index = [0]
        symbols_new = [symbols[0]]
        for i in range(1,len(npoint)):
            if npoint[i-1] == 1:
                continue
            symbol_index.append(symbol_index[-1] + npoint[i-1])
            if npoint[i] == 1 and i != len(npoint) - 1:
                symbols_new.append(symbols[i] + "|" + symbols[i+1])
            else:
                symbols_new.append(symbols[i])
                

        band_idx = [] # we need split the band data into different segments
        x_start = 0  # the start index in x axis
        data_start = 0 # the start index in band data
        points = [npoint[0]]
        for i in range(1, len(npoint)):
            if npoint[i-1] == 1 or i == len(npoint) - 1:
                band_idx.append([x_start, x_start + sum(points), data_start, data_start + sum(points)])
                x_start += sum(points) - 1  # need merge the two points into one
                data_start += sum(points)
                points = [npoint[i]]
            else:
                points.append(npoint[i])
        band_idx[-1][1] += 1
        band_idx[-1][3] += 1
        return band_idx, symbol_index, symbols_new
    
    def plot_band(self, band, kpt, filename, efermi = None):
        if kpt is None or kpt[1] != "line":
            band_idx = [[0, len(band[0]), 0, len(band[0])]]
            symbol_index = None
            symbols = None
        else:
            band_idx, symbol_index, symbols = self.rearrange_plotdata(band, kpt[0])
            print(band_idx, symbol_index, symbols)
        import matplotlib.pyplot as plt
        fontsize = 12
        fig, ax = plt.subplots(1,2,figsize=(8,4))
        for i, idx in enumerate(band_idx):
            for iband in band:
                ax[0].plot(range(idx[0], idx[1]), iband[idx[2]:idx[3]], "-", color="black",linewidth=0.8)
                ax[1].plot(range(idx[0], idx[1]), iband[idx[2]:idx[3]], "-", color="black",linewidth=0.8) 
        ax[0].set_xlim(0, band_idx[-1][1]-1)
        ax[1].set_xlim(0, band_idx[-1][1]-1)
        if symbols is not None:
            ax[0].set_xticks(symbol_index)
            ax[0].set_xticklabels(symbols, fontsize=fontsize)
            ax[1].set_xticks(symbol_index)
            ax[1].set_xticklabels(symbols, fontsize=fontsize)
        ax[0].set_xlabel("K points", fontsize=fontsize)
        ax[1].set_xlabel("K points", fontsize=fontsize)
        ax[0].set_ylabel("Energy (eV)", fontsize=fontsize)
        ax[1].set_ylabel("Energy (eV)", fontsize=fontsize)

        # plot the x = i for each symbol
        for i in range(1, len(symbol_index)-1):
            ax[0].axvline(x=symbol_index[i], color="black", linestyle="-", linewidth=0.3)
            ax[1].axvline(x=symbol_index[i], color="black", linestyle="-", linewidth=0.3)

        if efermi is not None:
            ax[0].axhline(y=efermi, color="red", linestyle="--", linewidth=0.8)
            ax[1].axhline(y=efermi, color="red", linestyle="--", linewidth=0.8)
            bg = cal_band_gap([np.array(band).T], efermi)
            #ax[0].text(band_idx[-1][1], efermi+0.3, f"E_f={efermi:.2f}eV (BG={bg:.2f}eV)", color="red", fontsize=fontsize-2)
            #ax[1].text(band_idx[-1][1], efermi+0.3, f"E_f={efermi:.2f}eV (BG={bg:.2f}eV)", color="red", fontsize=fontsize-2)
            ax[1].set_ylim(efermi+self.band_range[0],efermi+self.band_range[1])
            title = f"Band Structure (E_fermi={efermi:.2f}eV, BandGap={bg:.2f}eV)"
            print(f"E_fermi={efermi:.2f}eV, BandGap={bg:.2f}eV")
        else:  
            ax[1].set_ylim(self.band_range[0],self.band_range[1])
            title = "Band Structure"
        #ax[0].set_title(title, fontsize=fontsize)
        
        # set a total title
        plt.suptitle(title, fontsize=fontsize)
        
        plt.tight_layout()
        plt.savefig(filename,dpi=400)
        plt.close()
                
        
        
                        
        
        
        

    