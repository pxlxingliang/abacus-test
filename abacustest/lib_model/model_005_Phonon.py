import os,sys,glob,json
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,WriteInput,WriteKpt
from ..model import Model
from . import comm


class Phonon(Model):
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
        return "phonon"
    
    @staticmethod
    def description():
        '''
        Description of the model
        '''
        return "Prepare and postprocess the phonon calculation"
    
    @staticmethod
    def add_args(parser):
        '''
        Add arguments for the model
        The arguments can not be command, model, modelcommand 
        '''
    
    def run(self,params):
        '''
        Parse the parameters and run the model
        '''
        return

    @staticmethod
    def prepare_args(parser):
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand '''
        parser.description = "This script is used to preapre the calculation of phonon calculation by abacus and phonopy."
        parser.add_argument('--abacus_command', type=str,  default="OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log",help='the command to run abacus job')
        parser.add_argument('--abacus_image', type=str,  default="registry.dp.tech/deepmodeling/abacus-intel:latest",help='the image to run abacus job')
        parser.add_argument('--abacus_machine', type=str,  default="c32_m64_cpu",help='the machine to run abacus job')
        parser.add_argument('-p', '--phonopy', type=str,  default="phonopy",help='the path of phonopy executable, default is phonopy') 
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is the current path.')
        parser.add_argument('-s', '--setting',default="setting.conf",help='the setting file for phonopy, default is setting.conf')
        
    
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        abacus_command = params.abacus_command
        abacus_image = params.abacus_image
        abacus_machine = params.abacus_machine
        phonopy = params.phonopy
        jobs = params.job if len(params.job) == 1 else params.job[1:]
        print(params)
        
        setting = PreparePhono(jobs,abacus_command,abacus_image,abacus_machine,phonopy,params.setting)
        if setting:
            comm.doc_after_prepare("Phonon", setting["run_dft"]["example"], ["setting.json"],has_prepare=False)  
    
    @staticmethod
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand'''
        parser.description = "This script is used to post-process the calculation of phonon calculation by abacus and phonopy."
        parser.add_argument('-p', '--phonopy', type=str,  default="phonopy",help='the path of phonopy executable, default is phonopy') 
        parser.add_argument('-j', '--job',default=["."], action="extend",nargs="*" ,help='the path of abacus inputs, default is the current path.')
        parser.add_argument('-s', '--setting',default="setting.conf",help='the setting file for phonopy, default is setting.conf')

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        phonopy = params.phonopy
        jobs = params.job if len(params.job) == 1 else params.job[1:]
        print(params)
        PostprocessPhonon(jobs,phonopy,params.setting)
    
def PreparePhono(jobs,abacus_command,abacus_image,abacus_machine,phonopy,settingf="setting.conf") :
    allexamples = []
    all_fail = True
    cwd = os.getcwd()
    for ijob in jobs:
        os.chdir(ijob)
        allfiles = os.listdir()
        if not os.path.isfile("INPUT"):
            print("Error: INPUT not exist in",ijob)
            continue
        input_param = ReadInput("INPUT")
        input_param["calculation"] = "scf"  # force set to scf
        input_param["suffix"] = "ABACUS"    # force set to ABACUS
        input_param["cal_force"] = "1"
        if not os.path.isfile("setting.conf"):
            print("Error: setting.conf not exist in",ijob)
            continue
        return_code, out, err = comm.run_command(f"{phonopy} {settingf} --abacus -d",shell=True)
        if return_code != 0:
            print(out,"\n",err)
            print("Error: run phonopy prepare failed in",ijob)
            continue
        allstrus = glob.glob("STRU-*")
        for istru in allstrus:
            idx = istru.split("-")[-1]
            os.makedirs("ABACUS-"+idx)
            os.system(f"cp -r {' '.join(allfiles)} ABACUS-{idx}")
            os.system(f"cp -r {istru} ABACUS-{idx}/STRU ")
            WriteInput(input_param,os.path.join("ABACUS-"+idx,"INPUT"))
            if os.path.isfile("phonopy_disp.yaml"):
                os.system(f"cp phonopy_disp.yaml ABACUS-{idx}/")
            allexamples.append(os.path.join(ijob,f"ABACUS-{idx}"))
        all_fail = False
        os.chdir(cwd)
    
    if all_fail:
        print("Error: all jobs are failed.")
        return {}
    setting = {
        "save_path": ".",
        "bohrium_group_name": "Phonopy",
        "run_dft": {
            "example": allexamples,
            "command": abacus_command,
            "image": abacus_image,
            "bohrium": {
                "scass_type": abacus_machine,
                "job_type": "container",
                "platform": "ali",
            },
        },
    }
    comm.dump_setting(setting)  
    return setting

def PostprocessPhonon(jobs,phonon,settingf="setting.conf"):
    if len(jobs) == 0:
        print("No jobs are specified.")
        return
    
    cwd = os.getcwd()
    supermetrics = {}
    for ijob in jobs:
        if not os.path.isdir(ijob):
            continue
        os.chdir(ijob)
        if not os.path.isfile("phonopy_disp.yaml"):
            print("Error: phonopy_disp.yaml not exist in",ijob)
            continue
        
        return_code, out, err = comm.run_command(f"{phonon} -f */OUT.*/running_*.log",shell=True)
        if return_code != 0:
            print(out,"\n",err)
            print("Error: catch force failed in",ijob)
            continue
        
        return_code, out, err = comm.run_command(f"{phonon} -p {settingf} --abacus -s",shell=True)
        if return_code != 0:
            print(out,"\n",err)
            print("Error: run phonopy postprocess failed in",ijob)
            continue
        
        pngfile ="phonon_band_structure.png"
        if plot_phonon(pngfile):
            print("Error: plot phonon failed in",ijob)
            continue
        else:
            supermetrics[f"{ijob}-phonon"] = {"type":"image","file":os.path.join(ijob,pngfile)}
        os.chdir(cwd)
        
    if len(supermetrics) > 0:
        json.dump(supermetrics, open("supermetrics.json", "w"), indent=4)
    

def plot_phonon(savefile="phonon_band_structure.png"):
    import yaml
    import matplotlib.pyplot as plt
    import numpy as np

    # Load the band.yaml file
    if not os.path.isfile("band.yaml"):
        print("Error: band.yaml not exist.")
        return 1
    
    with open("band.yaml", 'r') as stream:
        data = yaml.safe_load(stream)

    npath = data['npath']
    segment_nqpoint = data['segment_nqpoint']

    # Get the Q-point path and the frequencies at each Q-point
    q_points = [q['q-position'] for q in data['phonon']]
    distances = [q['distance'] for q in data['phonon']]
    frequencies = [[band['frequency'] for band in q['band']] for q in data['phonon']]

    # Assume `distances` is a list with distances from Gamma point in the reciprocal space
    # so that we have a continuous X-axis.
    #distance_points = np.cumsum([0] + [np.linalg.norm(np.subtract(q_points[i], q_points[i + 1])) for i in range(len(q_points) - 1)])
    frequency_array = np.array(frequencies)

    # Plot the phonon band structure
    for band in range(frequency_array.shape[1]):
        plt.plot(distances, frequency_array[:, band], color='blue')

    # Adding labels and title
    plt.xlabel('Distance along q-point path')
    plt.ylabel('Frequency (THz)')
    plt.title('Phonon Band Structure')

    # x is the distance between the q-points
    # y is the frequency of the phonon mode
    # plot the figure
    plt.xlim(0, distances[-1])
    # find the ylim by using the min and max of the frequencies
    plt.ylim(np.min(frequency_array)-2, np.max(frequency_array)+2)
    plt.axhline(y=0, color='black', linewidth=0.5)
    plt.axvline(x=0, color='black', linewidth=0.5)
    plt.axvline(x=distances[-1], color='black', linewidth=0.5)
    for i in range(npath):
        num = sum(segment_nqpoint[:i])
        x = distances[num]
        plt.axvline(x, color='black', linewidth=0.5)
    
    plt.savefig(savefile)
    plt.close()