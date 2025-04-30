import os,sys,glob,json
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,WriteInput,WriteKpt
from ..model import Model
from . import comm
import numpy as np
import matplotlib.pyplot as plt
import traceback
from abacustest.constant import RECOMMAND_IMAGE

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
        parser.add_argument('--abacus_image', type=str,  default=RECOMMAND_IMAGE,help='the image to run abacus job')
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
        parser.add_argument('--onlyplot', action='store_true', help='only plot the phonon band structure, do not run the phonopy postprocess')

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        phonopy = params.phonopy
        jobs = params.job if len(params.job) == 1 else params.job[1:]
        print(params)
        PostprocessPhonon(jobs,phonopy,params.setting,params.onlyplot)
    
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

def PostprocessPhonon(jobs,phonon,settingf="setting.conf",only_plot=False):
    if len(jobs) == 0:
        print("No jobs are specified.")
        return
    
    cwd = os.getcwd()
    supermetrics = {}
    for ijob in jobs:
        if not os.path.isdir(ijob):
            continue
        os.chdir(ijob)
        
        if not only_plot:
            if not os.path.isfile("phonopy_disp.yaml"):
                print("Error: phonopy_disp.yaml not exist in",ijob)
                continue
            
            os.system(f"{phonon} -f */OUT.*/running_*.log")
            #return_code, out, err = comm.run_command(f"{phonon} -f */OUT.*/running_*.log",shell=True)
            #if return_code != 0:
            #    print(out,"\n",err)
            #    print("Error: catch force failed in",ijob)
            #    continue

            os.system(f"{phonon} -p {settingf} --abacus -s")
            #return_code, out, err = comm.run_command(f"{phonon} -p {settingf} --abacus -s",shell=True)
            #if return_code != 0:
            #    print(out,"\n",err)
            #    print("Error: run phonopy postprocess failed in",ijob)
            #    continue
        
        pngfile ="phonon_band_structure.png"
        if plot_phonon("band.yaml",save_file=pngfile) is None:
            print("Error: plot phonon failed in",ijob)
            continue
        else:
            supermetrics[f"{ijob}-phonon"] = {"type":"image","file":os.path.join(ijob,pngfile)}
            print(f"Plot phonon band structure in {ijob} saved in {os.path.join(ijob,pngfile)} successfully!!!")
        os.chdir(cwd)
        
    if len(supermetrics) > 0:
        json.dump(supermetrics, open("supermetrics.json", "w"), indent=4)
    
def rearrange_labels(distances, frequencies, labels, new_labels):
    new_distances = []
    new_frequencies = []
    end_distance = 0
    for label in new_labels:
        if label in labels:
            label_idx = labels.index(label)
            reverse = False
        elif label[::-1] in labels:
            label_idx = labels.index(label[::-1])
            reverse = True
        else:
            
            raise ValueError(f"Label {label} not found in labels.")
            
        idis = distances[label_idx]
        ifreq = frequencies[label_idx]
        if reverse:
            idis = [-1*i for i in idis[::-1] ]
            ifreq = ifreq[::-1]
        new_distances.append([i + end_distance - idis[0] for i in idis])
        new_frequencies.append(ifreq)
        end_distance += idis[-1] - idis[0]
    return new_distances, new_frequencies, new_labels

def read_band_yaml(filename):
    import yaml
    with open(filename, 'r') as file:
        data = yaml.safe_load(file)

    labels = []
    for path in data['labels']:
        labels.append((path[0], path[1]))
    
    segments = data['segment_nqpoint']
    qpoints = data['phonon']
    
    distances = []
    frequencies = []
    for q in qpoints:
        distances.append(q['distance'])
        freqs = [band['frequency'] for band in q['band']]
        frequencies.append(freqs)
    
    frequencies = np.array(frequencies) # shape (nq, nband)
    
    split_indices = np.cumsum(segments)[:-1]
    
    # split distances and frequencies
    new_distances = []
    new_frequencies = []
    begin_idx = 0
    for end_idx in split_indices:
        new_distances.append(distances[begin_idx:end_idx])
        new_frequencies.append(frequencies[begin_idx:end_idx])
        begin_idx = end_idx
    new_distances.append(distances[begin_idx:])
    new_frequencies.append(frequencies[begin_idx:])

    return new_distances, new_frequencies, labels

def plot_phonon(yaml_file="band.yaml", 
                type_label = None, 
                high_symmetry_path = None,
                save_file="band.png"):
    '''Plot the phonon band structure from yaml file
    
    Args:
        yaml_file (str|list): The yaml files to plot
        type_label (list): The label of the band structure
        high_symmetry_path (list of tuples): The high symmetry path, used to rearrange the labels
        savefile (str): The file name to save the plot
    
    Returns:
        datas (list): The data of the band structure, each element is a tuple of (distances, frequencies, q_labels)
        labels (list): The labels of the band structure
        
    If you want to plot multiple yaml files, type_label should be a list of the same length as yaml_file.
    '''

    if isinstance(yaml_file, str):
        yaml_file = [ yaml_file ]

    if type_label is not None and isinstance(type_label, str):
        type_label = [type_label]
    elif type_label is None and len(yaml_file) > 1:
        type_label = [i for i in yaml_file]
    
    if type_label is not None and len(type_label) != len(yaml_file):
        print("Error: The length of label and yaml_file should be the same.")
        return None

    # check files
    datas = []
    labels = []
    for i, file in enumerate(yaml_file):
        if not os.path.isfile(file):
            print(f"Error: {file} not exist.")
            continue
        try:
            distances, frequencies, q_labels = read_band_yaml(file)
            if high_symmetry_path is not None:
                distances, frequencies, q_labels = rearrange_labels(distances, frequencies, q_labels, high_symmetry_path)
        except:
            traceback.print_exc()
            print(f"Error: read {file} failed.")
            continue
        
        datas.append((distances, frequencies, q_labels))
        if type_label is not None:
            labels.append(type_label[i])
            
    if len(datas) == 0:
        print("Error: No valid yaml file.")
        return None 

    colors = [
        "#d62728",  # 红色
        "#1f77b4",  # 蓝色
        "#2ca02c",  # 绿色
        "#ff7f0e",  # 橙色
        "#9467bd",  # 紫色
        "#8c564b",  # 棕色
        "#e377c2",  # 粉色
        "#7f7f7f",  # 灰色
        "#bcbd22",  # 黄绿色
        "#17becf"   # 青色
    ]
    if len(datas) > len(colors):
        colors = plt.cm.viridis(np.linspace(0, 1, len(datas)))
    plt.figure(figsize=(8, 6))
    for i, (distances, frequencies, q_labels) in enumerate(datas):
        for j, distance in enumerate(distances):
            freq = np.array(frequencies[j])
            for k in range(freq.shape[1]):
                if j == 0 and k==0 and type_label is not None:
                    label = labels[i]
                else:
                    label = None
                plt.plot(distance, freq[:, k], "-", lw=1, label=label, color=colors[i])

    for x in distances[1:]:
        plt.axvline(x[0], color='black', linestyle='--', lw=0.8)

    print(q_labels)
    tick_pos = [distances[0][0]]
    tick_labels = [q_labels[0][0]]
    for idx, distance in enumerate(distances):
        if idx == 0:continue
        tick_pos.append(distance[0])
        if q_labels[idx-1][-1] == q_labels[idx][0]:
            label = q_labels[idx][0]
        else:
            label = f"{q_labels[idx-1][-1]}|{q_labels[idx][0]}"
        tick_labels.append(label)
    
    tick_pos.append(distances[-1][-1])
    tick_labels.append(q_labels[-1][1])

    print(tick_pos)
    print(tick_labels)
    plt.xticks(tick_pos, tick_labels, fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0, distances[-1][-1])
    plt.axhline(0, color='gray', linestyle='--', lw=0.8, alpha=0.5)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    y_max = np.max(frequencies)
    y_min = np.min(frequencies) - 0.1 * y_max
    if type_label is not None:
        # set the maximum Y to 1.3 * y_max
        plt.ylim(y_min, y_max * 1.3)
        if len(labels) <= 3:
            ncol = len(labels)
        else:
            ncol = int(len(labels)**0.5) + 1
        plt.legend(fontsize=12, frameon=False, loc='upper center', ncol=ncol)
    
    plt.xlabel("Q-Point", fontsize=14)
    plt.ylabel("Frequency (THz)", fontsize=14)
    
    plt.tight_layout()
    plt.savefig(save_file,dpi=300)
    plt.close()
    
    return (datas, labels)
    
    
    