# abacustest
Do the performance test of ABACUS \
Install:
`pip install .`
And then, you can use command `abacustest`. \
The sub_commands:
- `submit`
- `collectdata`
- `outresult`
- `status`
- `prepare`
- `download`

Please use `abacustest -h` to get the usages

## 1. submit
```
usage: abacustest submit [-h] [-p PARAM] [-s SAVE] [--override OVERRIDE] [--debug [DEBUG]] [--download DOWNLOAD]

This script is used to run a testing

optional arguments:
  -h, --help            show this help message and exit
  -p PARAM, --param PARAM
                        the job setting file, default is job.json
  -s SAVE, --save SAVE  the folder where the results will be put in, default: result/date_of_today (e.g. result/20230101)
  --override OVERRIDE   when the save folder exists, if override it. 0: no, 1: yes.
  --download DOWNLOAD   if wait the finish of flow and download the results, 0: no, 1: yes. Default is 1
  --debug [DEBUG]       debug mode for dflow
```
job.json (-p) is needed as input \
The installation of pydflow and its related related dependencies (https://github.com/deepmodeling/dflow#2--quick-start) are required for `submit`.

### 1.1 job.json
This file defines the detail of the jobs. \
An example is like:
```
{
    "bohrium_group_name": "abacustest",
    "config":{
        "bohrium_username":         "xxx",
        "bohrium_password":         "xxx",
        "bohrium_project_id":           111
    },
    "ABBREVIATION":{
            "ABACUS310_IMAGE": "registry.dp.tech/dptech/abacus:3.1.0",
            "PYTHON_IMAGE": "python:3.8"
    },

    "save_path":"result/abacus-pw",
    "max_parallel": 100,

    "run_dft":[
        {"ifrun": true,
         "sub_save_path": "",
         "image": "ABACUS310_IMAGE",
         "example":[["00[0-2]"],"00[3-5]"],
         "group_size" : 1,
         "bohrium": {"scass_type":"c8_m16_cpu","job_type":"container","platform":"ali"},
         "command": "mpirun -np 8 abacus > log",
         "extra_files":[],
         "outputs":[]
        },
        {"ifrun": true,
         "image": "ABACUS310_IMAGE",
         "example":[["00[6-7]"],"00[8-9]"],
         "group_size" : 1,
         "bohrium": {"scass_type":"c16_m32_cpu","job_type":"container","platform":"ali"},
         "command": "mpirun -np 8 abacus > log",
         "extra_files":[],
         "outputs":[]
        }
    ],

    "post_dft":{
                "ifrun": false,
                "command": "ls",
                "extra_files": [],
                "image":   "python:3.8",
                "metrics":{
                    "before_command":true,
			        "dft_type":"abacus",
			        "metrics_name": [],
			        "save_file": "metrics.json",
			        "newmethods": []
		        },
                "outputs": []
    }
}
```
- ``bohrium_group_name``: If use bohrium, you can define the group name at here.
- `ABBREVIATION`: define some abbreviation, and is only valid for `image`.
- `config`: The setting of config information. If you do not use bohrium, this key can be removed.
  - `bohrium_username`: the username to login Bohrium (the email or phone number that used to register Bohrium).
  - `bohrium_password`: the password of Bohrium.
  - `bohrium_project_id`: the project id of Bohrium. \
  Besides, you can set the bash environment variables (uppercase letter) to replace the value of `config`, such as: ```export BOHRIUM_USERNAME="xxx"```, and then you can remove the key `config` in job.json. If both are defined, the value in job.json will be used.
- `save_path`: define the path to save the results of this test. If this key is not defined or if is deined to be "" or None, the value will be replaced to be the path defined by "-s" (the param of abacustest, the default of "-s" is "result/date_of_today" like "result/20230101")
- `max_parallel`: define the max number of parallel jobs. Default is 100.
- `run_dft`: define the detail of the running of your jobs. The value is a list of dictionaries. You can set any number of dictionaries.
  - `ifrun`: if set it to be `false`, will skip this part. 
  - `sub_save_path`: the path to save the results of this part in `save_path`, which means the real save path will be "save_path/sub_save_path". If this key is not defined, or is defined to be "" ornull, the real save path will be "save_path".
  - `image`: define the image name. Also you can use the name defined in `ABBREVIATION`. Progrma will firstly check if the image name is defined in `ABBREVIATION`, if yes, the name will be replacedby the value in `ABBREVIATION`, and if no, the image name will be kept to be value of `image`.
  - `example`: The folder names of your jobs. Here assume that you have 5 jobs and the folder names are 000, 001, 002, ..., 004. You can write as ["000", "001", "002", "003", "004"], and also you can write as ["00[0-4]"] that can be recongnized by `glob.glob`. Besides, you can put some folders in a list, and they will be set as one group and run the jobs serially. Such as: [["00[0-1]"],"0[2-4]"], "000" and "001" is one group, each of "002", "003" and "004" is one group and total 4 groups.
  - `group_size`: define how many example groups to run on a single machine.
  - `bohrium`: if you want to submit the job to Bohrium, you need set this key, and if you do not want to use Bohrium, please remove this key or set it to be `false`.
  - `command`: the command to run the job. It is same for all jobs defined in `example`.
  - `extra_files`: if you need extra files (such as "collectdata-abacus.json"), you can defined them at here. Before run the command, these files will be copied to each example folder.
  - `metrics`: define the metrics information.
    - `before_command`: if generate the metrics before run the command. Default is True.
    - `dft_type`: the job type, can be one of: "abacus", "qe", "vasp", and also 0/1/2 is ok, which refers to "abacus"/"qe"/"vasp" respectively.
    - `metrics_name`: a list of your metric name. Can find all names by `abacustest collectdata --outparam -t 0`. If is a null list, will catch all registered metrics.
    - `save_file`: the file name to store the values of metrics, which is a json file type.
    - `newmethods`: if the self-defined methods is needed, please add the file name (generally is a .py file) in `extra_files` and the module name (generally is file name without .py) in here.Detail of newmethods please see [2.3 import self-defined methods](#2.3 import self-defined methods).
  - `outputs`: specify which files and folders will be downloaded. If set to be "[]", all of the files in the folder will be downloaded.\
  - `dispatcher`: if you need to use other platform (not Bohrium), you can delete key `bohrium` and add a key `dispatcher`, which is a dict, like:
    ```
    "dispatcher": {
                "machine_dict": {
                    "remote_root": "/home/username/abacustest_work",   # Here is a path where the username can access
                    "remote_profile": {
                        "hostname": "xxx.xx.xxx.xxx",
                        "username": "Username",         
                        "password": "password",
                        "port": 22
                    }
                },
                "resources_dict": {
                    "number_node": 1,
                    "cpu_per_node": 8,
                    "gpu_per_node": 1,
                    "queue_name": "Normal"
                }
            }
    ```
    If you do not need gpu, then you should delete key `gpu_per_node`.

- `post_dft`: define the detail of post processing, and now all examples will be put at one same place. The key are same as those in `run_dft`, but no need the definition of example. 


### 1.2 submit a test
```
abacustest submit -p job.json -s result/test
```


## 2. collectdata
```
usage: abacustest collectdata [-h] [-j [JOBS [JOBS ...]]] [-t {0,1,2}] [-p PARAM] [-o OUTPUT] [-m [MODULES [MODULES ...]]] [--newmethods [NEWMETHODS [NEWMETHODS ...]]] [--outparam [OUTPARAM]]

This script is used to collect some key values from the output of ABACUS/QE/VASP jobs

optional arguments:
  -h, --help            show this help message and exit
  -j [JOBS [JOBS ...]], --jobs [JOBS [JOBS ...]]
                        the path of jobs
  -t {0,1,2}, --type {0,1,2}
                        0:abacus, 1:qe, 2:vasp. Default: 0
  -p PARAM, --param PARAM
                        the parameter file, should be .json type
  -o OUTPUT, --output OUTPUT
                        the file name to store the output results, default is "metrics.json"
  -m [MODULES [MODULES ...]], --modules [MODULES [MODULES ...]]
                        add extra modules. Default only module 'job-type' will be loaded, such as: 'abacus' for abacus type. You can check all modules by --outparam
  --newmethods [NEWMETHODS [NEWMETHODS ...]]
                        the self-defined python modules, and shuold be format of import, such as "abc"(the file name is abc.py), "a.b.c" (the file is a/b/c.py)
  --outparam [OUTPARAM]
                        output the registed parameters, you can set the type by -t or --type to choose abacus/qe/vasp. 0: No, 1: yes
```
This function is used to get some key values of the ABACUS/QE/VASP jobs. \

You need one param.json file to define which key values you want to catch. \
An example is like:
```
{"PARAM":
        ["natom","kpt","ibzk","nelec","nbands","force","stress",
        "scf_steps","total_time","force_time","stress_time","energy_per_atom","band_gap",
        {"SELF_DEFINED_NAME": "An eval string to get the value"}]
}
```
Only one key "PARAM" is recongnized by `collectdata`, the value is a list of keys. \

If you wish to customize the key name, you can add a dictionary to the list. The key in this dictionary is the name you want to define, and the value is an eval string used to obtain the value that the key should hold in {}.\
For instance:
- To get the value of energy per atom, you can write: {"energy per atom": "{energy}/{natom}"}. Here, "energy" and "natom" are the keys, and the value of "energy per atom" is the result of "energy" divided by "natom".
- To rename the key name of energy to "energy (eV)", you can write: {"energy (eV)": "energy"}.



The results will be save as json type, and you can define the file name by `-o`.

### 2.1 Get the keys
You can set `--outparam` to print out all of the keys of one type. Such as: `collectdata --outparam 1 -t 0` will print out the keys of ABACUS:
```
Job type: abacus, all modules:  ['abacus']
             version:	abacus:Abacus.GetVersion()    	the version of ABACUS
               ncore:	abacus:Abacus.GetNcore()      	the mpi cores
          normal_end:	abacus:Abacus.GetNormalEnd()  	if the job is nromal ending
               INPUT:	abacus:Abacus.GetInputParameter()	a dict to store the setting in OUT.xxx/INPUT
                 kpt:	abacus:Abacus.GetKptParam()   	list, the K POINTS setting in KPT file
              nbands:	abacus:Abacus.GetLogParam()   	number of bands
            converge:	abacus:Abacus.GetLogParam()   	if the SCF is converged
...
```
The avial modules are listed in the first line, and you can add some of them if needed.

### 2.2 import collectdata.RESULT in your python script
You can also use this collectdata function in your script by `from abacustest.lib_collectdata.collectdata import RESULT`, \
`RESULT` is a function, the arguments is: 
```
def RESULT(fmt="abacus",outparam=False,**kwargs):
    #fmt: format type, abacus/qe/vasp is recongnized
    #outparam: if output all registered parameters
    #**kwargs: the other input commands for the init function of ResultAbacus/Qe/Vasp

#init function of ResultAbacus
class ResultAbacus(Result):
    def __init__(self,path = ".",output = None):
        #path： path of your ABACUS job
        #output: the filename of screen output of ABACUS job, 
        
#init function of ResultQe        
class ResultQe(Result):
    def __init__(self,path = ".",output = None):  
    #path： path of your qe job  
    #output: the filename of output of the qe job,  

# ResultVasp     
class ResultVasp(Result):
    def __init__(self,path = "."):    
    #path： path of your VASP job 
```    

Assume that you have an ABACUS job, whose path is `abacusjob`, and you can use the below codes to catch information of enery, force and stress:
```
from abacustest.lib_collectdata.collectdata import RESULT

abacusresult = RESULT(fmt="abacus",path=abacusjob)
energy = abacusresult["total_energy"] 
force = abacusresult["force"]
stress = abacusresult["stress"]
```

### 2.3 import self-defined methods
One can self define the methods by write a python script, and add it to collectdata system by below steps:
1. define a class and inherit the collect result class.
2. define a class function, and register it.
such as:
```
from abacustest.lib_collectdata.resultAbacus import ResultAbacus
class MyAbacus(ResultAbacus):
    @ResultAbacus.register(key_name="description of the key")
    def function_name(self):
        "your codes to get the value of key"
        self[key_name] = value
``` 
If you want to define a vasp method, you need change the import line and class name, like:
```
from abacustest.lib_collectdata.resultVasp import ResultVasp
class MyVasp(ResultVasp):
    @ResultVasp.register(key_name="description of the key")
    def function_name(self):
```
The class name `MyVasp` is self-defined, and only need to conform to the naming rule of class. \
The key_name is self-define string, such as: "natom", "scf_steps". \
Do not forget to save the value by `self[key_name] = value`. \

You can directly add the codes in you python script. \
Or you can save the codes in a file endwith ".py" (such as "mymethod.py"), and then add the parameters `--newmethods` in command, like:
```
abacustest collectdata --newmethods "mymethod"
``` 


## 3. outresult
```
abacustest outresult [-h] [-r [RESULT [RESULT ...]]]
```
Print out the results from collectdata.


## 4. status
```
usage: abacustest status [-h] paramf job_id
```
Not ready now.

## 5. prepare
```
usage: abacustest prepare [-h] -p PARAM [-s SAVE]

This script is used to prepare the INPUTS OF ABACUS JOB

optional arguments:
  -h, --help            show this help message and exit
  -p PARAM, --param PARAM
                        the parameter file, should be .json type
  -s SAVE, --save SAVE  where to store the inputs, default is abacustest
```
A template of param.json is:
```
{
  "prepare":{
      "example_template":["example1","example2"],
      "input_template":"INPUT",
      "kpt_template":"KPT",
      "strus": ["1.vasp","2.vasp"],
      "stru_format": "vasp",
      "mix_input":{
          "ecutwfc":[50,60,70],
          "kspacing":[0.1,0.12,0.13]
      },
      "mix_kpt":[],
      "pert_stru":{
            "pert_number": 0,
            "cell_pert_frac": null,
            "atom_pert_dist": null,
            "mag_rotate_angle": null,
            "mag_tilt_angle": null,
            "mag_norm_dist": null
            },
      "pp_dict":{},
      "orb_dict":{},
      "pp_path":,
      "orb_path":,
      "dpks_descriptor":,
      "extra_files":[],
      "link_example_template_extra_files":true,
      "abacus2qe": false,
      "qe_setting":{
        "version": 7.0,
        "electrons": {:},
        "system": {:},
      }
      "abacus2vasp": false,
      "potcar": , // the path of POTCARs, or a dict of POTCAR, such as: {"H":"H.psp8","O":"O.psp8"}
      "vasp_setting": {
        "emax_coef": 1.5  // the coefficient of ecutwfc, the real ENCUT = E_MAX * emax_coef, where E_MAX is the recommended value in POTCAR
      },
      "abacus2cp2k": false,
      "cp2k_setting":{},
      "link_example_template_extra_files":true
  }
}
```
Only key "prepare" is recongnized by `abacustest prepare`.    

- `example_template`: a list of examples, should at least has `STRU` file in each example folder. 
- `input_template`: the template of INPUT file. If is not null, all example will use this file as INPUT file. 
- `kpt_template`: the template of KPT file. If is not null, all example will use this file as KPT file. 
- `strus`/`stru_format`: Translate the structure file from `stru_format` to ABACUS stru. `strus` should be str or list of str which is the path that dpdata can read, and `stru_format` should be str that dpdata supportted. The common format is "poscar" for POSCAR, "deepmd/npy" for deepmd npy. Special format is "cif" for CIF file. If these two keys are setted, then `example_template` will be ignored, and the translated STRU will be put in a new folder named by 000000, 000001, 000002, ...
- `mix_input`: the mix of INPUT parameters. If is not null, will generate all combinations of the parameters for each example. Such as: {"ecutwfc":[50,60,70],"kspacing":[0.1,0.12,0.13]}, will generate 9 INPUTs. If one need combine the parameters, should use '|' to connect them, and the value should be also combined by '|'. Such as: {"ecutwfc|kspacing":["50|0.1","60|0.12","70|0.13"]}. 
- `mix_kpt`: the mix of KPT parameters. If is not null, will generate all combinations of the parameters for each example. There are three types to define the kpt parameters: 
    - One Int number defines the K points in a/b/c direction, and the shift in G space is (0 0 0). Such as: 4, means 4 4 4 0 0 0. 
    - Three Int number defines the K points in a/b/c direction, and the shift in G space is (0 0 0). Such as: [4,4,4], means 4 4 4 0 0 0. 
    - Three Int number defines the K points in a/b/c direction, and three Float defines the shift in G space. Such as: [4,4,4,1,1,1], means 4 4 4 1 1 1.  
    So, an example of mix_kpt can be: [2,[3,3,3],[4,4,4,1,1,1]] 
- `pert_stru`: the perturbation of structure. If is not null, will generate the perturbed structure for each example. 
    - `pert_number`: the number of perturbed structures. If is 0, will not generate the perturbed structure. This is the final number of examples, and each new structure will be perturbed based on below parameters.
    - `cell_pert_frac`: the maximum perturbation of cell. If is null, will not perturb the cell. Such as 0.01, means the cell vector will be perturbed by maximum 1%. 
    - `atom_pert_dist`: the maximum perturbation distance of atom. Unit in Angstrom.
    - `mag_rotate_angle`: the angle of rotation of magnetic moment. Unit in degree. If is null, will not rotate the magnetic moment. 
    - `mag_tilt_angle`: the angle of tilt of magnetic moment. Unit in degree. If is null, will not tilt the magnetic moment. 
    - `mag_norm_dist`: the distance of magnetic moment to the normal direction. If is null, will not move the magnetic moment. Unit in $\mu_B$.\
    
    Note1: the perturbation of magnetic moment is only valid for the spin-constrained atom that the "sc" of at leaset one magnetic component is 1.
    
    Note2: the value can also be a list of two values, which means the minimum and maximum value of perturbation. Such as for `atom_pert_dist`: [0.1,0.15], means the perturbation distance of atom is between 0.1 and 0.15 Angstrom, and if the two values are same, the perturbation distance is fixed.
- `pp_dict`: the pseudopotential dict. The key is the element name, and the value is the pseudopotential file name. Such as: {"H":"H.psp8","O":"O.psp8"}. 
- `orb_dict`: the orbital dict. The key is the element name, and the value is the orbital file name. Such as: {"H":"H.orb","O":"O.orb"}. 
- `pp_path`: the path of pseudopotential files. There should has an extra "element.json" file that defines the element name and the pseudopotential file name. Such as: {"H":"H.psp8","O":"O.psp8"}, or abacustest will read the first two letters of the pseudopotential file name as the element name. If one element has been defined in both `pp_dict` and `pp_path`, the value in `pp_dict` will be used. 
- `orb_path`: the path of orbital files. There should has an extra "element.json" file that defines the element name and the orbital file name. Such as: {"H":"H.orb","O":"O.orb"}, or abacustest will read the first two letters of the orbital file name as the element name. If one element has been defined in both `orb_dict` and `orb_path`, the value in `orb_dict` will be used. 
- `dpks_descriptor`: the descriptor of dpks. If is not null, will link the dpks file for each example. 
- `extra_files`: the extra files that will be lniked to each example folder. Such as: ["abc.py","def.json"].
- `link_example_template_extra_files`: if link the extra files in each example folder. Default is true (will link or copy all files in example folder). If is false, will only copy or link the required files for ABACUS job, such as: INPUT/KPT/STRU/.upf/.orb. 
- `abacus2qe`: if convert the ABACUS input to QE input. Default is false. Now support the convert of cell, coordinate, kpt, pp, and normal scf/relax/cell-relax calculation, and settings of symmetry/smearing/mixing/scf_thr/force_thr/stress_thr, and the magnetic setting of atom type.
- `qe_setting`: the setting of QE input. Here has three types settings: 
    - specify params in "system","control","electrons","ions","cell". Like: "system": {"ibrav":0},...
    - specify a block, and key is the title of block, and value is a list of all block lines. Like: "HUBBARD (ortho-atom)" : ["U Fe1-3d 5.3",..]
    - some special keys: 
        - "version": the version of QE, default is 7.0. In different version, the format of QE input is different.
- `abacus2vasp`: if convert the ABACUS input to VASP input. Default is false. Now support the convert of cell, coordinate, kpt, and normal scf/relax/cell-relax calculation, and settings of symmetry/smearing/mixing/scf_thr/force_thr/dft_plus_u/nupdown/lspinorb/noncolin, and the magnetic setting of atom type. Will set EDIFF based on scf_thr, and is scf_thr/1e-2 for PW and scf_thr/1e-1 for LCAO.
- `potcar`: the path of POTCARs, or a dict of POTCAR, such as: {"H":"H.psp8","O":"O.psp8"}. If is a path, there should has some subfolders, and the folder name is the element name, and there should has a POTCAR file in each subfolder. If is a dict, the key is the element name, and the value is the POTCAR file name.
- `vasp_setting`: two types of settings:
    - specify the value if INCAR. Like: "ENCUT": 500, "EDIFF": 1e-5, ...
    - some special settings:
        - "emax_coef": the coefficient of ecutwfc, the real ENCUT = E_MAX * emax_coef, where E_MAX is the recommended value in POTCAR. If it is not defined, then ENCUT=ecutwfc * Ry2eV
- `abacus2cp2k`: if convert the abacus input to cp2k input. Detault is false. Now support the convert of cell/coordinate/kpt/calculation/force_thr/stress_thr/smearing/mixing/scf_thr(for PW = abacus_value*1e3, for LCAO=abacus_value*1e2).
- `cp2k_setting`: some extra setting of cp2k. Should be a dict, and the key is the name of cp2k input, and the value is the value of the input. Such as: {"FORCE_EVAL": {"DFT": {"SCF": {"EPS_SCF": 1e-6}}}}.
- `link_example_template_extra_files`: if link the extra files in each example folder. Default is true (will link or copy all files in example folder). If is false, will only copy or link the required files for ABACUS job, such as: INPUT/KPT/STRU/.upf/.orb.

If there has more than two types of mixing, will put inputs in a subfolder named by 00000, 00001, 00002, ...

## 6. download
``` 
usage: abacustest download [-h] [-p PARAM] [-s SAVE] job_id

Download the results of the dflow job

positional arguments:
  job_id                the job id of dflow

optional arguments:
  -h, --help            show this help message and exit
  -p PARAM, --param PARAM
                        the file for bohrium account information, default is "job.json"
  -s SAVE, --save SAVE  the folder where the results will be put in, default: result
```
The job_id is the job id of dflow, and the param file is the same as the param file of `submit`. \

## 7. example
[examples](https://github.com/pxlxingliang/abacus-test/tree/develop/example)

