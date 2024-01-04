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

Please use `abacustest -h` to get the usages

## 1. submit
```
usage: abacustest submit [-h] [-p PARAM] [-s SAVE] [--override OVERRIDE] [--outinfo OUTINFO] [--debug [DEBUG]]

This script is used to run a testing

optional arguments:
  -h, --help            show this help message and exit
  -p PARAM, --param PARAM
                        the job setting file, default is job.json
  -s SAVE, --save SAVE  the folder where the results will be put in, default: result/date_of_today (e.g. result/20230101)
  --override OVERRIDE   when the save folder exists, if override it. 0: no, 1: yes.
  --outinfo OUTINFO     if output detail informations, 0: no, 1: yes
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

    "run_dft":[
        {"ifrun": true,
         "sub_save_path": "",
         "image": "ABACUS310_IMAGE",
         "example":[["00[0-2]"],"00[3-5]"],
         "group_size" : 1,
         "bohrium": {"scass_type":"c8_m16_cpu","job_type":"container","platform":"ali"},
         "command": "mpirun -np 8 abacus > log",
         "extra_files":[],
         "outputs":["log","result.json","OUT.*"]
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
                "command": "collectdata.py collectdata-abacus.json -o result.json -j 00*",
                "extra_files": [],
                "image":   "python:3.8",
                "metrics":{
                    "before_command":true,
			        "dft_type":"abacus",
			        "metrics_name": [],
			        "save_file": "result.json",
			        "newmethods": []
		        },
                "outputs": []
    }
}
```
- ``bohrium_group_name``: If use bohrium, you can define the group name at here.
- `ABBREVIATION`: define some abbreviation, and is only valid for `image`.
- `config`: The setting of config information.
- `save_path`: define the path to save the results of this test. If this key is not defined or if is deined to be "" or None, the value will be replaced to be the path defined by "-s" (the param of abacustest, the default of "-s" is "result/date_of_today" like "result/20230101")
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
                        the file name to store the output results, default is "result.json"
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
        {"INPUT":["ecutwfc","mixing_beta"]}]
}
```
Only one key "PARAM" is recongnized by `collectdata`, the value is a list of keys. \
If the value of one key (such as "INPUT") is a dictionary, you can write as {key: [key1_in_dict, key2_in_dict]}. \
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
        "your codes the get the value of key"
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
      "example_template":["example_path"],
      "input_template":"INPUT",
      "kpt_template":"KPT",
      "stru_template":"STRU",
      "mix_input":{
          "ecutwfc":[50,60,70],
          "kspacing":[0.1,0.12,0.13]
      },
      "mix_kpt":[],
      "mix_stru":[],
      "pp_dict":{},
      "orb_dict":{},
      "pp_path":,
      "orb_path":,
      "dpks_descriptor":,
      "extra_files":[],
      "abacus2qe": false
  }
}
```
Only key "prepare" is recongnized by `abacustest prepare`.    

- `example_template`: the template of example folder, such as: ["example_path"], ["example_path1","example_path2"]. 
- `input_template`: the template of INPUT file. If is not null, all example will use this file as INPUT file. 
- `kpt_template`: the template of KPT file. If is not null, all example will use this file as KPT file. 
- `stru_template`: the template of STRU file. If is not null, all example will use this file as STRU file. 
- `mix_input`: the mix of INPUT parameters. If is not null, will generate all combinations of the parameters for each example. Such as: {"ecutwfc":[50,60,70],"kspacing":[0.1,0.12,0.13]}, will generate 9 INPUTs. 
- `mix_kpt`: the mix of KPT parameters. If is not null, will generate all combinations of the parameters for each example. There are three types to define the kpt parameters: 
    - One Int number defines the K points in a/b/c direction, and the shift in G space is (0 0 0). Such as: 4, means 4 4 4 0 0 0. 
    - Three Int number defines the K points in a/b/c direction, and the shift in G space is (0 0 0). Such as: [4,4,4], means 4 4 4 0 0 0. 
    - Three Int number defines the K points in a/b/c direction, and three Float defines the shift in G space. Such as: [4,4,4,1,1,1], means 4 4 4 1 1 1.  
    So, an example of mix_kpt can be: [2,[3,3,3],[4,4,4,1,1,1]] 
- `mix_stru`: the mix of STRU file. If is not null, will generate all combinations for each example.
- `pp_dict`: the pseudopotential dict. The key is the element name, and the value is the pseudopotential file name. Such as: {"H":"H.psp8","O":"O.psp8"}. 
- `orb_dict`: the orbital dict. The key is the element name, and the value is the orbital file name. Such as: {"H":"H.orb","O":"O.orb"}. 
- `pp_path`: the path of pseudopotential files. There should has an extra "element.json" file that defines the element name and the pseudopotential file name. Such as: {"H":"H.psp8","O":"O.psp8"}, or abacustest will read the first two letters of the pseudopotential file name as the element name. If one element has been defined in both `pp_dict` and `pp_path`, the value in `pp_dict` will be used. 
- `orb_path`: the path of orbital files. There should has an extra "element.json" file that defines the element name and the orbital file name. Such as: {"H":"H.orb","O":"O.orb"}, or abacustest will read the first two letters of the orbital file name as the element name. If one element has been defined in both `orb_dict` and `orb_path`, the value in `orb_dict` will be used. 
- `dpks_descriptor`: the descriptor of dpks. If is not null, will link the dpks file for each example. 
- `extra_files`: the extra files that will be lniked to each example folder. Such as: ["abc.py","def.json"]. 
- `abacus2qe`: if convert the ABACUS input to QE input. Default is false.

If there has more than two types of mixing, will put inputs in a subfolder named by 00000, 00001, 00002, ...


## 6. example
some exmaples

