# abacustest
Do the performance test of ABACUS \
Install:
`pip install .`
And then, you can use command `abacustest`. \
There are 4 sub_commands:
- `submit`
- `collectdata`
- `outresult`
- `status`

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
job.json (-p) is needed as input

### 1.1 job.json
This file defines the detail of the jobs. \
An example is like:
```
{
    "if_run":{
        "abacus-pw": false,
        "abacus-pw1": false,
    },

    "config":{
        "lbg_username":         "xxx",
        "lbg_password":         "xxx",
        "project_id":           111
    },
    "ABBREVIATION":{
            "ABACUS310_IMAGE": "registry.dp.tech/dptech/abacus:3.1.0",
            "PYTHON_IMAGE": "python:3.8"
    },

    "param":{
        "abacus-pw":{
            "save_path":"result/abacus-pw",
            "run_dft":[
                {"ifrun": true,
                 "sub_save_path": "",
                 "image": "ABACUS310_IMAGE",
                 "example":[["00[0-2]"],"00[3-5]"],
                 "ngroup" : 0,
                 "bohrium": {"scass_type":"c8_m16_cpu","job_type":"container","platform":"ali"},
                 "command": "mpirun -np 8 abacus > log",
                 "extra_files":["collectdata-abacus.json"],
                 "outputs":["log","result.json","OUT.*"]
                },
                {"ifrun": true,
                 "image": "ABACUS310_IMAGE",
                 "example":[["00[6-7]"],"00[8-9]"],
                 "ngroup" : 3,
                 "bohrium": {"scass_type":"c16_m32_cpu","job_type":"container","platform":"ali"},
                 "command": "mpirun -np 8 abacus > log",
                 "extra_files":[],
                 "outputs":[]
                }
            ]
        },
        "abacus-pw1":{
            "save_path":"result/abacus-pw1",
            "run_dft":[
                {"ifrun": true,
                 "image": "registry.dp.tech/dptech/abacus:3.1.0",
                 "example":[["00[0-2]"],"00[3-5]"],
                 "ngroup" : 0,
                 "bohrium": {"scass_type":"c8_m16_cpu","job_type":"container","platform":"ali"},
                 "command": "mpirun -np 4 abacus > log && collectdata.py collectdata-abacus.json -o result.json",
                 "extra_files":["collectdata-abacus.json"],
                 "outputs":["log","result.json","OUT.*"]
                }
            ],
            "post_dft":{
                        "ifrun": false,
                        "command": "collectdata.py collectdata-abacus.json -o result.json -j 00*",
                        "extra_files": ["collectScf.py"],
                        "image":   "python:3.8",
                        "outputs": []
            }
        }
    }
}
```
- `if_run`: define if the test will be run.
- `ABBREVIATION`: define some abbreviation, and is only valid for `image`.
- `param`: the detail setting for each test.
  - test-name: the name of your test, and is same to that in `if_run`. If the name of one test is not defined in `if_run`, it will be ignored.
    - `save_path`: define the path to save the results of this test. If this key is not defined or if is deined to be "" or None, the value will be replaced to be the path defined by "-s" (the args of abacustest, the default of "-s" is "result/date_of_today" like "result/20230101")
    - `run_dft`: define the detail of the running of your jobs. The value is a list of dictionaries. You can set any number of dictionaries.
      - `ifrun`: if set it to be `false`, will skip this part. 
      - `sub_save_path`: the path to save the results of this part in `save_path`, which means the real save path will be "save_path/sub_save_path". If this key is not defined, or is defined to be "" or null, the real save path will be "save_path".
      - `iamge`: define the image name. Also you can use the name defined in `ABBREVIATION`. Progrma will firstly check if the image name is defined in `ABBREVIATION`, if yes, the name will be replaced by the value in `ABBREVIATION`, and if no, the iamge name will be kept to be value of `iamge`.
      - `example`: The folder names of your jobs. Here assume that you have 5 jobs and the folder names are 000, 001, 002, ..., 004. You can write as ["000", "001", "002", "003", "004"], and also you can write as ["00[0-4]"] that can be recongnized by `glob.glob`. Besides, you can put some folders in a list, and they will be set as one group and run the jobs serially. Such as: [["00[0-1]"],"00[2-4]"], "000" and "001" is one group, each of "002", "003" and "004" is one group and total 4 groups.
      - `ngroup`: split `example` to `ngroup` groups, and will apply `ngroup` machines to run the jobs.
      - `bohrium`: if you want to submit the job to Bohrium, you need set this key, and if you do not want to use Bohrium, please remove this key or set it to be `false`.
      - `command`: the command to run the job. It is same for all jobs defined in `example`.
      - `extra_files`: if you need extra files (such as "collectdata-abacus.json"), you can defined them at here. Before run the command, these files will be copied to each example folder.
      - `outputs`: specify which files and folders will be downloaded. If set to be "[]", all of the files in the folder will be downloaded.\
  
    - `post_dft`: define the detail of post processing, and now all examples will be put at one same place. The key are same as those in `run_dft`, but no need the definition of example. 

### 1.2 user.json
This file defines the detail of the account information to Bohrium. \
An example is like: 
```
{
        "lbg_username":         "xxx",
        "lbg_password":         "xxx",
        "project_id":           111
}
```
If you do not use Bohrium, this file and the keys is also needed, but the values are not matter.

### 1.3 submit a test
```
abacustest submit -p job.json -u user.json -s result/test
```


## 2. collectdata
```
usage: abacustest collectdata [-h] [-j [JOBS [JOBS ...]]] [-t {0,1,2}] [-p PARAM] [-o OUTPUT]
                              [--outparam [OUTPARAM]]

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
  --outparam [OUTPARAM]
                        output the registed parameters, you can set the type by -t or --type to choose
                        abacus/qe/vasp. 0: No, 1: yes
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
             version:	           Abacus.GetVersion()	the version of ABACUS
               ncore:	             Abacus.GetNcore()	the mpi cores
          normal_end:	         Abacus.GetNormalEnd()	if the job is nromal ending
               INPUT:	    Abacus.GetInputParameter()	a dict to store the setting in OUT.xxx/INPUT
...
```

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
from lib_collectdata.resultAbacus import ResultAbacus
class MyAbacus(ResultAbacus):
    @ResultAbacus.register(key_name="description of the key")
    def function_name(self):
        "your codes the get the value of key"
        self[key_name] = value
``` 
If you want to define a vasp method, you need change the import line and class name, like:
```
from lib_collectdata.resultVasp import ResultVasp
class MyVasp(ResultVasp):
    @ResultAbacus.register(key_name="description of the key")
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
usage: abacustest outresult [-h] [-p PARAM]
```
Out put some specified parameters of some jobs, and also can calculate some specified metrics.


## 4. status
```
usage: abacustest status [-h] paramf job_id
```
Not ready now.


## 5. example
some exmaples

