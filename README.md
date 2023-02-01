# abacustest
do the performance test of ABACUS
install:
`pip install .`

There are two commands:
- `abacustest`
- `collectdata`
- `outresult`

Please use `abacustest/collectdata/outresult -h` to get the usage

## abacustest
```
usage: abacustest [-h] [-p PARAM] [-u USER] [-s SAVE] [--override OVERRIDE] [--outinfo OUTINFO]. 
```
Two files are needed as input: job.json and user.json

### job.json
This file defines the detail of the jobs. \
An example is like:
```
{
    "if_run":{
        "abacus-pw": false,
        "abacus-pw1": false,
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
                 "image": "ABACUS310_IMAGE",
                 "example":[["00[0-2]"],"abacus-pw/00[3-5]"],
                 "ngroup" : 0,
                 "bohrium": {"scass_type":"c8_m16_cpu","job_type":"container","platform":"ali"},
                 "command": "mpirun -np 8 abacus > log",
                 "collectdata_script":["collectdata-abacus.json"],
                 "collectdata_command":"collectdata.py SCRIPT/collectdata-abacus.json -o result.json",
                 "outputs":["log","result.json","OUT.*"]
                },
                {"ifrun": true,
                 "image": "ABACUS310_IMAGE",
                 "example":[["00[6-7]"],"abacus-pw/00[8-9]"],
                 "ngroup" : 3,
                 "bohrium": {"scass_type":"c16_m32_cpu","job_type":"container","platform":"ali"},
                 "command": "mpirun -np 8 abacus > log",
                 "collectdata_script":[],
                 "collectdata_command":"",
                 "outputs":[]
                }
            ]
        },
        "abacus-pw1":{
            "save_path":"result/abacus-pw1",
            "run_dft":[
                {"ifrun": true,
                 "image": "registry.dp.tech/dptech/abacus:3.1.0",
                 "example":[["00[0-2]"],"abacus-pw/00[3-5]"],
                 "ngroup" : 0,
                 "bohrium": {"scass_type":"c8_m16_cpu","job_type":"container","platform":"ali"},
                 "command": "mpirun -np 4 abacus > log",
                 "collectdata_script":["collectdata-abacus.json"],
                 "collectdata_command":"collectdata.py SCRIPT/collectdata-abacus.json -o result.json",
                 "outputs":["log","result.json","OUT.*"]
                }
            ],
            "post_dft":{
                        "ifrun": false,
                        "script": ["collectScf.py"],
                        "image":   "python:3.8",
                        "command": "",
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
    - `save_path`: define the path to save the results of this test
    - `run_dft`: define the detail of the running of your jobs. The value is a list of dictionaries. You can set any number of dictionaries.
      - `ifrun`: if set it to be `false`, will skip this part. 
      - `iamge`: define the image name. Also you can use the name defined in `ABBREVIATION`. Progrma will firstly check if the image name is defined in `ABBREVIATION`, if yes, the name will be replaced by the value in `ABBREVIATION`, and if no, the iamge name will be kept to be value of `iamge`.
      - `example`: The folder names of your jobs. Here assume that you have 5 jobs and the folder names are 000, 001, 002, ..., 004. You can write as ["000", "001", "002", "003", "004"], and also you can write as ["00[0-4]"] that can be recongnized by `glob.glob`. Besides, you can put some folders in a list, and they will be set as one group and run the jobs serially. Such as: [["00[0-1]"],"00[2-4]"], "000" and "001" is one group, each of "002", "003" and "004" is one group and total 4 groups.
      - `ngroup`: split `example` to `ngroup` groups, and will apply `ngroup` machines to run the jobs.
      - `bohrium`: if you want to submit the job to Bohrium, you need set this key, and if you do not want to use Bohrium, please remove this key or set it to be `false`.
      - `command`: the command to run the job. It is same for all jobs defined in `example`.
      - `collectdata_script`: if you want to do some post processing, and it need extra files (such as "collectdata-abacus.json"), you can defined them at here.
      - `collectdata_command`: the command to run the post processing. BE NOTICED: the command will be executed in the folder defined in `example`. Such as "000", program will actually `cd 000` and then execute the command. Besides, program will create a new folder "SCRIPT" in each job floder, and copy the files defined in `collectdata_script` to SCRIPT. So if the command need specify the file in `collectdata_script`, you need add the "SCRIPT" in the path. Such as command: "collectdata.py SCRIPT/collectdata-abacus.json". If you want to use collectdata function of this program, you can directly use "collectdata.py".
      - `outputs`: specify which files and folders will be downloaded. If set to be "[]", all of the files in the folder will be downloaded.\
  
    - `post_dft`: define the detail of post processing, and now all examples will be put at one same place. The usage of key `ifurn`,`image`,`bohrium`,`command`,`outputs` are same as those in `run_dft`. Key `script` is to define the files used in this step.

### user.json
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

### submit a test
```
abacustest -p job.json -u user.json -s result/test
```


## collectdata
```
usage: collectdata [-h] [-j [JOBS [JOBS ...]]] [-t {0,1,2}] [-p PARAM] [-o OUTPUT] [--outparam [OUTPARAM]]
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
If the value of one key (such as "INPUT") is a dictionary, you can write as {key: [key1_in_dict, key2_in_dict]}. 



## example
some exmaples

