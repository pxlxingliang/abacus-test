# abacustest
abacustest is mainly designed to handle pre-processing and post-processing for calculations performed by ABACUS. Moreover, it supports high-throughput submission of computing tasks to the Bohrium cloud platform and supercomputers for calculation.

# Installation
- Install via pip:
```
pip install abacustest
```

- Or install from source code:
```
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

# Usage
After installation, you can use the command-line tool `abacustest`.​
Please use `abacustest -h` to get the usage instructions.​

**NOTE**: in abacustest, an **example** is a folder that contains the ABACUS input files, such as INPUT, KPT, STRU, etc. 

## 1. prepare
The `prepare` function of abacustest supports preparing input files for ABACUS. It requires a parameter configuration file in JSON format (e.g., `param.json`).
To initiate the input file preparation process, execute the command:
```
abacustest prepare -p param.json -s abacustest
```

`-s` specifies the directory where the input files will be stored. If not specified, it defaults to "abacustest". 

The structure of `param.json` is as follows:
```
{
  "prepare": {
    "example_template": ["example1", "example2"],
    "input_template": "INPUT",
    ...
  }
}
```
Only the key "prepare" is recognized by `abacustest prepare`, and its value is a dictionary that defines the parameters required for preparing the input files. The parameters containing the following sections:

### 1.1 Converting Structure Files to ABACUS STRU Format

`abacustest` supports converting structure files into ABACUS input files. This is useful for preparing the STRU files needed for ABACUS calculations. \

The `param.json` may be like:
```
{
  "prepare":{
    "strus": ["1.vasp", "2.vasp"],
    "stru_format": "poscar",
    "input_template": "INPUT",
    "kpt_template": "KPT",
    "pp_dict": {"H": "H.upf", "O": "O.upf"},
    "orb_dict": {"H": "H.orb", "O": "O.orb"},
    "pp_path": "path/to/pp",
    "orb_path": "path/to/orb"
  }
}
```

 The configuration file and its description are divided into the following sections:

*   `strus`/`stru_format`: Convert structure files from `stru_format` to ABACUS STRU format. `strus` should be a string or list of strings representing paths readable by dpdata, and `stru_format` should be a string representing a format supported by dpdata. Common formats include "poscar" for POSCAR files and "deepmd/npy" for deepmd npy files. A special format is "cif" for CIF files, which will be parsed by ASE. The converted STRU files will be placed in new folders named 000000, 000001, 000002, etc. 


*   `input_template`: The template for the INPUT file. If specified, all examples will use this file as their INPUT file.


*   `kpt_template`: The template for the KPT file. If specified, all examples will use this file as their KPT file.


*   `pp_dict`: A dictionary of pseudopotentials. Keys are element names, and values are pseudopotential filenames. For example: `{"H": "H.upf", "O": "O.upf"}`. The related files will be written into STRU file, and create a soft link to the pseudopotential file in each example folder.


*   `orb_dict`: A dictionary of orbitals. Keys are element names, and values are orbital filenames. For example: `{"H": "H.orb", "O": "O.orb"}`. The related files will be written into STRU file, and create a soft link to the orbital file in each example folder.


*   `pp_path`: The path to pseudopotential library. An additional "element.json" JSON file should be present, defining element names and corresponding pseudopotential filenames (e.g., `{"H": "H.psp8", "O": "O.psp8"}`). If not, abacustest will use the first two letters of pseudopotential filenames as element names. If an element is defined in both `pp_dict` and `pp_path`, the `pp_dict` value takes precedence.


*   `orb_path`: The path to orbital library. An additional "element.json" file should be present, defining element names and corresponding orbital filenames (e.g., `{"H": "H.orb", "O": "O.orb"}`). If not, abacustest will use the first two letters of orbital filenames as element names. If an element is defined in both `orb_dict` and `orb_path`, the `orb_dict` value takes precedence.

### 1.2 Preparing Input Files for ABACUS Parameter Testing 

`abacustest` also supports preparing input files for ABACUS parameter testing. This is useful for generating multiple input files with different parameters for high-throughput calculations. \

The `param.json` may be like:
```
{
  "prepare": {
    "example_template": ["example1", "example2"],
    "mix_input": {"ecutwfc": [50, 60, 70], 
    "kspacing": [0.1, 0.12]},
    "mix_kpt": [2, [3, 3, 3], [4, 4, 4, 1, 1, 1]]
  }
}
```

*   `example_template`: A list of examples, each of which must contain at least a `STRU` file.

*   `mix_input`: Combinations of INPUT parameters. If specified, all combinations of the parameters will be generated for each example. For instance, `{"ecutwfc": [50, 60, 70], "kspacing": [0.1, 0.12]}` will generate 3*2 = 6 INPUT files. To combine parameters, use '|' to connect them, with values also combined by '|'. For example: `{"ecutwfc|kspacing": ["50|0.1", "60|0.12", "70|0.13"]}`.


*   `mix_kpt`: Combinations of KPT parameters that automatically generate K points by Gamma-centered Monkhorst-Pack method. If specified, all combinations of the parameters will be generated for each example. There are three ways to define KPT parameters:

    Thus, an example of `mix_kpt` could be: `[2, [3, 3, 3], [4, 4, 4, 1, 1, 1]]`


    *   A single integer defines the number of K points in the a/b/c directions, with a shift of (0 0 0) in G space. For example, 4 represents 4 4 4 0 0 0.


    *   Three integers define the number of K points in the a/b/c directions, with a shift of (0 0 0) in G space. For example, `[4, 4, 4]` represents 4 4 4 0 0 0.


    *   Three integers define the number of K points in the a/b/c directions, and three floats define the shift in G space. For example, `[4, 4, 4, 1, 1, 1]` represents 4 4 4 1 1 1.

**NOTE**: if the mixed inputs are large than 1, then `abacustest` will create sub-folders named 000000, 000001, 000002, etc. in each example folder.

### 1.3 Perturbing Configurations

`abacustest` supports perturbing configurations to generate new structures for ABACUS calculations.

The `param.json` may be like:
```
{
  "prepare": {
    "example_template": ["example1", "example2"],
    "pert_stru": {
      "pert_number": 10,
      "cell_pert_frac": 0.01,
      "atom_pert_dist": 0.1,
      "mag_rotate_angle": 5,
      "mag_tilt_angle": 10,
      "mag_norm_dist": 0.05
    }
  }
}
```

*   `pert_stru`: Random structure perturbations. If specified, perturbed structures will be generated for each example.
 
    *   `pert_number`: The number of perturbed structures. If set to 0, no perturbed structures will be generated. This is the final number of examples, with each new structure perturbed based on the parameters below.


    *   `cell_pert_frac`: The maximum perturbation fraction for the unit cell. If null, the cell will not be perturbed. For example, 0.01 means the cell vectors will be perturbed by a maximum of 1%.


    *   `atom_pert_dist`: The maximum perturbation distance for atoms, in Angstrom.


    *   `mag_rotate_angle`: The rotation angle of magnetic moments for all magnetic atoms, in degrees. If null, magnetic moments will not be rotated. All magnetic atoms will be rotated by the same angle.


    *   `mag_tilt_angle`: The tilt angle for magnetic moments, in degrees. If null, magnetic moments will not be tilted. The tilt is applied to all magnetic atoms, but for each atom the tilt angle is randomly chosen.


    *   `mag_norm_dist`: The distance of magnetic moments from the normal direction, in $\mu_B$. If null, magnetic moments will not be displaced.

  Note 1: Magnetic moment perturbations are only valid for spin-constrained atoms where the "sc" value of at least one magnetic component is 1.

  Note 2: Values can also be a list of two numbers, representing the minimum and maximum perturbation values. For example, `atom_pert_dist: [0.1, 0.15]` means atomic perturbation distances range from 0.1 to 0.15 Angstrom. If both values are the same, the perturbation distance is fixed.

### 1.4 Converting ABACUS Input Files to VASP Input Files

`abacustest` supports converting ABACUS input files to VASP input files. This is useful for preparing VASP calculations based on existing ABACUS input files.

The `param.json` may be like:
```
{
  "prepare": {
    "abacus2vasp": true,
    "potcar": {"H": "H.pot", "O": "O.pot"},
    "vasp_setting": {
      "ENCUT": 500,
      "EDIFF": 1e-5,
      "emax_coef": 1.0
    }
  }
}
```

*   `abacus2vasp`: Whether to convert ABACUS inputs to VASP inputs. Defaults to false. Currently supports conversion of cell, coordinates, k-points, standard scf/relax/cell-relax calculations, symmetry, smearing, charge mixing, scf\_thr, force\_thr, dft\_plus\_u, nupdown, lspinorb, noncolin, and atomic magnetic settings. EDIFF is set based on scf\_thr: scf\_thr/1e-2 for PW and scf\_thr/1e-1 for LCAO.


*   `potcar`: The path to POTCAR files or a dictionary of POTCARs (e.g., `{"H": "H.psp8", "O": "O.psp8"}`). If a path, subfolders named after elements should contain POTCAR files, which the organizational form is consistent with the pseudopotential organizational form provided by VASP (e.g., vasp/PAW_PBE). If a dictionary, keys are element names and values are POTCAR filenames.


*   `vasp_setting`: Additional settings for VASP, and two types of settings:

    *   Specify values in INCAR. For example: `"ENCUT": 500, "EDIFF": 1e-5, ...`.

    *   Special settings:

        *   "emax\_coef": The coefficient for ecutwfc. The actual ENCUT = E\_MAX \* emax\_coef, where E\_MAX is the recommended value in the POTCAR. If not defined, ENCUT = ecutwfc \* Ry2eV. If ENCUT is also defined, it will be used directly.

### 1.5 Converting ABACUS Input Files to QE Input Files
`abacustest` supports converting ABACUS input files to QE input files. This is useful for preparing QE calculations based on existing ABACUS input files.

The `param.json` may be like:
```
{
  "prepare": {
    "abacus2qe": true,
    "qe_setting": {
      "version": "7.0",
      "system": {
        "ibrav": 0,
        "ecutwfc": 50
        }
    }
  }
}
```

*   `abacus2qe`: Whether to convert ABACUS inputs to QE inputs. Defaults to false. Currently supports conversion of cell, coordinates, k-points, standard scf/relax/cell-relax calculations, symmetry, smearing, mixing, scf\_thr, force\_thr, stress\_thr, atomic magnetic settings, and hubbard U settings.


*   `qe_setting`: Additional settings for QE inputs. There are three types of settings:



    *   Specify parameters in "system", "control", "electrons", "ions", "cell". For example: `"system": {"ibrav": 0}, ...`.


    *   Specify a block, where the key is the block title and the value is a list of all lines in the block. For example: `"HUBBARD (ortho-atom)": ["U Fe1-3d 5.3", ...]`.


    *   Special keys:

        *   "version": The QE version, defaulting to 7.0. QE input formats vary between versions.

### 1.6 Converting ABACUS Input Files to CP2K Input Files
`abacustest` supports converting ABACUS input files to CP2K input files. This is useful for preparing CP2K calculations based on existing ABACUS input files.

The `param.json` may be like:
```
{
  "prepare": {
    "abacus2cp2k": true,
    "cp2k_setting": {
      "FORCE_EVAL": {
        "DFT": {
          "SCF": {
            "EPS_SCF": 1e-6
          }
        }
      }
    }
  }
}
```


*   `abacus2cp2k`: Whether to convert ABACUS inputs to CP2K inputs. Defaults to false. Currently supports conversion of cell, coordinates, k-points, calculation type, force\_thr, stress\_thr, smearing, mixing, and scf\_thr (for PW = abacus\_value*1e3, for LCAO = abacus\_value*1e2).


*   `cp2k_setting`: Additional settings for CP2K inputs. Should be a dictionary where keys are CP2K input names and values are the corresponding input values. For example: `{"FORCE_EVAL": {"DFT": {"SCF": {"EPS_SCF": 1e-6}}}}`.


### 1.7 Other Related Parameters
Additional parameters that control file handling and auxiliary settings.

*   `link_example_template_extra_files`: Bool, Whether to link extra files in each example folder. Defaults to true (links or copies all files in the example folder). If false, only copies or links files required for ABACUS jobs (e.g., INPUT, KPT, STRU, .upf, .orb).


*   `extra_files`: List, Additional files to be linked to each generated example folder. For example: `["abc.py", "def.json"]`.


*   `dpks_descriptor`: Str, The descriptor for DeepKS calculation. If specified, the dpks file will be linked to each example.


**NOTE**: All prepare parameters can be used in combination. The actual running sequence of abacustest prepare is as follows:​
1. Convert the structure defined in `strus` to ABACUS input files 
2. Use the generated input files as `example_template`, and perform configuration perturbations to generate a series of sub-examples.
3. Do INPUT/KPT mixing for each sub-examples.
4. Convert all generated sub-examples from ABACUS to inputs of other software.​

In these processes, the following rules apply:
- When both `example_template` and `strus` are specified, `example_template` is ignored. 
- When `example_template` contains INPUT or KPT, and `input_template` or `kpt_template` is specified at the same time, the INPUT/KPT in example_template is ignored.



## 2. submit

`abacustest` supports submitting jobs to Bohrium and other supercomputers.​

Before preparing to submit calculations, you need to prepare the `example`s (each `example` is a folder that containing the inputs for one ABACUS job), as well as an additional configuration file (such as "param.json") that defines the calculation parameters and settings.​

**NOTE**: since tasks are transferred to remote computing resources (such as Bohrium or other supercomputer platforms), **do not use absolute paths** in INPUT or STRU files to define pseudopotential or orbital files.

### 2.1 Configuration File
A configuration file mainly contains the following content:
```
{
    "bohrium_group_name": "abacustest",
    "save_path":"results",
    "max_parallel": 100,
    "run_dft":[
        {
          "ifrun": true,
          "sub_save_path": "123",
          "image": "registry.dp.tech/dptech/abacus-stable:LTSv3.10",
          "example":[["00[0-2]"],"00[3-5]"],
          "group_size" : 1,
          "bohrium": {
            "scass_type":"c8_m16_cpu",
            "job_type":"container",
            "platform":"ali"
          },
          "command": "mpirun -np 8 abacus > log",
          "extra_files":[]
        }
    ]
}
```
Where:

*   `save_path`: Defines the path to save the results. If this key is not defined, or is set to "" or `None`, the value will be replaced by the path defined by "-s" (a parameter of `abacustest`; the default for "-s" is "result").
*   `max_parallel`: Defines the maximum number of parallel jobs. The default is 100.
*   `run_dft`: Defines details of job execution. The value is a list of dictionaries, and you can set any number of dictionaries.
    *   `ifrun`: If set to `false`, this part will be skipped.
    *   `sub_save_path`: The path within `save_path` where results of this part will be saved, meaning the actual save path will be "save\_path/sub\_save\_path". If this key is not defined, or is set to "" or `null`, the actual save path will be "save\_path".
    *   `image`: Defines the image name. You can also use names defined in `ABBREVIATION`. The program will first check if the image name is defined in `ABBREVIATION`; if so, the name will be replaced with the value from `ABBREVIATION`; otherwise, the image name remains as specified.
    *   `example`: The folder names of your jobs. Assuming you have 5 jobs with folder names 000, 001, 002, ..., 004, you can write this as \["000", "001", "002", "003", "004"] or use glob patterns like \["00\[0-4]"].
    *   `group_size`: Defines how many example groups to run on a single machine, default is 1.
    *   `bohrium`: If you want to submit jobs to Bohrium, this key must be set. If Bohrium is not used, please remove it.
    *   `command`: The command to run the job, which is the same for all jobs defined in `example`.
    *   `extra_files`: If additional files are needed (e.g., "collectdata-abacus.json"), they can be defined here. These files will be copied to each example folder before the command is executed.

The `bohrium` can be replaced by `dispatcher` (https://docs.deepmodeling.com/projects/dpdispatcher/en/latest/), which support submitting jobs to more platforms, such as supercomputers. The `dispatcher` is a dictionary, and the content is like:
```
    "dispatcher": {
        "machine_dict": 
        {
          "remote_root": "/home/username/work_path",
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

Key "upload_packages" can be added in the configuration file to specify some local packages that need to be uploaded to the remote platform, such as:
```
{
    ...
    "upload_packages": ["my_package1","/home/user/my_package2"],
    ...
}
```
where the value is a list of local python packages. The packages will be uploaded to the remote platform before the job is submitted. This is useful when the remote image lacks some python packages that are needed by your job.

### 2.2 Submit a Job
Before submitting a job, you need to set below environment variables:
```
export BOHRIUM_USERNAME="your_bohrium_username"
export BOHRIUM_PASSWORD="your_bohrium_password"
export BOHRIUM_PROJECT_ID="your_bohrium_project_id"
```
These variables are used to authenticate your Bohrium account.

Or, you can add below content in the param.json, like:
```
{
    "config": {
        "bohrium_username": "<your username>",
        "bohrium_password": "<your password>",
        "bohrium_project_id": "<your project id>"
    },
    "run_dft":[...] 
}
```

To submit a job, you can use the command:
```
abacustest submit -p param.json &
```
where "&" is used to run the command in the background, and `abacustest` will automatically download the results after the job is finished. 

After the job is submitted, you will see output similar to:
```
Workflow has been submitted (ID: abacustest-kbrb2, UID: f14c5d95-655c-47b4-a709-c9a8138a40cf)
Workflow link: https://workflows.deepmodeling.com/workflows/argo/abacustest-kbrb2
job ID: abacustest-kbrb2, UID: f14c5d95-655c-47b4-a709-c9a8138a40cf
You can track the flow by using your browser to access the URL:
 https://workflows.deepmodeling.com/workflows/argo/abacustest-kbrb2?tab=workflow
```
You can track the job status by visiting the URL provided in the output. If the terminal is closed after the job is submitted, you can download the results later by using the command:
```
abacustest download -p param.json <ID>
```
where `<ID>` is the job ID of the submitted job, such as `abacustest-kbrb2` in the above example. Notice that the job ID may be invalid after a certain period of time, then you can download the results by using the UID, such as `f14c5d95-655c-47b4-a709-c9a8138a40cf` in the above example. 


**NOTICE**: if configuration file defines "prepare" section, then the `abacustest` will prepare the input files before submitting the job, and in this case, the results of "prepare" will be saved in current directory. If "example" is not defined in "run_dft" section, the generated inputs by "prepare" will be used as "example", while if "example" is defined in "run_dft", jobs in "example" will be used.


## 3. collectdata

`abacustest` supports extracting key values from the output of ABACUS/QE/VASP jobs. This functionality is valuable for retrieving critical metrics from calculation results, including energy, forces, stress, and other relevant parameters.

You need to prepare a parameter configuration file (e.g., `param.json`) that specifies the key values to be collected, formatted as follows:
```
{"PARAM":
        ["natom","kpt","ibzk","nelec","nbands","force","stress",
        "scf_steps","total_time","force_time","stress_time","energy_per_atom","band_gap",
        {"SELF_DEFINED_NAME": "An eval string to get the value"}]
}
```
Only the key "PARAM" is recognized by `collectdata`, its value is a list of target keys. These keys can be strings or dictionaries: 
- If a string, it must correspond to a key predefined in the result class. 
- If a dictionary, the key represents a custom name, and the value is an evaluation string enclosed in {} to compute the desired value.

Examples:
- To calculate energy per atom: {"energy per atom": "{energy}/{natom}"}(where "energy" and "natom" are predefined keys, and the result is their quotient).
- To rename the "energy" key to "energy (eV)": {"energy (eV)": "energy"}.

Executing the command:
```
abacustest collectdata -p param.json -j job1 job2 job3
```
will collect the specified keys from the output of the jobs `job1`, `job2`, and `job3`. The results will be saved in a JSON file named `metrics.json` in the current directory. 

If you want to save the results in a different file, you can use the `-o` option. 

If you want to collect data for QE/VASP jobs, you can use the `-t` option to specify the job type, such as `-t 0` for ABACUS, `-t 1` for QE, and `-t 2` for VASP. If not specified, it defaults to ABACUS.

The structure of the `metrics.json` file is as follows:
```
{
    "job1": {
        "natom": 10,
        "nelec": 20,
        "nbands": 12,
        "band_gap": 1.5
    },
    ...
}
```
Where each key corresponds to a job name and contains the collected metrics.

### 3.1 Get the keys
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

### 3.2 import collectdata.RESULT in your python script
You can also use this collectdata function in your python script by `from abacustest.lib_collectdata.collectdata import RESULT`.

For example, if you want to collect the total energy, force, and stress from an ABACUS job, you can use the following code:
```
from abacustest.lib_collectdata.collectdata import RESULT

abacusresult = RESULT(fmt="abacus",path="abacusjob") 
# path is the path to the ABACUS job folder, such as "abacusjob"
# fmt is the format of the job, such as "abacus", "qe", "vasp".

energy = abacusresult["total_energy"] 
force = abacusresult["force"]
stress = abacusresult["stress"]
```

If the key is not supported by the RESULT class, the value will be `None`.

## 4. Models
`abacustest` supports automatically generating calculation input files and postprocessing results for specific ABACUS calculations using predefined models. The core workflow consists of three steps:
1. Use `abacustest model <modelname> prepare [args]` to prepare input files using provided inputs
2. Submit ABACUS calculation jobs. In most cases, the previous step will automatically generate a "setting. json" file. You can directly use 'abacustest submit - p setting. json' (of course, you need to configure your own Bohrium account first, please refer to Part 2) to submit the calculation task to Bohrium. If you are not using the Bohrium platform, you can also check the content of "run_dft/examples" in seting.json, submit all directories defined by this field to your computing platform for calculation, and then perform post-processing after the calculation is completed. 
3. uSE `abacustest model <modelname> post [args]` to extract results from ABACUS calculation and postprocess these results to obtain the target property.

### 4.1 Vacancy formation energy
Formation energy of non-charged vacancies is defined as:
$$ E_\text{f, vac} = E_\text{orig} - E_\text{vac} - \mu_\text{A} = E_\text{orig} - E_\text{vac} - \frac{E_\text{crys}}{A}{n} $$
Where $E_\text{f, vac}$ is the vacancy formation energy, $E_\text{orig}$ is the energy of structure with no vacancies,
and $E_\text{crys}(A){n}$ is energy of the most stable crystale of the element at the vacancy site, n is the number of atoms in the
most stable crystal.

Input structure files (supports CIF, POSCAR, or ABACUS STRU format) or pre-prepared ABACUS input files are required, along with vacancy-related parameters and ABACUS calculation parameters (if structure files are provided).
If structure file are used, use the command such as following:
```
abacustest model vacancy prepare -j stru.cif --ftype cif -i 1 -s 1 1 1 --lcao
```
If abacus input file directory is used, use the command such as following:
```
abacustest model vacancy prepare -j 000000 -i 1 -s 1 1 1
```
For more details, use `abacustest model vacancy prepare -h` to get help message.

After calculation using prepared input files finished, you can use command to extract results and postprocess them:
```
abacustest model vacancy post -j 000000
```
Then vacancy formation energy will be printed and saved to `metrics_vacancy.json`. A file named `ref_energy.txt` will be generated
and contains the reference energy of the vacancy atom and can be used in later calculation. Use `abacustest model vacancy post -h` for more details.

### 4.2 Born effective charge
Born effective charge ($Z_{s,ij}^{*}$, BEC) represents the effective charge response of lattice ions under an electric field, or the response of the system's polarization to ionic displacements, and is usually calculated by:

$Z_{s,ij}^{*}$ = $\frac{\Omega}{e}$ $\frac{\partial P_i}{\partial u_{s,j}}$,
where $\Omega$ is the volume of unit cell, $P_i$ is the polarization in direction $i$, and $u_{s,j}$ is the move of atom $s$ in direction $j$. The polarization can be calculated by berry phase method. 

Based on the above formula, we can calculate the BEC by applying a small displacement to an atom along one direction and then computing the change in the system's polarization before and after the displacement. The command `abacustest model bec prepare` can automatically generate displaced configurations and input files for Berry phase calculations based on the ABACUS "example" you provide. More details can be found by `abacustest model bec prepare -h`.

After the ABACUS jobs generated by `prepare` step are finished, you can use `abacustest model bec post` to postprocess the results. A summary of the BEC will be printed on screen and written in file "bec_summary.txt", like:
```
Born effective charge tensor for:
.:atom0_Ti:
        X       Y       Z
X  7.2792  0.0023 -0.0046
Y  0.0000  7.2835 -0.0020
Z  0.0000  0.0021  7.2830
```
The row and column indices represent the displacement direction and the polarization direction. Note that the polarization calculated via the Berry phase is usually along the direction of the cell vectors. Here, we have performed the calculation for each cell vector direction and then transformed the values onto the Cartesian XYZ axes.

Besides, two files are generated simultaneously: the file "metrics.json" compiles detailed information on each sub-job (unless specifically defined otherwise, the unit of length is angstrom, and the unit of energy is eV); the file "metrics_bec.json" compiles the results for each atom. Note that `p_vec_org` denotes the polarization of the original structure along the cell vectors, and `mod_org` is the modulus of the polarization along the cell vectors (since the polarization of a periodic structure is indeterminate, for a series of consecutive structures, we can calculate the change in polarization by adding a specific modulus to render their polarization relatively continuous). Additionally, `p_vec_disp` represents the change in polarization along the cell vectors caused by atomic displacement in the XYZ directions.


## 5. example
[examples](https://github.com/pxlxingliang/abacus-test/tree/develop/example)
