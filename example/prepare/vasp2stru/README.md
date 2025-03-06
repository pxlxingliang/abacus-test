This is an example to show how to use `abacustest preapre` to translate the VASP POSCAR to ABACUS STRU.

The `setting.json` file is as follows:
```json
{
        "prepare":{
                "strus": ["*.vasp"],
                "stru_format": "poscar",
                "input_template": "INPUT",
                "pp_path": "pp_lib"
        }
}
```
It defines the VASP POSCAR files by `"strus": ["*.vasp"]`, and the format of the input structure file is POSCAR by `"stru_format": "poscar"`. 
Only these two options can also do the translation, but will not set the pseudopotential in STRU. 

The `pp_path` is the path of the pseudopotential library, which is used to set the pseudopotential in STRU. `abacustest` will automatically find the pseudopotential in the library, and set the pseudopotential in STRU. The pseudopotential library should contain one file `element.json` that defines the pseudopotential of each element, or the pseudopotential file name should be started with the element symbol and followed by a non-alphabet character, like `C.pbe.UPF`, `Al_pbe.UPF`, etc.

The `input_template` is the template file of the input file of the calculation. The `abacustest` will link the INPUT file in each new generated directory to the `input_template` file.

By executing the following command, the VASP POSCAR files in the current directory will be translated to ABACUS STRU files.
```abacustest prepare -p setting.json```

`abacustest` actually do two following steps:
1. Translate the VASP POSCAR to ABACUS STRU, and put the STRU files to `000000`, `000001`, ... In this step, the STRU has only the structure information, and the pseudopotential is not set.
2. Prepare the ABACUS inputs for each new generated directory, and will generate the inputs in `abacustest/000000`, `abacustest/000001`, ... The `INPUT` file in each directory is linked to the `input_template` file, and the pseudopotential is set in the STRU file.
