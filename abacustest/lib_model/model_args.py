"""
MODEL_ARGS - Model registry for abacustest

This dictionary defines all available models in abacustest.
Each model is identified by a unique key (model_name) that is used as the subcommand.

Format of each entry:
    "model_name": {
        "description": "Brief description of the model (shown in help)",
        "file": "Filename of the model module (without .py extension)",
        "class_name": "Name of the Model subclass in the module"
    }

Adding a new model:
1. Create a new file in abacustest/lib_model/ (e.g., model_xxx_MyModel.py)
2. Define a class that inherits from Model (from .model import Model)
3. Implement the required methods:
   - add_args(parser): Add command-line arguments
   - run(params) / run_prepare(params) / run_postprocess(params): Implement functionality
   - Optional: HAS_PREPARE_POST_COMMAND = False if model doesn't have prepare/post subcommands
4. Add an entry to MODEL_ARGS with:
   - Key: The model name (used as subcommand)
   - description: Brief description
   - file: Module filename (without .py)
   - class_name: The class name (must match the class definition)

Note: 
- The model_name in the key must be unique and should be lowercase with hyphens if needed
- The file should follow the naming pattern: model_###_Name.py (e.g., model_022_newmodel.py)
- The class_name must exactly match the class defined in the module
- No need to implement model_name() or description() methods in the class - they are defined here
"""

MODEL_ARGS = {
    "committest": {
        "description": "Prepare the test on different abacus commit",
        "file": "model_000_CommitTest",
        "class_name": "CommitTest"
    },
    "conv": {
        "description": "Do a convergence test",
        "file": "model_001_Conv",
        "class_name": "ConvEcutwfc"
    },
    "aserelax": {
        "description": "Prepare the inputs for relax by ASE + ABACUS",
        "file": "model_003_AseRelax",
        "class_name": "AseRelax"
    },
    "eos": {
        "description": "Prepare and postprocess the EOS calculation",
        "file": "model_004_Eos",
        "class_name": "Eos"
    },
    "phonon": {
        "description": "Prepare and postprocess the phonon calculation",
        "file": "model_005_Phonon",
        "class_name": "Phonon"
    },
    "fdforce": {
        "description": "finite difference of force",
        "file": "model_006_FDForce",
        "class_name": "fdforce"
    },
    "fdstress": {
        "description": "finite difference of stress",
        "file": "model_007_FDStress",
        "class_name": "fdstress"
    },
    "comparem": {
        "description": "Compare two metrics.json files",
        "file": "model_008_CompareMetrics",
        "class_name": "CompareMetricsModel"
    },
    "fdmagforce": {
        "description": "finite difference of magenetic force",
        "file": "model_009_FDMagForce",
        "class_name": "fdmagforce"
    },
    "magj": {
        "description": "Calculate the magnetic exchange interactions of two atom",
        "file": "model_010_MagJ",
        "class_name": "magj"
    },
    "sptest": {
        "description": "Show the results of FP test",
        "file": "model_011_sptest",
        "class_name": "SPTestModel"
    },
    "band": {
        "description": "Prepare and postprocess the calcualtion of band structure",
        "file": "model_012_band",
        "class_name": "BandModel"
    },
    "inputs": {
        "description": "Prepare the ABACUS inputs of specified model",
        "file": "model_013_inputs",
        "class_name": "InputsModel"
    },
    "vasp2abacus": {
        "description": "Transform VASP input files to ABACUS input files.",
        "file": "model_014_vasp2abacus",
        "class_name": "Vasp2AbacusModel"
    },
    "elastic": {
        "description": "Prepare and postprocess the elastic",
        "file": "model_015_elastic",
        "class_name": "ElasticModel"
    },
    "bec": {
        "description": "Calculate the Born effective charge by finite difference method.",
        "file": "model_016_bec",
        "class_name": "BECModel"
    },
    "vacancy": {
        "description": "Calculate the vacancy formation energy for uncharged systems",
        "file": "model_017_vacancy",
        "class_name": "VacancyModel"
    },
    "supercell": {
        "description": "extend the unit cell to supercell",
        "file": "model_018_supercell",
        "class_name": "SuperCellModel"
    },
    "vibration": {
        "description": "Calculate the vibration frequency of selected atoms",
        "file": "model_019_vibration",
        "class_name": "VibrationModel"
    },
    "workfunc": {
        "description": "Prepare and postprocess the work function calculation",
        "file": "model_020_workfunc",
        "class_name": "WorkFuncModel"
    },
    "dos-pdos": {
        "description": "Postprocess the DOS and PDOS data from ABACUS calculations",
        "file": "model_021_dos_pdos",
        "class_name": "DOSPDOSModel"
    }
}