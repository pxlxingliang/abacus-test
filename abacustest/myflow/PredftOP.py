import os
from . import globV,comm
from dflow import (
    Step,
    upload_artifact,
    S3Artifact
)

from pathlib import Path
from typing import List, Optional

from dflow.python import (
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    BigParameter
)

from . import comm

class PreDFT(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "examples": Artifact(List[Path],optional=True),
                "command": str,
                "extra_files": Artifact(Path,optional=True),
                "predft_setting" : BigParameter(dict,default={}),  #the setting of predft
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "outputs": Artifact(List[Path]),
                "work_directories": [List[str]]  #a list of a list of work directories in rundft
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        
        print(op_in)
        log = ""
        work_directories = []
        outputs = []
        if not op_in["examples"]:
            all_examples = [os.getcwd()]
            example_root = None
        else:
            all_examples = op_in["examples"]
            example_root = op_in["examples"].art_root
        
        for iexample in all_examples:
            if example_root != None:
                #current_path = str(iexample).split("inputs/artifacts/examples/")[1]
                current_path = os.path.relpath(str(iexample),str(example_root))
                root_path = False
            else:
                current_path = ""
                root_path = True
            #work_directories.append([])
            log += "\nPrepare example %s\n" % current_path
            if op_in["extra_files"] != None:
                #extra_file_path = str(op_in["extra_files"]).split("/inputs/artifacts/extra_files")[0] + "/inputs/artifacts/extra_files"
                extra_file_path = op_in["extra_files"].art_root
                comm.CopyFiles(extra_file_path,iexample,move=False)

            os.chdir(iexample)
            log += "COMMAND: %s\n" % str(op_in["command"])
            if op_in["command"] != None and op_in["command"].strip() != "":
                cmd = "ulimit -c 0; " + str(op_in["command"])
                return_code, out, err = comm.run_command(cmd)
                log += out + err
                #log += os.popen("(%s) 2>&1" % cmd).read()

            work_directories_filename = op_in["predft_setting"].get("work_directories_filename") 
            if not work_directories_filename:
                work_directories_filename = "example.txt"
                print(f"Have not define the work directories filename! Try to find the work directories in file '{work_directories_filename}'")            
            
            if os.path.isfile(work_directories_filename):
                with open(work_directories_filename) as f:
                    lines = f.readlines()
                for i in lines:
                    if not os.path.isdir(i.strip()):
                        print("Can not find work directory %s" % i.strip())
                    else:
                        #work_directories[-1].append(os.path.join(current_path,i.strip()))
                        work_directories.append(os.path.join(current_path,os.path.relpath(i.strip(),os.getcwd())))
                        
                        if root_path:
                            outputs.append(Path(os.path.relpath(i.strip(),os.getcwd())))
                        else:
                            outputs.append(Path.resolve(Path(i.strip())))
            else:
                print(f"Can not find work_directories_filename '{work_directories_filename}'! Return current path: {iexample}")
                if root_path:
                    outputs.append(Path(os.path.relpath(iexample,os.getcwd())))
                else:
                    outputs.append(Path.resolve(Path(iexample)))
                #current_path = str(iexample).split("inputs/artifacts/dflow_examples")[1] #.split("/",1)[1]
                #current_path = "." if "/" not in current_path.strip().strip("/") else current_path.split("/",1)[1]
                #work_directories[-1].append(current_path)
                work_directories.append(current_path)
                
        
        print(outputs)
        print(work_directories)
        return OPIO({
            "outputs": outputs,
            "work_directories": work_directories
        })


def produce_predft(predft_set,stepname,example_path,gather_result=False):
    #predft
    #produce step
    allsteps = []
    allstepname = []
    all_save_path = []
    
    comm.printinfo("\nPreparing pre_dft")
    #assert example_path or predft_set.get("example") != None, "example in predft is not defined"
    if "example" not in predft_set:
        example_list = example_path if example_path else []
        example_source = None
        example_source_type = "local"
    else:
        example_list = predft_set["example"]
        example_source = predft_set.get("example_source")
        example_source_type = predft_set.get("example_source_type", "local")
    examples, examples_name = comm.transfer_source_to_artifact(
        example_list,
        source=example_source,
        source_type=example_source_type,
        only_folder=True,
        oneartifact=False)
    
    extrafiles, extrafiles_name = comm.transfer_source_to_artifact(
        predft_set.get("extra_files", []),
        source=predft_set.get("extra_files_source"),
        source_type=predft_set.get("extra_files_source_type", "local"),
        only_folder=False,
        oneartifact=True)
    
    if not examples:
        comm.printinfo("No examples found!")
    
    executor, bohrium_set = comm.ProduceExecutor(predft_set, group_name=stepname + "/predft")
    image = globV.get_value("ABBREVIATION").get(predft_set.get("image"), predft_set.get("image"))
    model_output_artifact = S3Artifact(key="{{workflow.name}}/predft" )
    
    #ngroup = predft_set.get("ngroup",1)
    #if ngroup != None and ngroup < 1:
    #    ngroup = len(examples)
    group_size = predft_set.get("group_size",len(examples_name))
    if group_size != None and group_size < 1:
        group_size = len(examples_name)
    
    igroup = 0
    istep = 0
    if examples:
        new_examples,new_examples_name = comm.SplitGroupSize(examples,examples_name,group_size)
    else:
        new_examples_name = [[]]
    for iexample_name in new_examples_name:
        igroup += 1
        dflow_stepname = f"predft-{igroup}"
        bohri_stepname = stepname+f"/predft-{igroup}"
        space = "\n" + (len(bohri_stepname)+2)*" "
        comm.printinfo("%s: %s" % (bohri_stepname,space.join(iexample_name)))
        pt = PythonOPTemplate(PreDFT,image=image,envs=comm.SetEnvs())
        artifacts={}
        if iexample_name:
            artifacts["examples"] = upload_artifact(iexample_name,archive=globV.get_value("COMPRESS"))

        #artifacts={"examples": upload_artifact(iexample_name,archive=None)}
        
        #get extrafiles
        if extrafiles:
            artifacts["extra_files"]=extrafiles[0][0]

        step = Step(name=dflow_stepname, template=pt,
                    parameters={
                        "command": predft_set.get("command",""),
                        "predft_setting" : {"work_directories_filename":predft_set.get("work_directories_filename")}
                    },
                    artifacts=artifacts)
        if executor != None:
            step.executor = executor
        if gather_result:
            step.outputs.artifacts["outputs"].save = [model_output_artifact]
            step.template.outputs.artifacts["outputs"].archive = None
        allsteps.append(step)
        allstepname.append(dflow_stepname)
        all_save_path.append("")
    
    '''
    If use slice, the examples can only be an artifact of List[Path], and if use group_size, the group size should divide len(examples)
    We try to support user specified parallel measurement, including serial calculation of specified examples on one machine
    '''
    comm.printinfo("image: %s" % image)
    comm.printinfo("set bohrium: %s" % str(bohrium_set))
    comm.printinfo("command: %s" % str(predft_set.get("command")))
    comm.printinfo("work_directories_filename: %s" % str(predft_set.get("work_directories_filename")))

    #allstepname = [stepname]
    #all_save_path = [[save_path, ""]]
    return allsteps, allstepname, all_save_path


