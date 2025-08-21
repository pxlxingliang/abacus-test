import os,sys,glob
from . import globV,comm
from dflow import (
    Step,
    argo_range,
    S3Artifact,
    argo_len,
    argo_sequence
)

from pathlib import Path
from typing import List

from dflow.python import (
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices
)

from . import comm


class RunDFT(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "examples": Artifact(List[Path]),
                "command": str,
                "sub_save_path": str,
                "extra_files": Artifact(Path,optional=True),
                "outputfiles":[str]
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "outputs": Artifact(List[Path])
            }
        )
                
    def get_cpu_info(self,example_path):
        # get cpuinfo to CPUINFO.log
        log = ""
        try:
            os.system(f"lscpu > {example_path}/CPUINFO.log")
        except:
            pass
        if os.path.isfile(f"{example_path}/CPUINFO.log"):
            with open(f"{example_path}/CPUINFO.log") as f1:
                log = "\n\nCPUINFO:\n" + f1.read()
        return log
                    
    def run_one_job(self,op_in):
        root_path = op_in["examples"].art_root
        print("root_path:",root_path,file=sys.stderr)
        cwd = os.getcwd()
        outpath = []
        for iexample in op_in["examples"]:
            example_path = os.path.relpath(str(iexample),str(root_path))
            if op_in["sub_save_path"] != None and str(op_in["sub_save_path"]).strip() != "":
                example_path = os.path.join(str(op_in["sub_save_path"]),example_path)
                            
            work_path = iexample
            print("work path:",work_path,file=sys.stderr)

            #copy extra_files to work path
            print("sub_save_path:",op_in["sub_save_path"],file=sys.stderr) 
            if op_in["extra_files"] != None:
                extra_file_path = op_in["extra_files"].art_root
                comm.CopyFiles(extra_file_path,work_path,move=False)

            #run command
            os.chdir(work_path)
            log = f"COMMAND: {op_in['command']}"
            if op_in["command"].strip() != "":
                cmd = "ulimit -c 0; " + str(op_in["command"])
                return_code, out, err = comm.run_command(cmd)
                log += out + err
               
            os.chdir(work_path)
            logfile_name = "STDOUTER.log"
            i = 1
            while os.path.isfile(logfile_name):
                logfile_name = "STDOUTER_%d.log" % i
                i += 1
            outfiles = []
            if len(op_in["outputfiles"]) == 0:
                outfiles = os.listdir(".")
            else:
                for i in op_in["outputfiles"]:
                    for j in glob.glob(i):
                        if os.path.exists(j):
                            outfiles.append(j)
                        else:
                            log += "\n%s is not exist" % j
            comm.LinkFiles(outfiles,os.path.join(cwd,example_path))
            
            os.chdir(cwd) 
            log += self.get_cpu_info(example_path)
            logfile = Path(os.path.join(example_path,logfile_name))  
            logfile.write_text(log)
            outpath.append(Path(example_path)) 
        return outpath
    
    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:

        print("op_in:",op_in,file=sys.stderr)
        outpath = self.run_one_job(op_in)
        print("outpath:",str(outpath),file=sys.stderr)
        
        op_out = OPIO(
            {
                "outputs": outpath
            }
        )
        return op_out

def produce_rundft(rundft_sets,predft_step,stepname,example_path,gather_result=False):
    # if predft_step is not None, then the input of rundft is the output of predft,
    # and the example of rundft is invalid, and slices is used as parallel computing scheme.
    # if predft_step is None, then the input of rundft is examples, 
    # and the settings of examples can be read to perform parallel computing flexibly.
    # that is, serial tasks can be placed in the settings of a list.
    #rundft_set is a list of dict, each dict is a rundft setting
    # the returned output_artifact is None if not gather_result, else is model_output_artifact
    comm.printinfo("\nPreparing run_dft")
    allsteps = []
    allstepname = []
    all_save_path = []
    
    rundft_idx = 0
    model_output_artifact = S3Artifact(key="{{workflow.name}}/rundft")
    output_artifact = None
    for rundft_set in rundft_sets:
        rundft_idx += 1
        dflow_stepname = f"rundft-{rundft_idx}"
        bohri_stepname = stepname + f"/rundft-{rundft_idx}"
        sub_savepath = comm.ParseSubSavePath(rundft_set.get("sub_save_path"))
        
        #get extra files
        extrafiles, extrafiles_name = comm.transfer_source_to_artifact(
            rundft_set.get("extra_files", []),
            source=rundft_set.get("extra_files_source"),
            source_type=rundft_set.get("extra_files_source_type", "local"),
            only_folder=False,
            oneartifact=True)
        
        executor, bohrium_set = comm.ProduceExecutor(rundft_set, group_name=bohri_stepname)
        image = globV.get_value("ABBREVIATION").get(rundft_set.get("image"), rundft_set.get("image"))
        group_size = int(rundft_set.get("group_size",1))

        # when use merge_sliced_step, the group_size is None
        if "dispatcher" in rundft_set and rundft_set["dispatcher"].get("merge_sliced_step",False):
            if "group_size" in rundft_set:
                comm.printinfo("group_size is not used when merge_sliced_step is True")
            group_size = None
        parameters = {
            "command": rundft_set.get("command", ""),
            "outputfiles": rundft_set.get("outputs", []),
            "sub_save_path": sub_savepath
        }
        
        if not predft_step:
            # if predft_step is none, then the input of rundft is examples,
            assert example_path or rundft_set.get("example") != None, "example in run_dft is not defined"
            if "example" not in rundft_set:
                example_list = example_path
                example_source = None
                example_source_type = "local"
            else:
                example_list = rundft_set["example"]
                example_source = rundft_set.get("example_source")
                example_source_type = rundft_set.get("example_source_type", "local")
            
            #examples, examples_name = comm.transfer_source_to_artifact(
            #    example_list,
            #    source=example_source,
            #    source_type=example_source_type,
            #    only_folder=True,
            #    oneartifact=False)
            #assert len(examples) > 0, "example in run_dft is not defined or the defined example is not exist!!!"
            '''
            new_examples,new_examples_name = comm.SplitGroupSize(examples,examples_name,group_size)
            istep = 0
            pt = PythonOPTemplate(RunDFT,image=image,envs=comm.SetEnvs())
            for iexample_name in new_examples_name:
                istep += 1
                dflow_stepname_istep = dflow_stepname + f"-{istep}"
                bohri_stepname_istep = bohri_stepname + f"-{istep}"
                space = "\n" + (len(bohri_stepname_istep)+2)*" "
                comm.printinfo("%s: %s" % (bohri_stepname_istep,space.join(iexample_name)))
                artifact_example = upload_artifact(iexample_name,archive=None)
                artifacts={"examples": artifact_example }
                if extrafiles:
                    artifacts["extra_files"]=extrafiles[0][0]

                step = Step(name=dflow_stepname_istep, template=pt,
                            parameters=parameters,
                            artifacts=artifacts,
                            continue_on_failed=True,
                            key=dflow_stepname_istep,
                            )
                if executor != None:
                    step.executor = executor
    
                if gather_result:
                    step.template.outputs.artifacts["outputs"].save = [model_output_artifact]
                    step.template.outputs.artifacts["outputs"].archive = None

                allstepname.append(dflow_stepname_istep)
                all_save_path.append(sub_savepath)
                allsteps.append(step) 
            '''
            examples, examples_name = comm.transfer_source_to_artifact(
                example_list,
                source=example_source,
                source_type=example_source_type,
                only_folder=True,
                oneartifact=True)
            comm.printinfo("Work path:",os.getcwd())
            comm.printinfo("example_name:",examples_name)
            assert len(examples) > 0, "example in run_dft is not defined or the defined example is not exist!!!"
            pt = PythonOPTemplate(RunDFT,image=image,envs=comm.SetEnvs(),
                    slices=Slices(
                        "{{item}}",
                        #sub_path=True,
                        input_artifact = ["examples"],
                        output_artifact = ["outputs"],
                        group_size=group_size,
                        )
                    )
            
            artifacts={"examples": examples[0][0]}
            if extrafiles:
                artifacts["extra_files"]=extrafiles[0][0]
            nexample = len(examples_name[0])
            step = Step(name=dflow_stepname, template=pt,
                        parameters=parameters,
                        artifacts=artifacts,
                        continue_on_failed=True,
                        with_param=argo_range(nexample),
                        key=dflow_stepname+"-{{item}}"
                        )
            if executor != None:
                step.executor = executor
                
            if gather_result:
                step.template.outputs.artifacts["outputs"].save = [model_output_artifact]
                step.template.outputs.artifacts["outputs"].archive = None    

            allstepname.append(dflow_stepname)
            all_save_path.append(sub_savepath)
            allsteps.append(step)  
            
        else:
            artifacts_example = predft_step.outputs.artifacts['outputs']
            #produce step
            pt = PythonOPTemplate(RunDFT,image=image,envs=comm.SetEnvs(),
                    slices=Slices(
                        "int('{{item}}')",
                        input_artifact = ["examples"],
                        output_artifact = ["outputs"],
                        group_size=group_size,
                        )
                    )
            artifacts={"examples": artifacts_example }
            if extrafiles:
                artifacts["extra_files"]=extrafiles[0][0]
            step = Step(name=dflow_stepname, template=pt,
                        parameters=parameters,
                        artifacts=artifacts,
                        continue_on_failed=True,
                        with_sequence=argo_sequence(argo_len(predft_step.outputs.parameters['work_directories']), format='%03d'),
                        key=dflow_stepname+"-{{item}}"
                        )
            if executor != None:
                step.executor = executor
                
            if gather_result:
                step.template.outputs.artifacts["outputs"].save = [model_output_artifact]
                step.template.outputs.artifacts["outputs"].archive = None    

            allstepname.append(dflow_stepname)
            all_save_path.append(sub_savepath)
            allsteps.append(step)
        
        comm.printinfo("image: %s" % image)
        comm.printinfo("set bohrium: %s" % str(bohrium_set))
        comm.printinfo("command: %s" % str(rundft_set.get("command")))   
    if gather_result:
        output_artifact = model_output_artifact
    return allsteps, allstepname, all_save_path, output_artifact


