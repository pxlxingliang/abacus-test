import os,sys,glob,time,shutil,argparse,json,traceback,copy,re
from . import globV,comm
from dflow import (
    Workflow,
    Step,
    Steps,
    Inputs,
    Outputs,
    argo_range,
    SlurmRemoteExecutor,
    upload_artifact,
    download_artifact,
    InputArtifact,
    InputParameter,
    OutputArtifact,
    OutputParameter,
    ShellOPTemplate,
    S3Artifact,
    argo_len,
    argo_sequence,
)

from pathlib import Path
from typing import List, Optional

from dflow.python import (
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices,
    BigParameter,
    Parameter
)

from . import comm,metrics,tracking


class RunDFT(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "examples": Artifact(List[Path]),
                #"examples_name":str,
                "command": str,
                "sub_save_path": str,
                "extra_files": Artifact(Path),
                "outputfiles":[str],
                "metrics": BigParameter(dict,default={}),
                "super_metrics": BigParameter(dict,default={}),
                "upload_tracking": BigParameter(dict,default={}),
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "outputs": Artifact(List[Path])
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:

        print("op_in:",op_in,file=sys.stderr)
        outpath = []
        cwd = os.getcwd()
        
        # if define sub_save_path, create sub_save_path in root and copy examples to work path
        # else work path is examples
        for iexample in op_in["examples"]:
            example_path = str(iexample).split("/inputs/artifacts/examples/")[1]
            work_path = example_path
            if op_in["sub_save_path"] != None and str(op_in["sub_save_path"]).strip() != "":
                work_path = os.path.join(str(op_in["sub_save_path"]),example_path)
            os.makedirs(work_path,exist_ok=True)
            comm.CopyFiles(str(iexample),work_path,move=False)
            work_path_name = work_path
            work_path = os.path.abspath(work_path)
            print("work path:",work_path,file=sys.stderr)

            #copy extra_files to work path
            print("sub_save_path:",op_in["extra_files"],file=sys.stderr) 
            if op_in["extra_files"] != None:
                extra_file_path = str(op_in["extra_files"]).split("/inputs/artifacts/extra_files")[0] + "/inputs/artifacts/extra_files"
                comm.CopyFiles(extra_file_path,work_path,move=False)

            #run command
            os.chdir(work_path)
            log = ""
            if op_in["command"].strip() != "":
                cmd = str(op_in["command"])
                log += os.popen("(%s) 2>&1" % cmd).read()

            #check if need to upload to tracking
            tracking_setting = op_in["upload_tracking"]
            metrics_setting = op_in["metrics"]
            super_metrics_setting = op_in["super_metrics"]
            do_upload_tracking = False
            if tracking_setting and tracking_setting.get("ifurn",True):
                do_upload_tracking = True

            #read metrics
            try:
                os.chdir(work_path)
                tracking_values = metrics.ReadMetrics(metrics.Metrics.TransferMetricsOPIO(metrics_setting),do_upload_tracking,["."])
            except:
                traceback.print_exc()
                tracking_values = None
                
            #calculate super_metrics
            try:
                os.chdir(work_path)
                tracking_summary,report = metrics.ReadSuperMetrics(metrics.Metrics.TransferMetricsOPIO(super_metrics_setting),do_upload_tracking)  
            except:
                traceback.print_exc()
                tracking_summary = None  

            #upload tracking
            if do_upload_tracking:
                try:
                    tracking.upload_to_tracking(tracking_setting,tracking_values,tracking_summary,AIM_ACCESS_TOKEN=None)
                except:
                    traceback.print_exc()   

            #collect outputs
            os.chdir(work_path)
            logfile_name = "STDOUTER.log"
            i = 1
            while os.path.isfile(logfile_name):
                logfile_name = "STDOUTER_%d.log" % i
                i += 1
            
            os.chdir(cwd)
            logfile = Path(os.path.join(work_path_name,logfile_name))
            if len(op_in["outputfiles"]) == 0:
                outpath.append(Path(work_path_name))
            else:
                for i in op_in["outputfiles"]:
                    for j in glob.glob(os.path.join(work_path_name,i)):
                        if os.path.exists(j):
                            outpath.append(Path(j))
                        else:
                            log += "\n%s is not exist" % j
                outpath.append(logfile)
            print("log:",log,file=sys.stderr)
            logfile.write_text(log)
            #os.chdir(cwd)

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
    model_output_artifact = S3Artifact(key="{{workflow.name}}/%s/rundft" % stepname)
    output_artifact = None
    for rundft_set in rundft_sets:
        rundft_idx += 1
        rundft_stepname = stepname + f"-rundft-group{rundft_idx}"
        sub_savepath = comm.ParseSubSavePath(rundft_set.get("sub_save_path"))
        
        #get extra files
        extrafiles, extrafiles_name = comm.transfer_source_to_artifact(
            rundft_set.get("extra_files", []),
            source=rundft_set.get("extra_files_source"),
            source_type=rundft_set.get("extra_files_source_type", "local"),
            only_folder=False,
            oneartifact=True)
        
        executor, bohrium_set = comm.ProduceExecutor(rundft_set, group_name=rundft_stepname)
        image = globV.get_value("ABBREVIATION").get(rundft_set.get("image"), rundft_set.get("image"))
        group_size = int(rundft_set.get("group_size",1))
        parameters = {
            "command": rundft_set.get("command", ""),
            "outputfiles": rundft_set.get("outputs", []),
            "sub_save_path": sub_savepath,
            "metrics": rundft_set.get("metrics", []),
            "super_metrics": rundft_set.get("super_metrics", {}),
            "upload_tracking": rundft_set.get("upload_tracking", {}),
        }
        
        if not predft_step:
            # if predft_step is none, then the input of rundft is examples,
            assert example_path or rundft_set.get("example") != None, "example in rundft is not defined"
            if "example" not in rundft_set:
                example_list = example_path
                example_source = None
                example_source_type = "local"
            else:
                example_list = rundft_set["example"]
                example_source = rundft_set.get("example_source")
                example_source_type = rundft_set.get("example_source_type", "local")
            examples, examples_name = comm.transfer_source_to_artifact(
                example_list,
                source=example_source,
                source_type=example_source_type,
                only_folder=True,
                oneartifact=False)
            new_examples,new_examples_name = comm.SplitGroupSize(examples,examples_name,group_size)
            istep = 0
            for iexample_name in new_examples_name:
                istep += 1
                rundft_stepname_istep = rundft_stepname + f"-{istep}"
                space = "\n" + (len(rundft_stepname_istep)+2)*" "
                comm.printinfo("%s: %s" % (rundft_stepname_istep,space.join(iexample_name)))
                artifact_example = upload_artifact(iexample_name,archive=None)
                pt = PythonOPTemplate(RunDFT,image=image)
                artifacts={"examples": artifact_example }
                if extrafiles:
                    artifacts["extra_files"]=extrafiles[0][0]
                else:
                    pt.inputs.artifacts["extra_files"].optional = True
                step = Step(name=rundft_stepname_istep, template=pt,
                            parameters=parameters,
                            artifacts=artifacts,
                            key=rundft_stepname_istep,
                            )
                if executor != None:
                    step.executor = executor
    
                if gather_result:
                    step.template.outputs.artifacts["outputs"].save = [model_output_artifact]
                    step.template.outputs.artifacts["outputs"].archive = None

                allstepname.append(rundft_stepname_istep)
                all_save_path.append(sub_savepath)
                allsteps.append(step)   
        else:
            artifacts_example = predft_step.outputs.artifacts['outputs']
            #produce step
            pt = PythonOPTemplate(RunDFT,image=image,slices=Slices(
                        "int('{{item}}')",
                        input_artifact = ["examples"],
                        output_artifact = ["outputs"],
                        group_size=group_size,
                        ))
            artifacts={"examples": artifacts_example }
            if extrafiles:
                artifacts["extra_files"]=extrafiles[0][0]
            else:
                pt.inputs.artifacts["extra_files"].optional = True
            step = Step(name=rundft_stepname, template=pt,
                        parameters=parameters,
                        artifacts=artifacts,
                        with_sequence=argo_sequence(argo_len(predft_step.outputs.parameters['work_directories']), format='%03d'),
                        key=rundft_stepname+"-{{item}}"
                        )
            if executor != None:
                step.executor = executor
                
            if gather_result:
                step.template.outputs.artifacts["outputs"].save = [model_output_artifact]
                step.template.outputs.artifacts["outputs"].archive = None    

            allstepname.append(rundft_stepname)
            all_save_path.append(sub_savepath)
            allsteps.append(step)
        
        comm.printinfo("image: %s" % image)
        comm.printinfo("set bohrium: %s" % str(bohrium_set))
        comm.printinfo("command: %s" % str(rundft_set.get("command")))   
    if gather_result:
        output_artifact = model_output_artifact
    return allsteps, allstepname, all_save_path, output_artifact

def produce_step(setting,
                 flowname,
                 stepname,
                 prestep_artifact,
                 doslice_parameter,
                 default_group_size=1,
                 example_only_folder=True,
                 example_oneartifact=False,
                 gather_result=False):
    '''
    {
        "examples",
        "group_size",
        "image",
        "command",
        "extra_files",
        "outputs",
        "metrics",
        "super_metrics",
        "upload_tracking",
    }
    '''
    # if pre_step is not None, then the input is the output of predft,
    # and the example invalid, and slices is used as parallel computing scheme.
    # if pr_step is None, then the input of rundft is examples, 
    # and the settings of examples can be read to perform parallel computing flexibly.
    # that is, serial tasks can be placed in the settings of a list.
    #rundft_set is a list of dict, each dict is a rundft setting
    # the returned output_artifact is None if not gather_result, else is model_output_artifact
    comm.printinfo(f"\nPreparing {stepname}")
    
    #replace _ to - in flowname and stepname
    flowname = flowname.replace("_","-")
    stepname = stepname.replace("_","-")
    
    allsteps = []
    allstepname = []
    all_sub_save_path = []
    
    idx = 0
    model_output_artifact = S3Artifact(key="{{workflow.name}}/%s/%s" % (flowname,stepname))
    output_artifact = None
    if isinstance(setting,dict):
        setting_tmp = [setting]
    else:
        setting_tmp = setting
    for iset in setting_tmp:
        idx += 1
        group_name = flowname + f"-{stepname}"
        if isinstance(setting,list):
            group_name = group_name + f"-group{idx}"
        sub_savepath = comm.ParseSubSavePath(iset.get("sub_save_path"))
        
        #get extrafiles
        extrafiles, extrafiles_name = comm.transfer_source_to_artifact(
            iset.get("extra_files", []),
            source=iset.get("extra_files_source"),
            source_type=iset.get("extra_files_source_type", "local"),
            only_folder=False,
            oneartifact=True)
        
        executor, bohrium_set = comm.ProduceExecutor(iset, group_name=group_name)
        image = globV.get_value("ABBREVIATION").get(iset.get("image"), iset.get("image"))
        group_size = int(iset.get("group_size",default_group_size))
        parameters = {
            "command": iset.get("command", ""),
            "outputfiles": iset.get("outputs", []),
            "sub_save_path": sub_savepath,
            "metrics": iset.get("metrics", []),
            "super_metrics": iset.get("super_metrics", {}),
            "upload_tracking": iset.get("upload_tracking", {}),
        }
        
        if prestep_artifact == None:
            #get examples
            # if predft_step is none, then the input of rundft is examples,
            assert iset.get("example") != None, "example in rundft is not defined"
            examples, examples_name = comm.transfer_source_to_artifact(
                iset["example"],
                source=iset.get("example_source"),
                source_type=iset.get("example_source_type", "local"),
                only_folder=example_only_folder,
                oneartifact=example_oneartifact)
            
            new_examples,new_examples_name = comm.SplitGroupSize(examples,examples_name,group_size)
            istep = 0
            for iexample_name in new_examples_name:
                istep += 1
                group_name_istep = group_name + f"-{istep}"
                space = "\n" + (len(group_name_istep)+2)*" "
                comm.printinfo("%s: %s" % (group_name_istep,space.join(iexample_name)))
                artifact_example = upload_artifact(iexample_name,archive=None)
                pt = PythonOPTemplate(RunDFT,image=image)
                artifacts={"examples": artifact_example }
                if extrafiles:
                    artifacts["extra_files"]=extrafiles[0][0]
                else:
                    pt.inputs.artifacts["extra_files"].optional = True
                parameters_tmp = parameters.copy()
                #parameters_tmp["examples_name"] = "{{item}}"
                step = Step(name=group_name_istep, template=pt,
                            parameters=parameters_tmp,
                            artifacts=artifacts,
                            key=group_name_istep,
                            )
                if executor != None:
                    step.executor = executor

                if gather_result:
                    step.template.outputs.artifacts["outputs"].save = [model_output_artifact]
                    step.template.outputs.artifacts["outputs"].archive = None
                allstepname.append(group_name_istep)
                all_sub_save_path.append(sub_savepath)
                allsteps.append(step)  
            
        elif doslice_parameter:
            #produce step
            pt = PythonOPTemplate(RunDFT,image=image,slices=Slices(
                        "int('{{item}}')",
                        input_artifact = ["examples"],
                        output_artifact = ["outputs"],
                        group_size=group_size
                        ))
            artifacts={"examples": prestep_artifact}
            if extrafiles:
                artifacts["extra_files"]=extrafiles[0]
            else:
                pt.inputs.artifacts["extra_files"].optional = True
            step = Step(name=group_name, template=pt,
                        parameters=parameters,
                        artifacts=artifacts,
                        with_sequence=argo_sequence(argo_len(doslice_parameter), format='%03d'),
                        key=group_name+"-{{item}}"
                        )
            if executor != None:
                step.executor = executor
                
            if gather_result:
                step.template.outputs.artifacts["outputs"].save = [model_output_artifact]
                step.template.outputs.artifacts["outputs"].archive = None    

            allstepname.append(group_name)
            all_sub_save_path.append(sub_savepath)
            allsteps.append(step)
        
        else:
            #produce step
            pt = PythonOPTemplate(RunDFT,image=image)
            artifacts={"examples": prestep_artifact}
            if extrafiles:
                artifacts["extra_files"]=extrafiles[0]
            else:
                pt.inputs.artifacts["extra_files"].optional = True
            step = Step(name=group_name, template=pt,
                        parameters=parameters,
                        artifacts=artifacts,
                        key=group_name
                        )
            if executor != None:
                step.executor = executor
                
            if gather_result:
                step.template.outputs.artifacts["outputs"].save = [model_output_artifact]
                step.template.outputs.artifacts["outputs"].archive = None    

            allstepname.append(group_name)
            all_sub_save_path.append(sub_savepath)
            allsteps.append(step)
            
        comm.printinfo(f"{stepname} group{idx}")
        comm.printinfo("image: %s" % image)
        comm.printinfo("set bohrium: %s" % str(bohrium_set))
        comm.printinfo("command: %s" % str(iset.get("command"))) 
            
    if gather_result:
        output_artifact = model_output_artifact
        
    return allsteps, allstepname, all_sub_save_path, output_artifact
