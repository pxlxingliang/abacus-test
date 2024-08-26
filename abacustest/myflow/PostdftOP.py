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

def execute_postdft(work_path,allexample_path,extra_files,extra_files_art_root,tracking_setting,metrics_setting,super_metrics_setting,command,outputfiles):
    """
    should execute this function in abacustest submited path.
    
    
    work_path: str, work path, should be a absolute path
    allexample_path: list, all example path, should be relative path to work_path
    extra_files: should be a list of extra files, and relative path to current path
    
    
    """
    
    #copy extra_files to work path 
    if extra_files:
        for iextra in extra_files:
            comm.CopyFiles(iextra,work_path,move=False)
    if extra_files != None:
        extra_file_path = extra_files_art_root
        comm.CopyFiles(extra_file_path,work_path,move=False)
    #check if need to upload to tracking
    do_upload_tracking = False
    if tracking_setting and tracking_setting.get("ifurn",True):
        do_upload_tracking = True
        
    #read metrics
    metrics_setting_list = metrics.Metrics.TransferMetricsOPIO(metrics_setting)
    # need to split metrics_setting_list to before_command and after_command
    metrics_setting_list_before_command = []
    metrics_setting_list_after_command = []
    for imetric in metrics_setting_list:
        if imetric.get("before_command",True):
            metrics_setting_list_before_command.append(imetric)
        else:
            metrics_setting_list_after_command.append(imetric)
    tracking_values_before = tracking_values_after = None
    if metrics_setting_list_before_command:
        try:
            os.chdir(work_path)
            tracking_values_before = metrics.ReadMetrics(metrics_setting_list_before_command,do_upload_tracking,allexample_path)
        except:
            traceback.print_exc()
    
    #execute command
    os.chdir(work_path)
    log = f"COMMAND: {command}"
    if command.strip() != "":
        cmd = str(command)
        return_code, out, err = comm.run_command(cmd)
        log += out + err
        #log += os.popen("(%s) 2>&1" % cmd).read()
    
    if metrics_setting_list_after_command:
        try:
            os.chdir(work_path)
            tracking_values_after = metrics.ReadMetrics(metrics_setting_list_after_command,do_upload_tracking,allexample_path)
        except:
            traceback.print_exc()              
    # merge tracking_values_before and tracking_values_after
    tracking_values = []
    if tracking_values_before:
        tracking_values.extend(tracking_values_before)
    if tracking_values_after:
        tracking_values.extend(tracking_values_after)
    if not tracking_values:
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
              
    # get cpuinfo to CPUINFO.log
    os.chdir(work_path)
    try:
        os.system("lscpu > CPUINFO.log")
    except:
        pass
    cpuinfo_log = None if not os.path.isfile("CPUINFO.log") else "CPUINFO.log" 
    if cpuinfo_log:
        with open(cpuinfo_log) as f1:
            log += "\n\nCPUINFO:\n" + f1.read() 
                            
    #collect outputs
    outpath = []
    os.chdir(work_path)
    logfile_name = "STDOUTER.log"
    i = 1
    while os.path.isfile(logfile_name):
        logfile_name = "STDOUTER_%d.log" % i
        i += 1
    logfile = Path(logfile_name)
    if len(outputfiles) == 0:
        for ifile in glob.glob("*"):
            if ifile.startswith("."):
                continue 
            outpath.append(Path(ifile))
    else:
        for i in outputfiles:
            for j in glob.glob(i):
                if os.path.exists(j):
                    outpath.append(Path(j))
                else:
                    log += "\n%s is not exist" % j
        outpath.append(logfile)
        if cpuinfo_log:
            outpath.append(Path(cpuinfo_log))
        if os.path.isfile("version.dat"):
            outpath.append(Path("version.dat"))
    #print("log:",log,file=sys.stderr)
    logfile.write_text(log)
    print("outpath:",str(outpath),file=sys.stderr)
    return outpath


class PostDFT(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "examples": Artifact(List[Path]),
                "command": str,
                "extra_files": Artifact(Path,optional=True),
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
        
        # if define sub_save_path, create sub_save_path in root and copy examples to sub_save_path
        #example_path = str(op_in["examples"][0]).split("/inputs/artifacts/examples/")[0] + "/inputs/artifacts/examples/"
        #work_path = example_path
        work_path = op_in["examples"].art_root
        
        allexample_path = []
        for i in op_in["examples"]:
            #allexample_path.append(str(i).split("/inputs/artifacts/examples/")[1])
            allexample_path.append(os.path.relpath(str(i),str(work_path)))

        print("work path:",work_path,file=sys.stderr)
        
        #copy extra_files to work path 
        if op_in["extra_files"] != None:
            #extra_file_path = str(op_in["extra_files"]).split("/inputs/artifacts/extra_files")[0] + "/inputs/artifacts/extra_files"
            extra_file_path = op_in["extra_files"].art_root
            comm.CopyFiles(extra_file_path,work_path,move=False)

        #check if need to upload to tracking
        tracking_setting = op_in["upload_tracking"]
        metrics_setting = op_in["metrics"]
        super_metrics_setting = op_in["super_metrics"]
        do_upload_tracking = False
        if tracking_setting and tracking_setting.get("ifurn",True):
            do_upload_tracking = True
            
        #read metrics
        metrics_setting_list = metrics.Metrics.TransferMetricsOPIO(metrics_setting)
        # need to split metrics_setting_list to before_command and after_command
        metrics_setting_list_before_command = []
        metrics_setting_list_after_command = []
        for imetric in metrics_setting_list:
            if imetric.get("before_command",True):
                metrics_setting_list_before_command.append(imetric)
            else:
                metrics_setting_list_after_command.append(imetric)
        tracking_values_before = tracking_values_after = None
        if metrics_setting_list_before_command:
            try:
                os.chdir(work_path)
                tracking_values_before = metrics.ReadMetrics(metrics_setting_list_before_command,do_upload_tracking,allexample_path)
            except:
                traceback.print_exc()
        
        #execute command
        os.chdir(work_path)
        log = f"COMMAND: {op_in['command']}"
        if op_in["command"].strip() != "":
            cmd = "ulimit -c 0; " + str(op_in["command"])
            return_code, out, err = comm.run_command(cmd)
            log += out + err
            #log += os.popen("(%s) 2>&1" % cmd).read()
        
        if metrics_setting_list_after_command:
            try:
                os.chdir(work_path)
                tracking_values_after = metrics.ReadMetrics(metrics_setting_list_after_command,do_upload_tracking,allexample_path)
            except:
                traceback.print_exc()              

        # merge tracking_values_before and tracking_values_after
        tracking_values = []
        if tracking_values_before:
            tracking_values.extend(tracking_values_before)
        if tracking_values_after:
            tracking_values.extend(tracking_values_after)
        if not tracking_values:
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
                  
        # get cpuinfo to CPUINFO.log
        os.chdir(work_path)
        try:
            os.system("lscpu > CPUINFO.log")
        except:
            pass
        cpuinfo_log = None if not os.path.isfile("CPUINFO.log") else "CPUINFO.log" 
        if cpuinfo_log:
            with open(cpuinfo_log) as f1:
                log += "\n\nCPUINFO:\n" + f1.read() 
                                
        #collect outputs
        os.chdir(work_path)
        logfile_name = "STDOUTER.log"
        i = 1
        while os.path.isfile(logfile_name):
            logfile_name = "STDOUTER_%d.log" % i
            i += 1
        logfile = Path(logfile_name)
        if len(op_in["outputfiles"]) == 0:
            for ifile in glob.glob("*"):
                if ifile.startswith("."):
                    continue 
                outpath.append(Path(ifile))
        else:
            for i in op_in["outputfiles"]:
                for j in glob.glob(i):
                    if os.path.exists(j):
                        outpath.append(Path(j))
                    else:
                        log += "\n%s is not exist" % j
            outpath.append(logfile)
            if cpuinfo_log:
                outpath.append(Path(cpuinfo_log))
            if os.path.isfile("version.dat"):
                outpath.append(Path("version.dat"))
        #print("log:",log,file=sys.stderr)
        logfile.write_text(log)

        print("outpath:",str(outpath),file=sys.stderr)
        op_out = OPIO(
            {
                "outputs": outpath
            }
        )
        return op_out
    
def produce_postdft(setting,prestep_output,flowname,example_path):
    # if rundft_output is not None, then the input of postdft is the output of rundft,
    # and the example of rundft is invalid.
    # if rundft_step is None, then the input of rundft is examples
    # postdft_set is a dict
    comm.printinfo("\nPreparing post_dft")
    allsteps = []
    allstepname = []
    all_sub_save_path = []
    
    dflow_stepname = "postdft"
    bohri_stepname = flowname + f"/postdft"
    
    extrafiles, extrafiles_name = comm.transfer_source_to_artifact(
        setting.get("extra_files", []),
        source=setting.get("extra_files_source"),
        source_type=setting.get("extra_files_source_type", "local"),
        only_folder=False,
        oneartifact=True)
    
    pre_script = ""
    if setting.get("upload_datahub", None):
        pre_script += "import os\nos.system('pip install --upgrade dp-metadata-sdk==3.0.1 -i https://repo.mlops.dp.tech/repository/pypi-group/simple')\n"
    if setting.get("upload_tracking", None):
        pre_script += "import os\nos.system('pip install dp-tracking-sdk==3.15.2.post23 -i https://repo.mlops.dp.tech/repository/pypi-group/simple')\n"
        
    if not prestep_output:
        assert example_path or setting.get("example") != None, "example is not defined in postdft"
        if "example" not in setting:
            example_list = example_path
            example_source = None
            example_source_type = "local"
        else:
            example_list = setting["example"]
            example_source = setting.get("example_source")
            example_source_type = setting.get("example_source_type", "local")
        examples, examples_name = comm.transfer_source_to_artifact(
            example_list,
            source=example_source,
            source_type=example_source_type,
            only_folder=False,
            oneartifact=True)
        if examples:
            artifact_example = examples[0][0]
        else:
            comm.printinfo("No example for post-dft")
            return None,None,None
        space = "\n" + (len(bohri_stepname)+2)*" "
        comm.printinfo("%s: %s" % (bohri_stepname,space.join(examples_name[0])))
        example_path_list = examples_name[0]
    else:
        artifact_example = prestep_output
        example_path_list = example_path     
        
    executor, bohrium_set = comm.ProduceExecutor(setting, group_name=bohri_stepname)
    image = globV.get_value("ABBREVIATION").get(setting.get("image"), setting.get("image"))
    parameters = {
        "command": setting.get("command", ""),
        "outputfiles": setting.get("outputs", []),
        "metrics": setting.get("metrics", []),
        "super_metrics": setting.get("super_metrics", {}),
        "upload_tracking": setting.get("upload_tracking", {})
    }
    
    pt = PythonOPTemplate(PostDFT,image=image,envs=comm.SetEnvs())
    artifacts={"examples": artifact_example }
    if extrafiles:
        artifacts["extra_files"]=extrafiles[0][0]
    step = Step(name=dflow_stepname, template=pt,
                parameters=parameters,
                artifacts=artifacts,
                key=dflow_stepname
                )
    if executor != None:
        step.executor = executor

    allstepname.append(dflow_stepname)
    all_sub_save_path.append("")
    allsteps.append(step)
    comm.printinfo("image: %s" % image)
    comm.printinfo("set bohrium: %s" % str(bohrium_set))
    comm.printinfo("command: %s" % str(setting.get("command")))
    return allsteps, allstepname, all_sub_save_path
