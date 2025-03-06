import os,sys,glob,time,shutil,argparse,json,traceback,copy,re
from . import globV,comm,PredftOP,PostdftOP,RundftOP
from .. import prepare
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
    S3Artifact
)

from pathlib import Path
from typing import List
from dflow.plugins.dispatcher import DispatcherExecutor
from dflow.python import upload_packages

from dflow.plugins.bohrium import BohriumContext, BohriumExecutor
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

from dflow import config, s3_config
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient,create_job_group

from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.collectdata import parse_value

def SetConfig(private_set,debug=False):  
    config["archive_mode"] = None
    if debug:
        config["mode"] = "debug"
        host = "LOCAL"
        client = "LOCAL"
    else:
        comm.printinfo("set config info...")
        config["save_keys_in_global_outputs"] = False
        if private_set.get("dflow_host","").strip() != "":
            config["host"] = private_set.get("dflow_host","").strip()
        else:
            config["host"] = "https://workflows.deepmodeling.com"
    
        host = config["host"] + "/workflows/argo"
        
        if private_set.get("dflow_s3_config_endpoint","").strip() != "": 
            s3_config["endpoint"] =  private_set.get("dflow_s3_config_endpoint","").strip()
            comm.printinfo("set s3_config_endpoint ... ")
        
        comm.printinfo("set k8s_api_server ... ")    
        if private_set.get("dflow_k8s_api_server","").strip() != "":
            config["k8s_api_server"] = private_set.get("dflow_k8s_api_server","").strip()
        else:
            config["k8s_api_server"] = "https://workflows.deepmodeling.com"
            
        if private_set.get("dflow_token","").strip() != "":
            config["token"] = private_set.get("dflow_token","").strip()
            comm.printinfo("set token ... ")
            
        if host == "https://lbg-workflow-deepmodeling.dp.tech/workflows/argo":
            config["token"] = ""
            config["namespace"] = "deepmodeling"
            host = "https://lbg-workflow-deepmodeling.dp.tech/workflows/deepmodeling"

        if "bohrium_username" in private_set:
            bohrium.config["username"] = private_set.get('bohrium_username')
            comm.printinfo("set bohrium_username ... ")
            
        if "bohrium_password" in private_set:
            bohrium.config["password"] = private_set.get('bohrium_password','')
            comm.printinfo("set bohrium_password ... ")
        
        if "bohrium_ticket" in private_set:
            bohrium.config["ticket"] = private_set.get('bohrium_ticket')
            comm.printinfo("set bohrium_ticket ... ")
                
        if "bohrium_project_id" in private_set:
            bohrium.config["project_id"] = private_set.get('bohrium_project_id','')
            comm.printinfo("set bohrium_project_id ... ")

        s3_config["repo_key"] = "oss-bohrium"
        s3_config["storage_client"] = TiefblueClient()
        client = s3_config["storage_client"]

        if globV.get_value("BOHRIUM_EXECUTOR"):
            globV.set_value("BRM_CONTEXT",BohriumContext(
                    executor="mixed",
                    extra={},
                    username=bohrium.config["username"],
                    password=bohrium.config["password"]))

        
    globV.set_value("HOST", host)
    globV.set_value("storage_client", client)

def PrepareExample(setting):
    '''
    {
        "ifrun": True,
        "example_template":["example_path"]
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
        "pp_path": str,
        "orb_path": str,
        "dpks_descriptor":"",
        "extra_files":[],
    }
    '''
    allpath = []
    for ipathlist in prepare.DoPrepare(setting,""):  #put prepared examples to current path
        for ipath in ipathlist:
            allpath.append(ipath)
    return allpath
    
def ParseProcess(param):
    """
    param:{
        "prepare":{},
        "pre_dft":{},
        "run_dft":[{}],
        "post_dft":{}
    }
    """
    prepare_example = pre_dft = run_dft = post_dft = True
    if "prepare" not in param or (not param["prepare"].get("ifrun",True)):
        prepare_example = False
    else:
        prepare_example = param["prepare"]
    
    if "pre_dft" not in param or ( not param["pre_dft"].get("ifrun",True)):
        pre_dft = None
    else:
        pre_dft = param["pre_dft"]
        
    if "run_dft" not in param:
        run_dft = None
    else:
        tmp = []
        run_dft = []
        if isinstance(param["run_dft"],dict):
            tmp = [param["run_dft"]]
        elif isinstance(param["run_dft"],list):
            tmp = param["run_dft"]
        else:
            comm.printinfo("run_dft should be dict or list")
            
        for itmp in tmp:
            if isinstance(itmp,dict):
                if itmp.get("ifrun",True):
                    run_dft.append(itmp)
        if len(run_dft) == 0:
            run_dft = None
        
    if "post_dft" not in param or (not param["post_dft"].get("ifrun",True)):
        post_dft = None
    else:
        post_dft = param["post_dft"]
        
    return prepare_example,pre_dft,run_dft,post_dft 

def ProduceOneSteps(stepname,param):
    prepare_example,pre_dft,run_dft,post_dft = ParseProcess(param)
    save_path = comm.ParseSavePath(param.get("save_path",None))
    steps = Steps(name="abacustesting-steps",
                  outputs=Outputs(artifacts={"outputs" : OutputArtifact()}))
    allstepname, all_save_path = [], []
    has_step = False
    final_step = None
    
    comm.printinfo(f"\n{stepname}")
    comm.printinfo(f"prepare, pre_dft, run_dft, post_dft: {bool(prepare_example)}, {bool(pre_dft)}, {bool(run_dft)}, {bool(post_dft)}")
    
    # if define prepare, then will prepare example firstly
    if prepare_example:
        #return of Doprepare is a list of dict, each dict is an example
        # each element of the dict is a sub example
        example_path = []
        for ii in prepare.DoPrepare(prepare_example,""):
            if ii != None:
                for jj in ii:
                   if prepare.CheckExample(jj,str(ii[jj])):
                       example_path.append(jj)  
        if example_path == []:
            comm.printinfo("WARNING: defined prepare, but no examples matched, skip it!")
            example_path = None 
        elif (not pre_dft) and (not run_dft) and (not post_dft):
            # copy the example to save_path
            for iexample in example_path:
                relative_path = os.path.relpath(iexample,os.getcwd())
                if os.path.isdir(os.path.join(save_path,relative_path)):
                    shutil.rmtree(os.path.join(save_path,relative_path))
                shutil.copytree(iexample,os.path.join(save_path,relative_path))
        #elif post_dft:
        #    # check the path of metrics, if no path, then use example_path as input
        #    if "metrics" in post_dft:
        #        if isinstance(post_dft["metrics"],dict):
        #            if "path" not in post_dft["metrics"]:
        #                post_dft["metrics"]["path"] = example_path
        #        elif isinstance(post_dft["metrics"],list):
        #            for imetric in range(len(post_dft["metrics"])):
        #                if isinstance(post_dft["metrics"][imetric],dict):
        #                    if "path" not in post_dft["metrics"][imetric]:
        #                        post_dft["metrics"][imetric]["path"] = example_path
    else:
        example_path = None  
    
    # in the predft/rundft/postdft step, in not define example, then will use example_path as input
    # if define example, then will use example as input
    predft_step = None
    if pre_dft:
        if not run_dft and not post_dft:
            gather_result = False
        else:
            gather_result = True
            
        step0, allstepname, all_save_path = PredftOP.produce_predft(pre_dft,stepname,example_path,gather_result)
        if step0:
            steps.add(step0)
            has_step = True
            final_step = step0[0]
            predft_step = step0[0]
    
    # post dft need the output of run dft
    # and if has post dft, the output of run dft will be saved in model_output_artifact
    if run_dft:
        if post_dft:
            gather_result = True
        else:
            gather_result = False

        rundft_step,allstepname,all_save_path,output_artifact = RundftOP.produce_rundft(run_dft,predft_step,stepname,example_path,gather_result)
            
        if len(rundft_step) > 0:
            steps.add(rundft_step)
            has_step = True
            final_step = rundft_step[0]
        else:
            comm.printinfo("WARNING: has defined run_dft, but no examples matched, skip it!")
            output_artifact = None
            
    else:
        output_artifact = None
    
    # if output_artifact is not None, it means that the output of run dft will be saved in model_output_artifact
    # and if output_artifact is None, post dft will use the examples defined by example
    if post_dft:
        #postdft_step, allstepname, all_save_path = produce_postdft(post_dft,stepname,output_artifact,save_path)
        postdft_step, allstepname, all_save_path = PostdftOP.produce_postdft(post_dft,output_artifact,stepname,example_path)
        if postdft_step is not None:
            steps.add(postdft_step)
            has_step = True
            final_step = postdft_step[0]

    if final_step:
        steps.outputs.artifacts['outputs']._from = final_step.outputs.artifacts["outputs"]
    
    step = None
    if has_step:
        step = Step(name="abacustesting",template=steps,key="abacustesting") 
    if all_save_path != None:
        all_save_path = [[save_path,i] for i in all_save_path]
    return step,allstepname,all_save_path

def ProduceAllStep(alljobs):
    allstep = []
    allstepname = []
    allsave_path = []

    stepname = alljobs.get("bohrium_group_name","abacustesting")
    step,stepnames,save_path = ProduceOneSteps(stepname,alljobs)
    if step != None:  
        allstep.append(step)
        allstepname.append(stepnames)
        allsave_path.append(save_path) 
        comm.printinfo("\nComplete the preparing for %s.\n" % (stepname))

    return allstep,allstepname,allsave_path



