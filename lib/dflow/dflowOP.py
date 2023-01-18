#!/usr/bin/env python3
import os,sys,glob,time,shutil,argparse,json
from . import globV
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
import subprocess, os, shutil, glob
from pathlib import Path
from typing import List
from dflow import config, s3_config
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient
from dflow.plugins.dispatcher import DispatcherExecutor

#from dflow.plugins.bohrium import BohriumContext, BohriumExecutor
from dflow.python import (
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices
)

class RunDFT(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "example": Artifact(Path),
                "command": str,
                "sub_path": [str],
                "collectdata_script": Artifact(Path),
                "collectdata_pythonlib": Artifact(Path),
                "collectdata_pythonlib_folders":[str],
                "collectdata_command": str,
                "outputfiles":[str]
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "outputs": Artifact(List[Path]),
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:

        def GetBakFile(sfile):
            n = 1
            bk = sfile + ".bak%d" % n
            while os.path.exists(bk):
                n += 1
                bk = sfile + ".bak%d" % n
            return bk

        outpath = []
        for sub_path in op_in["sub_path"]:
            work_path = os.path.join(op_in["example"],sub_path)
            print("work path:",work_path)
            log = ""
            os.chdir(work_path)
            if op_in["command"].strip() != "":
                log += "\n".join([i for i in os.popen(op_in["command"])])
        
            os.chdir(work_path)
            script_folder = "SCRIPT"
            if op_in["collectdata_script"] != None:
                if os.path.isdir(script_folder):
                    os.rename(script_folder,GetBakFile(script_folder))
                script_source_path = str(op_in["collectdata_script"]).split("collectdata_script")[0] + "collectdata_script"
                shutil.copytree(script_source_path,script_folder)

            cmd = ''
            if op_in["collectdata_pythonlib"] != None:
                if len(op_in["collectdata_pythonlib_folders"]) == 1:
                    tpath = os.path.split(str(op_in["collectdata_pythonlib"]).rstrip("/"))[0]
                    #cmd += "export PYTHONPATH=%s:$PYTHONPATH && " % tpath
                    cmd += "export PATH=%s:$PATH && " % tpath
                else:
                    for i in op_in["collectdata_pythonlib_folders"]: 
                        tpath = os.path.join(op_in["collectdata_pythonlib"],os.path.split(i.strip("/"))[0])
                        cmd += "export PYTHONPATH=%s:$PYTHONPATH && " % tpath

            if op_in["collectdata_command"] != "":
                cmd += op_in["collectdata_command"]
                log += "\ncommand:" + cmd + "\n"
                log += "\n".join([i for i in os.popen(cmd)])
            
            os.chdir(work_path)
            logfile_name = "STDOUTER"
            if os.path.isfile(logfile_name):
                logfile_name = GetBakFile(logfile_name)
            logfile = Path(logfile_name)

            if len(op_in["outputfiles"]) == 0:
                outpath.append(Path(work_path))
            else:
                for i in op_in["outputfiles"]:
                    for j in glob.glob(i):
                        if os.path.exists(j):
                            outpath.append(Path(os.path.join(os.getcwd(),j)))
                        else:
                            log += "\n%s is not exist" % j
                outpath.append(Path.resolve(logfile))

            print(log)
            logfile.write_text(log)

        print(str(outpath))
        op_out = OPIO(
            {
                "outputs": outpath 
            }
        )
        return op_out

class PostDFT(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "script": Artifact(Path),
                "data": Artifact(Path),
                "command": str,
                "outputfiles":[str]
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "outputs": Artifact(List[Path]),
                "log": Artifact(Path),
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
       
        outt = ''
        cwd = os.getcwd()
        
        input_data = str(op_in["data"])
        data_path = input_data.split('data')[0] + 'data'
        alldata = os.listdir(data_path)

        if op_in["script"] != None:
            input_script = str(op_in["script"])
            script_path = input_script.split('script')[0] + 'script'
            allscript = os.listdir(script_path)
        else:
            script_path = cwd
            allscript = []

        #move all scripts into datapath:
        for i in allscript:
            if i[0] == '.': continue
            shutil.move(os.path.join(script_path,i),data_path)
        
        #cd data path and run command
        os.chdir(data_path)
        if op_in["command"].strip() != "":
            outt += "execute command:\n%s\n" % op_in["command"]
            outt += "\n".join([i for i in os.popen(op_in["command"])])
        
        os.chdir(data_path)
        for i in allscript:
            if os.path.isdir(i):
                shutil.rmtree(i)
            elif os.path.isfile(i):
                os.unlink(i)
        
        logfile = "STDOUTER"
        n = 0
        while os.path.isfile(logfile):
            n += 1
            logfile = "STDOUTER%d" % n
        logfile = Path(logfile)
        logfile.write_text(outt)

        output = []
        if len(op_in["outputfiles"]) == 0:
            for i in os.listdir("."):
                output.append(Path.resolve(Path(i)))
        else:
            output.append(Path.resolve(logfile))
            for i in op_in["outputfiles"]:
                for j in glob.glob(i):
                    output.append(Path.resolve(Path(j)))

        op_out = OPIO(
            {
                "outputs": output,
                "log": logfile
            }
        )
        return op_out

def ProduceExecutor(param):
    if "bohrium" in param and param["bohrium"]:
        bohrium_set = {}
        for key in param["bohrium"]:
            if key == 'scassType':
                bohrium_set['scass_type'] = param["bohrium"][key]
            elif key == 'jobType':
                bohrium_set['job_type'] = param["bohrium"][key]
            else:
                bohrium_set[key] = param["bohrium"][key]

        if 'platform' not in bohrium_set:
            bohrium_set['platform'] = 'ali'

        dispatcher_executor = DispatcherExecutor(
            machine_dict={
                "batch_type": "Bohrium",
                "context_type": "Bohrium",
                "remote_profile": {"input_data": bohrium_set},
                },
        )
        printinfo("set bohrium: %s"%str(bohrium_set))
        return dispatcher_executor
    else:
        return None

def ProduceRunDFTStep(examples,param,name,DoSlices=True):
    collectdata_scripts = []
    for i in param.get("collectdata_script",[]):
        collectdata_scripts += glob.glob(i)
    
    '''
    collectdata_pythonlib = []
    for i in param.get("collectdata_pythonlib",[]):
        collectdata_pythonlib += glob.glob(i)
    '''
    abacustestpath = os.path.join(*(__file__.split('/')[:-3]))
    collectdata_pythonlib = [abacustestpath]
    printinfo("image: %s" % param['image'])
    if DoSlices:
        pt = PythonOPTemplate(RunDFT,image=param['image'],
                    slices=Slices(sub_path = True,
                                  input_artifact=["example"],
                                  output_artifact=["outputs"]))
        sub_path = [""]
    else:
        pt = PythonOPTemplate(RunDFT,image=param['image'])
        if len(examples) == 1:
            sub_path = [""]
        else:
            sub_path = list(examples)

    artifacts = {"example": upload_artifact(examples,archive=None)}
    if len(collectdata_scripts) > 0:
        artifacts["collectdata_script"] = upload_artifact(collectdata_scripts)
    else:
        pt.inputs.artifacts["collectdata_script"].optional = True
    if len(collectdata_pythonlib) > 0:
        artifacts["collectdata_pythonlib"] = upload_artifact(collectdata_pythonlib)
        if len(collectdata_pythonlib) == 1:
            collectdata_pythonlib_folders = [""]
        else:
            collectdata_pythonlib_folders = collectdata_pythonlib
    else:
        pt.inputs.artifacts["collectdata_pythonlib"].optional = True
        collectdata_pythonlib_folders = [""]

    step = Step(name=name,template=pt,
            parameters = {"command": param['command'], "sub_path":sub_path,
                          "collectdata_command":param.get("collectdata_command",""),
                          "collectdata_pythonlib_folders": collectdata_pythonlib_folders,
                          "outputfiles":param.get("outputs",[])},
            artifacts =  artifacts,continue_on_failed=True)

    executor = ProduceExecutor(param)
    if executor != None:
        step.executor = executor

    return step

def ProducePostDFTStep(param,data):
    '''
    param has key:
        image (str): container image name
        script (list): the scripts
        command (str): the command to run the scripts    
    '''
    pt = PythonOPTemplate(PostDFT,image=param['image'])
    
    artifacts = {"data": data}
    allscript = []
    for i in param.get('script',[]):
        allscript += glob.glob(i)
    if len(allscript) > 0:
        artifacts["script"] = upload_artifact(allscript)
    else:
        pt.inputs.artifacts["script"].optional = True

    step = Step(name="post-dft",template=pt,
            parameters = {"command": param.get('command',""),
                          "outputfiles":param.get("outputs",[])},
            artifacts = artifacts,continue_on_failed=True)

    executor = ProduceExecutor(param)
    if executor != None: step.executor = executor
    return step

def ProduceOneSteps(stepname,param):
    if "run_dft" not in param:
        return None
    rundft_set = param["run_dft"]

    steps = Steps(name=stepname+"-steps",
                  outputs=Outputs(artifacts={"outputs" : OutputArtifact()}))
    #rundft
    model_output_artifact = S3Artifact(key="{{workflow.name}}/%s"%stepname)
    def step1(example_tmp,rundft,istep,DoSlices):
        example_tmp.sort()
        print("\nrun-dft-%d:\n    %s" % (istep,"\n    ".join(example_tmp)))
        step1_tmp = ProduceRunDFTStep(example_tmp,rundft,"run-dft-%d"%istep,DoSlices=DoSlices)
        step1_tmp.template.outputs.artifacts["outputs"].save = [model_output_artifact]
        step1_tmp.template.outputs.artifacts["outputs"].archive = None
        return step1_tmp
    rundft_step = []
    istep = 0
    print("\n%s" % stepname)
    for rundft in rundft_set:
        if 'ifrun' in rundft and not rundft['ifrun']: continue
        examples = []
        for i in rundft['example']:
            if isinstance(i,(list,tuple)):
                example_tmp = []
                for j in i:
                    example_tmp += glob.glob(j)
                if len(example_tmp) > 0:
                    rundft_step.append(step1(example_tmp,rundft,istep,False))
                    istep += 1
            elif isinstance(i,str):
                examples += glob.glob(i)
            else:
                print(i,"element of 'example' should be a list, tuple or str")

        if len(examples) > 0:
            rundft_step.append(step1(examples,rundft,istep,True))
            istep += 1
    if len(rundft_step) == 0:
        print("No examples matched in %s, skip it!" % stepname)
        return None
    steps.add(rundft_step)
    
    #post dft
    if "post_dft" not in param or ("ifrun" in param["post_dft"] and not param["post_dft"]['ifrun']):
        step3 = rundft_step[0]
        if isinstance(step3,list):
            step3 = step3[0]
    else:
        step3 = ProducePostDFTStep(param["post_dft"],model_output_artifact)
        steps.add(step3)

    steps.outputs.artifacts['outputs']._from = step3.outputs.artifacts["outputs"]
    return steps

def ProduceAllStep(alljobs):
    allstep = []
    stepname = []
    for k in alljobs:
        steps = ProduceOneSteps(k,alljobs[k])
        if steps == None:
            continue
        step = Step(name=k,template=steps)
        allstep.append(step)
        stepname.append(k)
        print("Complete the preparing of %s\n" % k)

    globV.set_value("STEPNAME",stepname)
    return allstep

def ParamParser(param):
    alljobs = {}
    for k,v in param.get('if_run',{}).items():
        if v : alljobs[k] = {}
    for k in alljobs:
        if k not in param['param']:
            print("%s is not set in param,skip" % k)
            continue
        if "run_dft" not in param['param'][k]:
            print("run_dft is not set for param/%s" % k)
            sys.exit(1)
        if "post_dft" not in param['param'][k]:
            #print("post_dft is not set for param/%s" % k)
            param['param'][k]["post_dft"] = {"ifrun":False}

        alljobs[k]["run_dft"] = param['param'][k]["run_dft"]
        alljobs[k]["post_dft"] = param['param'][k]["post_dft"]
    return alljobs

def SetSaveFolder(storefolder=None):
    if storefolder == None:
        import datetime
        from time import strftime
        today = datetime.datetime.now()
        today = today.strftime("%Y%m%d")
        storefolder = os.path.join("result",today)
    globV.set_value("RESULT",storefolder)
    printinfo("set save floder: %s" % storefolder)

def MakeSaveFolder(storefolder=None):
    if storefolder == None:
        storefolder = globV.get_value("RESULT")
    if os.path.isdir(storefolder) and not globV.get_value("OVERRIDE"):
        n = 1
        bk = storefolder + ".bk%d" % n
        while os.path.isdir(bk):
            n += 1
            bk = storefolder + ".bk%d" % n
        os.rename(storefolder,bk)
        print("Folder %s is exist, rename to %s" % (storefolder,bk))
        #os.makedirs(storefolder)

def set_env(param):
    globV.set_value("OUTINFO", param.outinfo)
    globV.set_value("OVERRIDE", param.override)
    printinfo("Set enviroment ...")

    printinfo(param)
    if not os.path.isfile(param.user):
        print("ERROR: Can not find the bohrium account setting file '%s' " % param.user)
        sys.exit(1)

    if not os.path.isfile(param.param):
        print("ERROR: Can not find the test setting file '%s' " % param.param)
        sys.exit(1)

    SetSaveFolder(param.save)

    private_set = json.load(open(param.user))
    config["host"] = "https://workflows.deepmodeling.com"
    config["k8s_api_server"] = "https://workflows.deepmodeling.com"
    bohrium.config["username"] = private_set.get('lbg_username','')
    bohrium.config["password"] = private_set.get('lbg_password','')
    bohrium.config["project_id"] = private_set.get('project_id','')
    s3_config["repo_key"] = "oss-bohrium"
    s3_config["storage_client"] = TiefblueClient()

    globV.set_value("HOST", config["host"])

def waitrun(wf):
    wfid = wf.id
    print("\nResults will be downloaded to %s\n" % globV.get_value("RESULT"))
    hasmakefolder = False
    allstepnames = globV.get_value("STEPNAME")
    hasdownload = False
    while True:
        if len(allstepnames) == 0:
            return
        for i in allstepnames:
            step = wf.query_step(name = i)
            if len(step) > 0: step = step[0]
            else: continue
            if step.phase not in ["Pending","Running"]:
                allstepnames.remove(i)
            else:
                continue
            if step.phase != 'Succeeded':
                print("%s is not Succeeded, please check on: %s, workflow ID is: %s" %
                      (i,globV.get_value("HOST"),wfid))

            if True:
                if not hasmakefolder:
                    MakeSaveFolder()
                    hasmakefolder = True

                print("%s is finished, download the results!" % i)
                download_artifact(step.outputs.artifacts["outputs"],path=globV.get_value("RESULT"))
                hasdownload = True
        time.sleep(4)

    if hasdownload:
        print("Results are downloaded to %s\n" % globV.get_value("RESULT"))

def printinfo(istr):
    if globV.get_value("OUTINFO"):
        print(istr)

def Parser():
    print("Parse commands ...")
    parser = argparse.ArgumentParser(description="This script is used to run a testing")
    parser.add_argument('-p', '--param', type=str, default="param.json",help='the parameter setting file')
    parser.add_argument('-u', '--user', type=str, default="user.json",help='the file for bohrium account information, default is "user.json"')
    parser.add_argument('-s', '--save', type=str, default=None,help='the folder where the results will be put in, default: result/date_of_today (e.g. result/20230101)')
    parser.add_argument('--override', type=int, default=1,help="if the save folder is exist, if override it. 0: no, 1: yes. ")
    parser.add_argument('--outinfo', type=int, default=0,help='if output detail informations, 0: no, 1: yes')
    return parser.parse_args()

def RunJobs():
    param = Parser()
    set_env(param)
    alljobs = ParamParser(json.load(open(param.param)))
    allstep = ProduceAllStep(alljobs)
    
    wf = Workflow(name="abacus-test")
    wf.add(allstep)
    wf.submit()
    print("You can track the flow by using your browser to access the URL:\n %s" % globV.get_value("HOST"))

    waitrun(wf)

