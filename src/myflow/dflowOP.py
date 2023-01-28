#!/usr/bin/env python3
import os,sys,glob,time,shutil,argparse,json
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
    S3Artifact
)
import subprocess, os, shutil, glob
from pathlib import Path
from typing import List
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
                    if os.path.isfile(op_in["collectdata_pythonlib"]):
                        tpath = os.path.split(str(op_in["collectdata_pythonlib"]).rstrip("/"))[0]
                    elif os.path.isdir(op_in["collectdata_pythonlib"]):
                        tpath = str(op_in["collectdata_pythonlib"])
                    else:
                        print("str(op_in[collectdata_pythonlib]:",str(op_in["collectdata_pythonlib"]),"is either a file or dir")
                        tpath = str(op_in["collectdata_pythonlib"])
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
        comm.printinfo("set bohrium: %s"%str(bohrium_set))
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
    abacustestpath = "/".join(__file__.split('/')[:-3])
    collectdata_pythonlib = [abacustestpath]
    image = globV.get_value("ABBREVIATION").get(param['image'],param['image'])
    comm.printinfo("image: %s" % image)
    if DoSlices:
        pt = PythonOPTemplate(RunDFT,image=image,
                    slices=Slices(sub_path = True,
                                  input_artifact=["example"],
                                  output_artifact=["outputs"]))
        sub_path = [""]
    else:
        pt = PythonOPTemplate(RunDFT,image=image)
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
    image = globV.get_value("ABBREVIATION").get(param['image'],param['image'])
    pt = PythonOPTemplate(PostDFT,image=image)
    
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
            if len(examples) > 1:
                rundft_step.append(step1(examples,rundft,istep,True))
            else:
                rundft_step.append(step1(examples,rundft,istep,False))
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
