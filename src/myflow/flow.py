#!/usr/bin/env python3
import os,sys,glob,time,shutil,argparse,json

from dflow import (
    Workflow,
    download_artifact,
)
import  os, shutil, glob
from dflow import config, s3_config
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient
from . import globV,dflowOP,comm

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
        
    globV.set_value("ABBREVIATION",param.get('ABBREVIATION',{}))
        
    return alljobs

def SetSaveFolder(storefolder=None):
    if storefolder == None:
        import datetime
        from time import strftime
        today = datetime.datetime.now()
        today = today.strftime("%Y%m%d")
        storefolder = os.path.join("result",today)
    globV.set_value("RESULT",storefolder)
    comm.printinfo("set save floder: %s" % storefolder)

def MakeSaveFolder(storefolder=None):
    storefolder = globV.get_value("RESULT") if storefolder == None else storefolder
    if not os.path.isdir(storefolder):
        os.makedirs(storefolder)
    elif not globV.get_value("OVERRIDE"):
        n = 1
        bk = storefolder + ".bk%d" % n
        while os.path.isdir(bk):
            n += 1
            bk = storefolder + ".bk%d" % n
        os.rename(storefolder,bk)
        print("Folder %s is exist, rename to %s" % (storefolder,bk))
        os.makedirs(storefolder)
    
def WriteParamUserFile(storefolder=None,override=False):
    storefolder = globV.get_value("RESULT") if storefolder == None else storefolder
    paraf = os.path.join(storefolder,globV.get_value("PARAM_FNAME"))
    userf = os.path.join(storefolder,globV.get_value("USER_FNAME"))   
    
    if not override:
        def getfname(paraf):
            paraf_tmp = paraf
            n = 1
            while os.path.isfile(paraf_tmp):
                paraf_tmp = paraf + "%d" % n
                n += 1
            return paraf_tmp
            
        paraf = getfname(paraf)
        userf = getfname(userf)

    with open(paraf,'w') as f1: f1.write(globV.get_value("PARAM_CONTEXT")) 
    with open(userf,'w') as f1: f1.write(globV.get_value("USER_CONTEXT")) 
    
def set_env(param):
    globV.set_value("OUTINFO", param.outinfo)
    globV.set_value("OVERRIDE", param.override)
    comm.printinfo("Set enviroment ...")

    comm.printinfo(param)
    if not os.path.isfile(param.user):
        print("ERROR: Can not find the bohrium account setting file '%s' " % param.user)
        sys.exit(1)

    if not os.path.isfile(param.param):
        print("ERROR: Can not find the test setting file '%s' " % param.param)
        sys.exit(1)

    globV.set_value("PARAM_FNAME", os.path.split(param.param)[1])
    with open(param.param) as f1: 
        globV.set_value("PARAM_CONTEXT", f1.read())
    globV.set_value("USER_FNAME", os.path.split(param.user)[1])
    with open(param.user) as f1: 
        globV.set_value("USER_CONTEXT", f1.read())
    
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
                    WriteParamUserFile()
                    hasmakefolder = True

                print("%s is finished, download the results!" % i)
                download_artifact(step.outputs.artifacts["outputs"],path=globV.get_value("RESULT"))
                hasdownload = True
        time.sleep(4)

    if hasdownload:
        print("Results are downloaded to %s\n" % globV.get_value("RESULT"))

def Parser():
    print("Parse commands ...")
    parser = argparse.ArgumentParser(description="This script is used to run a testing")
    parser.add_argument('-p', '--param', type=str, default="param.json",help='the parameter setting file')
    parser.add_argument('-u', '--user', type=str, default="user.json",help='the file for bohrium account information, default is "user.json"')
    parser.add_argument('-s', '--save', type=str, default=None,help='the folder where the results will be put in, default: result/date_of_today (e.g. result/20230101)')
    parser.add_argument('--override', type=int, default=1,help="when the save folder exists, if override it. 0: no, 1: yes. ")
    parser.add_argument('--outinfo', type=int, default=1,help='if output detail informations, 0: no, 1: yes')
    return parser.parse_args()

def RunJobs():
    param = Parser()
    set_env(param)
    alljobs = ParamParser(json.load(open(param.param)))
    allstep = dflowOP.ProduceAllStep(alljobs)
    
    wf = Workflow(name="abacus-test")
    wf.add(allstep)
    wf.submit()
    print("You can track the flow by using your browser to access the URL:\n %s" % globV.get_value("HOST"))

    waitrun(wf)
        

