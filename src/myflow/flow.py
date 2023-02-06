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
            comm.printinfo("%s is not set in param,skip" % k)
            continue
        if "run_dft" not in param['param'][k]:
            comm.printinfo("run_dft is not set for param/%s" % k)
            sys.exit(1)
        if "post_dft" not in param['param'][k]:
            #comm.printinfo("post_dft is not set for param/%s" % k)
            param['param'][k]["post_dft"] = {"ifrun":False}

        alljobs[k]["run_dft"] = param['param'][k]["run_dft"]
        alljobs[k]["post_dft"] = param['param'][k]["post_dft"]
        alljobs[k]["save_path"] = param['param'][k].get("save_path",None)
        
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
        comm.printinfo("Folder %s is exist, rename to %s" % (storefolder,bk))
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
    comm.printinfo("\nSet enviroment ...")

    comm.printinfo(param)
    if not os.path.isfile(param.user):
        comm.printinfo("ERROR: Can not find the bohrium account setting file '%s' " % param.user)
        sys.exit(1)

    if not os.path.isfile(param.param):
        comm.printinfo("ERROR: Can not find the test setting file '%s' " % param.param)
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

def waitrun(wf,stepnames,allsave_path):
    wfid = wf.id
    #comm.printinfo("\nResults will be downloaded to %s\n" % globV.get_value("RESULT"))
    makedfolder = []
    allstepnames = set(stepnames)
    finishstep = []
    hasdownload = False
    #print(stepnames)
    while True:
        if len(finishstep) == len(stepnames):
            return
        for i,ia in enumerate(allstepnames):
            steps = wf.query_step(name = ia)
            #comm.printinfo("i:%d ia:%s step number: %d" % (i,ia,len(steps)))
            if len(steps) > 0: 
                idx = -1
                for i,step in enumerate(steps):
                    idx += stepnames[idx+1:].index(ia) + 1   #find the idx of ia in stepnames, to get the savepath
                    save_path = allsave_path[idx] if allsave_path[idx] != None else globV.get_value("RESULT")
                    if step.name in finishstep:
                        continue
                    if step.phase in ["Pending","Running"]:
                        continue
                    finishstep.append(step.name) 

                    if save_path not in makedfolder:
                        MakeSaveFolder(save_path)
                        WriteParamUserFile(storefolder=save_path)
                        makedfolder.append(save_path)
                    comm.printinfo("%4d/%4d: %s is finished, download the results to %s!" % (len(finishstep), len(stepnames),step.name,save_path))
                    if step.phase != 'Succeeded':
                        comm.printinfo("    This job is not Succeeded, please check on: %s, workflow ID is: %s" %
                              (globV.get_value("HOST"),wfid))
                    download_artifact(step.outputs.artifacts["outputs"],path=save_path)
                    hasdownload = True
        time.sleep(4)

    if hasdownload:
        comm.printinfo("Results are downloaded to %s\n" % globV.get_value("RESULT"))

def Parser():
    comm.printinfo("\nParse commands ...")
    parser = argparse.ArgumentParser(description="This script is used to run a testing")
    parser.add_argument('-p', '--param', type=str, default="job.json",help='the job setting file, default is job.json')
    parser.add_argument('-u', '--user', type=str, default="user.json",help='the file for bohrium account information, default is "user.json"')
    parser.add_argument('-s', '--save', type=str, default=None,help='the folder where the results will be put in, default: result/date_of_today (e.g. result/20230101)')
    parser.add_argument('--override', type=int, default=1,help="when the save folder exists, if override it. 0: no, 1: yes. ")
    parser.add_argument('--outinfo', type=int, default=1,help='if output detail informations, 0: no, 1: yes')
    return parser.parse_args()

def RunJobs():
    param = Parser()
    set_env(param)
    alljobs = ParamParser(json.load(open(param.param)))
    allstep,stepname,allsave_path = dflowOP.ProduceAllStep(alljobs)
    
    wf = Workflow(name="abacustest")
    wf.add(allstep)
    wf.submit()
    comm.printinfo("job ID: %s" % wf.id)
    comm.printinfo("You can track the flow by using your browser to access the URL:\n %s\n" % globV.get_value("HOST"))

    waitrun(wf,stepname,allsave_path)
        

