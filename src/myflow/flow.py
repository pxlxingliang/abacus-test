#!/usr/bin/env python3
import os,sys,glob,time,shutil,argparse,json,traceback
import numpy as np

from dflow import (
    Workflow,
    download_artifact,
)

import  os, shutil, glob
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
        alljobs[k]["upload_datahub"] = param['param'][k].get("upload_datahub",False)
    
    if "run_dft" in param:
        alljobs = {"testing":{}}
        alljobs["testing"]["save_path"] = param.get("save_path",None)
        alljobs["testing"]["run_dft"] = param.get("run_dft")
        alljobs["testing"]["post_dft"] = param.get("post_dft",{"ifrun":False})
        alljobs["testing"]["upload_datahub"] = param.get("upload_datahub",None)
    
    if "dataset_info" in param:
        globV.set_value("dataset_info",param['dataset_info']) 
    else:
        globV.set_value("dataset_info",None)   
        
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
    
    if not override:
        def getfname(paraf):
            paraf_tmp = paraf
            n = 1
            while os.path.isfile(paraf_tmp):
                paraf_tmp = paraf + "%d" % n
                n += 1
            return paraf_tmp
            
        paraf = getfname(paraf)


    with open(paraf,'w') as f1: f1.write(globV.get_value("PARAM_CONTEXT")) 
    
def set_env(param):
    globV.set_value("OUTINFO", param.outinfo)
    globV.set_value("OVERRIDE", param.override)
    comm.printinfo("\nSet enviroment ...")
    comm.printinfo(param)

    #read job.json
    if not os.path.isfile(param.param):
        comm.printinfo("ERROR: Can not find the test setting file '%s' " % param.param)
        sys.exit(1)
    comm.printinfo("Read parameter setting from %s" % param.param)
    globV.set_value("PARAM_FNAME", os.path.split(param.param)[1])
    with open(param.param) as f1: 
        globV.set_value("PARAM_CONTEXT", f1.read())       
    param_context = json.load(open(param.param))
    
    #read user config information
    bohrium_executor = False
    if "bohrium_executor" in param_context:
        bohrium_executor = bool(param_context["bohrium_executor"])
    globV.set_value("BOHRIUM_EXECUTOR",bohrium_executor)
    
    if "config" in param_context:
        user_context = param_context.get("config")
    else:
        comm.printinfo("ERROR: please set the config information by \"config\".")
        sys.exit(1)
    globV.set_value("PRIVATE_SET", user_context)
    dflowOP.SetBohrium(user_context,debug=param.debug) 

    #set save folder    
    SetSaveFolder(param.save)


def waitrun(wf,stepnames,allsave_path,postdft_local_jobs,test_name,upload_datahub):
    '''
    stepnames = [[test1_stepname1,test1_stepname2,...],[test2_stepname1,test2_stepname2,...],...]
    allsave_path = [[[save_path,sub_save_path],[save_path,sub_save_path],...],[]...] similar to stepnames
    postdft_local_jobs = [[],[save_path,param["post_dft"]],[],..], null list means no postdft_local,
    '''
    tmp_local_path = len(stepnames) * [False]
    finishtest = []
    makedfolder = []
    for i,istep in enumerate(stepnames):
        finishtest.append(len(istep)*[False])
        if len(postdft_local_jobs[i]) == 0:
            tmp_local_path.append(False)
        else:
            tmp_local_path.append(dflowOP.ProduceRadomPath(".tmp"))
            os.makedirs(tmp_local_path[-1])
    finishtest = np.array(finishtest)  
          
    wfid = wf.id

    while False in finishtest:
        for i,istep in enumerate(stepnames):
            for j,jfinish in enumerate(finishtest[i]):
                if jfinish:
                    continue
                steps = wf.query_step(name = istep[j])
                if len(steps) > 0: 
                    step = steps[0]
                    if step.phase in ["Pending","Running"]:
                        continue
                    finishtest[i][j] = True

                    if len(postdft_local_jobs[i]) > 0:
                        #do postdft on local 
                        #mkdir the tmp work path
                        if not tmp_local_path[i]:
                            tmp_local_path[i] = dflowOP.ProduceRadomPath(".")
                            os.makedirs(tmp_local_path[i])
                        save_path = tmp_local_path[i]
                    else:
                        save_path = allsave_path[i][j][0]
                        part_save_path = os.path.join(allsave_path[i][j][0],allsave_path[i][j][1])                      
                            
                        if part_save_path not in makedfolder:
                            MakeSaveFolder(part_save_path)
                            WriteParamUserFile(storefolder=part_save_path)
                            makedfolder.append(part_save_path)
                        comm.printinfo("%s is finished (remaining %d jobs for this test), download the results to %s!" % 
                                       (step.name,len(finishtest[i]) - np.sum(finishtest[i]),part_save_path))
                    if step.phase != 'Succeeded':
                        comm.printinfo("    This job is not Succeeded, please check on: %s, workflow ID is: %s" %
                              (globV.get_value("HOST"),wfid))
                        
                    try:
                        download_artifact(step.outputs.artifacts["outputs"],path=save_path)
                    except:
                        traceback.print_exc()
            if False not in finishtest[i] and len(postdft_local_jobs[i]) > 0:
                comm.printinfo("Test '%s' has finished the run_dft, do post_dft on local %s\n..." % (test_name[i], tmp_local_path[i]))
                dflowOP.RunPostDFTLocal(tmp_local_path[i],postdft_local_jobs[i][0],postdft_local_jobs[i][1])
                comm.printinfo("Test '%s' has finished the post_dft, move the outputs to %s" % (test_name[i],postdft_local_jobs[i][0]))
                part_save_path = postdft_local_jobs[i][0]
                if part_save_path not in makedfolder:
                    WriteParamUserFile(storefolder=part_save_path)
                    makedfolder.append(part_save_path)
                    
            if False not in finishtest[i] and upload_datahub[i]:
                steps = wf.query_step(key = test_name[i])
                if len(steps) > 0:
                    if steps[0].phase == "Succeeded":
                        comm.printinfo("upload the outputs of '%s' to datahub" % test_name[i])
                        dflowOP.Upload2Datahub(steps[0].outputs.artifacts['outputs'],upload_datahub[i])
                    else:
                        comm.printinfo("test '%s' is finished, but the phase is %s, will not upload to datahub" % (test_name[i],steps[0].phase))
                else:
                    comm.printinfo("WARNING: test '%s' is finished, but query_step can not find it, please check the dflow!" % test_name[i])
                    
                              
        time.sleep(4)

def RunJobs(param):
    set_env(param)
    alljobs = ParamParser(json.load(open(param.param)))
    allstep,stepname,allsave_path,postdft_local_jobs,test_name,upload_datahub = dflowOP.ProduceAllStep(alljobs)

    if len(allstep) == 0:
        comm.printinfo("No step is produced, exit!!!")
        sys.exit(1)
    
    if globV.get_value("BOHRIUM_EXECUTOR"):
        wf = Workflow(name="abacustest",context=globV.get_value("BRM_CONTEXT"))
    else:
        wf = Workflow(name="abacustest")

    wf.add(allstep)
    wf.submit()
    if param.command == 'mlops-submit':
        return
    comm.printinfo("job ID: %s, UID: %s" % (wf.id,wf.uid))
    comm.printinfo("You can track the flow by using your browser to access the URL:\n %s\n" % globV.get_value("HOST"))

    waitrun(wf,stepname,allsave_path,postdft_local_jobs,test_name,upload_datahub)
        
def CheckStatus(param):
    if os.path.isfile(param.param):
        private_set = json.load(open(param.param))
        if "USER" in private_set:
            private_set = private_set["USER"]
    else:
        print("config file is not found!\nUse the default setting!")
        private_set = {}
    
        
    dflowOP.SetBohrium(private_set,debug=False) 
    
    jobid = param.job_id
    
    wf = Workflow(id = jobid)
        
    try:
        return wf.query_status()
    except:
        comm.printinfo("Query status error")
        traceback.print_exc()
        return "not-running"
        
