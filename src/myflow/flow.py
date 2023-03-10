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
    """
    {
        "dataset_info":{
            "type": "local"/"datahub",
            "dataset_urn":,
            "download":false
        },
        "running_setting" : RUNNING_SETTING_FILE,

        "bohrium_goup_name": ,
        "ABBREVIATION":{},
        "save_path" : PATH_THE_FINAL_RESULT_WILL_BE_DOWNLOADED_TO,
        "run_dft" : {},
        "post_dft": {},
        "upload_datahub": {},
        "report":{}
    }
    """
    
    alljobs = {}

    alljobs["save_path"] = param.get("save_path",None)
    alljobs["run_dft"] = param.get("run_dft")
    alljobs["post_dft"] = param.get("post_dft",{"ifrun":False})
    alljobs["upload_datahub"] = param.get("upload_datahub",None)  #used to upload local files to datahub
    alljobs["upload_tracking"] = param.get("upload_tracking",None)  #used to upload tracking
    alljobs["report"] = param.get("report",None)
    alljobs["bohrium_group_name"] = param.get("bohrium_group_name","abacustesting")
    alljobs["ABBREVIATION"] = param.get("ABBREVIATION",{})
   
    dataset_info = param.get("dataset_info")
    globV.set_value("dataset_info",dataset_info)
    if dataset_info:
        if dataset_info.get("type") == "datahub" and dataset_info.get("download"):
            urn = dataset_info.get("dataset_urn","").strip().split()
            for iurn in urn:
                #download data
                data_path = comm.GetBakFile("data")
                uri,storage_client = dflowOP.GetURI(urn=iurn)
                dflowOP.DownloadURI(uri,path=data_path)
                comm.printinfo("download datahub data to %s" % data_path)

    #read the setting file if is setted
    setting_file = param.get("running_setting")
    if setting_file:
        if dataset_info and dataset_info.get("type") == "datahub" and not param.get("running_setting_from_local",False):
            urn = dataset_info.get("dataset_urn","").strip().split()
            if len(urn) == 1:
                iurn = urn[0]
            elif len(urn) > 1:
                iurn = urn[-1]
            uri,storage_client = dflowOP.GetURI(urn=iurn)
            uri += "/" + setting_file
            comm.printinfo("download  the setting file '%s' from datahub " % setting_file)
            dflowOP.DownloadURI(uri)    

        if not os.path.isfile(setting_file):
            comm.printinfo("Has specify the 'running_setting', but can not find file %s, skip to read it!" % setting_file)
        else:   
            comm.printinfo("Read setting from %s" % setting_file) 
            setting = json.load(open(setting_file))
            if "save_path" in setting:
                alljobs["save_path"] = setting["save_path"]
            if "run_dft" in setting:
                alljobs["run_dft"] = setting["run_dft"]
            if "post_dft" in setting:
                alljobs["post_dft"] = setting["post_dft"]
            if "upload_datahub" in setting:
                alljobs["upload_datahub"] = setting["upload_datahub"]
            if "upload_tracking" in setting:
                alljobs["upload_tracking"] = setting["upload_tracking"]
            if "report" in setting:
                alljobs["report"] = setting["report"]  
            if "bohrium_group_name" in setting:
                alljobs["bohrium_group_name"] = setting["bohrium_group_name"]
            if "ABBREVIATION" in setting:  
                alljobs["ABBREVIATION"] = setting["ABBREVIATION"]

    #if upload datahub or tracking is defined, then replace the definition in post_dft by current setting
    if alljobs["upload_datahub"] != None:
        if "upload_datahub" not in alljobs["post_dft"]:
            alljobs["post_dft"]["upload_datahub"] = {}
        for k,v in alljobs["upload_datahub"].items():
            alljobs["post_dft"]["upload_datahub"][k] = v

    if alljobs["upload_tracking"] != None:
        if "upload_tracking" not in alljobs["post_dft"]:
            alljobs["post_dft"]["upload_tracking"] = {}
        for k,v in alljobs["upload_tracking"].items():
            alljobs["post_dft"]["upload_tracking"][k] = v

    #print(alljobs)
    #sys.exit(1)
    globV.set_value("ABBREVIATION",alljobs.get('ABBREVIATION',{}))
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
    globV.set_value("PARAM", param_context)
    
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
    
    
    report = param_context.get("report",{})
    globV.set_value("REPORT", report)
    
    from datetime import datetime
    globV.set_value("BEGIN", str(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

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

def ReportMetrics():
    report_setting = globV.get_value("report",{})
    if not report_setting or not report_setting.get("ifrun",True):
        return

    from abacustest import outresult
    allresults = outresult.GetAllResults(report_setting)
    if not allresults:
        return
    
    split_example=None if len(allresults["type_name"]) == 1 else "----"
    cc_outparam,allparam_value = outresult.OutParam(allresults,split_example=split_example)
    cc_outmetrics,allmetric_value = outresult.OutMetrics(allresults,allparam_value)
    
    #report
    from datetime import datetime
    param_context = globV.get_value("PARAM")
    report  = """\n\n\t\tABACUS TESTING REPORT\n"""
    report += "testing begin: %s\n" % globV.get_value("BEGIN")
    report += "testing   end: %s\n" % str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    report += "run_dft setting:\n"
    for irun in param_context.get("run_dft",[]):
        image = globV.get_value("ABBREVIATION").get(irun.get("image"),irun.get("image"))
        bohrium = irun.get("bohrium")
        report += "\timage: %s\n\tbohrium:%s\n\texample:%s\n" % (image,str(bohrium),str(irun.get("example")))
        
    if len(allresults["metrics"]) > 0:
        report += cc_outmetrics
    report += cc_outparam

    report_file = report_setting.get("save_file")
    if report_file != None:
        with open(report_file,'w') as f1:
            f1.write(report)
    comm.printinfo(report)

def RunJobs(param):
    set_env(param)
    alljobs = ParamParser(json.load(open(param.param)))
    allstep,stepname,allsave_path,postdft_local_jobs,test_name,upload_datahub = dflowOP.ProduceAllStep(alljobs)

    if len(allstep) == 0:
        comm.printinfo("No step is produced, exit!!!")
    else:
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
    
    if globV.get_value("REPORT"):
        ReportMetrics()
        
def CheckStatus(param):
    if os.path.isfile(param.param):
        private_set = json.load(open(param.param))
        if "config" in private_set:
            private_set = private_set["config"]
        elif "USER" in private_set:
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
        
