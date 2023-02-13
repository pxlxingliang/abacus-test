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

from dflow import config, s3_config
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient

def SetBohrium(private_set,debug=False):  
    if debug:
        config["mode"] = "debug"
        host = "LOCAL"
        client = "LOCAL"
    else:
        #config["host"] = "https://workflows.deepmodeling.com"
        #config["k8s_api_server"] = "https://workflows.deepmodeling.com"
        comm.printinfo("set config info...")
        if private_set.get("config_host","").strip() != "":
            config["host"] = private_set.get("config_host","").strip()
            host = config["host"]
            comm.printinfo("set config['host']: %s" % host)
        
        if private_set.get("s3_config_endpoint","").strip() != "": 
            s3_config["endpoint"] =  private_set.get("s3_config_endpoint","").strip()
            comm.printinfo("set s3_config['endpoint']: %s" % s3_config["endpoint"])
        
        if private_set.get("config_k8s_api_server","").strip() != "":
            config["k8s_api_server"] = private_set.get("config_k8s_api_server","").strip()
            comm.printinfo("set config['k8s_api_server']: %s" % config["k8s_api_server"])
            
        if private_set.get("config_token","").strip() != "":
            config["token"] = private_set.get("config_token","").strip()
            comm.printinfo("set config['token']: ...")
   
        bohrium.config["username"] = private_set.get('lbg_username','')
        bohrium.config["password"] = private_set.get('lbg_password','')
        bohrium.config["project_id"] = private_set.get('project_id','')
        comm.printinfo("set bohrium.config['username']/['password']/['project_id']: %s/.../%s" 
                       % (bohrium.config["username"],bohrium.config["project_id"]))
        s3_config["repo_key"] = "oss-bohrium"
        s3_config["storage_client"] = TiefblueClient()
        client = s3_config["storage_client"]
    
    globV.set_value("HOST", host)
    globV.set_value("storage_client", client)

class RunDFT(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "abacustest_example": Artifact(Path),
                "command": str,
                "sub_path": [str],
                "sub_save_path": str,
                "collectdata_script": Artifact(Path),
                "collectdata_script_name": [str],
                "collectdata_pythonlib": Artifact(Path),
                "collectdata_pythonlib_folders":[str],
                "collectdata_lib": Artifact(Path),
                "collectdata_command": str,
                "outputfiles":[str],
                "sum_save_path_example_name": bool
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

        print(op_in["abacustest_example"])
        print(op_in["sub_path"])
        print(op_in["sub_save_path"])
        print(op_in["collectdata_script"])
        example_pathname = "dflow_abacustest_example_"
        outpath = []
        root_path_0 = str(op_in["abacustest_example"]).split(example_pathname)[0]
        sub_save_path = str(op_in['sub_save_path'])

        for isub,sub_path in enumerate(op_in["sub_path"]):
            root_path = os.path.join(root_path_0,example_pathname+str(isub))
            example_path = root_path
            if bool(op_in["sum_save_path_example_name"]):
                sub_save_path = os.path.join(op_in['sub_save_path'],sub_path)
            print("sub_save_path",sub_save_path)
            if sub_save_path != "":
                os.chdir(root_path)
                sub_save_path_list = sub_save_path.split('/')

                #if exist folder name that is same as sub_save_path, create a tmp sub_save_path
                #and then rename the tmp to sub_save_path
                hasexist = False
                if os.path.isdir(sub_save_path_list[0]): 
                    sub_save_path_list0 = sub_save_path_list[0]
                    sub_save_path_list[0] = GetBakFile(sub_save_path_list[0])
                    sub_save_path_tmp = "/".join(sub_save_path_list)
                    hasexist = True
                else:
                    sub_save_path_tmp = sub_save_path
                os.makedirs(sub_save_path_tmp)
                for idir in os.listdir("."):
                    if Path(idir) not in [Path(".dflow"),Path(sub_save_path_list[0])]: 
                        shutil.move(idir,sub_save_path_tmp)
                if hasexist:
                    os.replace(sub_save_path_list[0],sub_save_path_list0)
                
                example_path = os.path.join(root_path,op_in['sub_save_path'])

            work_path = os.path.join(example_path,sub_path)
                
            print("work path:",work_path)
            log = ""
            os.chdir(work_path)
            if op_in["command"].strip() != "":
                print("BEGIN THE DFT RUNNING...")
                log += os.popen("(%s) 2>&1" % op_in["command"]).read()
                print(log)
        
            os.chdir(work_path)
            script_folder = "SCRIPT"
            if op_in["collectdata_script"] != None:
                if os.path.isdir(script_folder):
                    os.rename(script_folder,GetBakFile(script_folder))
                os.makedirs(script_folder)
                for iii in range(len(op_in["collectdata_script_name"])):
                    script_source_path = str(op_in["collectdata_script"]).split("dflow_collectdata_script_%d"%iii)[0] + "dflow_collectdata_script_%d"%iii
                    if bool(op_in["sum_save_path_example_name"]):
                        dst_path = os.path.join(script_folder,op_in["collectdata_script_name"][iii])
                    else:
                        dst_path = script_folder
                    if os.path.isfile(script_source_path):
                        shutil.copy(script_source_path,dst_path)
                    elif  os.path.isdir(script_source_path):
                        shutil.copytree(script_source_path,dst_path,dirs_exist_ok=True)

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

            if op_in["collectdata_lib"] != None:
                tpath = str(op_in["collectdata_lib"])
                script = os.path.join(tpath,"collectdata.py")
                for chkcmd,headline in [("env python -h","#!/usr/bin/env python"),
                                        ("env python3 -h","#!/usr/bin/env python3"),
                                        ("/usr/bin/python -h","#!/usr/bin/python"),
                                        ("/usr/bin/python3 -h","#!/usr/bin/python3")]:
                    if not os.system(chkcmd+">/dev/null 2>&1"):
                        with open(script) as f1: lines = ["%s\n"%headline] + f1.readlines()
                        with open(script,'w') as f1: f1.writelines(lines)
                        break
                cmd += "export PATH=%s:$PATH && " % tpath
                os.chmod(script,0o777)
                shutil.copy(script,os.path.join(tpath,"collectdata"))
            
            if op_in["collectdata_command"] != "":
                cmd += op_in["collectdata_command"]
                log += "\ncommand:" + cmd + "\n"
                log += os.popen("(%s) 2>&1" % cmd).read()
                print(log)
            
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
            outt +=  os.popen("(%s) 2>&1" % op_in["command"]).read()
        
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
            image_pull_policy = "IfNotPresent"
        )
        #comm.printinfo("set bohrium: %s"%str(bohrium_set))
        return dispatcher_executor,bohrium_set
    else:
        return None,None

def ProduceRunDFTStep(step_name,
                      command,   #command of do dft running
                      example_artifact, 
                      image, 
                      collectdata_script = [],
                      collectdata_script_name = [] ,
                      collectdata_command = "", 
                      collectdata_pythonlib = [], 
                      outputs=[], 
                      executor = None, 
                      example_names =[""],
                      DoSlices=True,  
                      sub_save_path = "",
                      datahub = False):
    #define template
    if DoSlices:
        pt = PythonOPTemplate(RunDFT,image=image,
                    slices=Slices(sub_path = True,
                                  input_artifact=["abacustest_example"],
                                  output_artifact=["outputs"]))
        sub_path = [""]
    else:
        pt = PythonOPTemplate(RunDFT,image=image)
        sub_path = list(example_names)

    #define example artifacts
    artifacts = {"abacustest_example": example_artifact} 
    
    #collectdata_script
    if len(collectdata_script) > 0:
        artifacts["collectdata_script"] = collectdata_script
    else:
        pt.inputs.artifacts["collectdata_script"].optional = True
        
    #collectdata_pythonlib    
    if len(collectdata_pythonlib) > 0:
        artifacts["collectdata_pythonlib"] = upload_artifact(collectdata_pythonlib)
        if len(collectdata_pythonlib) == 1:
            collectdata_pythonlib_folders = [""]
        else:
            collectdata_pythonlib_folders = collectdata_pythonlib
    else:
        pt.inputs.artifacts["collectdata_pythonlib"].optional = True
        collectdata_pythonlib_folders = [""]
    
    #upload collectdata_lib    
    abacustestpath = "/".join(__file__.split('/')[:-2])
    artifacts["collectdata_lib"] = upload_artifact(abacustestpath)
    
    #produce step
    step = Step(name=step_name,template=pt,
            parameters = {"command": command, "sub_path":sub_path,
                          "sum_save_path_example_name": datahub,
                          "collectdata_script_name": collectdata_script_name,
                          "sub_save_path": sub_save_path,
                          "collectdata_command":collectdata_command,
                          "collectdata_pythonlib_folders": collectdata_pythonlib_folders,
                          "outputfiles":outputs},
            artifacts =  artifacts,continue_on_failed=True)

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
    comm.printinfo("image: %s" % image)
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

    executor,bohrium_set = ProduceExecutor(param)
    if executor != None: step.executor = executor
    return step

def GetBakFile(sfile):
    n = 1
    bk = sfile + ".bak%d" % n
    while os.path.exists(bk):
        n += 1
        bk = sfile + ".bak%d" % n
    return bk

def ProduceRadomPath(tmp_path):
    import random,string
    workpath = os.path.join(tmp_path,"." + "".join(random.sample(string.ascii_letters+string.digits,11)))
    while os.path.isdir(workpath):
        workpath = os.path.join(tmp_path,"." + "".join(random.sample(string.ascii_letters+string.digits,11)))
    return workpath
    
def RunPostDFTLocal(work_path, save_path, param):
    cwd = os.getcwd()
    command = param.get("command","")
    scripts = []
    for iscript in param.get("script",[]):
        scripts += glob.glob(iscript)

    bkf = False
    if len(scripts) > 0:
        SCRIPT = os.path.join(work_path,"SCRIPT")
        if os.path.isdir(SCRIPT):
            bkf = GetBakFile(SCRIPT)
            os.rename(SCRIPT,bkf)
        os.makedirs(SCRIPT)
        for iscript in scripts:
            if os.path.isfile(iscript):
                shutil.copy(iscript,SCRIPT,follow_symlinks=True)
            elif os.path.isdir(iscript):
                shutil.copytree(iscript,SCRIPT,symlinks=False,ignore_dangling_symlinks=True,dirs_exist_ok=True)

    if command == None or command == False or str(command).strip() == "":
        pass
    else:
        comm.printinfo(os.popen("(%s) 2>&1" % command).read()) 
        
    os.chdir(work_path)
    files_need_move = [] 
    outputs = param.get("outputs",[])
    if len(outputs) > 0:
        for ifile in outputs:
            for iifile in glob.glob(ifile):
                files_need_move.append(iifile)
    
    os.chdir(cwd)
    if len(files_need_move) == 0:
        shutil.copytree(work_path,save_path)
    else:
        for ifile in files_need_move:
            path,file_tmp = os.path.split(ifile)
            new_path = os.path.join(save_path,path)
            if not os.path.isdir(new_path):
                os.makedirs(new_path)
            shutil.move(os.path.join(work_path,ifile),new_path)
    shutil.rmtree(work_path)
    return

def FindLocalExamples(example):
    #use glob.glob find all examples, and transfer to artifact
    #example = [*[*],*]
    examples = []    
    examples_name = []  
    for i in example:
        if isinstance(i,list):
            example_tmp = []
            for j in i:
                example_tmp += glob.glob(j)
            if len(example_tmp) > 0:
                examples.append([upload_artifact(i,archive=None) for i in example_tmp])
                examples_name.append(example_tmp)
        elif isinstance(i,str):
            for ii in glob.glob(i):
                examples.append([upload_artifact(ii,archive=None)])
                examples_name.append([ii])
        else:
            comm.printinfo(i,"element of 'example' should be a list, or str")
            
    return examples,examples_name

def GetURI(urn,privateset=None):
    if privateset == None:
        privateset = globV.get_value("PRIVATE_SET")
    project = privateset.get("datahub_project")
    gms_url = privateset.get("datahub_gms_url")
    gms_token = privateset.get("datahub_gms_token")
    bohrium_username = privateset.get("datahub_username")
    bohrium_password = privateset.get("datahub_password")
    
    from dp.metadata import MetadataContext
    from dp.metadata.utils.storage import TiefblueStorageClient
    metadata_storage_client = TiefblueStorageClient(bohrium_username,bohrium_password)
    with MetadataContext(project=project,endpoint = gms_url,token = gms_token,storage_client=metadata_storage_client) as context:
        dataset = context.client.get_dataset(urn)
        if dataset == None:
            comm.printinfo("ERRO: can not catch the dataset for urn:'%s'. \nSkip it!!!\n" % urn)
            return None,None
        uri = str(dataset.uri)  
    storage_client =  context.storage_client #globV.get_value("storage_client")
    return uri,storage_client

def ProduceArtifact(storage_client,uri,name):
    tmp_artifact = []
    tmp_name = []
    #uri = "/12180/" + str(globV.get_value("PRIVATE_SET").get("project_id")) + "/" + uri
    #print(uri)
    res = storage_client.list(uri,recursive=False)
    if len(res) == 0:
        comm.printinfo("Can not find uri:%s" % uri)
    else:
        tmp_artifact.append(S3Artifact(key=res[0]))
        tmp_name.append(name)
        if len(res) > 1:
            comm.printinfo("WARNING: more then one result for uri: %s, use the first one" % uri)
        #comm.printinfo("example name: %s, uri: %s" % (name,uri) )
        #download_artifact(tmp_artifact[0],path='result-test/test')
        #sys.exit(1)
    return tmp_artifact,tmp_name

def FindDataHubExamples(example,scripts,urn):
    examples,examples_name,collectdata_script,collectdata_script_name = [],[],[],[]
    uri,storage_client = GetURI(urn)  
    if uri == None:
        return examples,examples_name
    for ie in example:
        if isinstance(ie,list):
            tmp_artifact = []
            tmp_name = []
            for j in ie:
                if j.strip() =="": continue
                uri_tmp = uri + "/" + j
                tmp_artifact1,tmp_name1 = ProduceArtifact(storage_client,uri_tmp,j)
                tmp_artifact += tmp_artifact1
                tmp_name += tmp_name1
        elif isinstance(ie,str):
            if ie.strip() == "": continue
            uri_tmp = uri + "/" + ie
            tmp_artifact,tmp_name = ProduceArtifact(storage_client,uri_tmp,ie)
        else:
            comm.printinfo(ie,"element of 'example' should be a list, or str")
        
        if len(tmp_artifact) > 0:
            examples.append(tmp_artifact)    
            examples_name.append(tmp_name)
    
    for iscript in scripts:
        if iscript.strip() == "": continue
        uri_tmp = uri + "/" + iscript
        tmp_script, script_name = ProduceArtifact(storage_client,uri_tmp,iscript)
        collectdata_script += tmp_script
        collectdata_script_name += script_name
        
    return examples,examples_name,collectdata_script,collectdata_script_name

def SplitGroup(examples,examples_name,ngroup): 
    #examples = [*[*],*]
    newexamples = [] 
    newexamples_name = []
    se = 0
    mod = len(examples) % ngroup
    for i in range(ngroup):
        example_tmp = []
        example_name_tmp = []
        add = 1 if mod > 0 else 0 
        ee = se + int(len(examples)/ngroup) + add
        for ie in range(se,ee):
            example_tmp += examples[ie]
            example_name_tmp += examples_name[ie]
        newexamples.append(example_tmp)
        newexamples_name.append(example_name_tmp)
        
        if mod > 0: mod -= 1
        se = ee
    return newexamples,newexamples_name

def ParseSavePath(save_path):
    if save_path != None:
        save_path = save_path.strip()
        if save_path == '':
            save_path = globV.get_value("RESULT")
    else:
        save_path = globV.get_value("RESULT")
    return save_path

def ParseSubSavePath(sub_save_path):
    if sub_save_path == None:
        sub_save_path = ""
    elif isinstance(sub_save_path,str):
        sub_save_path = sub_save_path.strip()
        if Path(sub_save_path) == Path("."):
            sub_save_path = ""
    else:
        comm.printinfo("the type of 'sub_save_path' should be 'str', but not '%s'. %s" % (type(sub_save_path),str(sub_save_path)))
        sub_save_path = ""
    return sub_save_path     

def ParsePostDft(param,stepname):
    post_dft = True
    post_dft_local = False
    model_output_artifact = None
    if "post_dft" not in param or ("ifrun" in param["post_dft"] and not param["post_dft"]['ifrun']):
        post_dft = False
    else:
        image = param["post_dft"].get("image")
        if image == None or str(image).strip() == "":
            post_dft_local = True
        else:
            model_output_artifact = S3Artifact(key="{{workflow.name}}/%s"%stepname) #if has post dft, save all dft results in one place
    
    return post_dft,post_dft_local,model_output_artifact
    
def ProduceOneSteps(stepname,param):
    if "run_dft" not in param:
        return None,None,None,None
    rundft_set = param["run_dft"]
    post_dft,post_dft_local,model_output_artifact = ParsePostDft(param,stepname)
    save_path = ParseSavePath(param.get("save_path",None))
    
    steps = Steps(name=stepname+"-steps",
                  outputs=Outputs(artifacts={"outputs" : OutputArtifact()})) 
    #rundft
    allstepname = []
    all_save_path = []
    rundft_step = []
    istep = 1
    comm.printinfo("\n%s\nPreparing run_dft" % stepname)
    for rundft in rundft_set:
        if 'ifrun' in rundft and not rundft['ifrun']: continue
        
        sub_save_path = ParseSubSavePath(rundft.get("sub_save_path",""))
        executor,bohrium_set = ProduceExecutor(rundft)                   
        #split the examples to ngroup and produce ngroup step
        if globV.get_value("dataset_info") == None:
            example_source = rundft.get("example_source","local").strip()
        else:
            example_source = globV.get_value("dataset_info").get("type")
            urn = globV.get_value("dataset_info").get("dataset_urn",None)
                
        if example_source == 'datahub':
            urn = urn if urn != None else rundft["urn"]
            examples,examples_name,collectdata_script,collectdata_script_name = \
                FindDataHubExamples(rundft['example'],rundft.get("collectdata_script",[]), urn)
            datahub = True
            comm.printinfo("Example_source is 'datahub', use the example and collectdata_script in 'datahub'")
        else:
            datahub = False
            examples,examples_name = FindLocalExamples(rundft['example'])
            collectdata_script = [upload_artifact(i) for i in rundft.get("collectdata_script",[])]
            collectdata_script_name = rundft.get("collectdata_script",[])
            if example_source != "local":
                comm.printinfo("example_source should be local or datahub, but not %s, will find example on local" % example_source)
        print(collectdata_script_name)        
        if len(examples) > 0:
            image = globV.get_value("ABBREVIATION").get(rundft.get("image"),rundft.get("image"))
            ngroup = rundft.get("ngroup",0)
            if ngroup == None or ngroup < 1:
                ngroup = len(examples)

            for example_tmp,example_name_tmp in zip(*SplitGroup(examples,examples_name,ngroup)):
                stepname_tmp = "%s-run-dft-%d"%(stepname,istep)
                space = "\n" + (len(stepname_tmp)+2)*" "
                comm.printinfo("%s: %s" % (stepname_tmp,space.join(example_name_tmp)))
                step1_tmp = ProduceRunDFTStep(stepname_tmp,
                      command = rundft.get("command"),   #command of do dft running
                      example_artifact = example_tmp, 
                      image = image, 
                      collectdata_script = collectdata_script,
                      collectdata_script_name = collectdata_script_name, 
                      collectdata_command = rundft.get("collectdata_command",""), 
                      collectdata_pythonlib = rundft.get("collectdata_pythonlib",[]), 
                      outputs=rundft.get("outputs",[]), 
                      executor = executor, 
                      example_names = example_name_tmp,
                      DoSlices=False,  
                      sub_save_path = sub_save_path,
                      datahub = datahub)
                
                if post_dft and not post_dft_local:
                    step1_tmp.template.outputs.artifacts["outputs"].save = [model_output_artifact]
                    step1_tmp.template.outputs.artifacts["outputs"].archive = None
                else:
                    allstepname.append(stepname_tmp)
                    all_save_path.append([save_path,sub_save_path])
            
                
                rundft_step.append(step1_tmp)
                istep += 1
                
            comm.printinfo("image: %s" % image) 
            comm.printinfo("set bohrium: %s"%str(bohrium_set))   
            comm.printinfo("results will be download to %s\n" % os.path.join(save_path,sub_save_path))

    if len(rundft_step) == 0:
        comm.printinfo("No examples matched in %s, skip it!" % stepname)
        return None,None,None,None
    steps.add(rundft_step)
    
    #post dft
    if post_dft: comm.printinfo("\nPreparing post_dft")
    post_dft_local_jobs = []
    if not post_dft or post_dft_local:
        step3 = rundft_step[0]
        if isinstance(step3,list):
            step3 = step3[0]
        if post_dft_local:
            comm.printinfo("image is '%s', will run this step on local" % str(param["post_dft"].get("image","")))
            post_dft_local_jobs = [save_path,param["post_dft"]]
    else:
        step3 = ProducePostDFTStep(param["post_dft"],model_output_artifact)
        steps.add(step3)
        allstepname = [stepname]
        all_save_path = [[save_path,""]]

    steps.outputs.artifacts['outputs']._from = step3.outputs.artifacts["outputs"]
    step = Step(name=stepname,template=steps)
    return step,allstepname,all_save_path,post_dft_local_jobs

def ProduceAllStep(alljobs):
    allstep = []
    allstepname = []
    allsave_path = []
    postdft_local_jobs = []
    test_name = []
    for k in alljobs:
        step,stepname,save_path,postdft_local_job = ProduceOneSteps(k,alljobs[k])
        if step == None:
            continue  
        allstep.append(step)
        allstepname.append(stepname)
        allsave_path.append(save_path) 
        postdft_local_jobs.append(postdft_local_job) 
        comm.printinfo("\nComplete the preparing for %s.\n" % (k))
        test_name.append(k)

    return allstep,allstepname,allsave_path,postdft_local_jobs,test_name
