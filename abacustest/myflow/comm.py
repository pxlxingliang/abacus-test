from . import globV
import os,shutil,glob,json
from pathlib import Path

def printinfo(istr,*args):
    LOGFILE = "abacustest.log"
    output = " ".join([str(istr)]+[str(i) for i in args])
    with open(LOGFILE,'a+') as f1:
        f1.write(output + "\n")
    print(output,flush=True)
        
def GetBakFile(sfile):
    while sfile[-1] == '/':
        sfile = sfile[:-1]
    n = 1
    bk = sfile + ".bak%d" % n
    while os.path.exists(bk):
        n += 1
        bk = sfile + ".bak%d" % n
    return bk

def CopyFiles(path1,path2,move = False):
    '''copy the files in path1 to path2'''
    abspath1 = os.path.abspath(path1)
    abspath2 = os.path.abspath(path2)
    
    if not os.path.isdir(abspath2):
        os.makedirs(abspath2)
    if not os.path.isdir(abspath1):
        print(f"{abspath1} is not a exist path")
        return
    
    def CopyFile(old_path,new_path):
        for ifile in os.listdir(old_path):
            iold_path = os.path.join(old_path,ifile)
            if os.path.isfile(iold_path):
                shutil.copy(iold_path,new_path)
            else:
                if not os.path.isdir(os.path.join(new_path,ifile)):
                    os.makedirs(os.path.join(new_path,ifile))
                CopyFile(iold_path,os.path.join(new_path,ifile))

    if abspath2.startswith(abspath1):
        '''
        If path2 is a son path of path1,
        we will firstly create a tmp path, and copy/move
        files in path1 to the tmp path, and then replace tmp path
        to path2.
        '''
        tmp_path = GetBakFile(os.path.split(abspath1)[0])
        os.makedirs(tmp_path)
        if move:
            for i in os.listdir(abspath1):
                shutil.move(os.path.join(abspath1,i),tmp_path)
        else:
            #shutil.copytree(abspath1,tmp_path,dirs_exist_ok=True)
            CopyFile(abspath1,tmp_path)
               
        #shutil.copytree(tmp_path,abspath2,dirs_exist_ok=True)
        CopyFile(tmp_path,abspath2)
        shutil.rmtree(tmp_path)
        
    else:
        if os.path.isfile(abspath1):
            if move:
                shutil.move(abspath1,abspath2)
            else:
                shutil.copy(abspath1,abspath2)
            return
        elif os.path.isdir(abspath1):
            if move:
                for i in os.listdir(abspath1):
                    shutil.move(os.path.join(abspath1,i),abspath2)
            else:
                #shutil.copytree(abspath1,abspath2,dirs_exist_ok=True)
                CopyFile(abspath1,abspath2)

def LinkFiles(src_list, dst_path):
    # the src_list should be a list of relative path, and the dst_path should be a absolute path
    # will create the soft link in dst_path with the same path tree as src_list
    for src in src_list:
        src = src.rstrip("/")
        src_root, src_name = os.path.split(src)   
        if os.path.isfile(os.path.join(dst_path,src)):
            os.remove(os.path.join(dst_path,src))
        elif os.path.isdir(os.path.join(dst_path,src)):
            shutil.rmtree(os.path.join(dst_path,src))
                        
        if not os.path.exists(os.path.join(dst_path,src_root)):
            os.makedirs(os.path.join(dst_path,src_root))
        os.symlink(os.path.abspath(src),os.path.join(dst_path,src))           

def CollectFileName(paths):
    #Recursively find all files under the paths
    allfiles = []
    for ipath in glob.glob(os.path.join(paths,"*")):
        if os.path.isfile(ipath):
            allfiles.append(ipath)
        elif os.path.isdir(ipath):
            allfiles += CollectFileName(ipath)
    return allfiles            

def ProduceExecutor(param,group_name="abacustesting"):
    from dflow.plugins.bohrium import BohriumContext, BohriumExecutor
    from dflow.plugins.dispatcher import DispatcherExecutor
    from dflow.plugins.bohrium import TiefblueClient,create_job_group
    
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

        bohrium_set["bohr_job_group_id"] = create_job_group(group_name)
        
        if not globV.get_value("BOHRIUM_EXECUTOR"):    
            dispatcher_executor = DispatcherExecutor(
                machine_dict={
                    "batch_type": "Bohrium",
                    "context_type": "Bohrium",
                    "remote_profile": {"input_data": bohrium_set},
                    },
                image_pull_policy = "Always",
                retry_on_submission_error=3,
                #image=param["bohrium"].get("dispatcher_image","registry.dp.tech/public/dptechnology/dpdispatcher:v0.6.0"),
            )
            if "dispatcher_image" in param["bohrium"]:
                dispatcher_executor.image = param["bohrium"]["dispatcher_image"]
            #comm.printinfo("set bohrium: %s"%str(bohrium_set))
            return dispatcher_executor,bohrium_set
        else:
            executor = BohriumExecutor(
                executor= "bohrium_v2",
                extra={
                    "scassType": bohrium_set['scass_type'],
                    "platform": bohrium_set['platform'] ,
                    "projectId": globV.get_value("PRIVATE_SET").get("bohrium_project_id"),
                    "jobType":  bohrium_set['job_type']
                }
            )
            return executor,bohrium_set
    elif "dispatcher" in param and param["dispatcher"]:
        '''
        host: remote host
        queue_name: queue name
        port: SSH port
        username: username
        private_key_file: private key file for SSH
        image: image for dispatcher
        command: command for dispatcher
        remote_command: command for running the script remotely
        map_tmp_dir: map /tmp to ./tmp
        machine_dict: machine config for dispatcher
        resources_dict: resources config for dispatcher
        task_dict: task config for dispatcher
        json_file: JSON file containing machine and resources config


        dispatcher_executor = DispatcherExecutor(
            host = param["dispatcher"].get("host",None),
            queue_name = param["dispatcher"].get("queue_name",None),
            port = param["dispatcher"].get("port",22),
            username = param["dispatcher"].get("username","root"),
            private_key_file = param["dispatcher"].get("private_key_file",None),
        #    image = image,
        #    command = command,
        #    remote_command = remote_command,
        #    map_tmp_dir = map_tmp_dir,
            machine_dict = param["dispatcher"].get("machine_dict",None),
            resources_dict = param["dispatcher"].get("resources_dict",None),
            task_dict = param["dispatcher"].get("task_dict",None),
            json_file = param["dispatcher"].get("json_file",None)
            )
        '''
        dispatcher_executor = DispatcherExecutor(**param["dispatcher"])
        import copy
        tmp_param = copy.deepcopy(param["dispatcher"])
        hide_config_in_dispatcher(tmp_param)
        return dispatcher_executor,tmp_param
    else:
        return None,None

def hide_config_in_dispatcher(tmp_param):
    # tmp_param is the value of dispatcher
    if "host" in tmp_param:
        tmp_param["host"] = "******"
    if "username" in tmp_param:
        tmp_param["username"] = "******"
    if "port" in tmp_param:
        tmp_param["port"] = "******"
    if "private_key_file" in tmp_param:
        tmp_param["private_key_file"] = "******"
    if "machine_dict" in tmp_param and "remote_profile" in tmp_param["machine_dict"]:
        if "hostname" in tmp_param["machine_dict"]["remote_profile"]:
            tmp_param["machine_dict"]["remote_profile"]["hostname"] = "******"
        if "username" in tmp_param["machine_dict"]["remote_profile"]:
            tmp_param["machine_dict"]["remote_profile"]["username"] = "******"
        if "password" in tmp_param["machine_dict"]["remote_profile"]:
            tmp_param["machine_dict"]["remote_profile"]["password"] = "******"
        if "port" in tmp_param["machine_dict"]["remote_profile"]:
            tmp_param["machine_dict"]["remote_profile"]["port"] = "******"

def FindLocalExamples_new(example,only_folder=False,oneartifact=False):
    from dflow import upload_artifact
    #use glob.glob find all examples, and transfer to artifact
    #example = [[*],*]    
    examples_name = []  
    examples_name1 = []
    for i in example:
        if isinstance(i,list):
            example_tmp = []
            for j in i:
                example_tmp += glob.glob(j)
            if len(example_tmp) > 0:
                example_tmp.sort()
                tmp = []
                for ii in example_tmp:
                    # if ii is __MACOSX, skip it
                    if os.path.split(ii)[-1].startswith("__MACOSX"):
                        continue
                    if only_folder and not os.path.isdir(ii):
                        continue
                    tmp.append(ii)
                    examples_name1.append(ii)
                if len(tmp) > 0:
                    tmp.sort()
                    examples_name.append(tmp)
        elif isinstance(i,str):
            for ii in glob.glob(i):
                if os.path.split(ii)[-1].startswith("__MACOSX"):
                    continue
                if only_folder and not os.path.isdir(ii):
                        continue
                examples_name.append([ii])
                examples_name1.append(ii)
        else:
            printinfo(i,"element of 'example' should be a list, or str")
    
    examples = None
    if len(examples_name) > 0:  
        if oneartifact:
            examples_name = list(set(examples_name1))
            examples_name.sort()
            examples = [[upload_artifact(examples_name,archive=globV.get_value("COMPRESS"))]]
            examples_name = [examples_name]
        else:
            examples = []
            examples_name.sort()
            for ii in examples_name:
                examples.append([upload_artifact(iii,archive=globV.get_value("COMPRESS")) for iii in ii])    
    
    # if oneartifact is True, examples = [[artifact]], examples_name = [[example1,example2,...]]
    # else examples = [[artifact1,artifact2,...],[artifact3,artifact4,...],...], examples_name = [[example1,example2,...],[example3,example4,...],...]
    return examples,examples_name

def transfer_source_to_artifact(example,source=None,source_type="local",only_folder=True,oneartifact=False):
    # example specify the example name
    # if onearitfact is True, then the source will be transfered to one artifact
    # else, the source will be transfered to a list of artifacts (for rundft)
    # for rundft, need set only_folder = True
    if isinstance(example,str):
        example = [example]
    
    examples,examples_name = [],[]   
     
    if source_type == "local":
        examples,examples_name = FindLocalExamples_new(example,only_folder=only_folder,oneartifact=oneartifact)
    
    # examples = [[artifact1,artifact2,...],[artifact3,artifact4,...],...]
    # examples_name = [[example1,example2,...],[example3,example4,...],...]
    return examples,examples_name

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

def SplitGroupSize(examples,examples_name,group_size):
    #examples = [*[*],*]
    newexamples = [[]] 
    newexamples_name = [[]]
    
    i = 0
    j = 0
    while i < len(examples):
        newexamples[-1] += examples[i]
        newexamples_name[-1] += examples_name[i]
        j += 1
        i += 1
        if i >= len(examples):
            break
        if j >= group_size:
            j = 0
            newexamples.append([])
            newexamples_name.append([])
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
        printinfo("the type of 'sub_save_path' should be 'str', but not '%s'. %s" % (type(sub_save_path),str(sub_save_path)))
        sub_save_path = ""
    return sub_save_path 

def SetEnvs():
    if globV.get_value("PRIVATE_SET",None) == None:
        return None
    from dflow import Secret
    envs = {}
    for k,v in globV.get_value("PRIVATE_SET").items():
        if isinstance(v,dict):
            envs[k.upper()] = Secret("'" + json.dumps(v) + "'")
        else:
            envs[k.upper()] = Secret(str(v))
    return envs

def run_command(
        cmd,
        shell = True
):
    import subprocess,select
    process = subprocess.Popen(
        cmd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        shell=shell,
        executable='/bin/bash'
    )
    out = ""
    err = ""
    while True:
        # 监视stdout和stderr文件描述符的可读状态
        readable, _, _ = select.select([process.stdout, process.stderr], [], [])

        # 读取已经准备好的输出
        for fd in readable:
            if fd == process.stdout:
                line = process.stdout.readline()
                print(line.decode()[:-1])
                out += line.decode()
            elif fd == process.stderr:
                line = process.stderr.readline()
                if line.decode()[:-1].strip() != "":
                    print("STDERR:", line.decode()[:-1])
                err += line.decode()

        # 如果子进程已经结束，则退出循环
        return_code = process.poll()
        if return_code is not None:
            break
    return return_code, out, err

