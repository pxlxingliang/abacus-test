import os,shutil,copy,argparse,sys,json,time

def gen_special_value():
    sv = {}
    date = time.strftime("%Y%m%d", time.localtime())
    sv["DATE"] = date
    return sv

def prepare_files(settings):
    # prepare the files, and then push to remote like github
    # 1. create the folder remote
    # 2. copy the files/folders to the folder
    # 3. support some special format:
    #   - {DATE}: the date of today, like 20210801
    '''
    "remote":{
        "files":{
            "file1_local": "file1_remote",
            "file2_local": "file2_remote",
            "file3_local": "/bda/{DATE}/fiel3_remote"
        }
    }
    '''
    if not settings:
        return None
    
    remote_folder = "remote"
    n = 1
    while os.path.isdir(remote_folder):
        remote_folder += str(n)
        n += 1
    os.mkdir(remote_folder)
    
    files_setting = settings.get("files",{})
    SV = gen_special_value()
    for ifile,jfile in files_setting.items():
        jfile = jfile.format(**SV)
        if not os.path.exists(ifile):
            print(f"Error: {ifile} does not exist!")
            continue
        if os.path.isfile(ifile):
            jdir = os.path.join(remote_folder,os.path.dirname(jfile))
            if not os.path.isdir(jdir):
                os.makedirs(jdir)
            shutil.copy(ifile,os.path.join(remote_folder,jfile))
        elif os.path.isdir(ifile):
            shutil.copytree(ifile,os.path.join(remote_folder,jfile),dirs_exist_ok=True)
    return remote_folder

def push_to_github(remote_path,github_setting,config_setting,extra_setting={}):
    # github_setting = 
    '''
        "github":{
            "repo": "",
            "branch": "",
            "commit_msg": ""
        }
        "add_date": "bda/dates"
    '''
    repo = github_setting.get("repo")
    branch = github_setting.get("branch")
    username = config_setting.get("github_username")
    email = config_setting.get("github_email")
    token = config_setting.get("github_token")
    commit_msg = github_setting.get("commit_msg","update from abacustest")
    
    if not repo:
        print("Error: the repo is not specified!")
        return 1
    if not branch:
        print("Error: the branch is not specified!")
        return 1
    if not username and not email and not token:
        print("Error: the username/email/token is not specified in config!")
        return 1

    remote_path = os.path.abspath(remote_path)
    
    github_path = "github"
    n = 1
    while os.path.isdir(github_path):
        github_path = "github" + str(n)
        n += 1
    os.mkdir(github_path)
    
    # add date to dates
    add_date_path = extra_setting.get("add_date")
    add_date_command = ""
    if add_date_path:
        date = time.strftime("%Y%m%d", time.localtime())
        add_date_command = f"""
# check if dates file exists
if [ -f "{add_date_path}" ]; then
    # check if the date exists
    if ! grep -q "{date}" "{add_date_path}"; then
        echo "{date}" >> "{add_date_path}"
    fi
else
    echo "{date}" > "{add_date_path}"
fi
        """
    
    # gen bash script
    repo_url_with_token=f"https://oauth2:{token}@github.com/{repo}.git"    
    bash_script = f'''#!/bin/bash
wget --no-proxy https://dp-public.oss-cn-beijing.aliyuncs.com/k8s/speedUp.sh -O speedUp.sh && bash speedUp.sh >> start.log
cd {github_path}
git init
git remote add origin "{repo_url_with_token}"

# Configure git user
git config user.name "{username}"
git config user.email "{email}"

# configure git debug mode
set GIT_TRACE=true
set GIT_CURL_VERBOSE=true
set GIT_SSH_COMMAND=ssh -vvv

# Checkout the repository
# try 10 times
fetch=0
for i in {{1..10}}
do
    git fetch --depth 1 origin {branch}
    if [ $? -eq 0 ]; then
        fetch=1
        break
    fi
    sleep 10
done
if [ $fetch -eq 0 ]; then
    echo "Error: can not fetch the branch {branch}!"
    cd ..
    exit 1
fi
git checkout {branch}

# copy the files
cp -r {remote_path}/* .

# check if add date
{add_date_command}

# commit and push
git add .
git commit -m "{commit_msg}"

# try 10 times
push=0
for i in {{1..10}}
do

    git push origin {branch}
    if [ $? -eq 0 ]; then
        push=1
        break
    fi  
    sleep 10
done
if [ $push -eq 0 ]; then
    echo "Error: can not push the branch {branch}!"
    cd ..
    exit 1
fi
cd ..
'''
    pushf = "push.sh"
    n=1
    while os.path.isfile(pushf):
        pushf = f"push{n}.sh"
        n+=1
    open(pushf,"w").write(bash_script)
    os.system(f"bash {pushf}")
    shutil.rmtree(github_path)
    os.remove(pushf)
    return 0

def push_to_remote(remote_path,setting,config_setting): 
    if setting.get("github"):
        push_to_github(remote_path,setting.get("github"),config_setting,setting)
    else:
        print("Error: the remote setting is not supported!")
        return 1

def ReportArgs(parser):  
    parser.description = '''push the local files to remote, like github
    Need to specify the parameter file, which should be .json type, and the content should be like:
    "config":{
        "github_username": "",
        "github_email": "",
        "github_token": ""
    },
    "add_date": "bda/dates",
    "remote":{
        "files":{
            "file1_local": "file1_remote",
            "file2_local": "file2_remote",
        },
        "github":{
            "repo": "",
            "branch": "",
            "commit_msg": ""
        }
    }
    the key of "files" is the local file/folder, and the value is the remote file/folder.
    '''
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type',required=True)
    return parser

def Remote(param):
    param_file = param.param
    if not os.path.exists(param_file):
        print(f"Error: {param_file} does not exist!")
        sys.exit(1)
    settings = json.load(open(param_file))
    config_setting = settings.get("config",{})
    remote_setting = settings.get("remote",{})

    if remote_setting:
        remote_path = prepare_files(remote_setting)
        if remote_path:
            push_to_remote(remote_path,remote_setting,config_setting) 

def main():
    parser = argparse.ArgumentParser()
    param = ReportArgs(parser).parse_args()
    Remote(param)
    
if __name__ == "__main__":
    main()