import os,json
import subprocess

BOHRIUM_DES = '''If you use Bohrium to accelerate the calculation, you need to set below environment variables:
    export BOHRIUM_USERNAME=<your username> BOHRIUM_PASSWORD=<your password> BOHRIUM_PROJECT_ID=<your project id> 
Or you can add below block in your setting.json:
    "config": {
        "bohrium_username": "<your username>",
        "bohrium_password": "<your password>",
        "bohrium_project_id": "<your project id>"
    } 
Or you can use dispatcher to submit the calculation to the remote server,
and please refer to the dispatcher section for more details in https://github.com/pxlxingliang/abacus-test/tree/develop.'''

def doc_after_prepare(model_name, jobs, files,has_prepare=True):
    print(f"\nPrepare the {model_name} successfully!!!")
    print("\nYou can choose one of below methods to do the calculation:")
    print("1. run with abacustest locally")
    print(BOHRIUM_DES,)
    print("After finishing the setting, you can run below command to submit the calculation:")
    print("        abacustest submit -p setting.json &")
    print("\n2. run with bohrium app abacustest:")
    print("Compress the inputs and submit to https://app.bohrium.dp.tech/abacustest/?request=GET%3A%2Fapplications%2Fabacustest with Reuse Submodel.")
    print("Compress command:")
    print("    tar -zcvf %s.tar.gz %s" % (model_name, " ".join(jobs+files)))
    if has_prepare:
        print("\n3. run with your own script:")
        print("You can execute below command and abacustest will generate all inputs in folder abacustest/, and then you can submit the calculation by yourself.")
        outfolder = "" if len(jobs) > 1 else f"-s {jobs[0]}"
        print("    abacustest prepare -p setting.json " + outfolder + "\n\n")
    
def dump_setting(setting):
    if os.path.isfile("setting.json"):
        i = 1
        settingbakf = f"setting.json.bak{i}"
        while os.path.isfile(settingbakf):
            i += 1
            settingbakf = f"setting.json.bak{i}"
        print("Warning: setting.json exists. Will bakup it to", settingbakf)
        os.system(f"mv setting.json {settingbakf}")
    json.dump(setting, open("setting.json", "w"), indent=4)

def get_physical_cores():
    cmd = "lscpu | grep 'Core(s) per socket:' | awk '{print $4}'"
    cores = subprocess.check_output(cmd, shell=True).decode().strip()
    return int(cores)