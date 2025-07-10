import json
import os
import subprocess
import glob
import shutil
import select
import copy
#from dp.metadata import MetadataContext
#from dp.metadata.utils.storage import TiefblueStorageClient
from dflow import download_artifact, S3Artifact, config, s3_config
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient


def create_path(output_path):
    work_path = "abacustest"  # os.path.join(output_path,"abacustest")
    # os.path.join(output_path,"abacustest_download")
    download_path = "abacustest_download"
    for ipath in [output_path, work_path, download_path]:
        os.makedirs(ipath, exist_ok=True)
    return {"output_path": os.path.abspath(output_path),
            "work_path": os.path.abspath(work_path),
            "download_path": os.path.abspath(download_path)}


def register_dflow(private_set):
    config["host"] = private_set.get("dflow_host", "").strip()
    s3_config["endpoint"] = private_set.get(
        "dflow_s3_config_endpoint", "").strip()
    config["k8s_api_server"] = private_set.get(
        "dflow_k8s_api_server", "").strip()
    config["token"] = private_set.get("dflow_token", "").strip()
    bohrium.config["username"] = private_set.get('bohrium_username', '')
   # bohrium.config["password"] = private_set.get('bohrium_password','')
    bohrium.config["ticket"] = private_set.get('bohrium_ticket', '')
    bohrium.config["project_id"] = private_set.get('bohrium_project_id', '')
    s3_config["repo_key"] = "oss-bohrium"
    s3_config["storage_client"] = TiefblueClient()


def read_config(opts):
    # parse config
    CONFIG_KEYS = ["bohrium_username", "bohrium_password", "bohrium_ticket", "bohrium_project_id",
                   "dflow_host", "dflow_s3_config_endpoint", "dflow_k8s_api_server", "dflow_token",
                   "aim_access_token", "dflow_labels", "github_username", "github_email", "github_token"]
    my_config = {}
    for ikey in CONFIG_KEYS:
        config_key = "Config_" + ikey
        if hasattr(opts, config_key):
            if ikey == "dflow_labels":
                my_config[ikey] = getattr(opts, config_key)
            elif getattr(opts, config_key).strip():
                my_config[ikey] = getattr(opts, config_key).strip()
    register_dflow(my_config)
    return my_config

def run_command(
        cmd,
        shell=True
):
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
        readable, _, _ = select.select(
            [process.stdout, process.stderr], [], [])

        # 读取已经准备好的输出
        for fd in readable:
            if fd == process.stdout:
                line = process.stdout.readline()
                print(line.decode()[:-1])
                out += line.decode()
            elif fd == process.stderr:
                line = process.stderr.readline()
                print("STDERR:", line.decode()[:-1])
                err += line.decode()

        # 如果子进程已经结束，则退出循环
        return_code = process.poll()
        if return_code is not None:
            break
    return return_code, out, err


def exec_abacustest(allparams, work_path, command="abacustest submit -p param.json"):
    # work_path: to run abacustest
    # ouput_path: the write the report files
    # write param.json
    params = copy.deepcopy(allparams)
    '''
    if "config" in params:
        for ik,iv in params["config"].items():
            os.environ[ik.upper()] = str(iv)
        del params["config"]
    '''
    json.dump(params,
              open(os.path.join(work_path, "param.json"), 'w'),
              indent=4)

    cwd = os.getcwd()
    os.chdir(work_path)

    print("Execute abacustest ...")
    '''
    return_code, stdout, stderr = run_command(command)
    if return_code != 0:
        print("Error Encountered!")
        #print(stdout)
        print(stderr)
    else:
        #print(stdout)
        print(stderr)
    '''
    import abacustest.abacustest as abatest
    import argparse
    from abacustest.main import print_head
    print_head()
    parser = argparse.ArgumentParser(description="abacustest")
    subparser = parser.add_subparsers(dest="command")
    abatest.AbacusTestArgs(subparser.add_parser("submit"))
    parser = parser.parse_args(["submit", "-p", "param.json"])
    abatest.abacustest(parser)

    os.chdir(cwd)
    stdout, stderr = "", ""
    return stdout, stderr


def produce_html_table(table):
    content = "<table border=\"2px\">"
    for i, it in enumerate(table):
        if i == 0:
            content += "<thead>"
        else:
            content += "<tbody>"

        content += "<tr>"
        if len(it) == 1:
            content += f"<td colspan=\"{len(table[0])}\">" + \
                str(it[0]) + "</td>"
        else:
            for j in it:
                content += "<td>" + str(j) + "</td>"
        content += "</tr>"

        if i == 0:
            content += "</thead>"
        else:
            content += "</tbody>"
    content += "</table>"
    return content

'''
def get_datahub_dataset(bohrium_username, bohrium_password, bohrium_project, urn, download_path=None):
    metadata_storage_client = TiefblueStorageClient(
        bohrium_username, bohrium_password, bohrium_project)
    with MetadataContext(storage_client=metadata_storage_client) as context:
        dataset = context.client.get_dataset(urn)
        # print(dataset,dataset.uri)
        if dataset != None and download_path != None:
            # context.client.download_dataset(dataset,download_path)
            artifact = S3Artifact(key=dataset.uri)
            download_artifact(artifact, path=download_path)
        return dataset
'''

def unpack(filepath, output_path, filetype=None, get_support_filetype=False):
    # if file type is not supportted, just copy the file to output_path
    if get_support_filetype:
        return ["zip", "tar", "gz", "bz2", "tgz"]

    import zipfile
    import tarfile
    import gzip
    import bz2

    if not os.path.isfile(filepath):
        raise Exception("File %s not exists!" % filepath)

    if not os.path.isdir(output_path):
        os.makedirs(output_path, exist_ok=True)

    if filetype == None:
        if filepath.endswith(".tar.gz"):
            filetype = "tgz"
        else:
            filetype = os.path.splitext(filepath)[1][1:]

    if filetype == "zip":
        with zipfile.ZipFile(filepath, 'r') as zip_ref:
            zip_ref.extractall(output_path)
    elif filetype == "tar":
        with tarfile.open(filepath, "r") as tar_ref:
            tar_ref.extractall(output_path)
    elif filetype == "gz":
        with gzip.open(filepath, 'rb') as gz_ref:
            with open(output_path, 'wb') as out_ref:
                out_ref.write(gz_ref.read())
    elif filetype == "bz2":
        with bz2.open(filepath, 'rb') as bz2_ref:
            with open(output_path, 'wb') as out_ref:
                out_ref.write(bz2_ref.read())
    elif filetype == "tgz":
        with tarfile.open(filepath, "r:gz") as tar_ref:
            tar_ref.extractall(output_path)
    else:
        print("File type \'%s\' is not supported!" % filetype)
        filename = os.path.split(filepath)[-1]
        if os.path.isfile(os.path.join(output_path, filename)):
            os.remove(os.path.join(output_path, filename))
        shutil.copyfile(filepath, os.path.join(output_path, filename))

    print("Unpack %s to %s" % (filepath, output_path))
    return output_path


def pack(packfile_list, packfile_name, pack_type="zip"):
    if pack_type == "zip":
        import zipfile
        with zipfile.ZipFile(packfile_name, 'w') as zip_ref:
            for ifile in packfile_list:
                if os.path.isfile(ifile):
                    zip_ref.write(ifile)
                else:
                    for root, __, ifile in os.walk(ifile):
                        for i in ifile:
                            zip_ref.write(os.path.join(root, i))
    elif pack_type == "tar":
        import tarfile
        with tarfile.open(packfile_name, "w") as tar_ref:
            for ifile in packfile_list:
                tar_ref.add(ifile)
    elif pack_type == "gz":
        import gzip
        with gzip.open(packfile_name, 'wb') as gz_ref:
            for ifile in packfile_list:
                with open(ifile, 'rb') as in_ref:
                    gz_ref.write(in_ref.read())
    elif pack_type == "bz2":
        import bz2
        with bz2.open(packfile_name, 'wb') as bz2_ref:
            for ifile in packfile_list:
                with open(ifile, 'rb') as in_ref:
                    bz2_ref.write(in_ref.read())
    elif pack_type == "tgz":
        import tarfile
        with tarfile.open(packfile_name, "w:gz") as tar_ref:
            for ifile in packfile_list:
                tar_ref.add(ifile)
    else:
        print("Pack type \'%s\' is not supported!" % pack_type)
        return None
    return packfile_name


def download_url(url, output_path="./"):
    import requests
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        if not os.path.isdir(output_path):
            os.makedirs(output_path, exist_ok=True)
        filename = os.path.join(output_path, url.split("/")[-1])
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        return filename
    else:
        print(f"dwonload ({url}) failed, status code:", response.status_code)
        return None


def clean_dictorys(ipath):
    for ifile in glob.glob(os.path.join(ipath, "*")):
        if os.path.isdir(ifile):
            shutil.rmtree(ifile)
        else:
            os.remove(ifile)


def move_results_to_output(work_path, output_path, result_folder):
    target_result_folder = os.path.join(output_path, result_folder)
    if os.path.abspath(target_result_folder) != os.path.abspath(output_path):
        n = 1
        while os.path.isdir(target_result_folder):
            target_result_folder = os.path.join(
                output_path, result_folder + "_" + str(n))
            n += 1
    
    src_folder = os.path.join(work_path, result_folder)
    if os.path.abspath(src_folder) != os.path.abspath(work_path):
        shutil.move(os.path.join(work_path, result_folder), target_result_folder)
    else:
        # move all files in work_path to target_result_folder
        for ifile in glob.glob(os.path.join(work_path, "*")):
            shutil.move(ifile, target_result_folder)

def pack_results(output_path,result_path):
    cwd = os.getcwd()
    os.chdir(os.path.join(output_path,result_path))
    allfiles = os.listdir(".")
    alldirs = []
    allfile = []
    for f in allfiles:
        if os.path.isdir(f):
            alldirs.append(f)
        elif os.path.isfile(f) and not f.endswith(".zip") and not f.startswith("."):
            allfile.append(f)
    packed_file_name = "results.zip"
    pack(alldirs,packed_file_name,"zip")
    pack(allfile,"files.zip","zip")
    os.chdir(cwd)

def gen_dir(dir1):
    dir1 = dir1.rstrip("/")
    n = 0
    while os.path.isdir(dir1):
        dir1 = dir1 + "_" + str(n)
        n += 1
    return dir1

def capture_output(func, *args, **kwargs):
    """ Capture the output of a function and return it as a string.
    Args:
        func (callable): The function whose output is to be captured.
        *args: Positional arguments to pass to the function.
        **kwargs: Keyword arguments to pass to the function.
    Returns:
        str: The captured output of the function.
    """
    import sys
    from io import StringIO
    if not callable(func):
        raise ValueError("The provided func must be callable.")
    
    old_stdout = sys.stdout
    redirected_output = sys.stdout = StringIO()
    
    try:
        return_values = func(*args, **kwargs)
        output = redirected_output.getvalue()
    finally:
        sys.stdout = old_stdout
    
    return return_values, output
    