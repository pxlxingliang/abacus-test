import json,os,sys,subprocess,glob,shutil
import traceback
from . import comm_echarts
from dp.metadata import MetadataContext
from dp.metadata.utils.storage import TiefblueStorageClient
from dflow import download_artifact,S3Artifact,config,s3_config
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient


def create_path(output_path):
    work_path = "abacustest" #os.path.join(output_path,"abacustest")
    download_path = "abacustest_download" #os.path.join(output_path,"abacustest_download")
    for ipath in [output_path,work_path,download_path]:
        os.makedirs(ipath,exist_ok=True)
    return {"output_path": os.path.abspath(output_path),
            "work_path": os.path.abspath(work_path),
            "download_path": os.path.abspath(download_path)}

def register_dflow(private_set):
    config["host"] = private_set.get("config_host","").strip()  
    s3_config["endpoint"] =  private_set.get("s3_config_endpoint","").strip()
    config["k8s_api_server"] = private_set.get("config_k8s_api_server","").strip()
    config["token"] = private_set.get("config_token","").strip()
    bohrium.config["username"] = private_set.get('lbg_username','')
    bohrium.config["password"] = private_set.get('lbg_password','')
    bohrium.config["project_id"] = private_set.get('project_id','')
    s3_config["repo_key"] = "oss-bohrium"
    s3_config["storage_client"] = TiefblueClient()

def read_config(opts):
    #parse config
    CONFIG_KEYS=["lbg_username","lbg_password","project_id",
                 "config_host","s3_config_endpoint","config_k8s_api_server","config_token",
                 "datahub_project","datahub_gms_token","datahub_gms_url","AIM_ACCESS_TOKEN",
                 "dflow_labels"]
    my_config = {}
    for ikey in CONFIG_KEYS:
        config_key = "Config_" + ikey
        if hasattr(opts,config_key):
            if ikey == "dflow_labels":
                my_config[ikey] = getattr(opts,config_key)
            elif getattr(opts,config_key).strip():
                my_config[ikey] = getattr(opts,config_key).strip()
    register_dflow(my_config)
    return my_config

def run_command(
        cmd,
        shell = True
):
    pp = subprocess.Popen(
        cmd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        shell=shell,
        executable='/bin/bash'
    )
    out = ""
    while True:
        output = pp.stdout.readline()
        return_code = pp.poll()
        if output == b'' and return_code is not None:
            break
        if output:
            print(output.decode()[:-1])
            out += output.decode()
    err = pp.stderr.read().decode()
    return return_code, out, err

def exec_abacustest(allparams,work_path,command = "abacustest submit -p param.json"):
    #work_path: to run abacustest
    #ouput_path: the write the report files
    #write param.json
    json.dump(allparams,
              open(os.path.join(work_path, "param.json"), 'w'),
              indent=4)
    
    cwd = os.getcwd()
    os.chdir(work_path)

    print("Execute abacustest ...")
    return_code, stdout, stderr = run_command(command)
    if return_code != 0:
        print("Error Encountered!")
        #print(stdout)
        print(stderr)
    else:
        #print(stdout)
        print(stderr)
    os.chdir(cwd)
    return stdout,stderr

def produce_metrics_superMetrics_reports(allparams,work_path,output_path):
    from abacustest import outresult
    import pandas as pd
    from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement

    if "save_path" not in allparams:
        allparams["save_path"] = "result"
    reports = []
    chart_section = []
    metric_filename = allparams.get("post_dft",{}).get("metrics",{}).get("save_file","metrics.json")
    metric_file = os.path.join(work_path,allparams["save_path"],metric_filename)
    print(metric_file)
    if os.path.isfile(metric_file):
        try:
            allresults = json.load(open(metric_file))
            csv_filename = os.path.splitext(metric_filename)[0] + ".csv"
            outresult.pandas_out(allresults,os.path.join(output_path,csv_filename))
            reports.append(ReportSection(title="metrics",elements=[AutoReportElement(title='metrics', path=csv_filename,description="")]))

            #produce the Echarts option for each metric
            pddata = pd.DataFrame.from_dict(allresults)
            example_name = pddata.columns.to_list() #get the column name
            metric_name = pddata.index.to_list()  #get the row/index name
            type_set = (int,float,type(None),bool)
            for imetric in metric_name:
                ivalue = pddata.loc[imetric,:].to_list()
                if False not in [isinstance(i,type_set) for i in ivalue]:
                    options = comm_echarts.get_bar_option(imetric,example_name,ivalue)
                    options["xAxis"][0]["axisLabel"] = {
                        "rotate": 15,
                        "interval": int(len(example_name)/15)
                    }
                    chart_section.append(ChartReportElement(options=options,title=imetric))
        except:
            traceback.print_exc()

    customized_metrics_filename = allparams.get("post_dft",{}).get("metrics",{}).get("value_from_file",None)
    if customized_metrics_filename:
        customized_metrics_file = os.path.join(work_path,allparams["save_path"],customized_metrics_filename)
        if os.path.isfile(customized_metrics_file):
            try:
                customized_metrics = json.load(open(customized_metrics_file))
                pddata = pd.DataFrame(customized_metrics,index=[0])
                print(pddata)
                csv_filename = os.path.splitext(customized_metrics_filename)[0] + ".csv"
                pddata.to_csv(os.path.join(output_path,csv_filename),index=False)
                reports.append(ReportSection(title="customized_metrics",elements=[AutoReportElement(title='customized_metrics', path=csv_filename,description="")]))
            except:
                traceback.print_exc()


    super_metric_filename = allparams.get("post_dft",{}).get("super_metrics",[{}])[0].get("save_file","superMetrics.json")
    super_metric_file = os.path.join(work_path,allparams["save_path"],super_metric_filename)
    if os.path.isfile(super_metric_file):
        try:
            supermetrics = json.load(open(super_metric_file))
            pddata = pd.DataFrame(supermetrics,index=[0])
            print(pddata)
            pddata.to_csv(os.path.join(output_path,"superMetrics.csv"),index=False)
            reports.append(ReportSection(title="super_metrics",
                                         elements=[
                                             AutoReportElement(title='super_metrics', path="superMetrics.csv", description="")]))
        except:
            traceback.print_exc()

    
    customized_super_metrics_filename = allparams.get("post_dft",{}).get("super_metrics",[{}])[0].get("value_from_file",None)
    if customized_super_metrics_filename:
        customized_super_metrics_file = os.path.join(work_path,allparams["save_path"],customized_super_metrics_filename)
        if os.path.isfile(customized_super_metrics_file):
            try:
                customized_super_metrics = json.load(open(customized_super_metrics_file))
                pddata = pd.DataFrame(customized_super_metrics,index=[0])
                print(pddata)
                csv_filename = os.path.splitext(customized_super_metrics_filename)[0] + ".csv"
                pddata.to_csv(os.path.join(output_path,csv_filename),index=False)
                reports.append(ReportSection(title="customized_super_metrics",elements=[AutoReportElement(title='customized_super_metrics', path=csv_filename,description="")]))
            except:
                traceback.print_exc()

    if chart_section:
        reports.append(ReportSection(title="metrics chart",elements=chart_section,ncols=2))
    return reports

def produce_html_table(table):
    content = "<table border=\"2px\">"
    for i,it in enumerate(table):
        if i == 0:
            content += "<thead>"
        else:
            content += "<tbody>"

        content += "<tr>"
        for j in it:
            content += "<td>" + str(j) + "</td>"
        content += "</tr>"

        if i == 0:
            content += "</thead>"
        else:
            content += "</tbody>"
    content += "</table>"
    return content

def get_datahub_dataset(bohrium_username,bohrium_password,bohrium_project,urn,download_path=None):
    metadata_storage_client = TiefblueStorageClient(bohrium_username,bohrium_password,bohrium_project)
    with MetadataContext(storage_client=metadata_storage_client) as context:
        dataset = context.client.get_dataset(urn)
        #print(dataset,dataset.uri)
        if dataset != None and download_path != None:
            #context.client.download_dataset(dataset,download_path)
            artifact = S3Artifact(key=dataset.uri)
            download_artifact(artifact,path=download_path)
        return dataset

def unpack(filepath, output_path, filetype = None, get_support_filetype = False):
    if get_support_filetype:
        return ["zip","tar","gz","bz2","tgz"]
    
    import zipfile,tarfile,gzip,bz2

    if not os.path.isfile(filepath):
        raise Exception("File %s not exists!" % filepath)
    
    if not os.path.isdir(output_path):
        os.makedirs(output_path,exist_ok=True)

    if filetype == None:
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
        raise Exception("File type %s not supported!" % filetype)
    
    print("Unpack %s to %s" % (filepath,output_path))
    return output_path

def download_url(url, output_path="./"):
    import requests
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        if not os.path.isdir(output_path):
            os.makedirs(output_path,exist_ok=True)
        filename = os.path.join(output_path,url.split("/")[-1])
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        return filename
    else:
        print(f"dwonload ({url}) failed, status code:", response.status_code)
        return None

def clean_dictorys(ipath):
    for ifile in glob.glob(os.path.join(ipath,"*")):
        if os.path.isdir(ifile):
            shutil.rmtree(ifile)
        else:
            os.remove(ifile)

def move_results_to_output(work_path,output_path,result_folder):
    target_result_folder = os.path.join(output_path,result_folder)
    n = 1
    while os.path.isdir(target_result_folder):
        target_result_folder = os.path.join(output_path,result_folder + "_" + str(n))
        n += 1
    shutil.move(os.path.join(work_path,result_folder),target_result_folder)
