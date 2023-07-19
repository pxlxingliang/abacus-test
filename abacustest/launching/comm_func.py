import json,os,sys,subprocess,glob,shutil,select
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
   # bohrium.config["password"] = private_set.get('lbg_password','')
    bohrium.config["ticket"] = private_set.get('bohrium_ticket','')
    bohrium.config["project_id"] = private_set.get('project_id','')
    s3_config["repo_key"] = "oss-bohrium"
    s3_config["storage_client"] = TiefblueClient()

def read_config(opts):
    #parse config
    CONFIG_KEYS=["lbg_username",
                 #"lbg_password",
                 "bohrium_ticket",
                 "project_id",
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
                print("STDERR:", line.decode()[:-1])
                err += line.decode()

        # 如果子进程已经结束，则退出循环
        return_code = process.poll()
        if return_code is not None:
            break
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
    parser = argparse.ArgumentParser(description="abacustest")
    subparser = parser.add_subparsers(dest="command")
    abatest.AbacusTestArgs(subparser.add_parser("submit"))
    parser = parser.parse_args(["submit","-p","param.json"])
    abatest.abacustest(parser)
    
    os.chdir(cwd)
    stdout,stderr = "",""     
    return stdout,stderr

def produce_metrics(metric_file,output_path,report_titile="metrics"):
    from abacustest import outresult
    import pandas as pd
    from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement
    
    metric_filename = os.path.split(metric_file)[-1]
    allresults = json.load(open(metric_file))
    csv_filename = os.path.splitext(metric_filename)[0] + ".csv"
    _, _, savefile_names = outresult.pandas_out(allresults,os.path.join(output_path,csv_filename))
    report_elements = []
    for ifilename in savefile_names:
        ifilename = os.path.split(ifilename)[-1]
        report_elements.append(AutoReportElement(title=os.path.splitext(ifilename)[0], path=ifilename,description=""))
    
    chart_elements = []
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
            chart_elements.append(ChartReportElement(options=options,title=imetric))
    
    return report_elements,chart_elements

def produce_supermetrics(supermetric_file,output_path,work_path, save_path,report_titile="supermetrics"):
    import pandas as pd
    from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement
    
    supermetrics_filename = os.path.split(supermetric_file)[-1]
    super_metrics = json.load(open(supermetric_file))
    normal_metrics = {}
    special_metrics = {}
    for ikey,ivalue in super_metrics.items():
        if isinstance(ivalue,dict):
            special_metrics[ikey] = ivalue
        else:
            normal_metrics[ikey] = ivalue
    
    report_element = None
    if normal_metrics:
        pddata = pd.DataFrame(normal_metrics,index=[0])
        print(pddata)
        csv_filename = os.path.splitext(supermetrics_filename)[0] + ".csv"
        pddata.to_csv(os.path.join(output_path,csv_filename),index=False)
        report_element = AutoReportElement(title=report_titile, path=csv_filename,description="")
    
    special_section = None   
    special_elements = [] 
    if special_metrics:
        for ikey,ivalue in special_metrics.items():
            file_name = ivalue.get("file",None)
            if file_name:
                ifile_name = os.path.join(work_path,save_path,file_name)
                if os.path.isfile(ifile_name):
                    special_elements.append(AutoReportElement(title=ikey, path=os.path.join(save_path,file_name),description=""))
                    
    if special_elements:
        special_section = ReportSection(title=report_titile,elements=special_elements)
        
    return report_element,special_section

def produce_metrics_superMetrics_reports(allparams,work_path,output_path):
    from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement

    save_path = allparams.get("save_path","result")
    reports = []
    allmetrics_files = []
    allsupermetrics_files = []
    
    #produce the metrics reports
    # 1. metrics from the custome defined metrics file
    customized_metrics_filename = allparams.get("post_dft",{}).get("metrics",{}).get("value_from_file",None)
    if customized_metrics_filename:
        allmetrics_files.append(os.path.join(work_path,save_path,customized_metrics_filename)) 
    
    # 2. metrics from the metrics file
    metric_filename = allparams.get("post_dft",{}).get("metrics",{}).get("save_file")
    if metric_filename:
        allmetrics_files.append(os.path.join(work_path,save_path,metric_filename))
    
    # 3. metrics from the undefined metrics file
    allmetrics_files += glob.glob(os.path.join(work_path,save_path,"metric*.json"))
    metrics_report = []
    metrics_chart_elements = []
    for metric_file in list(set(allmetrics_files)): 
        try:
            metric_filename = os.path.split(metric_file)[-1]
            tmp_report_elements,tmp_chart_elements = produce_metrics(metric_file,output_path,report_titile=metric_filename)
            if tmp_report_elements:
                metrics_report += tmp_report_elements
            if tmp_chart_elements:
                metrics_chart_elements += tmp_chart_elements
        except:
            traceback.print_exc()
            print(f"Error: report metrics from file \"{metric_file}\" failed!") 
    
    #produce the supermetrics reports
    # 1. supermetrics from the custome defined supermetrics file        
    super_metric_filename = allparams.get("post_dft",{}).get("super_metrics",[{}])[0].get("save_file")
    if super_metric_filename:
        allsupermetrics_files.append(os.path.join(work_path,save_path,super_metric_filename))
    
    # 2. supermetrics from the supermetrics file
    customized_super_metrics_filename = allparams.get("post_dft",{}).get("super_metrics",[{}])[0].get("value_from_file",None)
    if customized_super_metrics_filename:
        allsupermetrics_files.append(os.path.join(work_path,save_path,customized_super_metrics_filename))   
    
    # 3. supermetrics from the undefined supermetrics file
    allsupermetrics_files += glob.glob(os.path.join(work_path,save_path,"superMetric*.json")) + \
                             glob.glob(os.path.join(work_path,save_path,"super_metric*.json")) + \
                             glob.glob(os.path.join(work_path,save_path,"supermetric*.json"))
    supermetrics_report = []
    supermetrics_special_section = []
    for super_metric_file in list(set(allsupermetrics_files)):
        try:
            super_metric_filename = os.path.split(super_metric_file)[-1]
            tmp_report_element,tmp_special_section = produce_supermetrics(super_metric_file,output_path,work_path,save_path,report_titile=super_metric_filename)
            if tmp_report_element:
                supermetrics_report.append(tmp_report_element)
            if tmp_special_section:
                supermetrics_special_section.append(tmp_special_section)
        except:
            traceback.print_exc()
            print(f"Error: report supermetrics from file \"{super_metric_file}\" failed!")
    
    #produce the report
    # 1. metrics report
    if metrics_report:
        reports.append(ReportSection(title="metrics",elements=metrics_report,ncols=1))
    
    # 2. supermetrics report
    if supermetrics_report:
        reports.append(ReportSection(title="supermetrics",elements=supermetrics_report,ncols=1))
        
    # 3. metrics chart
    if metrics_chart_elements:
        reports.append(ReportSection(title="metrics chart",elements=metrics_chart_elements,ncols=2))
        
    # 4. supermetrics special section
    if supermetrics_special_section:
        for i in supermetrics_special_section:
            reports.append(i)
 
    return reports

def produce_html_table(table):
    content = "<table border=\"2px\">"
    for i,it in enumerate(table):
        if i == 0:
            content += "<thead>"
        else:
            content += "<tbody>"

        content += "<tr>"
        if len(it) == 1:
            content += f"<td colspan=\"{len(table[0])}\">" + str(it[0]) + "</td>"
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
    # if file type is not supportted, just copy the file to output_path
    if get_support_filetype:
        return ["zip","tar","gz","bz2","tgz"]
    
    import zipfile,tarfile,gzip,bz2

    if not os.path.isfile(filepath):
        raise Exception("File %s not exists!" % filepath)
    
    if not os.path.isdir(output_path):
        os.makedirs(output_path,exist_ok=True)

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
        if os.path.isfile(os.path.join(output_path,filename)):
            os.remove(os.path.join(output_path,filename))
        shutil.copyfile(filepath,os.path.join(output_path,filename))
    
    print("Unpack %s to %s" % (filepath,output_path))
    return output_path

def pack(packfile_list,packfile_name,pack_type="zip"):
    if pack_type == "zip":
        import zipfile
        with zipfile.ZipFile(packfile_name, 'w') as zip_ref:
            for ifile in packfile_list:
                for root,__,ifile in os.walk(ifile):
                    for i in ifile:
                        zip_ref.write(os.path.join(root,i))
    elif pack_type == "tar":
        import tarfile
        with tarfile.open(packfile_name, "w") as tar_ref:
            for ifile in packfile_list:
                tar_ref.add(ifile)
    elif pack_type == "gz":
        import gzip
        with open(packfile_name, 'wb') as out_ref:
            with gzip.open(packfile_name, 'wb') as gz_ref:
                gz_ref.write(out_ref.read())
    elif pack_type == "bz2":
        import bz2
        with open(packfile_name, 'wb') as out_ref:
            with bz2.open(packfile_name, 'wb') as bz2_ref:
                bz2_ref.write(out_ref.read())
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
