
import json,os,sys,subprocess

def create_path(output_path):
    work_path = os.path.join(output_path,"abacustest")
    download_path = os.path.join(output_path,"abacustest_download")
    for ipath in [output_path,work_path,download_path]:
        os.makedirs(ipath,exist_ok=True)
    return {"output_path": os.path.abspath(output_path),
            "work_path": os.path.abspath(work_path),
            "download_path": os.path.abspath(download_path)}

def read_config(opts):
    #parse config
    print("read config setting")
    CONFIG_KEYS=["lbg_username","lbg_password","project_id",
                 "config_host","s3_config_endpoint","config_k8s_api_server","config_token",
                 "datahub_project","datahub_gms_token","datahub_gms_url","AIM_ACCESS_TOKEN"]
    config = {}
    for ikey in CONFIG_KEYS:
        config_key = "Config_" + ikey
        if hasattr(opts,config_key) and getattr(opts,config_key).strip():
                config[ikey] = getattr(opts,config_key).strip()
    return config

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

    return_code, stdout, stderr = run_command(command)
    if return_code != 0:
        print("Error Encountered!")
        #print(stdout)
        print(stderr)
    else:
        #print(stdout)
        print(stderr)
    os.chdir(cwd)

def produce_metrics_superMetrics_reports(allparams,work_path,output_path):
    from abacustest import outresult
    import pandas as pd
    from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement

    reports = []
    metric_filename = allparams.get("post_dft",{}).get("metrics",{}).get("save_file","result.json")
    metric_file = os.path.join(work_path,allparams["save_path"],metric_filename)
    if os.path.isfile(metric_file):
        allresults = json.load(open(metric_file))
        outresult.pandas_out(allresults,os.path.join(output_path,"metrics.csv"))
        reports.append(ReportSection(title="metrics",elements=[AutoReportElement(title='metrics', path="metrics.csv",description="")]))
        print("create reportsection")
    
    super_metric_filename = allparams.get("post_dft",{}).get("super_metric",[{}])[0].get("save_file","superMetric.json")
    super_metric_file = os.path.join(work_path,allparams["save_path"],super_metric_filename)
    if os.path.isfile(super_metric_file):
        supermetrics = json.load(open(super_metric_file))
        pddata = pd.DataFrame.from_dict(supermetrics)
        pddata.to_csv(os.path.join(output_path,"superMetrics.csv"),index=False)
        reports.append(ReportSection(title="super_metrics",
                                     elements=[
                                         AutoReportElement(title='super_metrics', path="superMetrics.csv", description="")]))
    return reports