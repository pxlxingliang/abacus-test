import json,os,sys,subprocess
import copy 
from . import comm_echarts

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

def convert_metrics(metrics_list):
    #the metrics in launching may have the format of KEY1:KEY2
    #this should be transfer to KEY1:{[KEY2]} in abacustest metrics
    #and transfer to KEY1/KEY2 in abacustest supermetrics
    new_metrics = []
    dict_tmp = {}
    for imetric in metrics_list:
        if ":" in imetric:
            allkeys = imetric.split(":")
            if allkeys[0] not in dict_tmp:
                dict_tmp[allkeys[0]] = []
            dict_tmp[allkeys[0]].append(allkeys[1])
        else:
            new_metrics.append(imetric)
    if dict_tmp:
        new_metrics.append(copy.deepcopy(dict_tmp))
    return new_metrics

def convert_supermetrics_metrics_name(imetric):
    if ":" in imetric:
        allkeys = imetric.split(":")
        return "%s/%s" % (allkeys[0].strip(),allkeys[1].strip())
    else:
        return imetric.strip()


def produce_metrics_superMetrics_reports(allparams,work_path,output_path):
    from abacustest import outresult
    import pandas as pd
    from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement

    reports = []
    chart_section = []
    metric_filename = allparams.get("post_dft",{}).get("metrics",{}).get("save_file","result.json")
    metric_file = os.path.join(work_path,allparams["save_path"],metric_filename)
    print(metric_file)
    if os.path.isfile(metric_file):
        allresults = json.load(open(metric_file))
        outresult.pandas_out(allresults,os.path.join(output_path,"metrics.csv"))
        reports.append(ReportSection(title="metrics",elements=[AutoReportElement(title='metrics', path="metrics.csv",description="")]))

        #produce the Echarts option for each metric
        pddata = pd.DataFrame.from_dict(allresults)
        example_name = pddata.columns.to_list() #get the column name
        metric_name = pddata.index.to_list()  #get the row/index name
        type_set = (int,float,type(None),bool)
        for imetric in metric_name:
            ivalue = pddata.loc[imetric,:].to_list()
            if False not in [isinstance(i,type_set) for i in ivalue]:
                print(imetric,example_name,ivalue)
                options = comm_echarts.get_bar_option(imetric,example_name,ivalue)
                options["xAxis"][0]["axisLabel"] = {
                    "rotate": 15,
                    "interval": int(len(example_name)/15)
                }
                chart_section.append(ChartReportElement(options=options,title=imetric))

    super_metric_filename = allparams.get("post_dft",{}).get("super_metrics",[{}])[0].get("save_file","superMetrics.json")
    super_metric_file = os.path.join(work_path,allparams["save_path"],super_metric_filename)
    if os.path.isfile(super_metric_file):
        supermetrics = json.load(open(super_metric_file))
        pddata = pd.DataFrame.from_dict(supermetrics)
        print(pddata)
        pddata.to_csv(os.path.join(output_path,"superMetrics.csv"),index=False)
        reports.append(ReportSection(title="super_metrics",
                                     elements=[
                                         AutoReportElement(title='super_metrics', path="superMetrics.csv", description="")]))
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
