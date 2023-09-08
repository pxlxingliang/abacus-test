import json
import os
import sys
import subprocess
import glob
import shutil
import select
import traceback
import copy
from . import comm_echarts
from dp.metadata import MetadataContext
from dp.metadata.utils.storage import TiefblueStorageClient
from dflow import download_artifact, S3Artifact, config, s3_config
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient
from dp.launching.report import Report, AutoReportElement, ReportSection, ChartReportElement

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
                   "aim_access_token", "dflow_labels"]
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
    parser = argparse.ArgumentParser(description="abacustest")
    subparser = parser.add_subparsers(dest="command")
    abatest.AbacusTestArgs(subparser.add_parser("submit"))
    parser = parser.parse_args(["submit", "-p", "param.json"])
    abatest.abacustest(parser)

    os.chdir(cwd)
    stdout, stderr = "", ""
    return stdout, stderr

def sort_lists(lists):
    # lists is a list of list
    # lists = [[x1,x2,...],[y1,y2,...],...]
    # return the sorted lists
    # the first list will be sorted, and the other lists will be sorted according to the first list
    # if the length of lists is not equal, return the original lists
    if len(lists) in [0,1]:
        return lists
    else:
        # check if all elements are list and the length of each list is equal
        if False in [isinstance(i,list) for i in lists]:
            return lists
        if len(set([len(i) for i in lists])) != 1:
            return lists
        
        lists_tmp = copy.deepcopy(lists)
        new_lists = [[] for i in range(len(lists))]
        for ii in sorted(zip(*lists_tmp)):
            for ij,jj in enumerate(list(ii)):
                new_lists[ij].append(jj)
        
        return new_lists

def add_ref(allresults, ref_data):
    # allresults: the metrics data
    # ref_data: the reference data
    # return the metrics data with reference data
    # allresults = {"example1":{"metric1":value1,"metric2":value2,...},"example2":{"metric1":value1,"metric2":value2,...},...}
    # ref_data = {"ref1":{"example1":{"metric1":value1,"metric2":value2,...},"example2":{"metric1":value1,"metric2":value2,...},...},"ref2":{"example1":{"metric1":value1,"metric2":value2,...},"example2":{"metric1":value1,"metric2":value2,...},...},...}
    # 1. the reference metrics wiil be named as "<metrics>_ref_<name of reference>"
    # 2. if example/metrics in ref_data is not in allresults, the ref will not be added
    # 3. if example in ref_data is a father path of example in allresults, the ref will be added

    # get the example name
    example_name = list(allresults.keys())
    
    # get the metrics name
    metric_name = []
    for ivalue in allresults.values():
        metric_name += list(ivalue.keys())
    metric_name = list(set(metric_name))
    
    # get the reference name
    ref_name = list(ref_data.keys())
    # get the metric_name in ref_data
    ref_metric_name = []
    for iref in ref_name:
        for imetric in ref_data[iref].values():
            ref_metric_name += list(imetric.keys())
    #only metrics in metric_name will be added
    #print(ref_name,ref_metric_name,metric_name)
    ref_metric_name = list(set(ref_metric_name) & set(metric_name))
    
    # if ref_metric_name is empty, return
    if len(ref_metric_name) == 0:
        return allresults,[]
    
    new_results = {}
    for ikey,ivalue in allresults.items():
        new_results[ikey] = ivalue.copy()
        
        for iref in ref_name:
            # check if ikey in ref_data[iref]
            if ikey in ref_data[iref]:
                for imetric in ref_metric_name:
                    new_results[ikey][f"{imetric}_ref_{iref}"] = ref_data[iref][ikey].get(imetric,None)
            else:
                # check if path in iref is a father path of ikey
                for iref_example in ref_data[iref]:
                    if ikey.startswith(iref_example):
                        for imetric in ref_metric_name:
                            new_results[ikey][f"{imetric}_ref_{iref}"] = ref_data[iref][iref_example].get(imetric,None)
                        break
    return new_results,[f"{imetric}_ref_{iref}" for imetric in ref_metric_name for iref in ref_name]   

def plot_delta_Y(y_list, legend_list, example_name, imetric):
    # y_list is a list of list
    # y_list = [[y1,y2,...],[y1,y2,...],...]
    # return the delta y_list and percentage delta y_list
    # the first list is the value of this job, and the other lists are the values of different reference
    # length of y_list should be at least 2, and equal to length of legend_list
    chart_elements = []
    
    def flat_list(list1):
        new_list = []
        for i in list1:
            if isinstance(i,list):
                new_list += flat_list(i)
            else:
                new_list.append(i)
        return new_list
    
    # need to check if y_list[0] is a list of list
    #print(imetric)
    if True in [isinstance(i,list) for i in y_list[0]]:
        # if y_list[0] is a list of list, we need to plot the delta Y
        #print(y_list)
        new_y_list = copy.deepcopy(y_list)
        # need to transfer each list to one dimension list
        for iy in range(len(new_y_list)):
            for iiy in range(len(new_y_list[iy])):
                if isinstance(new_y_list[iy][iiy],list):
                    new_y_list[iy][iiy] = flat_list(new_y_list[iy][iiy])
        
        #print("new_y_list",new_y_list)
        #print(len(new_y_list),len(new_y_list[0]))
        delta_y = []
        for iy in range(1,len(new_y_list)): # loop for each reference
            delta_y.append([])
            for jy in range(len(new_y_list[0])): # loop for each example
                if not isinstance(new_y_list[0][jy],list) or not isinstance(new_y_list[iy][jy],list) or len(new_y_list[0][jy]) != len(new_y_list[iy][jy]):
                    #print(type(new_y_list[0][jy]),type(new_y_list[iy][jy]))
                    delta_y[-1].append(None)
                    continue
                try:
                    ivalue = [new_y_list[0][jy][i] - new_y_list[iy][jy][i] for i in range(len(new_y_list[0][jy]))]
                    delta_y[-1].append(ivalue)
                except:
                    delta_y[-1].append(None)
        # now we need to calculate some statistics
        # 1. the max absolute value
        # 2. the norm of delta_y
        
        max_abs = []
        norm = []
        legend_abs = []
        legend_norm = []
        for iy in range(len(delta_y)):# loop for each reference
            iabs = []
            inorm = []
            for jy in range(len(delta_y[iy])):# loop for each example
                if delta_y[iy][jy] == None:
                    iabs.append(None)
                    inorm.append(None)
                    continue
                try:
                    ivalue = max([abs(i) for i in delta_y[iy][jy]])
                except:
                    ivalue = None
                iabs.append(ivalue)
                try:
                    ivalue = sum([i**2 for i in delta_y[iy][jy]])**0.5
                except:
                    ivalue = None
                inorm.append(ivalue)
            if set(iabs) != {None}:
                max_abs.append(iabs)
                legend_abs.append(legend_list[iy+1])
            if set(inorm) != {None}:
                norm.append(inorm)
                legend_norm.append(legend_list[iy+1])
        
        # plot max_abs and norm
        
        if len(max_abs) > 0:
            tmp_ = sort_lists([example_name]+max_abs)
            options = comm_echarts.produce_multiple_y(f"{imetric}(max(abs(this job - reference)))", tmp_[0], tmp_[1:], legend_abs, x_type="category", y_type="value")
            options["xAxis"][0]["axisLabel"] = {
                "rotate": 15,
                "interval": int(len(example_name)/15)
                } 
            chart_elements.append(ChartReportElement(
                    options=options, title=f"{imetric}(max(abs(this job - reference)))"))
        if len(norm) > 0:
            tmp_ = sort_lists([example_name]+norm)
            options = comm_echarts.produce_multiple_y(f"{imetric}(Norm(this job - reference))", tmp_[0], tmp_[1:], legend_norm, x_type="category", y_type="value")
            options["xAxis"][0]["axisLabel"] = {
                "rotate": 15,
                "interval": int(len(example_name)/15)
                } 
            chart_elements.append(ChartReportElement(
                    options=options, title=f"{imetric}(Norm(this job - reference))"))      
    else:
        delta_y_list = []
        percentage_delta_y_list = []
        delta_legend_list1 = []
        delta_legend_list2 = []

        for iy in range(1,len(y_list)):
            # need to check if the two value can do minus
            values1 = []
            values2 = []
            for i in range(len(y_list[iy])):
                if isinstance(y_list[0][i],(int,float)) and isinstance(y_list[iy][i],(int,float)):
                    ivalue = y_list[0][i] - y_list[iy][i]                     
                else:
                    ivalue = None
                values1.append(ivalue)

                if isinstance(y_list[0][i],(int,float)) and isinstance(y_list[iy][i],(int,float)) and y_list[iy][i] != 0:
                    ivalue = (y_list[0][i] - y_list[iy][i])/y_list[iy][i]
                else:
                    ivalue = None
                values2.append(ivalue)
            if set(values1) != {None}:
                delta_y_list.append(values1)
                delta_legend_list1.append(legend_list[iy])
            if set(values2) != {None}:
                percentage_delta_y_list.append(values2)
                delta_legend_list2.append(legend_list[iy])
        if len(delta_y_list) > 0:  
            tmp_ = sort_lists([example_name]+delta_y_list)
            options = comm_echarts.produce_multiple_y(f"{imetric}(Delta = this job - reference)", tmp_[0], tmp_[1:], delta_legend_list1, x_type="category", y_type="value")
            options["xAxis"][0]["axisLabel"] = {
                "rotate": 15,
                "interval": int(len(example_name)/15)
                } 
            chart_elements.append(ChartReportElement(
                    options=options, title=f"{imetric}(Delta = this job - reference)"))
        if len(percentage_delta_y_list) > 0:
            tmp_ = sort_lists([example_name]+percentage_delta_y_list)
            options = comm_echarts.produce_multiple_y(f"{imetric}(Delta/Reference)", tmp_[0], tmp_[1:], delta_legend_list2, x_type="category", y_type="value")
            options["xAxis"][0]["axisLabel"] = {
                "rotate": 15,
                "interval": int(len(example_name)/15)
                } 
            chart_elements.append(ChartReportElement(
                    options=options, title=f"{imetric}(Delta/Reference)"))
    return chart_elements

def plot_two_metrics(pddata,
                     x_name, y_name,
                     example_name,
                     ref_metric_name=[],
                     x_type="category", y_type="value", shift_type=0):
    '''
    pddata: the pandas data
    ref_metric_name: the reference metric name
    x_name: the name of x
    y_name: the name of y
    example_name: a list of example name
    x_type: the type of x in echart  # now only support category
    y_type: the type of y in echart
    shift_type: the type to shift the y value
        -1: has reference, will minus the reference value
        0: do not shift
        1: shift to make the min value is 0
        2: shift to make the last value is 0
    if has reference, all value will minus the reference value
    if has multiple reference, will minus each reference, and set shift_type to -1
    '''
    x_type = "category"

    if x_name not in pddata.index.to_list() or y_name not in pddata.index.to_list():
        return []

    # check if has reference, only chech for y_name
    ref_names = []
    ref_values = []
    for imetric in ref_metric_name:
        if imetric.startswith(f"{y_name}_ref_") and imetric in pddata.index.to_list():
            ref_names.append(imetric.split("_ref_")[-1])
            ref_values.append(pddata.loc[imetric, :].to_list())
            shift_type = -1
    
    #print("x_name,y_name:", x_name, y_name)
    all_x = pddata.loc[x_name, :].to_list()
    all_y = pddata.loc[y_name, :].to_list()
    chart_elements = []
    if len(set(all_x)) <= 1:
        return []

    # the example name may be a/00000, a/00001, ..., b/00000, b/00001
    # we need plot each chart for a, b, ...
    # so we need to split the example name to get the prefix
    # and then plot the chart for each prefix
    # we need to plot the chart for all example that do not match the prefix/00000 format
    # we will seperate the example to several parts: prefix/00000 and others
    all_xy = {"": [[], [], f"{x_name} VS {y_name}"]}
    for i in range(len(ref_names)): 
        all_xy[""].append([])  # add a new list for each reference

    # remove / in the end of the example_name
    example_name = [i.rstrip("/") for i in example_name]
    for i, ie in enumerate(example_name):
        prefix = os.path.dirname(ie)
        basename = os.path.basename(ie)
        if basename.isdigit():
            if prefix not in all_xy:
                all_xy[prefix] = [[], [], f"{x_name} vs {y_name} ({prefix})"]
                for iref in range(len(ref_names)):
                    all_xy[prefix].append([])
            all_xy[prefix][0].append(all_x[i])
            all_xy[prefix][1].append(all_y[i])
            for iref in range(len(ref_names)):
                all_xy[prefix][iref+3].append(ref_values[iref][i])
        else:
            all_xy[""][0].append(f"{all_x[i]}({prefix})")
            all_xy[""][1].append(all_y[i])
            for iref in range(len(ref_names)):
                all_xy[""][iref+3].append(ref_values[iref][i])

    # if one prefix only has one example, add the example to ""
    for iprefix, ivalue in all_xy.items():
        if len(ivalue[0]) == 1:
            all_xy[""][0].append(f"{ivalue[0][0]}({iprefix})")
            all_xy[""][1].append(ivalue[1][0])
            for iref in range(len(ref_names)):
                all_xy[""][iref+3].append(ivalue[iref+3][0])
            del all_xy[iprefix]

    # plot the chart for each prefix
    for iprefix, ivalue in all_xy.items():
        # if ivalues is empty or only one value, continue
        if len(ivalue[0]) <= 1:
            continue

        x = ivalue[0]
        y = [ivalue[1]]
        for iref in range(len(ref_names)):
            y.append(ivalue[iref+3])
        title_full = title = ivalue[2]

        # need shift the y
        y_real = [i for i in y[0] if i != None]  # get the none-none value of y[0]
        if len(y_real) == 0:
            continue
        legend = [y_name]
        if shift_type == 1:
            title_full = title_full + " (minus the minimum)"
            min_e = min(y_real)
            y = [[i if i == None else i-min_e for i in y[0]]]
        elif shift_type == 2:
            title_full = title_full + " (minus the last value)"
            min_e = y_real[-1]
            y = [[i if i == None else i-min_e for i in y[0]]]
        elif shift_type == -1:
            title_full = title_full + f" (minus reference)"
            title += f" (minus reference)"
            new_y = [[] for i in range(len(y)-1)]
            for iy in range(len(y[0])):
                for iy2 in range(len(y)-1):
                    if y[0][iy] != None and y[iy2+1][iy] != None:
                        new_y[iy2].append(y[0][iy]-y[iy2+1][iy])
                    else:
                        new_y[iy2].append(None) 
            y = new_y
            legend = [f"ref_{i}" for i in ref_names]
        if y_type == "log":
            title_full = title_full + " (abs(delta_y))"
        
        tmp_ = sort_lists([x]+y)     
        options =comm_echarts.produce_multiple_y(
            title, tmp_[0], tmp_[1:], legend, x_type=x_type, y_type=y_type)
        
        chart_elements.append(ChartReportElement(options=options, title=title_full))
    return chart_elements


def plot_drho(drho_input, example_input):
    # drho_input: a list of drho of all examples
    # example_input: a list of example name
    # return a list of ChartReportElement

    chart_report = []
    # only plot drho if drho is not None
    drho_list = []
    example_name = []
    for i in range(len(drho_input)):
        if isinstance(drho_input[i], list):
            drho_list.append(drho_input[i])
            example_name.append(example_input[i])
    if len(drho_list) == 0:
        return []

    # sort the example_name and drho_list
    example_name, drho_list = zip(*sorted(zip(example_name, drho_list)))
    max_example_in_one_chart = 10
    all_data = []
    for i in range(0, len(example_name), max_example_in_one_chart):
        end = i+max_example_in_one_chart if i + \
            max_example_in_one_chart < len(example_name) else len(example_name)
        all_data.append([example_name[i:end], drho_list[i:end]])
    idrho = 0
    for legend, drho in all_data:
        x = [i+1 for i in range(max([len(j) for j in drho]))]
        # print(x,drho,legend)
        options = comm_echarts.produce_multiple_y(
            f"drho{idrho}", x, drho, legend, x_type="category", y_type="log")
        chart_report.append(ChartReportElement(
            options=options, title=f"drho{idrho}"))
        idrho += 1
    return chart_report

def fill_none(allresults):
    # if some metrics is not in some example, we need to add None to the example
    allkeys = []
    for ivalue in allresults.values():
        allkeys += list(ivalue.keys())
    allkeys = list(set(allkeys))
    for ikey in allkeys:
        for iexample in allresults:
            if ikey not in allresults[iexample]:
                allresults[iexample][ikey] = None

def produce_metrics(metric_file, output_path, ref_data={}, report_titile="metrics"):
    from abacustest import outresult
    import pandas as pd

    metric_filename = os.path.split(metric_file)[-1]
    allresults = json.load(open(metric_file))
    allresults,ref_metric_name = add_ref(allresults, ref_data)
    csv_filename = os.path.splitext(metric_filename)[0] + ".csv"
    _, _, savefile_names = outresult.pandas_out(
        allresults, os.path.join(output_path, csv_filename),print_result=False)
    report_elements = []
    for ifilename in savefile_names:
        ifilename = os.path.split(ifilename)[-1]
        report_elements.append(AutoReportElement(
            title=os.path.splitext(ifilename)[0], path=ifilename, description=""))

    chart_elements = []
    # produce the Echarts option for each metric
    fill_none(allresults)
    pddata = pd.DataFrame.from_dict(allresults)
    example_name = pddata.columns.to_list()  # get the column name
    metric_name = pddata.index.to_list()  # get the row/index name
    type_set = (int, float, type(None), bool,list)
    print("ref_metric_name",ref_metric_name)
    for imetric in metric_name:
        if imetric in ref_metric_name:
            continue
        ivalue = pddata.loc[imetric, :].to_list()
        if False not in [isinstance(i, type_set) for i in ivalue]:
            #check if imetric has refence
            ref_type = []
            for iref_metric in ref_metric_name:
                if iref_metric.startswith(f"{imetric}_ref_"):
                    ref_type.append(iref_metric.split("_ref_")[-1])
                    
            # if imetric has refence, we need to plot the imetric and its reference together
            # we need to split the imetric to imetric and its reference
            # comm_echarts.produce_multiple_y(produce_multiple_y(title,x,y_list,legend_list,x_type="category",y_type="value")
            y_list = [ivalue]
            legend_list = ["This_Job"]
            for iref in ref_type:
                if f"{imetric}_ref_{iref}" in metric_name:
                    y_list.append(pddata.loc[f"{imetric}_ref_{iref}", :].to_list())
                    legend_list.append(f"ref_{iref}")
                    
            if True not in [isinstance(i,list) for i in ivalue]:
                # do not plot a list
                #print(y_list,legend_list)
                tmp_ = sort_lists([example_name]+y_list)
                options = comm_echarts.produce_multiple_y(imetric, tmp_[0], tmp_[1:], legend_list, x_type="category", y_type="value")
                options["xAxis"][0]["axisLabel"] = {
                        "rotate": 15,
                        "interval": int(len(example_name)/15)
                    } 
                chart_elements.append(ChartReportElement(
                        options=options, title=imetric))
            
            
            if len(ref_type) > 0:
                # plot the delta Y and percentage delta Y
                try:
                    delta_y_elements = plot_delta_Y(y_list, legend_list, example_name, imetric)
                    chart_elements += delta_y_elements
                except:
                    traceback.print_exc()
                    print("Error: plot delta Y failed!")
                    

    # produce some special case chart
    # 1. if metrics has ecutwfc/kspacing and energy_per_atom, produce the ecutwfc vs energy_per_atom chart
    print(f"plot {metric_file} metric name:", metric_name)
    for x_name, y_name, y_type, shift_type in [
        ["INPUT/ecutwfc", "energy_per_atom","value",1],
        ["INPUT/kspacing", "energy_per_atom","value",1],
        ["INPUT/lcao_ecut", "energy_per_atom","log",2],
        ["INPUT/lcao_ecut", "band_gap","value",0]]:
        '''
        shift_type: the type to shift the y value
        0: do not shift
        1: shift to make the min value is 0
        2: shift to make the last value is 0
        '''
        try:
            chart_elements += plot_two_metrics(pddata,
                                               x_name, y_name,
                                               example_name,
                                               ref_metric_name,
                                               x_type="category", y_type=y_type, shift_type=shift_type)
        except:
            traceback.print_exc()
            print("Error: plot two metrics failed!", x_name, y_name)

    # if drho in metric_name, plot the drho
    # drho should be a list
    if "drho" in metric_name:
        try:
            drho_list = pddata.loc["drho", :].to_list()
            chart_elements += plot_drho(drho_list, example_name)
        except:
            traceback.print_exc()
            print("Error: plot drho failed!")

    return report_elements, chart_elements


def produce_supermetrics(supermetric_file, output_path, work_path, save_path, report_titile="supermetrics"):
    import pandas as pd

    supermetrics_filename = os.path.split(supermetric_file)[-1]
    super_metrics = json.load(open(supermetric_file))
    normal_metrics = {}
    special_metrics = {}
    for ikey, ivalue in super_metrics.items():
        if isinstance(ivalue, dict):
            special_metrics[ikey] = ivalue
        else:
            if isinstance(ivalue,list) and len(ivalue) == 1:
                normal_metrics[ikey] = ivalue[0]
            else:
                normal_metrics[ikey] = ivalue

    report_element = None
    if normal_metrics:
        pddata = pd.DataFrame(normal_metrics, index=[0],)
        t_pddata = pddata.transpose()
        print(t_pddata)
        csv_filename = os.path.splitext(supermetrics_filename)[0] + ".csv"
        t_pddata.to_csv(os.path.join(output_path, csv_filename), index=True, header=True)
        report_element = AutoReportElement(
            title=report_titile, path=csv_filename, description="")

    special_section = None
    special_elements = []
    if special_metrics:
        for ikey, ivalue in special_metrics.items():
            file_name = ivalue.get("file", None)
            if file_name:
                ifile_name = os.path.join(work_path, save_path, file_name)
                if os.path.isfile(ifile_name):
                    special_elements.append(AutoReportElement(
                        title=ikey, path=os.path.join(save_path, file_name), description=""))

    if special_elements:
        special_section = ReportSection(
            title=report_titile, elements=special_elements)

    return report_element, special_section, normal_metrics

def gen_sm_tag(allparams):
    # generate the tag for supermetrics to save on launching
    # the tag include below information:
    # 1. lbg_username from allparams["config"]
    # 2. schedule name from allparams["config"]["dflow_labels"]["launching-schedule"]
    # 3. job name from allparams["config"]["dflow_labels"]["launching-job"]
    return [
        f"user_{allparams.get('config',{}).get('bohrium_username','')}",
        f"schedule_{allparams.get('config',{}).get('dflow_labels',{}).get('launching-schedule','')}",
        f"job_{allparams.get('config',{}).get('dflow_labels',{}).get('launching-job','')}"    
    ]

def judge_sm(x,criteria):
    if criteria == None:
        return None
    if isinstance(criteria,str):
        try:
            print("SuperMetric value:",x, "\ncriteria:",criteria)
            sm_pass = eval(criteria)
            return bool(sm_pass)
        except:
            print("Error: the supermetrics criteria string is not correct!")
            return False
    else:
        print("Error: the supermetrics  criteria string is not correct!")
        return False
            

def produce_metrics_superMetrics_reports(allparams, work_path, output_path):

    save_path = allparams.get("save_path", "result")
    reports = []
    allmetrics_files = []
    allsupermetrics_files = []

    # produce the metrics reports 
    metrics_block = allparams.get("post_dft", {}).get("metrics", [])
    if not isinstance(metrics_block, list):
        metrics_block = [metrics_block]
    for imetrics_block in metrics_block:
        if isinstance(imetrics_block,dict):
            # 1. metrics from the custome defined metrics file
            if imetrics_block.get("value_from_file",None):
                allmetrics_files.append(os.path.join(
                    work_path, save_path, str(imetrics_block.get("value_from_file"))))
                
            # 2. metrics from the metrics file
            if imetrics_block.get("save_file",None):
                allmetrics_files.append(os.path.join(
                    work_path, save_path, str(imetrics_block.get("save_file"))))

    # 3. metrics from the undefined metrics file
    allmetrics_files += glob.glob(os.path.join(work_path,
                                  save_path, "metric*.json"))
    allmetrics_files = list(set(allmetrics_files))
    # 4. support the compare with reference, the reference file name shuold be "metrics_ref.json"
    # we need to find the reference file and get the reference metrics.
    # the format of the reference file is the same as the metrics file
    # or the key is the name of reference, and the value is the metrics
    # the reference metrics wiil be named as "<metrics>_ref_<name of reference>"
    # we first read the reference file and get the reference metrics, and save to a dict
    # if the ref file is metrics format, then the refence name will be ""
    ref_data = {}
    ref_file = os.path.join(work_path, save_path, "metrics_ref.json")
    if os.path.isfile(ref_file):
        ref_metrics = json.load(open(ref_file))
        # need to check if the reference metrics is a dict of dict of dict
        format_ok = True
        Three_layer = True
        if isinstance(ref_metrics, dict):
            for ikey, ivalue in ref_metrics.items():
                if not isinstance(ivalue, dict):
                    format_ok = False
                    break
                if Three_layer:
                    for jkey, jvalue in ivalue.items():
                        if not isinstance(jvalue, dict):
                            Three_layer = False
                            break
        else:
            format_ok = False
        if format_ok:
            if Three_layer:
                ref_data = ref_metrics
            else:
                ref_data = {"": ref_metrics}
                
        if ref_file in allmetrics_files:
            allmetrics_files.remove(ref_file)
        
    metrics_report = []
    metrics_chart_elements = []
    for metric_file in list(set(allmetrics_files)):
        print("PRODUCE REPORT for metric file:", metric_file)
        if not os.path.isfile(metric_file):
            print(f"\tError: metric file \"{metric_file}\" does not exist! Skip it!")
            continue
        try:
            metric_filename = os.path.split(metric_file)[-1]
            tmp_report_elements, tmp_chart_elements = produce_metrics(
                metric_file, output_path, ref_data=ref_data,report_titile=metric_filename)
            if tmp_report_elements:
                metrics_report += tmp_report_elements
            if tmp_chart_elements:
                metrics_chart_elements += tmp_chart_elements
        except:
            traceback.print_exc()
            print(f"Error: report metrics from file \"{metric_file}\" failed!")

    # produce the supermetrics reports
    # 1. supermetrics from the custome defined supermetrics file
    super_metric_filename = allparams.get("post_dft", {}).get(
        "super_metrics", [{}])[0].get("save_file")
    if super_metric_filename:
        allsupermetrics_files.append(os.path.join(
            work_path, save_path, super_metric_filename))

    # 2. supermetrics from the supermetrics file
    customized_super_metrics_filename = allparams.get("post_dft", {}).get(
        "super_metrics", [{}])[0].get("value_from_file", None)
    if customized_super_metrics_filename:
        allsupermetrics_files.append(os.path.join(
            work_path, save_path, customized_super_metrics_filename))

    # 3. supermetrics from the undefined supermetrics file
    allsupermetrics_files += glob.glob(os.path.join(work_path, save_path, "superMetric*.json")) + \
        glob.glob(os.path.join(work_path, save_path, "super_metric*.json")) + \
        glob.glob(os.path.join(work_path, save_path, "supermetric*.json"))
    
    # 4. read the criteria.json, which defined criteria for supermetrics
    if os.path.isfile(os.path.join(work_path, save_path, "criteria.json")):
        criterias = json.load(open(os.path.join(work_path, save_path, "criteria.json")))
    else:
        criterias = {}
    print("criteria:",criterias)
    
    supermetrics_report = []
    supermetrics_special_section = []
    supermetrics_save = {}
    # the summary of supermetrics, this table will be tranlated to a html table
    sm_summary = [["Super Metric","Value","Criteria"]]   
    for super_metric_file in list(set(allsupermetrics_files)):
        print("PRODUCE REPORT for supermetric file:", super_metric_file)
        if not os.path.isfile(super_metric_file):
            print(
                f"\tError: supermetric file \"{super_metric_file}\" does not exist! Skip it!")
            continue
        try:
            super_metric_filename = os.path.split(super_metric_file)[-1]
            tmp_report_element, tmp_special_section, tmp_smetrics = produce_supermetrics(
                super_metric_file, output_path, work_path, save_path, report_titile=super_metric_filename)
            if tmp_report_element:
                supermetrics_report.append(tmp_report_element)
            if tmp_special_section:
                supermetrics_special_section.append(tmp_special_section)
            #print("tmp_smetrics:",tmp_smetrics)
            if tmp_smetrics:
                for ik,iv in tmp_smetrics.items():
                    sm_summary.append([ik,iv,criterias.get(ik,None)])
                    cri_tmp = criterias.get(ik,None)
                    sm_pass = judge_sm(iv,cri_tmp)
                    if sm_pass == True:
                        # if True, set the iv to green
                        sm_summary[-1][1] = "<font color=\"green\">"+str(iv)+"</font>"
                    elif sm_pass == False:
                        # if False, set the iv to red
                        sm_summary[-1][1] = "<font color=\"red\">"+str(iv)+"</font>"
                            
                    if not isinstance(iv,(int,float,bool)):
                        continue
                    sname = ik
                    if ik in supermetrics_save:
                        print(f"\tWarning: supermetric \"{ik}\" is already in the supermetrics_save, will be replaced by the new one!")
                        sname = ik + "_" + os.path.splitext(os.path.basename(super_metric_file))[0]
                    if isinstance(iv,bool):
                        supermetrics_save[sname+":int"] = 1 if iv else 0
                    else:
                        supermetrics_save[sname+":float"] = iv
        except:
            traceback.print_exc()
            print(
                f"Error: report supermetrics from file \"{super_metric_file}\" failed!")

    # produce the summary of supermetrics
    if len(sm_summary) > 1:
        sm_summary = produce_html_table(sm_summary)
        criterias_comment = criterias.get("__COMMENT","")
        # write the summary to a html file
        with open(os.path.join(output_path,"supermetrics_summary.html"),"w") as f:
            # firstly add the head of html
            # set the table width to 100%, and set content to center
            # add a title to the table, "supermetrics summary"
            # add the criterias_comment to the end of the table
            f.write("<html><head><meta charset=\"utf-8\"><style>table {width:100%;text-align:center;}</style></head><body>")
            f.write("<table border=\"2px\"><thead><tr><td>Super Metrics Summary</td></tr></thead>")
            f.write(sm_summary)
            if criterias_comment:
                f.write("<tr><td>Note: "+criterias_comment+"</td></tr>")
            f.write("</body></html>")
        # add the summary to the report at the first position
        supermetrics_report = [] # only keep the summary
        supermetrics_report.insert(0,AutoReportElement(title="supermetrics summary",path="supermetrics_summary.html",description=""))
        

    # 2. supermetrics report
    #print("supermetrics_save:",supermetrics_save)
    if supermetrics_report:
        reports.append(ReportSection(title="supermetrics",
                       elements=supermetrics_report, ncols=1,
                       metrics=supermetrics_save,
                       tags=gen_sm_tag(allparams)))

    # 4. supermetrics special section
    if supermetrics_special_section:
        for i in supermetrics_special_section:
            reports.append(i)
            
    # produce the report
    # 1. metrics report
    if metrics_report:
        reports.append(ReportSection(title="metrics",
                       elements=metrics_report, ncols=1))
        
    # 3. metrics chart
    if metrics_chart_elements:
        reports.append(ReportSection(title="metrics chart",
                       elements=metrics_chart_elements, ncols=2))

    return reports


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
    for f in allfiles:
        if os.path.isdir(f):
            alldirs.append(f)
    packed_file_name = "results.zip"
    pack(alldirs,packed_file_name,"zip")
    os.chdir(cwd)

def gen_dir(dir1):
    dir1 = dir1.rstrip("/")
    n = 0
    while os.path.isdir(dir1):
        dir1 = dir1 + "_" + str(n)
        n += 1
    return dir1
    