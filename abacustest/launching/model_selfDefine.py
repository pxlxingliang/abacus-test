from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.typing import InputFilePath, OutputDirectory
from dp.launching.typing import (
    Field,
    DataSet,
)
from dp.launching.report import Report,AutoReportElement,ReportSection

from . import comm_class,comm_func,comm_class_exampleSource,comm_pmetrics
import json,traceback,os

import dp.launching.typing.addon.ui as ui
from dp.launching.typing.addon.sysmbol import Equal,NotEqual,Exists,NotExists

class SelfDefine(BaseModel):
    setting_upload: InputFilePath = Field(default=None,
                                          title="Upload setting file",
                                          st_kwargs_type=["json"],
                                          description="Use setting file by uploading. (priority: setting string > setting upload > setting dataset)",
                                          description_type="markdown")

    setting_dataset: DataSet = Field(title="setting dataset",
                                     default=None,
                                     description="Use setting file from dataset. (priority: setting string > setting upload > setting dataset)",)

    setting_string: String = Field(
        title="setting string",
        default="",
        description="Use setting by enter string. (priority: setting string > setting upload > setting dataset)",)
  

class SelfDefineModel(comm_class.TrackingSet,
                      SelfDefine,
                      comm_class_exampleSource.ExampleSet,
                      comm_class_exampleSource.ExampleSourceSet,
                      comm_class.ConfigSet,
                      comm_class.OutputSet,
                      BaseModel):
    ...

class SelfDefineDatasetsModel(comm_class.TrackingSet,
                      SelfDefine,
                      comm_class_exampleSource.ExampleSet,
                      comm_class_exampleSource.DatasetSet,
                      comm_class.ConfigSet,
                      comm_class.OutputSet,
                      BaseModel):
    ...
    
def SelfDefineModelRunner(opts):
    logs = comm_class.myLog()  
    paths = comm_func.create_path(str(opts.IO_output_path))
    output_path = paths["output_path"]
    work_path = paths["work_path"]
    download_path = paths["download_path"]

    logs.iprint("read source setting ...")
    datas = comm_class_exampleSource.read_source(opts,work_path,download_path,logs.iprint)
    if datas == None or not datas.get("all_files"):
        logs.iprint("Error: download examples or rundft_extrafiles or postdft_extrafiles failed!")
        return 1

    # parse inputs
    # setting inputs is prefered
    # and then the uploaded file
    # and then the file in dataset or examples
    if opts.setting_string.strip() != "":
        try:
            setting = json.loads(opts.setting_string)
        except:
            traceback.print_exc()
            return 1
    elif opts.setting_upload != None:
        try:
            setting = json.load(open(opts.setting_upload.get_path()))
        except:
            traceback.print_exc()
            return 1
    elif opts.setting_dataset:
        try:
            setting = json.load(open(opts.setting_dataset.get_full_path()))
        except:
            traceback.print_exc()
            return 1
    else:
        print("Please supply the setting information")
        return 1

    #read setting
    allparams = {"config": comm_func.read_config(opts)}
    for k,v in setting.items():
        if k == "config":
            for ik,iv in v.items():
                allparams["config"][ik] = iv
        #elif k in ["ABBREVIATION","save_path","run_dft","post_dft","report","dataset_info","upload_datahub","upload_tracking"]:
        else:
            allparams[k] = v
            
    tracking_set = comm_class.TrackingSet.parse_obj(opts)
    if tracking_set:
        if "aim_access_token" not in allparams["config"]:
            allparams["config"]["aim_access_token"] = str(tracking_set.get("token"))
        if "post_dft" not in allparams:
            allparams["post_dft"] = {}
        if "upload_tracking" not in allparams["post_dft"]:
            allparams["post_dft"]["upload_tracking"] = {}
        if "name" not in allparams["post_dft"]["upload_tracking"]:
            allparams["post_dft"]["upload_tracking"]["name"] = tracking_set.get("name")
        if "experiment" not in allparams["post_dft"]["upload_tracking"]:
            allparams["post_dft"]["upload_tracking"]["experiment"] = tracking_set.get("experiment")
        if "tags" not in allparams["post_dft"]["upload_tracking"]:
            allparams["post_dft"]["upload_tracking"]["tags"] = tracking_set.get("tags")
            
    #execut
    stdout,stderr = comm_func.exec_abacustest(allparams,work_path)
    logs.iprint(f"{stdout}\n{stderr}\nrun abacustest over!\n")
    reports = comm_pmetrics.produce_metrics_superMetrics_reports(allparams,work_path,output_path)
    
    logfname = "output.log"
    logs.write(os.path.join(str(opts.IO_output_path),logfname))
    log_section = ReportSection(title="",
                              elements=[AutoReportElement(title='', path=logfname, description="")])
    reports.append(log_section)

    if reports:
        report = Report(title="abacus test report",
                        sections=reports,
                        description="a report of abacustest")
        report.save(output_path)
        
    #move results to output_path
    comm_func.move_results_to_output(work_path,output_path,allparams.get("save_path","results"))
    comm_func.pack_results(output_path,allparams.get("save_path","results"))
    
    return 0
