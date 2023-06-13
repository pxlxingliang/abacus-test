from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.typing import InputFilePath, OutputDirectory
from dp.launching.typing import (
    Field
)
from dp.launching.report import Report,AutoReportElement,ReportSection

from . import comm_class,comm_func,comm_class_exampleSource
import json,traceback,os

class SelfDefine(BaseModel):
    IO_input_path:InputFilePath = Field(default = None,
                                        title="Upload setting file",
                                        st_kwargs_type = ["json"], 
                                        description="Please upload the setting file or enter the setting information in latter 'setting' section.",
                                        description_type="markdown")
    setting: String = Field(default = "",
                            description="Please enter the setting information in json format. If you fill in this field, the previous file will be ignored!")
    

class SelfDefineModel(comm_class.TrackingSet,
                      SelfDefine,
                      comm_class_exampleSource.ExampleSourceSet,
                      comm_class.ConfigSet,
                      comm_class.OutputSet,
                      BaseModel):
    ...
    
def SelfDefineModelRunner(opts:SelfDefineModel):
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

    #parse inputs
    if opts.setting.strip() != "":
        try:
            setting = json.loads(opts.setting)
        except:
            traceback.print_exc()
            return 1
    elif opts.IO_input_path != None:
        try:
            setting = json.load(open(opts.IO_input_path.get_path()))
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
        if "AIM_ACCESS_TOKEN" not in allparams["config"]:
            allparams["config"]["AIM_ACCESS_TOKEN"] = tracking_set.get("token")
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
    reports = comm_func.produce_metrics_superMetrics_reports(allparams,work_path,output_path)
    
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
    
    return 0