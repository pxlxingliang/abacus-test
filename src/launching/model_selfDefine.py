from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.typing import InputFilePath, OutputDirectory
from dp.launching.typing import (
    Field
)
from dp.launching.report import Report

from . import comm_class,comm_func
import json,traceback

class SelfDefine(BaseModel):
    IO_input_path:InputFilePath = Field(default = None,
                                        title="Upload setting file",
                                        st_kwargs_type = ["json"], 
                                        description="Please upload the setting file or enter the setting information in latter 'setting' section.",
                                        description_type="markdown")
    IO_output_path: OutputDirectory = Field(default="./output")
    setting: String = Field(default = "",
                            description="Please enter the setting information in json format. If you fill in this field, the previous file will be ignored!")
    

class SelfDefineModel(SelfDefine,comm_class.ConfigSet,BaseModel):
    ...
    
def SelfDefineModelRunner(opts:SelfDefineModel):
    paths = comm_func.create_path(str(opts.IO_output_path))
    output_path = paths["output_path"]
    work_path = paths["work_path"]
    download_path = paths["download_path"]

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
        elif k in ["ABBREVIATION","save_path","run_dft","post_dft","report","dataset_info","upload_datahub","upload_tracking"]:
            allparams[k] = v

    #execut
    comm_func.exec_abacustest(allparams,work_path)
    reports = comm_func.produce_metrics_superMetrics_reports(allparams,work_path,output_path)

    if reports:
        report = Report(title="abacus test report",
                        sections=reports,
                        description="a report of abacustest")
        report.save(output_path)
    
    return 0