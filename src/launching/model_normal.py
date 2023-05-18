from enum import Enum
from typing import Literal

from sqlalchemy import desc

from dp.launching.cli import to_runner,SubParser,run_sp_and_exit
from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.cli import to_runner, default_minimal_exception_handler
from dp.launching.typing import InputFilePath, OutputDirectory
from dp.launching.typing import (
    BaseModel,
    Set,
    Boolean,
    Field,
    DflowAccessToken,
    DflowArgoAPIServer,
    DflowK8sAPIServer,
    DflowStorageEndpoint,
    DflowStorageRepository,
    BohriumMachineType,
    BohriumImage,
    BohriumPlatform,
    BohriumJobType,
    BohriumUsername,
    BohriumPassword,
    BohriumProjectId
)
import os,zipfile,shutil,json

from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement

from . import (comm_class,
               comm_func,
               comm_class_exampleSource,
               comm_class_rundft,
               comm_class_postdft,
               comm_class_metrics)


class IOSet(BaseModel):
    IO_output_path: OutputDirectory = Field(default="./output")

class RunSet(BaseModel):
    
    ngroup: Int = Field(default=0,description="Number of groups to run in parallel. If set to 0, all examples will be run in parallel.",ge=0)
    
    rundft_command: String = Field(default="OMP_NUM_THREADS=1 mpirun -np 16 abacus > log",
                          title = "Command to run each example. Please note that the program will first enter each folder before executing this command",
                          description = "",)

class NormalModel(comm_class.TrackingSet,
                  comm_class_metrics.SuperMetricsSet,
                  comm_class_metrics.MetricsSet,
                  comm_class_postdft.PostdftImageSet,
                  comm_class_rundft.RundftImageSet,
                  RunSet,
                  comm_class_exampleSource.ExampleSourceSet,
                  comm_class.ConfigSet,
                  IOSet,
                  BaseModel):
    ...  

def ReadSetting(logs:comm_class.myLog,opts:NormalModel,work_path,download_path):
    """
    {
        "config":{},
        "run_dft":{
            "image":
            "bohrium":{"scass_type","job_type","platform"},
            "example_source": "datahub",
            "urn":
            "example":[],
            "command":
        },
        "post_dft":{
            "iamge":
            "bohrium":
            "metrics":{"path":,"metric_name":},
            "super_metric":,
            "upload_datahub":,
            "upload_tracking":
        }
    }
    """
    logs.iprint("read config setting ...")
    config = comm_func.read_config(opts)

    #parse rundft
    run_dft = [{}]

    #parse example
    logs.iprint("read example setting ...")

    #download all examples to work path
    if comm_class_exampleSource.parse_example_source(example_source_set=opts,
                                                     download_path=work_path,
                                                     configs=opts,
                                                     logs=logs.iprint) == None:
        return None

    #collect the examples
    for ifile in os.listdir(work_path):
        if os.path.isdir(os.path.join(work_path,ifile)):
            run_dft[-1]["example"].append(ifile)

    #read rundft image and command
    logs.iprint("read run dft setting ...")
    logs.iprint("\tcommand:",opts.rundft_command)
    run_dft[-1]["command"] = opts.rundft_command
    if opts.ngroup > 0:
        run_dft[-1]["ngroup"] = opts.ngroup
    
    #read rundft image
    logs.iprint("\timage:",opts.rundft_image_set.image)
    for k,v in comm_class_rundft.parse_image_set(opts.rundft_image_set).items():
        run_dft[-1][k] = v
    if "bohrium" in run_dft[-1]:   
        logs.iprint("\tbohrium:",run_dft[-1]["bohrium"])

    #read postdft image
    logs.iprint("read post dft setting ...")
    logs.iprint("\timage:",opts.postdft_image_set.image)
    need_post_dft = False
    post_dft = comm_class_postdft.parse_image_set(opts.postdft_image_set)
    if "bohrium" in post_dft:
        logs.iprint("\tbohrium:",post_dft["bohrium"])

    #read metrics setting
    logs.iprint("read metrics setting ...")
    metrics_set = comm_class_metrics.parse_metrics_set(opts,opts)
    if "metrics" in metrics_set:
        post_dft["metrics"] = metrics_set["metrics"]
        post_dft["metrics"]["path"] = run_dft[-1]["example"]
        need_post_dft = True
    if "super_metric" in metrics_set:
        post_dft["super_metric"] = metrics_set["super_metric"]
        need_post_dft = True

    #read tracking setting
    if need_post_dft:
        tracking_set = comm_class.TrackingSet.parse_obj(opts)
        if tracking_set:
            config["AIM_ACCESS_TOKEN"] = tracking_set.get("token")
            post_dft["upload_tracking"] = {
                "tags": tracking_set.get("tags"),
                "name": tracking_set.get("name"),
                "experiment": tracking_set.get("experiment")
            }

    allparams = {"config": config,
            "run_dft": run_dft,
            "save_path": "."}
    if need_post_dft:
        allparams["post_dft"] = post_dft
    logs.iprint("read setting over!\n")
    return allparams

def NormalModelRunner(opts: NormalModel) -> int:
    logs = comm_class.myLog()

    paths = comm_func.create_path(str(opts.IO_output_path))
    output_path = paths["output_path"]
    work_path = paths["work_path"]
    download_path = paths["download_path"]

    allparams = ReadSetting(logs,opts,work_path,download_path)
    if allparams == None:
        return 1

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

    return 0
