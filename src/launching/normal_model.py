from ast import Dict
from enum import Enum
from typing import Literal

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


from . import comm_class,comm_func


io_input_path_description = """A zip file contains all example folders. For each folder is one example, and containing all the required files. 

If you want to use the examples from datahub, please refer to the later 'Example Datahub Urn' section and there is no need to upload files here.
"""

class IOSet(BaseModel):
    IO_input_path:InputFilePath = Field(default = None,
                                        title="Upload examples locally",
                                        st_kwargs_type = ["zip"], 
                                        description=io_input_path_description,
                                        description_type="markdown")
    IO_output_path: OutputDirectory = Field(default="./output")

class SuperMetricsSet(BaseModel):
    #name : String
    param_name: comm_class.AbacusMetricEnum 
    method: comm_class.SuperMetricMethodEnum
    #normalization: Boolean

class RundftBohriumImage(BaseModel):
    type: Literal["set Bohrium"]

    image: String = Field(default=BohriumMachineType("registry.dp.tech/deepmodeling/abacus-intel:latest"),
                                             title="Bohrium Image",
                                             description = "",)
    bohrium_machine_type: BohriumMachineType = Field(default=BohriumMachineType("c32_m64_cpu"))
    bohrium_job_type : BohriumJobType = Field(default=BohriumJobType.CONTAINER)
    bohrium_plat_form : BohriumPlatform = Field(default=BohriumPlatform.ALI)

class PostdftBohriumImage(BaseModel):
    type: Literal["set Bohrium"]

    image: String = Field(default=BohriumMachineType("registry.dp.tech/dptech/abacustest:latest"),
                                             title="Bohrium Image",
                                             description = "",)
    bohrium_machine_type: BohriumMachineType = Field(default=BohriumMachineType("c2_m4_cpu"))
    bohrium_job_type : BohriumJobType = Field(default=BohriumJobType.CONTAINER)
    bohrium_plat_form : BohriumPlatform = Field(default=BohriumPlatform.ALI)

class RundftImage(BaseModel):
    type: Literal["no extra setting in rundft"]
    image: String = Field(default="registry.dp.tech/deepmodeling/abacus-intel:latest",
                          title = "",
                          description = "",) 

class PostdftImage(BaseModel):
    type: Literal["no extra setting in postdft"]
    image: String = Field(default="registry.dp.tech/dptech/abacustest:latest",
                          title = "",
                          description = "",)

class NoneSet(BaseModel):
    type: Literal["no extra setting"]

class UplaodTrackingSet(BaseModel):
    type: Literal["upload the metrics to tracking"]
    AIM_ACCESS_TOKEN: String = Field(description="Token to access tracking")
    test_name: String
    experiment_name: String
    tags: List[String] =  Field(default=[],description="")

example_datahub_urn_description = """If you want to use a datahub example, please enter the urn of the example. 
Please note that if you fill in this field, the previously uploaded examples will be ignored.
"""
class RunSet(BaseModel):

    example_datahub_urn: String = Field(default="",
                                        title = "Datahub URN of examples",
                                        description = example_datahub_urn_description)
    
    rundft_command: String = Field(default="OMP_NUM_THREADS=1 mpirun -np 16 abacus > log",
                          title = "Command to run each example",
                          description = "",)

    rundft_image_set: Union[RundftBohriumImage,RundftImage] = Field(discriminator="type",
                                                                    description = "set the image used in RUN DFT step")
    
    

    postdft_image_set: Union[PostdftImage,PostdftBohriumImage] = Field(discriminator="type",
                                                                      description = "set the image used in POST DFT step")
    
    postdft_metrics: Set[comm_class.AbacusMetricEnum]

    postdft_super_metrics: List[SuperMetricsSet] = Field(default=[])

    tracking: Union[NoneSet,UplaodTrackingSet] = Field(discriminator="type",
                                                       title = "Upload metrics to TRACKING?",
                                                       description = "")

class NormalModel(IOSet,comm_class.ConfigSet,RunSet,BaseModel):
    ...  

def ReadSetting(opts:NormalModel,work_path,download_path,hasdatahub=False):
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
    print("read config setting ...")
    config = comm_func.read_config(opts)

    #parse rundft
    run_dft = [{}]

    #parse example
    print("read example setting ...")
    example_datahub = opts.example_datahub_urn.strip()
    example_local = opts.IO_input_path

    print("\texample_datahub_urn:",example_datahub)
    print("\texample_local:",example_local)

    if example_datahub == "" and example_local == None:
        print("Please upload the examplee locally or supply datahub urn")
        return None
    elif example_datahub != "":
        run_dft[-1]["example_source"] = "datahub"
        run_dft[-1]["urn"] = example_datahub
        run_dft[-1]["example"] = ["*"]
    else:
        run_dft[-1]["example"] = []
        with zipfile.ZipFile(example_local.get_path(), "r") as zip_ref:
            zip_ref.extractall(download_path)
        for ifile in os.listdir(download_path):
            if os.path.isdir(os.path.join(download_path,ifile)):
                run_dft[-1]["example"].append(ifile)
            shutil.move(os.path.join(download_path,ifile),os.path.join(work_path,ifile))

    #read rundft image and command
    print("read run dft image command setting ...")
    print("\timage:",opts.rundft_image_set.image)
    print("\tcommand:",opts.rundft_command)


    run_dft[-1]["image"] = opts.rundft_image_set.image
    run_dft[-1]["command"] = opts.rundft_command
    if isinstance(opts.rundft_image_set, RundftBohriumImage):
        run_dft[-1]["bohrium"] = {
            "scass_type": opts.rundft_image_set.bohrium_machine_type,
            "job_type": opts.rundft_image_set.bohrium_job_type,
            "platform": opts.rundft_image_set.bohrium_plat_form
        }
        print("\tbohrium:",run_dft[-1]["bohrium"])

    #read postdft image
    print("read post dft image command setting ...")
    print("\timage:",opts.postdft_image_set.image)
    need_post_dft = False
    post_dft = {"image":opts.postdft_image_set.image}
    if isinstance(opts.postdft_image_set, PostdftBohriumImage):
        post_dft["bohrium"] = {
            "scass_type": opts.postdft_image_set.bohrium_machine_type,
            "job_type": opts.postdft_image_set.bohrium_job_type,
            "platform": opts.postdft_image_set.bohrium_plat_form
        }
        print("\tbohrium:",post_dft["bohrium"])

    
    #read metrics setting
    print("read metrics setting ...")
    allexamplepath = run_dft[-1]["example"]
    post_dft["metrics"] = {
        "path": allexamplepath,
        "dft_type": "abacus",
        "metrics_name": list(opts.postdft_metrics),
        "save_file": "metrics.json"
    }
    
    if len(opts.postdft_super_metrics) > 0:
        print("read super metrics setting ...")
        post_dft["super_metrics"] = [{
            "save_file": "superMetrics.json",
            "result_file": ["metrics.json"],
            "metrics":[{"name": f"{i.method}({i.param_name})" ,
                        "param_name": i.param_name,
                        "method": i.method,
                        "normalization": False}
                       for i in opts.postdft_super_metrics]
        }]
        for i in opts.postdft_super_metrics:
            if i.param_name not in post_dft["metrics"]["metrics_name"]:
                post_dft["metrics"]["metrics_name"].append(i.param_name)
        post_dft["super_metrics"][-1]["outparams"] = [[j, [j], -1] for j in post_dft["metrics"]["metrics_name"]]

    if len(post_dft["metrics"]["metrics_name"]) > 0:
        need_post_dft = True

    #read tracking setting
    if need_post_dft:
        if isinstance(opts.tracking,UplaodTrackingSet):
            print("read tracking setting ...")
            config["AIM_ACCESS_TOKEN"] = opts.tracking.AIM_ACCESS_TOKEN.strip()
            post_dft["upload_tracking"] = {
                "tags": opts.tracking.tags,
                "name": opts.tracking.test_name,
                "experiment": opts.tracking.experiment_name
            }

    allparams = {"config": config,
            "run_dft": run_dft,
            "save_path": "results"}
    if need_post_dft:
        allparams["post_dft"] = post_dft
    print("read setting over!\n")
    return allparams

def NormalModelRunner(opts: NormalModel) -> int:
    paths = comm_func.create_path(str(opts.IO_output_path))
    output_path = paths["output_path"]
    work_path = paths["work_path"]
    download_path = paths["download_path"]

    allparams = ReadSetting(opts,work_path,download_path,hasdatahub=False)
    if allparams == None:
        return 1

    comm_func.exec_abacustest(allparams,work_path)
    reports = comm_func.produce_metrics_superMetrics_reports(allparams,work_path,output_path)

    if reports:
        report = Report(title="abacus test report",
                        sections=reports,
                        description="a report of abacustest")
        report.save(output_path)

    return 0