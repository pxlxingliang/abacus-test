from enum import Enum
from genericpath import isdir
from pdb import run
from typing import Literal

from requests import post
from sqlalchemy import exists
from zmq import has

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
import os,zipfile,shutil,json,glob

from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement


from . import comm_class,comm_func


io_input_path_description = """A compressed file contains all example folders. For each folder is one example, and containing all the required files. 

If you want to use the examples from datahub, please refer to the later 'Example Datahub Urn' section and there is no need to upload files here.
"""

class IOSet(BaseModel):
    IO_input_path:InputFilePath = Field(default = None,
                                        title="Upload examples locally (optional)",
                                        st_kwargs_type = comm_func.unpack(None,None,get_support_filetype=True), 
                                        description=io_input_path_description,
                                        description_type="markdown")
    IO_input_path_extrafiles:InputFilePath = Field(default = None,
                                        title="Extra files (optional)",
                                        st_kwargs_type = comm_func.unpack(None,None,get_support_filetype=True), 
                                        description="Some others files needed in the examples. For example: customized script files. \
Please package the files required for rundft and postdft in one file.",
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

example_datahub_urn_description = """
"""
class RunSet(BaseModel):
    example_datahub_urn: String = Field(default="",
                                        title = "Example URN",
                                        description = "If you want to use a datahub example, please enter the urn of the example. \
Please note that if you fill in this field, the previously uploaded examples will be ignored.")

    example: String = Field(default="*",description = "You can choose to run only partial examples, and separate each example with space. \
Tips: you can use regex to select examples. For example: example_00[1-5]* example_[6,7,8]*. If you want to run all examples, please type '*'.")
    
    rundft_extrafiles_urn: String = Field(default=None,
                          title = "Rundft extra file URN",
                          description = "If you want to use the local files uploaded in Extra files, please not fill in anything here",)
    
    rundft_extrafiles: String = Field(default=None,
                                      title = "Rundft extra file",
                          description = "Before executing the rundft_command, these files will be copied to each example directory. \
Please separate each file with space. Tips: you can use regex to select examples. \
If you need all files, please type '*', and if you do not need any files, please do not fill in anything here.",)
    
    rundft_command: String = Field(default="OMP_NUM_THREADS=1 mpirun -np 16 abacus > log",
                          description = "Command to run each example. Please note that the program will first enter each folder before executing this command. \
During runtime, there will be no interaction between different examples",)
    
    rundft_image_set: Union[RundftBohriumImage,RundftImage] = Field(discriminator="type",
                                                                    description = "set the image used in RUN DFT step")
    
    postdft_extrafiles_urn: String = Field(default=None,
                                           title = "Postdft extra file URN",
                          description = "If you want to use the files uploaded in Extra files, please not fill in anything here",)
    
    postdft_extrafiles: String = Field(default=None,
                                       title = "Postdft extra file",
                          description = "Before executing the postdft_command, these files will be copied to the work directory. \
Please separate each file with space. Tips: you can use regex to select files. \
If you need all files, please type '*', and if you do not need any files, please do not fill in anything here.",)
    
    postdft_command: String = Field(default="",
                          description = "If you need to execute some custom scripts or commands, please enter the bash command here. \
Usually used to generate some custom metrics. At this step, the program will first collect the results of all examples in rundft, and then execute the command. \
The working directory is the same level directory as the outer layer of all examples.",)
    
    postdft_image_set: Union[PostdftImage,PostdftBohriumImage] = Field(discriminator="type",
                                                                      description = "set the image used in POST DFT step")
    
    postdft_metrics: Set[comm_class.AbacusMetricEnum] = Field(default=None,
                                                              description = "Choose the metrics to be collected.")

    postdft_super_metrics: List[SuperMetricsSet] = Field(default=[],description = "Choose the super metrics to be calculated.")

    postdft_metrics_savefile: String = Field(default=None,
                                             description = "If you need to display or upload the metrics generated by executing the postdft command, please enter the file name here. \
Should be a json file, and the key is the name of the metric, and the value is the value of the metric.")
    
    postdft_super_metrics_savefile: String = Field(default=None,
                                             description = "If you need to display or upload the super metrics generated by executing the postdft command, please enter the file name here. \
Should be a json file, and the key is the name of the metric, and the value is the value of the metric.")

    tracking: Union[NoneSet,UplaodTrackingSet] = Field(discriminator="type",
                                                       title = "Upload metrics to TRACKING? ",
                                                       description = "If uploading, it can be compared with historical data")

class AdvancedModel(IOSet,comm_class.ConfigSet,RunSet,BaseModel):
    ...  

def clean_dictorys(ipath):
    for ifile in glob.glob(os.path.join(ipath,"*")):
        if os.path.isdir(ifile):
            shutil.rmtree(ifile)
        else:
            os.remove(ifile)

def copy_download_to_work(download_path,work_path,needed_files):
    cwd = os.getcwd()
    os.chdir(download_path)
    alldictories = []
    allfiles = []

    for ifile in needed_files:
        for iifile in glob.glob(ifile):
            ipath,iiifile = os.path.split(iifile)
            if os.path.isdir(os.path.join(work_path,ipath)) == False:
                os.makedirs(os.path.join(work_path,ipath),exist_ok=True)
            if os.path.exists(os.path.join(work_path,iifile)):
                if os.path.isdir(iifile):
                    shutil.rmtree(os.path.join(work_path,iifile))
                else:
                    os.remove(os.path.join(work_path,iifile))
            if os.path.isdir(iifile):
                alldictories.append(iifile)
                shutil.copytree(iifile,os.path.join(work_path,iifile))
            else:
                shutil.copy(iifile,os.path.join(work_path,iifile))
            allfiles.append(iifile)

    os.chdir(cwd)
    if alldictories == []:
        alldictories = None
    else:
        alldictories = list(set(alldictories))
        alldictories.sort()
    if allfiles == []:
        allfiles = None
    else:   
        allfiles = list(set(allfiles))
        allfiles.sort()
    return alldictories,allfiles

def prepare_files(logs,needed_files,urn,io_input_path,
                  lbg_username,lbg_password,lbg_project_id,
                  download_path,work_path):
    if needed_files == None or len(needed_files) == 0:
        return None,None
    if urn != None and urn.strip() != "":
        dataset = comm_func.get_datahub_dataset(lbg_username, 
                                         lbg_password, 
                                         lbg_project_id, 
                                         urn,download_path)
        if dataset == None:
            logs.iprint(f"ERROR: The datahub urn ({urn}) is not valid!")
            logs.iprint(f"\tPlease check the datahub urn, and ensure that your Bohrium project ID has permission to access this data!")
            return None,None
        else:
            logs.iprint(f"download the dataset to local ...")
    elif io_input_path != None:
        logs.iprint(f"unpack the uploaded dataset ...")
        try:
            comm_func.unpack(io_input_path.get_path(),download_path)
        except:
            logs.iprint(f"ERROR: Unpack the uploaded dataset failed!")
            logs.iprint(f"\tPlease check the uploaded dataset!")
            return None,None
    else:
        return None,None
    #move all files to work path
    logs.iprint(f"copy the needed files to work path ...")
    alldictories,allfiles = copy_download_to_work(download_path,work_path,needed_files)
    clean_dictorys(download_path)
    return alldictories,allfiles


def ReadSetting(logs:comm_class.myLog,opts:AdvancedModel,work_path,download_path):
    """
    {
        "config":{},
        "run_dft":[{
            "image":
            "bohrium":{"scass_type","job_type","platform"},
            "example":[],
            "extra_files":[],
            "command":
        }],
        "post_dft":{
            "command":,
            "extra_files":,
            "iamge":,
            "bohrium":,
            "metrics":{
                "path":,
                "metric_name":
                "value_from_file":,},
            "super_metric":{
                "value_from_file"},
            "upload_datahub":,
            "upload_tracking":
        }
    }
    """
    logs.iprint("read config setting ...")
    config = comm_func.read_config(opts)

    #parse rundft
    run_dft = [{}]

    #download examples
    logs.iprint("read example setting ...")
    example_datahub = opts.example_datahub_urn.strip()
    example_local = opts.IO_input_path

    #get examples
    if opts.example == None or opts.example.strip() == "":
        logs.iprint("Please input the examples you want to run!")
        return None
    else:
        needed_example = opts.example.strip().split()
    logs.iprint("\texample:",needed_example)
    logs.iprint("\texample_datahub_urn:" + str(example_datahub))
    logs.iprint("\texample_local:" + str(example_local))

    examples_name,tmp = prepare_files(logs,needed_example,example_datahub,example_local,
                                      opts.Config_lbg_username,opts.Config_lbg_password,opts.Config_project_id,
                                      download_path,work_path)
    if examples_name == None:
        logs.iprint("Please upload the examples locally or supply datahub urn!")
        return None
    run_dft[-1]["example"] = examples_name
    logs.iprint("\texample:",examples_name)

    #get rundft extra files
    rundft_extrafiles_needed = None
    if opts.rundft_extrafiles != None and opts.rundft_extrafiles.strip() != "":
        rundft_extrafiles_needed = opts.rundft_extrafiles.strip().split()
    tmp,rundft_extrafiles_name = prepare_files(logs,rundft_extrafiles_needed,opts.rundft_extrafiles_urn,opts.IO_input_path_extrafiles,
                                               opts.Config_lbg_username,opts.Config_lbg_password,opts.Config_project_id,
                                               download_path,work_path)
    if rundft_extrafiles_name:
        run_dft[-1]["extra_files"] = rundft_extrafiles_name
        logs.iprint("\trundft_extrafiles:",rundft_extrafiles_name)

    #read rundft image and command
    logs.iprint("read run dft image command setting ...")
    logs.iprint("\timage:",opts.rundft_image_set.image)
    logs.iprint("\tcommand:",opts.rundft_command)
    run_dft[-1]["image"] = opts.rundft_image_set.image
    run_dft[-1]["command"] = opts.rundft_command
    if isinstance(opts.rundft_image_set, RundftBohriumImage):
        run_dft[-1]["bohrium"] = {
            "scass_type": opts.rundft_image_set.bohrium_machine_type,
            "job_type": opts.rundft_image_set.bohrium_job_type,
            "platform": opts.rundft_image_set.bohrium_plat_form
        }
        logs.iprint("\tbohrium:",run_dft[-1]["bohrium"])

    #read postdft image
    logs.iprint("read post dft image command setting ...")
    logs.iprint("\timage:",opts.postdft_image_set.image)
    need_post_dft = False
    post_dft = {"image":opts.postdft_image_set.image}
    if isinstance(opts.postdft_image_set, PostdftBohriumImage):
        post_dft["bohrium"] = {
            "scass_type": opts.postdft_image_set.bohrium_machine_type,
            "job_type": opts.postdft_image_set.bohrium_job_type,
            "platform": opts.postdft_image_set.bohrium_plat_form
        }
        logs.iprint("\tbohrium:",post_dft["bohrium"])
    if opts.postdft_command != None and opts.postdft_command.strip() != "":
        post_dft["command"] = opts.postdft_command.strip()
        need_post_dft = True
        logs.iprint("\tcommand:",post_dft["command"])

    #get postdft extra files
    postdft_extrafiles_needed = None
    if opts.postdft_extrafiles != None and opts.postdft_extrafiles.strip() != "":
        postdft_extrafiles_needed = opts.postdft_extrafiles.strip().split()
    tmp,postdft_extrafiles_name = prepare_files(logs,postdft_extrafiles_needed,opts.postdft_extrafiles_urn,opts.IO_input_path_extrafiles,
                                                opts.Config_lbg_username,opts.Config_lbg_password,opts.Config_project_id,
                                                download_path,work_path)
    if postdft_extrafiles_name:
        post_dft["extra_files"] = postdft_extrafiles_name
        logs.iprint("\tpostdft_extrafiles:",postdft_extrafiles_name)

    #read metrics setting
    need_metrcis = False
    need_super_metrics = False
    logs.iprint("read metrics setting ...")
    allexamplepath = run_dft[-1]["example"]
    metrics_name = []
    if opts.postdft_metrics != None and len(opts.postdft_metrics) > 0:
        need_metrcis = True
        metrics_name += list(opts.postdft_metrics)
    if len(opts.postdft_super_metrics) > 0:
        need_super_metrics = True
        need_metrcis = True
        logs.iprint("read super metrics setting ...")
        #complete metrics
        for i in opts.postdft_super_metrics:
            if i.param_name not in metrics_name:
                metrics_name.append(i.param_name)

    #convert some special metrics (such as: KEY1:KEY2,..)
    metrics_name= comm_func.convert_metrics(metrics_name)

    #set metrics
    metrics = {}
    if len(metrics_name) > 0:
        need_metrcis = True
        metrics["path"] = allexamplepath
        metrics["dft_type"] = "abacus"
        metrics["metrics_name"] = metrics_name
        metrics["save_file"] = "metrics.json"

    if opts.postdft_metrics_savefile != None and opts.postdft_metrics_savefile.strip() != "":
        need_metrcis = True
        metrics["value_from_file"] = opts.postdft_metrics_savefile.strip()
        metrics["save_file"] = "metrics.json" if metrics["value_from_file"] != "metrics.json" else "metrics_1.json"
    
    #set super metrics
    super_metrics = {}
    if len(opts.postdft_super_metrics) > 0:
        need_super_metrics = True
        super_metrics["save_file"] = "superMetrics.json"
        super_metrics["result_file"] = [metrics["save_file"]]
        super_metrics["metrics"] = []
        super_metrics["outparams"] = []

        for i in opts.postdft_super_metrics:
            metric_name = comm_func.convert_supermetrics_metrics_name(i.param_name)
            super_metrics["metrics"].append({
                "name": f"{i.method}({metric_name})" ,
                "param_name": metric_name,
                "method": i.method,
                "normalization": False})
            super_metrics["outparams"].append([metric_name, [metric_name], -1])

    if opts.postdft_super_metrics_savefile != None and opts.postdft_super_metrics_savefile.strip() != "":
        need_super_metrics = True
        super_metrics["value_from_file"] = opts.postdft_super_metrics_savefile.strip()
        super_metrics["save_file"] = "superMetrics.json" if super_metrics["value_from_file"] != "superMetrics.json" else "superMetrics_1.json"

    if need_metrcis:
        post_dft["metrics"] = metrics
        need_post_dft = True
    if need_super_metrics:
        post_dft["super_metrics"] = [super_metrics]
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
    logs.iprint("read setting over!\n")
    return allparams

def AdvancedModelRunner(opts: AdvancedModel) -> int:
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
