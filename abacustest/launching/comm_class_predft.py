from dp.launching.typing.basic import BaseModel, Int, String, Float, List, Optional, Union, Dict
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
    BohriumProjectId,
    BenchmarkLabels,
    BenchmarkTags
)
from enum import Enum
from typing import Literal
import re


class Image(BaseModel):
    type: Literal["Local Image"]
    image: String = Field(default="registry.dp.tech/dptech/abacustest:latest",
                          title="",
                          description="",)


class ImageBohrium(BaseModel):
    type: Literal["Use Bohrium Image"]

    image: String = Field(default=BohriumMachineType("registry.dp.tech/dptech/abacustest:latest"),
                          title="Bohrium Image Address",
                          description="",)
    bohrium_machine_type: BohriumMachineType = Field(
        default=BohriumMachineType("c2_m4_cpu"))
    bohrium_job_type: BohriumJobType = Field(default=BohriumJobType.CONTAINER)
    bohrium_plat_form: BohriumPlatform = Field(default=BohriumPlatform.ALI)
    
    #bohrium_machine_type: String = Field(default="c2_m4_cpu")
    #bohrium_job_type: String = Field(default="container")
    #bohrium_plat_form: String = Field(default="ali")


class PredftImageSet(BaseModel):
    predft_image_set: Union[Image,
                     ImageBohrium] = Field(discriminator="type",
                                           description="IMAGE and MACHINE used in PRE DFT")

class PredftGroupSizeSet(BaseModel):
    predft_group_size: Int = Field(
        default=0, description="Number of examples running in each parallel machine. If not set, all examples will be run by serial in one machine.", ge=0)


class PredftCommandSet(BaseModel):
    predft_command: String = Field(default="",
                                   description="Command of predft in each example. Please note that the program will first enter each folder before executing this command.\
If executing the command will generate some new example directories, please write these directories to a file with one line per directory. And enter the file in below \'predft_work_directories_filename\'",)
    
    predft_work_directories_filename: String = Field(default="",
                                                   description="If executing the predft_command will generate some new example directories, please write these directories to a file with one line per directory. And enter the file name here",)

def parse_image_set(image_set):
    if isinstance(image_set, ImageBohrium):
        return {
            "image": image_set.image,
            "bohrium": {
                "scass_type": image_set.bohrium_machine_type,
                "job_type": image_set.bohrium_job_type,
                "platform": image_set.bohrium_plat_form
            }
        }
    elif isinstance(image_set, Image):
        return {
            "image": image_set.image
        }
    else:
        raise ValueError(f"Unknown image set type: {type(image_set)}")

def construct_input(datas,opts,logs):
    # datas is a dict of examples, created by comm_class_exampleSource.read_source
    #read predft
    need_predft = False
    logs.iprint("read pre dft setting ...")
    pre_dft = {}
    
    #read predft example
    if datas.get("predft_example"):
        pre_dft["example"] = datas.get("predft_example")    
        logs.iprint("\texample:",pre_dft["example"])
    
    #read predft extra files
    if datas.get("predft_extrafile"):
        pre_dft["extra_files"] = datas.get("predft_extrafile")
        logs.iprint("\rundft_extrafile:",pre_dft["extra_files"])
        
    #read rundft command
    if hasattr(opts,"predft_command"):
        need_predft = True
        logs.iprint("\tcommand:",opts.predft_command)
        pre_dft["command"] = opts.predft_command
        pre_dft["work_directories_filename"] = opts.predft_work_directories_filename
    
    #read predft ngroup        
    #if hasattr(opts,"predft_group_size") and opts.predft_group_size > 0:
    #    logs.iprint("\tgroup size:",opts.predft_group_size)
    #    pre_dft["group_size"] = opts.predft_group_size
    
    #read predft image
    if hasattr(opts,"predft_image_set"):
        logs.iprint("\timage:",opts.predft_image_set.image)
        for k,v in comm_class_rundft.parse_image_set(opts.predft_image_set).items():
            pre_dft[k] = v
        if "bohrium" in pre_dft[-1]:   
            logs.iprint("\tbohrium:",pre_dft[-1]["bohrium"])
            
    return need_predft,pre_dft