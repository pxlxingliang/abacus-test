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
    image: String = Field(default="registry.dp.tech/deepmodeling/abacus-intel:latest",
                          title="",
                          description="",)


class ImageBohrium(BaseModel):
    type: Literal["Use Bohrium Image"]

    image: String = Field(default=BohriumMachineType("registry.dp.tech/deepmodeling/abacus-intel:latest"),
                          title="Bohrium Image Address",
                          description="",)
    bohrium_machine_type: BohriumMachineType = Field(default=BohriumMachineType("c32_m64_cpu"))
    bohrium_job_type: BohriumJobType = Field(default=BohriumJobType.CONTAINER)
    bohrium_plat_form: BohriumPlatform = Field(default=BohriumPlatform.ALI)
    #bohrium_machine_type: String = Field(default="c32_m64_cpu")
    #bohrium_job_type: String = Field(default="container")
    #bohrium_plat_form: String = Field(default="ali")


class RundftImageSet(BaseModel):
    rundft_image_set: Union[ImageBohrium,
                     Image] = Field(discriminator="type",
                                    description="IMAGE and MACHINE used in RUN DFT")

class RundftGroupSizeSet(BaseModel):
    rundft_group_size: Int = Field(
        default=1, description="Number of examples running in each parallel machine. If not set, all examples will be run in parallel.", ge=0)

class RundftCommandSet(BaseModel):
    rundft_command: String = Field(default="OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log",
                                   description="Command to run each example. Please note that the program will first enter each folder before executing this command. \
During runtime, there will be no interaction between different examples",)

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
    #read rundft
    need_rundft = False
    logs.iprint("read run dft setting ...")
    run_dft = [{}]
    
    #read rundft example
    if datas.get("rundft_example"):
        run_dft[-1]["example"] = datas.get("rundft_example")    
        logs.iprint("\texample:",run_dft[-1]["example"])
    
    #read rundft extra files
    if datas.get("rundft_extrafile"):
        run_dft[-1]["extra_files"] = datas.get("rundft_extrafile")
        logs.iprint("\rundft_extrafile:",run_dft[-1]["extra_files"])
        
    #read rundft command
    if hasattr(opts,"rundft_command"):
        need_rundft = True
        logs.iprint("\tcommand:",opts.rundft_command)
        run_dft[-1]["command"] = opts.rundft_command
    
    #read rundft ngroup        
    if hasattr(opts,"rundft_group_size") and opts.rundft_group_size > 0:
        logs.iprint("\tgroup size:",opts.rundft_group_size)
        run_dft[-1]["group_size"] = opts.rundft_group_size
    
    #read rundft image
    if hasattr(opts,"rundft_image_set"):
        logs.iprint("\timage:",opts.rundft_image_set.image)
        for k,v in parse_image_set(opts.rundft_image_set).items():
            run_dft[-1][k] = v
        if "bohrium" in run_dft[-1]:   
            logs.iprint("\tbohrium:",run_dft[-1]["bohrium"])

    return need_rundft, run_dft