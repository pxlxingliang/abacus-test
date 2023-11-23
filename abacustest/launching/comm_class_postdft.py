from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
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

from . import comm_class

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
    on_demand: Boolean = Field(default=False)
    
    #bohrium_machine_type: String = Field(default="c2_m4_cpu")
    #bohrium_job_type: String = Field(default="container")
    #bohrium_plat_form: String = Field(default="ali")


class PostdftImageSet(BaseModel):
    postdft_image_set: Union[Image,
                     ImageBohrium,
                     comm_class.ImageDispatcher] = Field(discriminator="type",
                                           description="IMAGE and MACHINE used in POST DFT")

class PostdftCommandSet(BaseModel):
    postdft_command: String = Field(default="",
                                    description="If you need to execute some custom scripts or commands, please enter the bash command here. \
Usually used to generate some custom metrics. At this step, the program will first collect the results of all examples in rundft, and then execute the command. \
The working directory is the same level directory as the outer layer of all examples.",)


def parse_image_set(image_set):
    if isinstance(image_set, ImageBohrium):
        return {
            "image": image_set.image,
            "bohrium": {
                "scass_type": image_set.bohrium_machine_type,
                "job_type": image_set.bohrium_job_type,
                "platform": image_set.bohrium_plat_form,
                "on_demand": 1 if image_set.on_demand else 0
            }
        }
    elif isinstance(image_set, Image):
        return {
            "image": image_set.image
        }
    elif isinstance(image_set, comm_class.ImageDispatcher):
        return comm_class.ImageDispatcher.construct_dispatcher(image_set)
    else:
        raise ValueError(f"Unknown image set type: {type(image_set)}")

def construct_input(datas,opts,logs):
    # datas is a dict of examples, created by comm_class_exampleSource.read_source
    #read postdft
    logs.iprint("read post dft setting ...")
    need_postdft = False
    post_dft = {}
    
    #read postdft example
    if datas.get("postdft_example"):
        post_dft["example"] = datas.get("postdft_example")    
        logs.iprint("\tpostdft example:",post_dft["example"])

    #read rundft extra files
    if datas.get("postdft_extrafile"):
        post_dft["extra_files"] = datas.get("postdft_extrafile")
        logs.iprint("\tpostdft_extrafile:",post_dft["extra_files"])
        
    #read postdft image
    if hasattr(opts,"postdft_image_set"):
        image_set = parse_image_set(opts.postdft_image_set)   
        logs.iprint("\timage:",image_set.get("image",""))
        for k,v in image_set.items():
            post_dft[k] = v
        if "bohrium" in post_dft:   
            logs.iprint("\tbohrium:",post_dft["bohrium"])
        elif "dispatcher" in post_dft:
            logs.iprint("\tdispatcher:",post_dft["dispatcher"]["resources_dict"])
    
    #read postdft command
    if hasattr(opts,"postdft_command"): 
        if opts.postdft_command != None and opts.postdft_command.strip() != "":
            post_dft["command"] = opts.postdft_command.strip()
            need_postdft = True
            logs.iprint("\tcommand:",post_dft["command"])
        
    #read postdft extra files
    if datas.get("postdft_extrafiles"):
        post_dft["extra_files"] = datas.get("postdft_extrafiles")
        logs.iprint("\tpostdft_extrafiles:",post_dft["extra_files"])
    
    return need_postdft, post_dft  