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
    #bohrium_machine_type: BohriumMachineType = Field(
    #    default=BohriumMachineType("c2_m4_cpu"))
    #bohrium_job_type: BohriumJobType = Field(default=BohriumJobType.CONTAINER)
    #bohrium_plat_form: BohriumPlatform = Field(default=BohriumPlatform.ALI)
    
    bohrium_machine_type: String = Field(default="c2_m4_cpu")
    bohrium_job_type: String = Field(default="container")
    bohrium_plat_form: String = Field(default="ali")


class PostdftImageSet(BaseModel):
    postdft_image_set: Union[Image,
                     ImageBohrium] = Field(discriminator="type",
                                           description="IMAGE and MACHINE used in POST DFT")

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