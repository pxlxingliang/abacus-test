import traceback
from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.typing import (
    BaseModel,
    Field,
    InputFilePath
)
from enum import Enum
from typing import Literal
from . import comm_func,comm_class
import os,shutil


class DataSetsEnum(String, Enum):
    dataset1 = "dataset1-pw-v1.0"
    dataset2 = "dataset2-lcao-v1.0"

    @classmethod
    def GetAddress(cls, dataset):
        # find index of last -, and the string before it is the dataset name
        # the string after it is the version
        dataset_name = "-".join(dataset.split("-")[:-1])
        return f"https://launching.mlops.dp.tech/artifacts/datasets/{dataset_name}.abacustest/packages/{dataset}.tar.gz"


class ExampleFromPreUpload(BaseModel):
    type: Literal["from pre-upload examples"]


class ExampleFromDatahub(BaseModel):
    type: Literal["from datahub examples"]
    urn: String = Field(title="Datahub URN of examples",
                        description="Please enter the urn of the example.")


class ExampleFromDatasets(BaseModel):
    type: Literal["from datasets"]
    dataset: DataSetsEnum = Field(title="datasets",
                                  description="Please choose the datasets.")


class ExampleSourceSet(BaseModel):
    ExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all example folders. For each folder is one example, and containing all the required files. \
If you want to use the examples from datahub, please refer to the later 'Example Datahub Urn' section and there is no need to upload files here.
""",
                                         description_type="markdown")
    ExampleSource: Union[ExampleFromPreUpload,
                         ExampleFromDatahub,
                         ExampleFromDatasets] = Field(title="Example source",
                                                      discriminator="type",
                                                      description="Please choose the example source.")


def parse_example_source(example_source_set:ExampleSourceSet, download_path, configs: comm_class.ConfigSet, logs=None):
    if logs == None:
        logs = print

    example_source = example_source_set.ExampleSource
    upload_path = example_source_set.ExampleSource_local
    if isinstance(example_source, ExampleFromPreUpload):
        if upload_path == None or upload_path.strip() == "":
            logs("Please upload the example file.")
            return None
        else:
            try:
                comm_func.unpack(upload_path.get_path(), download_path)
            except:
                traceback.print_exc()
                logs(f"ERROR: The example file ({upload_path.get_path()}) is not valid!")
                logs(f"\tPlease check the example file!")
                return None
    elif isinstance(example_source, ExampleFromDatahub):
        try:
            dataset = comm_func.get_datahub_dataset(configs.Config_lbg_username,
                                                    configs.Config_lbg_password,
                                                    configs.Config_project_id,
                                                    example_source.urn)
            if dataset == None:
                logs(f"ERROR: The datahub urn ({example_source.urn}) is not valid!")
                logs(f"\tPlease check the datahub urn, and ensure that your Bohrium project ID has permission to access this data!")
                return None
        except:
            traceback.print_exc()
            return None
    elif isinstance(example_source, ExampleFromDatasets):
        url = DataSetsEnum.GetAddress(example_source.dataset)
        try:
            package = comm_func.download_url(url, download_path)
            if package == None:
                logs(f"ERROR: download dataset ({example_source.dataset}) failed!")
                logs(f"\tPlease check the dataset!")
                return None
            else:
                comm_func.unpack(package, download_path,filetype="tgz")
                #remove the package
                os.remove(package)
        except:
            traceback.print_exc()
            return None
    return download_path
