from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.typing import InputFilePath, OutputDirectory
from dp.launching.typing import (
    Boolean,
    Field
)
from dp.launching.report import Report

from . import comm_class,comm_func
import json,traceback,os,shutil,datetime,zipfile
from dp.metadata import MetadataContext, Dataset
from dp.metadata.utils.storage import TiefblueStorageClient
from abacustest.myflow.comm import CollectFileName
from dp.launching.report import ChartReportElement,ReportSection,AutoReportElement

datahub_project = "abacustest_dataset"

io_input_path_description = f"""A compressed file contains all example folders. 

For each folder is one example, and containing all the required files.

The dataset will be uploaded to https://datahub.mlops.dp.tech/browse/dataset/corp/tiefblue/{datahub_project}
"""

class UplaodDataset(BaseModel):
    IO_input_path:InputFilePath = Field(title="Upload dataset",
                                        st_kwargs_type = comm_func.unpack(None,None,get_support_filetype=True), 
                                        description=io_input_path_description,
                                        description_type="markdown")
    IO_output_path: OutputDirectory = Field(default="./output")
    name: String = Field(title="Dataset Name",regex="^\\s*[a-zA-Z0-9_-.]+\\s*$",description="Can only contains letters, numbers, dot(.), _ and -. (regex is: [a-zA-Z0-9_-.])")
    overwrite: Boolean = Field(description="If overwrite when the dataset already exists? Only owner can overwrite it.")
    description: String = Field(default="")

class UplaodDatasetModel(UplaodDataset,comm_class.ConfigSet,BaseModel):
    ...

def upload(download_path, bohrium_username, bohrium_password, bohrium_project, dataset_name, new_dataset_name_tail, tags, properties, description, overwrite, logs):
    '''
    dataset.description   dataset.entity_type   dataset.group         dataset.properties    dataset.uri           
    dataset.display_name  dataset.gen_urn(      dataset.owners        dataset.tags          dataset.urn
    '''
    cwd = os.getcwd()
    os.chdir(download_path)
    mess = description + "\\\n" + "\\\n".join(CollectFileName("."))
    os.chdir(cwd)

    metadata_storage_client = TiefblueStorageClient(
        bohrium_username, bohrium_password, bohrium_project)
    with MetadataContext(storage_client=metadata_storage_client, project=datahub_project) as context:
        def produce_unique_urn(dataset_name, bohrium_username=None, overwrite=False):
            org_name = dataset_name
            n = 1
            while True:
                urn = Dataset.gen_urn(context, platform="tiefblue",
                                      name=dataset_name, auto_suffix=False)
                dataset = client.get_dataset(urn)
                if (dataset == None) or ():
                    return dataset_name, urn
                elif overwrite and bohrium_username in dataset.tags:
                    logs.iprint(
                        f"\nDataset '{dataset_name}' already exists!\n\turn: {urn}\n\ttags: {dataset.tags}\nOverwrite it!")
                    return dataset_name, urn

                logs.iprint(
                    f"\nDataset '{dataset_name}' already exists!\n\turn: {urn}\n\ttags: {dataset.tags}")
                dataset_name = org_name + f"_{n}"
                n += 1
                logs.iprint(f"Rename dataset to '{dataset_name}'")

        client = context.client
        urn = Dataset.gen_urn(context, platform="tiefblue",
                              name=dataset_name, auto_suffix=False)
        old_dataset = client.get_dataset(urn)
        if old_dataset != None and (bohrium_username not in old_dataset.tags) and (not dataset_name.endswith(new_dataset_name_tail)):
            logs.iprint(
                f"\nDataset '{dataset_name}' already exists!\n\turn: {urn}\n\ttags: {old_dataset.tags}")
            dataset_name += new_dataset_name_tail
            logs.iprint(f"Rename the dataset to '{dataset_name}'")
        if overwrite:
            dataset_name, urn = produce_unique_urn(
                dataset_name, bohrium_username, True)
        else:
            dataset_name, urn = produce_unique_urn(dataset_name)

        uri = client.upload_artifact(None, None, download_path)
        dataset = Dataset(
            urn=urn,
            uri=uri,
            display_name=dataset_name,
            tags=tags,
            description=mess,
            properties=properties
        )
        urn1 = client.create_dataset(dataset)
        dataset = client.get_dataset(urn1)
        if dataset == None:
            logs.iprint("\nUpload to datahub failed!\nPlease contact developer to deal with it.")
            return False
        else:
            website = f"https://datahub.mlops.dp.tech/dataset/{urn1}/Documentation?is_lineage_mode=false"
            logs.iprint(f"\nUpload dataset to datahub successfully.")
            logs.iprint(f"\tname: {dataset.display_name}")
            logs.iprint(f"\turn : {dataset.urn}")
            logs.iprint(f"\ttags: {dataset.tags}")
            logs.iprint(f"\twebsite: {website}")
            return True

def UplaodDatasetModelRunner(opts:UplaodDatasetModel):
    download_path = "download_data"
    if os.path.isdir(download_path):
        shutil.rmtree(download_path)
    os.makedirs(download_path,exist_ok=True)
    os.makedirs(str(opts.IO_output_path),exist_ok=True)

    with zipfile.ZipFile(opts.IO_input_path.get_path(), "r") as zip_ref:
        zip_ref.extractall(download_path)

    today = datetime.datetime.now().strftime('%Y/%m/%d-%H:%M')
    new_dataset_name_tail = "_" + opts.Config_bohrium_username.split("@")[0]
    overwrite = opts.overwrite
    tags = [opts.Config_bohrium_username,opts.Config_bohrium_project_id,today.split("-")[0],today]
    properties = {"bohrium_project":opts.Config_bohrium_project_id}

    logs = comm_class.myLog()
    doupload = upload(download_path,
           opts.Config_bohrium_username, 
           opts.Config_bohrium_password, 
           opts.Config_bohrium_project_id, 
           opts.name.strip(), 
           new_dataset_name_tail, 
           tags, 
           properties, 
           opts.description,
           overwrite, 
           logs)

    logfname = "output.log"
    logs.write(os.path.join(str(opts.IO_output_path),logfname))

    report = Report(title="",
                    sections=[ReportSection(title="",
                              elements=[AutoReportElement(title='', path=logfname, description="")])],
                    description="")
    report.save(str(opts.IO_output_path))
    if doupload:
        return 0
    else:
        return 1




