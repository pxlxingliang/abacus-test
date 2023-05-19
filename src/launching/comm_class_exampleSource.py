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
import os,shutil,glob


class DataSetsEnum(String, Enum):
    dataset1 = "dataset1-pw-v1.0"
    dataset2 = "dataset2-lcao-v1.0"

    @classmethod
    def GetAddress(cls, dataset):
        # find index of last -, and the string before it is the dataset name
        # the string after it is the version
        dataset_name = "-".join(dataset.split("-")[:-1])
        return f"https://launching.mlops.dp.tech/artifacts/datasets/{dataset_name}.abacustest/packages/{dataset}.tar.gz"

class NotRquired(BaseModel):
    type: Literal["not required"]

class FromPreUpload(BaseModel):
    type: Literal["from pre-upload examples"]


class FromDatahub(BaseModel):
    type: Literal["from datahub examples"]
    urn: String = Field(title="Datahub URN of examples",
                        description="Please enter the urn of the example.")


class FromDatasets(BaseModel):
    type: Literal["from datasets"]
    dataset: DataSetsEnum = Field(title="datasets",
                                  description="Please choose the datasets.")


class ExampleSourceSet(BaseModel):
    ExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all example folders. For each folder is one example, and containing all the required files. \
If you want to use the examples from datahub or datasets, please refer to the later 'ExampleSource' section.""",
                                         description_type="markdown")
    ExampleSource: Union[FromPreUpload,
                         FromDatahub,
                         FromDatasets] = Field(title="Example source",
                                                      discriminator="type",
                                                      description="Please choose the example source.")
    Example: String = Field(default="*",title="Examples",description = "You can choose to run only partial examples, and separate each example with space. \
Tips: you can use regex to select examples. For example: example_00[1-5]* example_[6,7,8]*. If you want to run all examples, please type '*'.")                    


class RundftExtraFileSet(BaseModel):
    RundftExtraFile_local: InputFilePath = Field(default=None,
                                                 title="Upload rundft extra file locally",
                                                 st_kwargs_type=comm_func.unpack(
                                                     None, None, get_support_filetype=True),
                                                 description="""If you need some extra files that may be used in rundft, you can upload the compressed file here.""")
    RundftExtraFile: Union[
        NotRquired,
        FromPreUpload,
        FromDatahub] = Field(title="Rundft Extra File Source",
                             discriminator="type",
                             description="Please choose the extra file source.")
    RundftExtraFile_needed_files: String = Field(default=None,
                                                 title="Rundft extra files",
                                                 description="Before executing the rundft_command, these files will be copied to each example directory. \
Please separate each file with space. Tips: you can use regex to select examples. \
If you need all files, please type '*', and if you do not need any files, please do not fill in anything here.")


class PostdftExtraFileSet(BaseModel):
    PostdftExtraFile_local: InputFilePath = Field(default=None,
                                                  title="Upload postdft extra file locally",
                                                  st_kwargs_type=comm_func.unpack(
                                                      None, None, get_support_filetype=True),
                                                  description="""If you need some extra files that may be used in postdft, you can upload the compressed file here.""")
    PostdftExtraFile: Union[
        NotRquired,
        FromPreUpload,
        FromDatahub] = Field(title="Postdft Extra File Source",
                             discriminator="type",
                             description="Please choose the extra file source.")
        
    PostdftExtraFile_needed_files: String = Field(default=None,
                                                  title="Postdft extra files",
                                                  description="Before executing the postdft_command, these files will be copied to the work directory. \
Please separate each file with space. Tips: you can use regex to select files. \
If you need all files, please type '*', and if you do not need any files, please do not fill in anything here.")
    

def parse_source(example_source,upload_path,download_path,configs: comm_class.ConfigSet, logs=None):
    if logs == None:
        logs = print
        
    if isinstance(example_source, FromPreUpload):
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
    elif isinstance(example_source, FromDatahub):
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
    elif isinstance(example_source, FromDatasets):
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
    elif isinstance(example_source, NotRquired):
        return False
    else:
        logs(f"ERROR: The example source ({example_source}) is not valid!")
        return None
    return download_path
    
def parse_example_source(example_source_set:ExampleSourceSet, download_path, configs: comm_class.ConfigSet, logs=None):
    example_source = example_source_set.ExampleSource
    upload_path = example_source_set.ExampleSource_local
    return parse_source(example_source,upload_path,download_path,configs,logs)

def parse_rundft_extrafile_source(extrafile_source_set:RundftExtraFileSet, download_path, configs: comm_class.ConfigSet, logs=None):
    example_source = extrafile_source_set.RundftExtraFile
    upload_path = extrafile_source_set.RundftExtraFile_local
    return parse_source(example_source,upload_path,download_path,configs,logs)

def parse_postdft_extrafile_source(extrafile_source_set:PostdftExtraFileSet, download_path, configs: comm_class.ConfigSet, logs=None):
    example_source = extrafile_source_set.PostdftExtraFile
    upload_path = extrafile_source_set.PostdftExtraFile_local
    return parse_source(example_source,upload_path,download_path,configs,logs)

def copy_download_to_work(download_path,work_path,needed_files):
    #needed_files is a string, each file is separated by space
    cwd = os.getcwd()
    os.chdir(download_path)
    alldictories = []
    allfiles = []

    for ifile in needed_files.split():
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

def read_source(opts,work_path,download_path,logs=None):
    # read setting in opts, and dowload example/rundft_extrafile/postdft_extrafile to downlaod path
    # and copy selected files to work path
    # return a dictory, which contains the list of example/rundft_extrafile/postdft_extrafile
    outdict = {}
    
    if logs == None:
        logs = print
    
    if hasattr(opts,"ExampleSource"):
        tmp = parse_example_source(opts,download_path,opts,logs)
        if tmp == None:
            return None
        else:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.Example)
            outdict["example"] = all_directories
        comm_func.clean_dictorys(download_path)    

    #read rundft extra files
    if hasattr(opts,"RundftExtraFile_needed_files"):
        tmp = parse_rundft_extrafile_source(opts,download_path,opts,logs.iprint)
        if tmp == None:
            return None
        if tmp:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.RundftExtraFile_needed_files)
            outdict["rundft_extrafile"] = all_files
        comm_func.clean_dictorys(download_path)
    
    #read postdft extra files
    if hasattr(opts,"PostdftExtraFile_needed_files"):
        tmp = parse_postdft_extrafile_source(opts,download_path,opts,logs.iprint)
        if tmp == None:
            return None
        if tmp:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.PostdftExtraFile_needed_files)
            outdict["postdft_extrafile"] = all_files
        comm_func.clean_dictorys(download_path)
    return outdict
