import traceback
from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.typing import (
    BaseModel,
    Field,
    InputFilePath,
    DataSet
)
from enum import Enum
from typing import Literal
from . import comm_func,comm_class
import os,shutil,glob

def GetDatasetAddress(package,dataset=None):
    # find index of last -, and the string before it is the dataset name
    # the string after it is the version
    if dataset == None:
        dataset = "-".join(package.split("-")[:-1])
    return f"https://launching.mlops.dp.tech/artifacts/datasets/{dataset}.abacustest/packages/{package}.tar.gz"
    

class DataSetsEnum(String, Enum):
    dataset1 = "dataset1-pw-v1.0"
    dataset2 = "dataset2-lcao-v1.0"
    dataset3 = "dataset3-exx-v1.0"

    @classmethod
    def GetAddress(cls, package,dataset=None):
        return GetDatasetAddress(package,dataset)

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
    dataset_unrecorded: String = Field(default=None,
                                       description="If the dataset you want is not in the above list, please enter download link of your datasets.")


class ExampleSourceSet(BaseModel):
    ExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload files locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all required files. \
If you want to use the examples from datahub or datasets, please refer to the later 'RundftExampleSource' section.""",
                                         description_type="markdown")
    ExampleSource: Union[FromPreUpload,
                         FromDatahub,
                         FromDatasets] = Field(title="Example source",
                                                      discriminator="type",
                                                      description="Please choose the example source.")
    Example: String = Field(default="*",title="Examples",description = "You can choose to use only partial files, and seperate each example with space. \
Tips: you can use regex to select files. For example: example_00[1-5]* example_[6,7,8]*. If you want to use all files, please type '*'.")  


class PrepareExampleSourceSet(BaseModel):
    PrepareExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload prepare-examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all example folders. For each folder is one example, and containing all the required files. \
If you want to use the examples from datahub or datasets, please refer to the later 'PrepareExampleSource' section.""",
                                         description_type="markdown")
    PrepareExampleSource: Union[FromPreUpload,
                                FromDatahub,
                         FromDatasets] = Field(title="Prepare Example source",
                                                      discriminator="type",
                                                      description="Please choose the example source.")
    PrepareExample: String = Field(default="*",title="Prepare Examples",description = "You can choose to run only partial examples of PrepareExampleSource, and separate each example with space. \
Tips: you can use regex to select examples. For example: example_00[1-5]* example_[6,7,8]*. If you want to run all examples, please type '*'.")     

class PredftExampleSourceSet(BaseModel):
    PredftExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload predft examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all example folders. For each folder is one example, and containing all the required files. \
If you want to use the examples from datahub or datasets, please refer to the later 'PredftExampleSource' section.""",
                                         description_type="markdown")
    PredftExampleSource: Union[FromPreUpload,
                         FromDatahub,
                         FromDatasets] = Field(title="Predft Example source",
                                                      discriminator="type",
                                                      description="Please choose the example source.")
    PredftExample: String = Field(default="*",title="Predft Examples",description = "You can choose to run only partial examples of PredftExampleSource, and separate each example with space. \
Tips: you can use regex to select examples. For example: example_00[1-5]* example_[6,7,8]*. If you want to run all examples, please type '*'.")                    


class RundftExampleSourceSet(BaseModel):
    RundftExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload rundft examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all example folders. For each folder is one example, and containing all the required files. \
If you want to use the examples from datahub or datasets, please refer to the later 'RundftExampleSource' section.""",
                                         description_type="markdown")
    RundftExampleSource: Union[FromPreUpload,
                         FromDatahub,
                         FromDatasets] = Field(title="Rundft Example source",
                                                      discriminator="type",
                                                      description="Please choose the example source.")
    RundftExample: String = Field(default="*",title="Rundft Examples",description = "You can choose to run only partial examples of RundftExampleSource, and separate each example with space. \
Tips: you can use regex to select examples. For example: example_00[1-5]* example_[6,7,8]*. If you want to run all examples, please type '*'.")                    


class PostdftExampleSourceSet(BaseModel):
    PostdftExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload postdft examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all required files in post dft. \
If you want to use the examples from datahub or datasets, please refer to the later 'PostdftExampleSource' section.""",
                                         description_type="markdown")
    PostdftExampleSource: Union[FromPreUpload,
                         FromDatahub,
                         FromDatasets] = Field(title="Postdft example source",
                                                      discriminator="type",
                                                      description="Please choose the example source.")
    PostdftExample: String = Field(default="*",title="Postdft examples",description = "You can choose to use only partial files of PostdftExampleSource, and separate each file or folder with space. \
Tips: you can use regex to select files. For example: example_00[1-5]* example_[6,7,8]*. If you want to use all files, please type '*'.")                    


class PrepareExtraFileSet(BaseModel):
    PrepareExtraFile_local: InputFilePath = Field(default=None,
                                                  title="Upload extra file for prepare locally",
                                                  st_kwargs_type=comm_func.unpack(
                                                      None, None, get_support_filetype=True),
                                                  description="""If some extra files that will be used for each newly created examples, you can upload the compressed file here.""")
    PrepareExtraFile: Union[
        NotRquired,
        FromPreUpload,
        FromDatahub,
        FromDatasets] = Field(title="Prepare Extra File Source",
                              discriminator="type",
                              description="Please choose the extra file source.")
    PrepareExtraFile_needed_files: String = Field(default=None,
                                                  title="Prepare extra files",
                                                  description="These files will be copied to each example directory. \
Please separate each file with space. Tips: you can use regex to select examples. \
If you need all files, please type '*', and if you do not need any files, please do not fill in anything here.")

class PredftExtraFileSet(BaseModel):
    PredftExtraFile_local: InputFilePath = Field(default=None,
                                                 title="Upload Predft extra file locally",
                                                 st_kwargs_type=comm_func.unpack(
                                                     None, None, get_support_filetype=True),
                                                 description="""If you need some extra files that may be used in predft, you can upload the compressed file here.""")
    PredftExtraFile: Union[
        NotRquired,
        FromPreUpload,
        FromDatahub,
        FromDatasets] = Field(title="Predft Extra File Source",
                             discriminator="type",
                             description="Please choose the extra file source.")
    PredftExtraFile_needed_files: String = Field(default=None,
                                                 title="Predft extra files",
                                                 description="Before executing the predft_command, these files will be copied to each example directory. \
Please separate each file with space. Tips: you can use regex to select examples. \
If you need all files, please type '*', and if you do not need any files, please do not fill in anything here.")


class RundftExtraFileSet(BaseModel):
    RundftExtraFile_local: InputFilePath = Field(default=None,
                                                 title="Upload rundft extra file locally",
                                                 st_kwargs_type=comm_func.unpack(
                                                     None, None, get_support_filetype=True),
                                                 description="""If you need some extra files that may be used in rundft, you can upload the compressed file here.""")
    RundftExtraFile: Union[
        NotRquired,
        FromPreUpload,
        FromDatahub,
        FromDatasets] = Field(title="Rundft Extra File Source",
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
        FromDatahub,
        FromDatasets] = Field(title="Postdft Extra File Source",
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
        if example_source.dataset_unrecorded != None and example_source.dataset_unrecorded.strip() != "":
            url = example_source.dataset_unrecorded.strip()
        else:
            url = DataSetsEnum.GetAddress(example_source.dataset)
        
        try:
            package = comm_func.download_url(url, download_path)
            if package == None:
                logs(f"ERROR: download dataset failed!\n\turl:{url}")
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

def copy_download_to_work(download_path,work_path,needed_files):
    #needed_files is a string, each file is separated by space
    cwd = os.getcwd()
    os.chdir(download_path)
    alldirectories = []
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
                alldirectories.append(iifile)
                shutil.copytree(iifile,os.path.join(work_path,iifile))
            else:
                shutil.copy(iifile,os.path.join(work_path,iifile))
            allfiles.append(iifile)

    os.chdir(cwd)
    if alldirectories == []:
        alldirectories = None
    else:
        alldirectories = list(set(alldirectories))
        alldirectories.sort()
    if allfiles == []:
        allfiles = None
    else:   
        allfiles = list(set(allfiles))
        allfiles.sort()
    return alldirectories,allfiles

def read_source(opts,work_path,download_path,logs=None):
    # read setting in opts, and dowload example/rundft_extrafile/postdft_extrafile to downlaod path
    # and copy selected files to work path
    # return a dictory, which contains the list of example/rundft_extrafile/postdft_extrafile
    outdict = {}
    
    if logs == None:
        logs = print
    
    if hasattr(opts,"ExampleSource"):
        tmp = parse_source(opts.ExampleSource,opts.ExampleSource_local,download_path,opts,logs)
        if tmp == None:
            return None
        else:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.Example)
            outdict["all_files"] = all_files
        comm_func.clean_dictorys(download_path)
        
    if hasattr(opts,"PrepareExampleSource"):
        tmp = parse_source(opts.PrepareExampleSource,opts.PrepareExampleSource_local,download_path,opts,logs)
        if tmp == None:
            return None
        else:
            #put all examples to example_template
            all_directories, all_files = copy_download_to_work(
                download_path, os.path.join(work_path,"example_template"), opts.PrepareExample)
            if all_directories:
                outdict["prepare_example"] = [os.path.join("example_template",i) for i in all_directories]
            else:
                outdict["prepare_example"] = None
        comm_func.clean_dictorys(download_path) 
        
    if hasattr(opts,"PredftExampleSource"):
        tmp = parse_source(opts.PredftExampleSource,opts.PredftExampleSource_local,download_path,opts,logs)
        if tmp == None:
            return None
        else:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.PredftExample)
            outdict["predft_example"] = all_directories
        comm_func.clean_dictorys(download_path)  
    
    if hasattr(opts,"RundftExampleSource"):
        tmp = parse_source(opts.RundftExampleSource,opts.RundftExampleSource_local,download_path,opts,logs)
        if tmp == None:
            return None
        else:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.RundftExample)
            outdict["rundft_example"] = all_directories
        comm_func.clean_dictorys(download_path)    
    
    #read postdft example source
    if hasattr(opts,"PostdftExampleSource"):
        tmp = parse_source(opts.PostdftExampleSource,opts.PostdftExampleSource_local,download_path,opts,logs)
        if tmp == None:
            return None
        else:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.PostdftExample)
            outdict["postdft_example"] = all_directories
        comm_func.clean_dictorys(download_path)  

    #read prepare extra files
    if hasattr(opts,"PrepareExtraFile_needed_files"):
        tmp = parse_source(opts.PrepareExtraFile,opts.PrepareExtraFile_local,download_path,opts,logs)
        if tmp == None:
            return None
        if tmp:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.PrepareExtraFile_needed_files)
            outdict["prepare_extrafile"] = all_files
        comm_func.clean_dictorys(download_path)
        
    #read predft extra files
    if hasattr(opts,"PredftExtraFile_needed_files"):
        tmp = parse_source(opts.PredftExtraFile,opts.PredftExtraFile_local,download_path,opts,logs)
        if tmp == None:
            return None
        if tmp:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.PredftExtraFile_needed_files)
            outdict["predft_extrafile"] = all_files
        comm_func.clean_dictorys(download_path)
        
    #read rundft extra files
    if hasattr(opts,"RundftExtraFile_needed_files"):
        tmp = parse_source(opts.RundftExtraFile,opts.RundftExtraFile_local,download_path,opts,logs)
        if tmp == None:
            return None
        if tmp:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.RundftExtraFile_needed_files)
            outdict["rundft_extrafile"] = all_files
        comm_func.clean_dictorys(download_path)
    
    #read postdft extra files
    if hasattr(opts,"PostdftExtraFile_needed_files"):
        tmp = parse_source(opts.PostdftExtraFile,opts.PostdftExtraFile_local,download_path,opts,logs)
        if tmp == None:
            return None
        if tmp:
            all_directories, all_files = copy_download_to_work(
                download_path, work_path, opts.PostdftExtraFile_needed_files)
            outdict["postdft_extrafile"] = all_files
        comm_func.clean_dictorys(download_path)
    return outdict