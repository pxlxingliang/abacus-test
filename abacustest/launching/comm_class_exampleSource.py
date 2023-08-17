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

class DatasetSet(BaseModel):
    dataset: DataSet = Field(title=None,
                             default=None,
                            description="Please enter your dataset in launching.")

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
        if package == "dataset1-pw-v1.0":
            return "https://launching.mlops.dp.tech/download/artifacts/datasets/abacusdailytest.abacustest/packages/v1_0-dataset1-pw.tar.gz"
        elif package == "dataset2-lcao-v1.0":
            return "https://launching.mlops.dp.tech/download/artifacts/datasets/abacusdailytest.abacustest/packages/v1_0-dataset2-lcao.tar.gz"
        elif package == "dataset3-exx-v1.0":
            return "https://launching.mlops.dp.tech/download/artifacts/datasets/abacusdailytest.abacustest/packages/v1_0-dataset3-exx.tar.gz"
        else:
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

class ExampleSourceSet(BaseModel):
    ExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload files locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all required files.""",
                                         description_type="markdown")
    
class ExampleSet(BaseModel):
    Example: String = Field(default="*",title="Examples",description = "You can choose to use only partial files, and seperate each example with space. \
Tips: you can use regex to select files. For example: example_00[1-5]* example_[6,7,8]*. If you want to use all files, please type '*'.")  


class PrepareExampleSourceSet(BaseModel):
    PrepareExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload prepare-examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all example folders. For each folder is one example, and containing all the required files. \
If you want to use the examples from datasets, please refer to the later 'PrepareExampleSource' section.""",
                                         description_type="markdown")

class PrepareExampleSet(BaseModel):                         
    PrepareExample: String = Field(default="*",title="Prepare Examples",description = "You can choose to run only partial examples of PrepareExampleSource, and separate each example with space. \
Tips: you can use regex to select examples. For example: example_00[1-5]* example_[6,7,8]*. If you want to run all examples, please type '*'.")     

class PredftExampleSourceSet(BaseModel):
    PredftExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload predft examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all example folders. For each folder is one example, and containing all the required files. \
If you want to use the examples from datasets, please refer to the later 'PredftExampleSource' section.""",
                                         description_type="markdown")
                         
class PredftExampleSet(BaseModel):                         
    PredftExample: String = Field(default="*",title="Predft Examples",description = "You can choose to run only partial examples of PredftExampleSource, and separate each example with space. \
Tips: you can use regex to select examples. For example: example_00[1-5]* example_[6,7,8]*. If you want to run all examples, please type '*'.")                    


class RundftExampleSourceSet(BaseModel):
    RundftExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload rundft examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all example folders. For each folder is one example, and containing all the required files. \
If you want to use the examples from datasets, please refer to the later 'RundftExampleSource' section.""",
                                         description_type="markdown")
                         
class RundftExampleSet(BaseModel):                         
    RundftExample: String = Field(default="*",title="Rundft Examples",description = "You can choose to run only partial examples of RundftExampleSource, and separate each example with space. \
Tips: you can use regex to select examples. For example: example_00[1-5]* example_[6,7,8]*. If you want to run all examples, please type '*'.")                    


class PostdftExampleSourceSet(BaseModel):
    PostdftExampleSource_local: InputFilePath = Field(default=None,
                                         title="Upload postdft examples locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all required files in post dft. \
If you want to use the examples from datasets, please refer to the later 'PostdftExampleSource' section.""",
                                         description_type="markdown")

class PostdftExampleSet(BaseModel):                         
    PostdftExample: String = Field(default="*",title="Postdft examples",description = "You can choose to use only partial files of PostdftExampleSource, and separate each file or folder with space. \
Tips: you can use regex to select files. For example: example_00[1-5]* example_[6,7,8]*. If you want to use all files, please type '*'.")                    


class PrepareExtraFileSet(BaseModel):
    PrepareExtraFile_local: InputFilePath = Field(default=None,
                                                  title="Upload extra file for prepare locally",
                                                  st_kwargs_type=comm_func.unpack(
                                                      None, None, get_support_filetype=True),
                                                  description="""If some extra files that will be used for each newly created examples, you can upload the compressed file here.""")
        
class PrepareExtraFileNeededSet(BaseModel):        
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

class PredftExtraFileNeededSet(BaseModel):        
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
        
class RundftExtraFileNeededSet(BaseModel):        
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

class PostdftExtraFileNeededSet(BaseModel):        
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
    elif isinstance(example_source, FromDatasets):
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

def download_source(opts,
                    example_source_local_name,
                    example_name,
                    work_path,
                    download_path,
                    dataset_work_path,
                    logs,):
    '''
    example_source_name: choose where to get the example source
    example_source_local_name: the local path of the example source (in IO class to upload source from local)
    example_name: a string to define which files will be used.
    dataset_work_path: if use dataset as source
    
    Only support two types:
        1. Not use dataset, and only need example_source_local_name.
            If example_source_local_name is None, return None, else unpack the example_source_local_name to download_path,
            and then find all_files and all_directories in download_path. Then copy allfiles to work_path.
        2. Use dataset, and only need example_name.
            now will check dataset_work_path, if it is not None, will copy the example_name to work_path
    '''
    all_files = all_directories = None
    
    logs(f"read {example_name} setting ...")
    if example_source_local_name and hasattr(opts,example_source_local_name) and getattr(opts,example_source_local_name) != None:
        local_file_path = getattr(opts,example_source_local_name).get_path()
        logs(f"\t{example_source_local_name}:",local_file_path)
        comm_func.unpack(local_file_path, download_path)
        all_directories, all_files = copy_download_to_work(
            download_path, work_path,"*")
        comm_func.clean_dictorys(download_path)
        return all_directories,all_files
    elif dataset_work_path and example_name != None and hasattr(opts,example_name):
        need_files = getattr(opts,example_name)
        if need_files != None and need_files.strip() != "":
            logs(f"\t{example_name}:",need_files)
            all_directories, all_files = copy_download_to_work(
                dataset_work_path, work_path,need_files.strip())
            return all_directories,all_files
    return all_directories,all_files

def get_dataset_work_path(opts):
    if hasattr(opts,"dataset"):
        if opts.dataset.get_full_path():
            dataset_work_path = opts.dataset.get_full_path()
            if os.path.isdir(dataset_work_path):
                return dataset_work_path
            elif os.path.isfile(dataset_work_path):
                return os.path.dirname(dataset_work_path)
    return None
    
def read_source(opts,work_path,download_path,logs=None):
    # read setting in opts, and dowload example/rundft_extrafile/postdft_extrafile to downlaod path
    # and copy selected files to work path
    # return a dictory, which contains the list of example/rundft_extrafile/postdft_extrafile
    outdict = {}
    
    if logs == None:
        logs = print
    
    dataset_work_path = get_dataset_work_path(opts)
    logs("dataset_work_path:",dataset_work_path)
    logs("\t",os.listdir(dataset_work_path))
    
    #read example source
    all_directories, all_files = download_source(opts,
                                                 "ExampleSource_local",
                                                 "Example",
                                                 work_path,
                                                 download_path,
                                                 dataset_work_path,
                                                 logs)
    outdict["all_files"] = all_files
    
    #read prepare example source
    all_directories, all_files = download_source(opts,
                                                 "PrepareExampleSource_local",
                                                 "PrepareExample",
                                                 work_path,
                                                 download_path,
                                                 dataset_work_path,
                                                 logs)    
    outdict["prepare_example"] = all_files
            
    #read predft example source
    all_directories, all_files = download_source(opts,
                                                 "PredftExampleSource_local",
                                                 "PredftExample",
                                                 work_path,
                                                 download_path,
                                                 dataset_work_path,
                                                 logs) 
    outdict["predft_example"] = all_directories
           
    #read rundft example source
    all_directories, all_files = download_source(opts,
                                                 "RundftExampleSource_local",
                                                 "RundftExample",
                                                 work_path,
                                                 download_path,
                                                 dataset_work_path,
                                                 logs) 
    outdict["rundft_example"] = all_directories
   
    #read postdft example source
    all_directories, all_files = download_source(opts,
                                                 "PostdftExampleSource_local",
                                                 "PostdftExample",
                                                 work_path,
                                                 download_path,
                                                 dataset_work_path,
                                                 logs) 
    outdict["postdft_example"] = all_directories

    #read prepare extra files
    all_directories, all_files = download_source(opts,
                                                 "PrepareExtraFile_local",
                                                 "PrepareExtraFile_needed_files",
                                                 work_path,
                                                 download_path,
                                                 dataset_work_path,
                                                 logs) 
    outdict["prepare_extrafile"] = all_files
        
    #read predft extra files
    all_directories, all_files = download_source(opts,
                                                 "PredftExtraFile_local",
                                                 "PredftExtraFile_needed_files",
                                                 work_path,
                                                 download_path,
                                                 dataset_work_path,
                                                 logs) 
    outdict["predft_extrafile"] = all_files
        
    #read rundft extra files
    all_directories, all_files = download_source(opts,
                                                 "RundftExtraFile_local",
                                                 "RundftExtraFile_needed_files",
                                                 work_path,
                                                 download_path,
                                                 dataset_work_path,
                                                 logs) 
    outdict["rundft_extrafile"] = all_files
    
    #read postdft extra files
    all_directories, all_files = download_source(opts,
                                                 "PostdftExtraFile_local",
                                                 "PostdftExtraFile_needed_files",
                                                 work_path,
                                                 download_path,
                                                 dataset_work_path,
                                                 logs) 
    outdict["postdft_extrafile"] = all_files
    return outdict
