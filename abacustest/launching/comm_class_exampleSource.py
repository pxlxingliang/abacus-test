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

# DatasetSet/ExampleSourceSet/ExampleSet define the files of examples
class DatasetSet(BaseModel):
    dataset: DataSet = Field(title=None,
                             default=None,
                            description="Please enter your dataset in launching.")

class ExampleSourceSet(BaseModel):
    ExampleSource_local: List[InputFilePath] = Field(
                                        #default=None,
                                         title="Examples and settings",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""Upload one or some compressed files contianing all examples and the setting file.""",
                                         description_type="markdown")
    
class ExampleSet(BaseModel):
    Example: String = Field(default="*",title="Examples",description = "You can choose to use only partial files, and seperate each example with space. \
Tips: you can use regex to select files. For example: example_00[1-5]* example_[6,7,8]*. If you want to use all files, please type '*'.")  

# script dataset is similar to example dataset, but it is used for scripts, and there has an attribute to define the unneeded files
class ScriptDatasetSet(BaseModel):
    script_dataset: DataSet = Field(title=None,
                                default=None,
                                description="Please enter the script dataset. Like: \"launching+datasets://reuse.abacustest@v1.0/stress_difference_test\". You can choose use the scripts from datasets or upload from local at below.")

class ScriptSourceSet(BaseModel):
    ScriptSource_local: InputFilePath = Field(default=None,
                                         title="Upload the script files locally",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""A compressed file contains all required files. Please make sure there has a setting.json file in the root directory, which controls the running of workflow""",
                                         description_type="markdown")

class ScriptExampleSet(BaseModel):
    Script_Unneeded: String = Field(default="",
                            title="Unneeded Script",
                            description = "If not all the files in the script dataset are needed, please enter the unneeded files or directories. Like: \"*.py,*.sh\". If all the files are needed, please leave it blank.")  

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
Tips: you can use regex to select examples. For example: example_00[1-5]* example_[6,7,8]*. If you want to run all examples, please enter '*'.")     

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

def copy_download_to_work(download_path,work_path,needed_files,reverse=False):
    #needed_files is a string, each file is separated by space
    # if reverse is True, will copy all files except needed_files
    if reverse:
        cwd = os.getcwd()
        os.chdir(download_path)
        unneeded_files = [] 
        for ifile in needed_files.split():
            for iifile in glob.glob(os.path.join(download_path,ifile)):
                unneeded_files.append(iifile)
        os.chdir(cwd)
        if unneeded_files == []:
            return copy_download_to_work(download_path,work_path,"*",reverse=False)
        else:
            # we firstly copy all files to a tmp directory, and then remove the unneeded files
            # then move the tmp directory to work_path
            tmp_path = "tmp"
            n = 0
            while os.path.exists(os.path.join(download_path,tmp_path)) and os.path.exists(os.path.join(work_path,tmp_path)):
                tmp_path = f"tmp{n}"
                n += 1
            tmp_path = os.path.join(work_path,tmp_path)
            copy_download_to_work(download_path,tmp_path,"*",reverse=False)
            # remove unneeded files
            os.chdir(tmp_path)
            for iifile in unneeded_files:
                if os.path.isdir(iifile):
                    shutil.rmtree(iifile)
                else:
                    os.remove(iifile)
            os.chdir(cwd)
            alldirectories,allfiles = copy_download_to_work(tmp_path,work_path,"*",reverse=False)
            shutil.rmtree(tmp_path)
            return alldirectories,allfiles
          
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
                    logs,
                    example_name_is_unneeded=False):
    '''
    example_source_name: choose where to get the example source
    example_source_local_name: the local path of the example source (in IO class to upload source from local)
    example_name: a string to define which files will be used or not be used (example_name_is_unneeded=True).
    dataset_work_path: if use dataset as source
    
    Support both two types:
        1. Not use dataset, and only need example_source_local_name.
            If example_source_local_name is None, return None, else unpack the example_source_local_name to download_path,
            and then find all_files and all_directories in download_path. Then copy allfiles to work_path.
        2. Use dataset, and only need example_name.
            now will check dataset_work_path, if it is not None, will copy the example_name to work_path
    Will firstly copy files from dataset, and the copy files from example_source_local_name
    '''
    all_files = all_directories = None
    
    need_files = None
    if example_name != None and hasattr(opts,example_name):
        need_files = getattr(opts,example_name)
    if example_name_is_unneeded and need_files == None:
        need_files = ""
        # be noticed, if example_name_is_unneeded is True, then need_files actually is unneeded_files, 
        # and if need_files is None, it means all files are needed, and so we set need_files to ""
    
    logs(f"read {example_name} setting ...")
    if dataset_work_path:
        if not example_name_is_unneeded:
            if need_files != None and need_files.strip() != "":
                logs(f"\t{example_name}:",need_files)
                all_directories, all_files = copy_download_to_work(
                    dataset_work_path, work_path,need_files.strip())
        else:
            logs(f"\t{example_name}:",need_files)
            all_directories, all_files = copy_download_to_work(
                dataset_work_path, work_path,need_files,reverse=True)
    
    if example_source_local_name and hasattr(opts,example_source_local_name) and getattr(opts,example_source_local_name) != None:
        if need_files == None: need_files = "*" 
        local_files = getattr(opts,example_source_local_name)
        
        if isinstance(local_files,list):
            local_file_path = [i.get_path() for i in local_files]
        else:
            local_file_path = [local_files.get_path()]
            
        logs(f"\t{example_source_local_name}:",local_file_path) 
        for i_local_file_path in local_file_path:
            comm_func.unpack(i_local_file_path, download_path)
        
        all_directories_tmp, all_files_tmp = copy_download_to_work(
                download_path, work_path,need_files,reverse=example_name_is_unneeded)
            
        comm_func.clean_dictorys(download_path)
        if all_directories_tmp:
            if all_directories == None:
                all_directories = all_directories_tmp
            else:
                all_directories.extend(list(set(all_directories_tmp)-set(all_directories)))
        if all_files_tmp:
            if all_files == None:
                all_files = all_files_tmp
            else:
                all_files.extend(list(set(all_files_tmp)-set(all_files)))
    
    return all_directories,all_files

def get_dataset_work_path(opts,dataset_name="dataset"):
    # dataset_name should be "dataset" or "script_dataset"
    if hasattr(opts,dataset_name):
        if getattr(opts,dataset_name):
            try:
                dataset_work_path = getattr(opts,dataset_name).get_full_path()
                if os.path.isdir(dataset_work_path):
                    return dataset_work_path
                elif os.path.isfile(dataset_work_path):
                    return os.path.dirname(dataset_work_path)
            except:
                return None
    # for reuse my_model, if dataset_name == "dataset" and has "my_model" in opts, then use "my_model" as dataset
    # if the origin dataset has enetered infoamtion, then it will be returned in previous step.
    # here, we only judge if my_model is not " "
    if dataset_name=="dataset" and hasattr(opts,"my_model") and hasattr(opts,"reuse_dataset"):
        my_model = getattr(opts,"my_model")
        print("my_model:",my_model)
        if my_model and my_model.strip() != "":
            print("reuse_dataset:",getattr(opts,"reuse_dataset"))
            dataset_work_path = getattr(opts,"reuse_dataset").get_full_path()
            dataset_work_path = os.path.join(dataset_work_path,my_model)
            if os.path.isdir(dataset_work_path):
                return dataset_work_path
            else:
                print(f"ERROR: The dataset ({dataset_work_path}) is not valid!")
                return None 
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
    
    #read script source
    # For original design, only use one dataset, but here we may need another script dataset
    # So, the dataset_work_path is not used here, but use the script_dataset
    script_dataset_work_path = get_dataset_work_path(opts,"script_dataset")
    all_directories, all_files = download_source(opts,
                                                 "ScriptSource_local",
                                                 "Script_Unneeded",
                                                 work_path,
                                                 download_path,
                                                 script_dataset_work_path,
                                                 logs,
                                                 example_name_is_unneeded=True)
    outdict["all_scripts"] = all_files
    return outdict
