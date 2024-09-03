import traceback
from dp.launching.typing.basic import BaseModel, Int, String, Float, List, Optional, Union, Dict
from dp.launching.typing import (
    BaseModel,
    Set,
    Boolean,
    Field,
    InputFilePath
)
from enum import Enum
from typing import Literal
import re,os,shutil
from . import comm_class_exampleSource,comm_func


class PPLibEnum(String, Enum):
    dataset1 = "SG15-ONCV-v1_0"

class OrbLibEnum(String, Enum):
    dataset1 = "SG15-Version1_0__StandardOrbitals-Version2_0"

class PPFromDatasets(BaseModel):
    type: Literal["from datasets"]
    dataset: PPLibEnum = Field(title="PP (pseudo potential) Lib datasets",
                        description="Please choose the PP Lib datasets.")
    '''
    dataset_unrecorded: String = Field(default=None,
                                       description="By default, the pp lib \"[SG15-ONCV-v1_0](https://launching.mlops.dp.tech/?request=GET%3A%2Fapplications%2Fabacustest%2Fdatasets%2Fpporb)\" will be used. \
If you want to use other pplib dataset, please enter download link here.")
    '''

class OrbFromDatasets(BaseModel):
    type: Literal["from datasets"]
    dataset: OrbLibEnum = Field(title="Orb Lib datasets",
                        description="Please choose the Orb Lib datasets.")
    '''
    dataset_unrecorded: String = Field(default=None,
                                       description="By default, the Orb lib \"[SG15-Version1_0__StandardOrbitals-Version2_0](https://launching.mlops.dp.tech/?request=GET%3A%2Fapplications%2Fabacustest%2Fdatasets%2Fpporb)\" will be used. \
If you want to use other pplib dataset, please enter download link here.")
    '''

#prepare INPUT template
class PrepareInputTemplateSet(BaseModel):
    PrepareInputTemplate_local: InputFilePath = Field(default=None,
                                                  title="Upload INPUT template locally",
                                                  description="""You can upload an INPUT file, and then all examples will use this INPUT""")

# prepare INPUT template path, for dataset case
class PrepareInputTemplatePathSet(BaseModel):
    PrepareInputTemplatePath: String = Field(default=None,title="INPUT template",description = "If you want to use an INPUT template for all examples, please enter the path of INPUT template here.")             

#prepare STRU template
class PrepareStruTemplateSet(BaseModel):
    PrepareStruTemplate_local: InputFilePath = Field(default=None,
                                                  title="Upload STRU template locally",
                                                  description="""You can upload a STRU file, and then all examples will use this STRU""")

#prepare STRU template path, for dataset case
class PrepareStruTemplatePathSet(BaseModel):
    PrepareStruTemplatePath: String = Field(default=None,title="STRU template",description = "If you want to use an STRU template for all examples, please enter the path of STRU template here.")
            
#prepare KPT template
class PrepareKptTemplateSet(BaseModel):
    PrepareKptTemplate_local: InputFilePath = Field(default=None,
                                                  title="Upload KPT template locally",
                                                  description="""You can upload a KPT file, and then all examples will use this KPT""")

#prepare KPT template path, for dataset case
class PrepareKptTemplatePathSet(BaseModel):
    PrepareKptTemplatePath: String = Field(default=None,title="KPT template",description = "If you want to use an KPT template for all examples, please enter the path of KPT template here.")

#prepare dpks descriptor file
class PrepareDPKSDescriptorSet(BaseModel):
    PrepareDPKSDescriptor_local: InputFilePath = Field(default=None,
                                                  title="Upload DeeP-KS Descriptor locally",
                                                  description="""If you want to use Deep-KS, you can upload the descriptor file here, or prepare the descriptor file in each example directory.""")

# prepare dpks descriptor path, for dataset case
class PrepareDPKSDescriptorPathSet(BaseModel):
    PrepareDPKSDescriptorPath: String = Field(default=None,title="DeeP-KS Descriptor",description = "If you want to use an DeeP-KS Descriptor for all examples, please enter the path of DeeP-KS Descriptor here.")
    
#prepare pp lib
class PreparePPLibSet(BaseModel):
    PreparePPLib_local: InputFilePath = Field(default=None,
                                                  title="Upload pseudopotential library locally",
                                                  description="""If you have a pseudopotential library and you want to use it for all elements, please upload it here. Please also prepare an \"element.json\" file in the library directory, and the key is name of element and value is the name of pseudopotential file.""")
    '''
    PreparePPLib: Union[
        comm_class_exampleSource.NotRquired,
        comm_class_exampleSource.FromPreUpload,
    #    comm_class_exampleSource.FromDatahub,
        PPFromDatasets] = Field(title="Prepare pseudopotential lib source",
                              discriminator="type",
                              description="Please choose the PP Lib source.")   
    '''
# prepare pp lib path, for dataset case
class PreparePPLibPathSet(BaseModel):
    PreparePPLibPath: String = Field(default=None,title="PP Lib",description = "If you want to use pseudopotential library, please enter the path of PP Lib here.")
    
#prepare ORB lib
class PrepareOrbLibSet(BaseModel):
    PrepareOrbLib_local: InputFilePath = Field(default=None,
                                                  title="Upload Orbital library locally",
                                                  description="""If you have an Orbital library, you can upload it here. Please also prepare a \"element.json\" file in the library directory, and the key is name of element and value is the name of Orbital file.""")
    '''
    PrepareOrbLib: Union[
        comm_class_exampleSource.NotRquired,
        comm_class_exampleSource.FromPreUpload,
    #    comm_class_exampleSource.FromDatahub,
        OrbFromDatasets] = Field(title="Prepare Orb Lib source",
                              discriminator="type",
                              description="Please choose the Orb Lib source.")    
    '''
# prepare orb lib path, for dataset case
class PrepareOrbLibPathSet(BaseModel):
    PrepareOrbLibPath: String = Field(default=None,title="Orb Lib",description = "If you want to use Orbital library, please enter the path of Orb Lib here.")

class PrepareSet(BaseModel):
    
    prepare_mix_kpt: String = Field(default=None,title="Additional KPT settings",description="You can set additional KPT SETTING for each examples. \
Please enter 1/3/6 values seperated by space (such as \"2\" means [2,2,2,0,0,0], which is 3 K ponits and the shift of mesh in K space; \"1 2 3\" means [1,2,3,0,0,0]; \"1 2 3 1 0 0\" means [1,2,3,1,0,0]).\
You can set multiple values (separated by comma), which will generate a set of ABACUS inputs for each value. (Scuch as: \"3, 1 2 3\")")

    prepare_mix_input: Optional[Dict[String,String]] = Field(default=None,title="Additional INPUT settings",description="You can set additional INPUT parameters for each examples. \
You can set multiple values (separated by comma) for each parameter, which will generate a set of ABACUS inputs for each value. If you want to combine several parameters, you can connect the parameters with |, and also conbine values with |. Such as: key is ecutwfc|kspacing, and value is 50|0.1, 60|0.2, 70|0.3. Commonly used parameters: \
\"calculation\", \"ecutwfc\", \"scf_thr\", \"scf_nmax\", \"basis_type\", \"smearing_method\", \"smearing_sigma\", \"mixing_type\", \"mixing_beta\",\"ks_solver\"")
    
def parse_prepare(prepare_set,work_path,download_path,logs):
    prepare = {}
    
    dataset_work_path = comm_class_exampleSource.get_dataset_work_path(prepare_set)
            
    # parse INPUT template
    if hasattr(prepare_set,"PrepareInputTemplate_local") and getattr(prepare_set,"PrepareInputTemplate_local") \
        and os.path.isfile(prepare_set.PrepareInputTemplate_local.get_path()):
            prepare["input_template"] = prepare_set.PrepareInputTemplate_local.get_path()
    elif dataset_work_path and hasattr(prepare_set,"PrepareInputTemplatePath") and getattr(prepare_set,"PrepareInputTemplatePath") and \
        os.path.isfile(os.path.join(dataset_work_path,prepare_set.PrepareInputTemplatePath)):
            prepare["input_template"] = os.path.join(dataset_work_path,prepare_set.PrepareInputTemplatePath)  
    
    # parse STRU template
    if hasattr(prepare_set,"PrepareStruTemplate_local") and getattr(prepare_set,"PrepareStruTemplate_local") and \
        os.path.isfile(prepare_set.PrepareStruTemplate_local.get_path()):
            prepare["stru_template"] = prepare_set.PrepareStruTemplate_local.get_path()
    elif dataset_work_path and hasattr(prepare_set,"PrepareStruTemplatePath") and getattr(prepare_set,"PrepareStruTemplatePath") and \
        os.path.isfile(os.path.join(dataset_work_path,prepare_set.PrepareStruTemplatePath)):
            prepare["stru_template"] = os.path.join(dataset_work_path,prepare_set.PrepareStruTemplatePath)
    
    # parse KPT template
    if hasattr(prepare_set,"PrepareKptTemplate_local") and getattr(prepare_set,"PrepareKptTemplate_local") and \
        os.path.isfile(prepare_set.PrepareKptTemplate_local.get_path()):
            prepare["kpt_template"] = prepare_set.PrepareKptTemplate_local.get_path()
    elif dataset_work_path and hasattr(prepare_set,"PrepareKptTemplatePath") and getattr(prepare_set,"PrepareKptTemplatePath") and \
        os.path.isfile(os.path.join(dataset_work_path,prepare_set.PrepareKptTemplatePath)):
            prepare["kpt_template"] = os.path.join(dataset_work_path,prepare_set.PrepareKptTemplatePath)
    
    # parse DPKS descriptor
    if hasattr(prepare_set,"PrepareDPKSDescriptor_local") and getattr(prepare_set,"PrepareDPKSDescriptor_local") and \
        os.path.isfile(prepare_set.PrepareDPKSDescriptor_local.get_path()):
            prepare["dpks_descriptor"] = prepare_set.PrepareDPKSDescriptor_local.get_path()
    elif dataset_work_path and hasattr(prepare_set,"PrepareDPKSDescriptorPath") and getattr(prepare_set,"PrepareDPKSDescriptorPath") and \
        os.path.isfile(os.path.join(dataset_work_path,prepare_set.PrepareDPKSDescriptorPath)):
            prepare["dpks_descriptor"] = os.path.join(dataset_work_path,prepare_set.PrepareDPKSDescriptorPath)
            
    # parse mix input
    if hasattr(prepare_set,"prepare_mix_input") and prepare_set.prepare_mix_input:
        prepare["mix_input"] = {}
        for key,value in prepare_set.prepare_mix_input.items():
            ivalue = value.split(",")
            prepare["mix_input"][key] = ivalue

    if hasattr(prepare_set,"prepare_mix_kpt") and prepare_set.prepare_mix_kpt:
        kpt_tmp = []
        for ivalue in prepare_set.prepare_mix_kpt.split(","): 
            iivalue = ivalue.split()
            if len(iivalue) == 1:
                try:
                    kpt_tmp.append(int(iivalue[0]))
                except:
                    traceback.print_exc()
                    logs(f"ERROR: The {ivalue} is not valid for KPT!")
            elif len(iivalue) == 3:
                try:
                    kpt_tmp.append([int(iivalue[0]),int(iivalue[1]),int(iivalue[2]),0,0,0])
                except:
                    traceback.print_exc()
                    logs(f"ERROR: The {ivalue} is not valid for KPT!")
            elif len(iivalue) == 6:
                try:
                    kpt_tmp.append([int(iivalue[0]),int(iivalue[1]),int(iivalue[2]),float(iivalue[3]),float(iivalue[4]),float(iivalue[5])])
                except:
                    traceback.print_exc()
                    logs(f"ERROR: The {ivalue} is not valid for KPT!")
            else:
                logs(f"ERROR: The {ivalue} is not valid for KPT!")
                        
        if kpt_tmp:    
            prepare["mix_kpt"] = kpt_tmp
        else:
            logs(f"ERROR: The {prepare_set.prepare_mix_kpt} is not valid for KPT!")
            
    # parse pp lib
    if hasattr(prepare_set,"PreparePPLib"):
        if isinstance(prepare_set.PreparePPLib,PPFromDatasets):
            tmp = None
            try:
                #if prepare_set.PreparePPLib.dataset_unrecorded != None and prepare_set.PreparePPLib.dataset_unrecorded.strip() != "":
                #    url = prepare_set.PreparePPLib.dataset_unrecorded.strip()
                #else:
                url = comm_class_exampleSource.GetDatasetAddress(prepare_set.PreparePPLib.dataset,dataset="pporb")
                package = comm_func.download_url(url, download_path)
                if package == None:
                    logs(f"ERROR: download dataset failed!\n\turl:{url}")
                    logs(f"\tPlease check the dataset!")
                else:
                    comm_func.unpack(package, download_path,filetype="tgz")
                    #remove the package
                    os.remove(package)
                    tmp = True
                    #move the package to work_path
            except:
                traceback.print_exc()
        else:
            tmp = comm_class_exampleSource.parse_source(prepare_set.PreparePPLib,prepare_set.PreparePPLib_local,download_path,prepare_set,logs)
        if tmp:
            pp_path = os.path.join(work_path,"pplib")
            if not os.path.isdir(pp_path):
                os.makedirs(pp_path)
            for ifile in os.listdir(download_path):
                shutil.move(os.path.join(download_path,ifile),os.path.join(pp_path,ifile))
            prepare["pp_path"] = "pplib"
        comm_func.clean_dictorys(download_path)
    elif dataset_work_path and hasattr(prepare_set,"PreparePPLibPath") and getattr(prepare_set,"PreparePPLibPath") and \
        os.path.isdir(os.path.join(dataset_work_path,prepare_set.PreparePPLibPath)):
            prepare["pp_path"] = os.path.join(dataset_work_path,prepare_set.PreparePPLibPath)
    
    #prepare orb lib
    if hasattr(prepare_set,"PreparePPLib"):
        if isinstance(prepare_set.PrepareOrbLib,OrbFromDatasets):
            tmp = None
            try:
                #if prepare_set.PrepareOrbLib.dataset_unrecorded != None and prepare_set.PrepareOrbLib.dataset_unrecorded.strip() != "":
                #    url = prepare_set.PrepareOrbLib.dataset_unrecorded.strip()
                #else:
                url = comm_class_exampleSource.GetDatasetAddress(prepare_set.PrepareOrbLib.dataset,dataset="pporb")
                package = comm_func.download_url(url, download_path)
                if package == None:
                    logs(f"ERROR: download dataset failed!\n\turl:{url}")
                    logs(f"\tPlease check the dataset!")
                else:
                    comm_func.unpack(package, download_path,filetype="tgz")
                    #remove the package
                    os.remove(package)
                    tmp = True
            except:
                traceback.print_exc()
        else:
            tmp = comm_class_exampleSource.parse_source(prepare_set.PrepareOrbLib,prepare_set.PrepareOrbLib_local,download_path,prepare_set,logs)
        if tmp:
            #move the package to work_path
            orb_path = os.path.join(work_path,"orblib")
            if not os.path.isdir(orb_path):
                os.makedirs(orb_path)
            for ifile in os.listdir(download_path):
                shutil.move(os.path.join(download_path,ifile),os.path.join(orb_path,ifile))
            prepare["orb_path"] = "orblib"
        comm_func.clean_dictorys(download_path)
    elif dataset_work_path and hasattr(prepare_set,"PrepareOrbLibPath") and getattr(prepare_set,"PrepareOrbLibPath") and \
        os.path.isdir(os.path.join(dataset_work_path,prepare_set.PrepareOrbLibPath)):
            prepare["orb_path"] = os.path.join(dataset_work_path,prepare_set.PrepareOrbLibPath)
                                                          
    return prepare    
    
def construct_input(datas,opts,work_path,download_path,logs):
    # datas is a dict of examples, created by comm_class_exampleSource.read_source
    #read predft
    need_prepare = False
    logs.iprint("read prepare setting ...")
    prepare = {}

    other_sets = parse_prepare(opts,work_path,download_path,logs)
    if other_sets:
        prepare = other_sets
        need_prepare = True
    
    if datas.get("prepare_example"):
        prepare["example_template"] = datas.get("prepare_example")    
        logs.iprint("\texample:",prepare["example_template"])
        need_prepare = True
    
    if datas.get("prepare_extrafile"):
        prepare["extra_files"] = datas.get("prepare_extrafile")
        logs.iprint("\tprepare_extrafile:",prepare["extra_files"])
        need_prepare = True
    
    return need_prepare, prepare
    
