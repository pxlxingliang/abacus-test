import traceback,json,sys
from dp.launching.typing.basic import BaseModel,String
from dp.launching.typing import Field
import os



from . import (comm_class,
               comm_func,
               comm_class_exampleSource,
               comm_report) 

'''
class Models(String, Enum):
    model00 = " "
    model0 = "000-abacus-branch"
    model1 = "001-abacus-dailytest"
    #model3 = 
    model4 = "004-lcao004_unstable_elpa"
    model4_1 = "004-lcao004_unstable_scalapack"
    model5 = "005-finite_diff_stress"
    model5_1 = "005-finite_diff_force"  
    model6 = "006_eos"
    model7 = "007_bandstru"
    model7_1 = "007_bandstru_abacusVSvasp"
    model8 = "008_dos"
    model10 = "010_LcaoVSPw"
    model11 = "011_abacusVSvasp_eos"
    model11_1 = "011_LcaoVSPw_eos"
    
class MyModel(BaseModel):
    my_model: Models = Field(
        title="Model",
        default="",description="The model you want to use. If you have entered the dataset information in the previous step, here will be ignored. DOCs: https://dptechnology.feishu.cn/docx/I7NKdmURHosJSMxn7NCc5RsEnxb")
'''
class NewSetting(BaseModel):
    predft_command: String = Field(default="",
                            title="New Predft Command",
                            description="If you don't want to use the command from the previously uploaded setting file, you can fill in the new command here",)
    rundft_image: String = Field(default="",
                          title="New Rundft Image",
                          description="If you don't want to use the image from the previously uploaded setting file, you can fill in the new command here.",)
    rundft_command: String = Field(default="",
                            title="New Rundft Command",
                            description="If you don't want to use the command from the previously uploaded setting file, you can fill in the new command here",)
    rundft_machine: String = Field(default="",
                            title="New Rundft Machine",
                            description="If you don't want to use the machine from the previously uploaded setting file, you can fill in the new command here.",)
    postdft_command: String = Field(default="",
                            title="New Postdft Command",
                            description="If you don't want to use the command from the previously uploaded setting file, you can fill in the new command here.",)

'''
group1 = ui.Group("if_use_reuse_dataset","test")

@group1
@ui.Visible(MyModel,("my_model"),Equal,("do not show this group"))
class ReuseDataset(BaseModel):
    reuse_dataset: DataSet = Field(title=None,
                             default="launching+datasets://reuse.abacustest@latest",
                            description="the reuse dataset")
'''                            
        
class ReuseModel(
    #comm_class.ConfigSetGithub,
    #comm_class.TrackingSet,
    NewSetting,
    #comm_class_exampleSource.ScriptExampleSet,
    #comm_class_exampleSource.ScriptSourceSet,
    #comm_class_exampleSource.ScriptDatasetSet,
    #comm_class_exampleSource.ExampleSet,
    comm_class_exampleSource.ExampleSourceSet,
    #comm_class_exampleSource.DatasetSet,
    #MyModel,
    #ReuseDataset,
                    comm_class.OutputSet,
                    comm_class.ConfigSet,
                    BaseModel):
    ...  

def ReuseModelRunner(opts:ReuseModel) -> int:
    try:
        logs = comm_class.myLog()

        paths = comm_func.create_path(str(opts.IO_output_path))
        output_path = paths["output_path"]
        work_path = paths["work_path"]
        download_path = paths["download_path"]

        logs.iprint("read source setting ...")
        datas = comm_class_exampleSource.read_source(opts,work_path,download_path,logs.iprint)
        
        setting_file = os.path.join(work_path,"setting.json")
        if not os.path.exists(setting_file):
            logs.iprint("\nERROR: Please supply the setting information!!!!")
            return 1
        try:
            setting = json.load(open(setting_file))
        except:
            traceback.print_exc()
            return 1
        
        #read setting
        allparams = {"config": comm_func.read_config(opts)}
        for k,v in setting.items():
            if k == "config":
                for ik,iv in v.items():
                    allparams["config"][ik] = iv
            #elif k in ["ABBREVIATION","save_path","run_dft","post_dft","report","dataset_info","upload_datahub","upload_tracking"]:
            else:
                allparams[k] = v
                
        # modify the command and image
        if opts.predft_command != None and opts.predft_command.strip() != "":
            new_predft_command = opts.predft_command.strip()
        else:
            new_predft_command = None
        if opts.rundft_image != None and opts.rundft_image.strip() != "":
            new_rundft_image = opts.rundft_image.strip()
        else:
            new_rundft_image = None
        if opts.rundft_command != None and opts.rundft_command.strip() != "":
            new_rundft_command = opts.rundft_command.strip()
        else:
            new_rundft_command = None
        if opts.rundft_machine != None and opts.rundft_machine.strip() != "":
            new_rundft_machine = opts.rundft_machine.strip()
        else:
            new_rundft_machine = None
        if opts.postdft_command != None and opts.postdft_command.strip() != "":
            new_postdft_command = opts.postdft_command.strip()
        else:
            new_postdft_command = None
        
        if "pre_dft" in allparams:
            if new_predft_command:
                allparams["pre_dft"]["command"] = new_predft_command
        
        if "run_dft" in allparams:
            if isinstance(allparams["run_dft"],list):
                for idx in range(len(allparams["run_dft"])):
                    if new_rundft_image:
                        allparams["run_dft"][idx]["image"] = new_rundft_image
                    if new_rundft_command:
                        allparams["run_dft"][idx]["command"] = new_rundft_command
                    if new_rundft_machine and "bohrium" in allparams["run_dft"][idx] and isinstance(allparams["run_dft"][idx]["bohrium"],dict):
                        allparams["run_dft"][idx]["bohrium"]["scass_type"] = new_rundft_machine
            elif isinstance(allparams["run_dft"],dict):
                if new_rundft_image:
                    allparams["run_dft"]["image"] = new_rundft_image
                if new_rundft_command:
                    allparams["run_dft"]["command"] = new_rundft_command
                if new_rundft_machine and "bohrium" in allparams["run_dft"] and isinstance(allparams["run_dft"]["bohrium"],dict):
                    allparams["run_dft"]["bohrium"]["scass_type"] = new_rundft_machine
        
        if "post_dft" in allparams:
            if new_postdft_command:
                allparams["post_dft"]["command"] = new_postdft_command
     
        #execut
        stdout,stderr = comm_func.exec_abacustest(allparams,work_path)
        logs.iprint(f"{stdout}\n{stderr}\nrun abacustest over!\n")
        
        comm_report.gen_report(opts,logs,work_path,output_path,allparams)
        
    except:
        traceback.print_exc()
        return 1

    return 0
