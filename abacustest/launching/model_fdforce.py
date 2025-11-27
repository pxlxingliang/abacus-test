import traceback,json,sys,os
from dp.launching.typing.basic import BaseModel,String,Float,Int
from dp.launching.typing import Field
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

from . import (comm_class,
               comm_func,
               comm_report,
               comm_class_exampleSource) 
from abacustest.lib_model import model_006_FDForce as FDForce


class NewSetting(BaseModel):
    fd_step: Float = Field(default=0.01,
                            title="FD Step (Angstrom)",
                            description="The finite difference step size. Unit in angstrom",)
    
    fd_number: Int = Field(default=5,
                           titile= "FD Number",
                            description="The number of finite difference steps. Will calculate extra 2*Fd_Number structures.",)
    
    abacus_image: String = Field(default=RECOMMAND_IMAGE,
                          title="Abacus Image",
                          description="The image to run abaucs.",)
    abacus_command: String = Field(default=RECOMMAND_COMMAND,
                            title="Abacus Command",
                            description="The command to execute abacus",)
    bohrium_machine: String = Field(default=RECOMMAND_MACHINE,
                            title="Bohrium Machine",
                            description="The bohrium machine type to run abacus",)


class FDForceModel(
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

def FDForceModelRunner(opts:FDForceModel) -> int:
    try:
        logs = comm_class.myLog()

        paths = comm_func.create_path(str(opts.IO_output_path))
        output_path = paths["output_path"]
        work_path = paths["work_path"]
        download_path = paths["download_path"]

        logs.iprint("read source setting ...")
        datas = comm_class_exampleSource.read_source(opts,work_path,download_path,logs.iprint)

        pwd = os.getcwd()
        # prepare the fdforce
        os.chdir(work_path)
        allfiles = datas["all_files"]
        allexamples = [i for i in allfiles if os.path.isdir(i) and os.path.isfile(os.path.join(i,"info.txt"))]
        if len(allexamples) == 0:
            logs.iprint("no example found, exit! Please prepare the inputs of each example in each folder with a info.txt file.")
            return 1
        subfolders = FDForce.PrepareFDForce("info.txt",allexamples,opts.fd_step,opts.fd_number).run()
        if len(subfolders) == 0:
            logs.iprint("no example found, exit! Please prepare the inputs of each example in each folder with a info.txt file.")
            return 1

        allparams = {
            "config": comm_func.read_config(opts),
            "save_path": ".",
            "bohrium_group_name": "FDForce",
            "run_dft": {
                "example": subfolders,
                "command": opts.abacus_command,
                "image": opts.abacus_image,
                "bohrium": {
                    "scass_type": opts.bohrium_machine,
                    "job_type": "container",
                    "platform": "ali",
                },
            },
        }
        os.chdir(pwd)

        # execut
        stdout,stderr = comm_func.exec_abacustest(allparams,work_path)
        os.chdir(work_path)
        FDForce.PostProcessFDForce(allexamples).run()
        os.chdir(pwd)

        
        
        comm_report.gen_report(opts,logs,work_path,output_path,allparams)
    except:
        traceback.print_exc()
        return 1

    return 0
