import traceback,json,sys,glob
from dp.launching.typing.basic import BaseModel,String
from dp.launching.typing import Field
import os

from abacustest.lib_model.model_005_Phonon import PreparePhono,PostprocessPhonon
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

from . import (comm_class,
               comm_func,
               comm_report,
               comm_class_exampleSource,
               readsetting) 


class Images(BaseModel):
    abacus_image: String = Field(default=RECOMMAND_IMAGE,
                          title="abacus Image",
                          description="If you do not want to use the default image, please enter the new image here.",)
    abacus_command: String = Field(default=RECOMMAND_COMMAND,
                            title="abacus Command",
                            description="The command to run abacus",)
    abacus_machine: String = Field(default=RECOMMAND_MACHINE,
                                   title="abacus Machine",
                                    description="The machine to run abacus",)
    phonopy_setting: String = Field(default="setting.conf",
                                   title="Phonopy Setting",
                                    description="the phonopy setting file name, and each example should have this file and containing DIM/BAND",)

class PhononModel(
    Images,
    comm_class_exampleSource.ExampleSourceSet,
                    comm_class.OutputSet,
                    comm_class.ConfigSet,
                    BaseModel):
    ...      
    
def PhononModelRunner(opts:PhononModel) -> int:   
    try:
        logs = comm_class.myLog()
        paths = comm_func.create_path(str(opts.IO_output_path))
        output_path = paths["output_path"]
        work_path = paths["work_path"]
        download_path = paths["download_path"]
        cwd = os.getcwd()
        logs.iprint("read source setting ...")
        datas = comm_class_exampleSource.read_source(opts,work_path,download_path,logs.iprint)
        allparams = {"config": comm_func.read_config(opts),
                     "save_path": ".",
                     "bohrium_group_name": "Phonopy"}
        
        examples = datas.get("all_files",[])
        if len(examples) == 0:
            return 1
        os.chdir(work_path)
        setting = PreparePhono(examples,opts.abacus_command,opts.abacus_image,opts.abacus_machine,"phonopy",opts.phonopy_setting)
        if len(setting) == 0:
            return 1
        
        # merge setting and config
        allparams.update(setting)

        #exec abacustest
        stdout,stderr = comm_func.exec_abacustest(allparams,work_path)
        
        #postprocess
        os.chdir(work_path)
        PostprocessPhonon(examples,"phonopy",opts.phonopy_setting)
        os.chdir(cwd)

        logs.iprint(f"{stdout}\n{stderr}\nrun abacustest over!\n")
        
        # produce metrics and superMetrics reports
        comm_report.gen_report(opts,logs,work_path,output_path,allparams)
    except:
        traceback.print_exc()
        return 1

    return 0