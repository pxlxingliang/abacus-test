import traceback,json,sys,os
from dp.launching.typing.basic import BaseModel,String,Float,Int,Boolean,Set, List, Optional
from dp.launching.typing import Field, InputFilePath
from enum import Enum
from abacustest.lib_model.model_014_vasp2abacus import Vasp2Abacus


from . import (comm_class,
               comm_func,
               comm_report,
               ) 

class VASPSourceSet(BaseModel):
    VaspSource: List[InputFilePath] = Field(
                                         title="VASP任务目录压缩包（支持多个任务）[ENCUT/EDIFF等参数不会进行转换]",
                                         st_kwargs_type=comm_func.unpack(
                                             None, None, get_support_filetype=True),
                                         description="""上传VASP任务目录的压缩包，可对INCAR/POSCAR/KPOINTS转换成ABACUS相应的输入文件""",
                                         format="markdown")
    
    InputTemplate: Optional[String] = Field( title="ABACUS额外输入参数设置",
                                default="basis_type lcao\necutwfc 100",
                                description="""ABACUS的额外输入参数设置，如果不需要额外设置，请留空""",
                                format="multi-line")

def unpack_files(vaspsource, work_path):
    """
    Link the files from the VASP source set to the work path.
    """
    for source in vaspsource:
        source = source.get_full_path()
        if os.path.isfile(source):
            comm_func.unpack(source, work_path)
    
    dirs = [i for i in os.listdir(work_path) if os.path.isdir(os.path.join(work_path, i))]
    if len(dirs) == 0:
        dirs = ["."]
    
    return dirs

def write_input_template(input_template, work_path):
    """
    Write the input template to a file in the work path.
    """
    if input_template.strip():
        input_file = os.path.join(work_path, "INPUT.template")
        with open(input_file, "w") as f:
            f.write(input_template)
        return input_file
    return None        

class Vasp2AbacusModel(
    VASPSourceSet,
    comm_class.OutputSet,
    BaseModel):
    ...  

def Vasp2AbacusRunner(opts:Vasp2AbacusModel) -> int:
    try:
        logs = comm_class.myLog()

        paths = comm_func.create_path(str(opts.IO_output_path))
        work_path = paths["work_path"]
        output_path = paths["output_path"]
        
        logs.iprint("read source setting ...")
        
        vaspjobs = unpack_files(opts.VaspSource, work_path) # only folder names
        input_template = write_input_template(opts.InputTemplate, work_path)

        pwd = os.getcwd()
        os.chdir(work_path)
        
        vasp2abacus = Vasp2Abacus(vaspjobs, "/root/pporb/ABACUS-V1/pp", "/root/pporb/ABACUS-V1/orb", input_template)
        jobs, running_output = comm_func.capture_output(vasp2abacus.run)
        
        logs += running_output
        
        # pack jobs
        pack_filename = "ABACUS.zip"
        comm_func.pack(jobs, pack_filename)
        
        logs.iprint(f"VASP jobs converted to ABACUS jobs, and are packed in {pack_filename}")

        os.chdir(pwd)
        comm_report.gen_report(opts,logs,work_path,output_path,None)
        
    except:
        traceback.print_exc()
        return 1

    return 0
