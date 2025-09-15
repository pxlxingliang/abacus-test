import traceback,json,sys,os
from dp.launching.typing.basic import BaseModel,String,Float,Int,Boolean,Set
from dp.launching.typing import Field
from enum import Enum
from abacustest.lib_model.model_013_inputs import PrepInput
from abacustest.lib_model.model_012_band import PostBand

from abacustest.constant import RECOMMAND_COMMAND, RECOMMAND_IMAGE, RECOMMAND_MACHINE

from . import (comm_class,
               comm_func,
               comm_report,
               comm_class_exampleSource) 

class StruTypeEnum(String, Enum):
    cif = "cif"
    poscar = "poscar"
    stru = "abacus/stru"

class PPTypeEnum(String, Enum):
    apnsv1 = "apns-V1"
    
def parse_pp_type(pp_type):
    if pp_type == "apns-V1":
        return "/root/pporb/ABACUS-V1/pp","/root/pporb/ABACUS-V1/orb"
    else:
        raise ValueError(f"Unknown pp type: {pp_type}")

class ABACUSImage(String, Enum):
    ltsv310 = "registry.dp.tech/dptech/abacus-stable:LTSv3.10"
    latest_intel = "registry.dp.tech/deepmodeling/abacus-intel:latest"
    latest_gnu = "registry.dp.tech/deepmodeling/abacus-gnu:latest"
    latest_cuda = "registry.dp.tech/deepmodeling/abacus-cuda:latest"

class JobTypeEnum(String, Enum):
    scf = "单点计算"
    relax = "原子弛豫"
    cellrelax = "晶格弛豫"
    band = "能带计算"

def parse_job_type(job_type):
    if job_type == "单点计算":
        return "scf"
    elif job_type == "原子弛豫":
        return "relax"
    elif job_type == "晶格弛豫":
        return "cell-relax"
    elif job_type == "能带计算":
        return "band"
    else:
        raise ValueError(f"Unknown job type: {job_type}")

class BasisTypeEnum(String, Enum):
    pw = "pw"
    lcao = "lcao"

class NewSetting(BaseModel):
    stru_type: StruTypeEnum = Field(default="cif",
                            title="Structure format",
                            description="The format of the structures",)
    
    pp_type: PPTypeEnum = Field(default="apns-V1",
                            title="Pseudopotential and orbital type",
                            description="Please choose one of the pseudopotential and orbital type",)
    
    job_type: JobTypeEnum = Field(default="单点计算",
                            title="Job type",
                            description="Please choose one of the job type",)
    
    basis_type: BasisTypeEnum = Field(default="lcao",
                            title="Basis type",
                            description="Please choose one of the basis type",)
    
    abacus_image: ABACUSImage = Field(default=RECOMMAND_IMAGE,
                          title="Abacus Image",
                          description="The image to run abaucs.",)
    
    abacus_command: String = Field(default=RECOMMAND_COMMAND,
                            title="Abacus Command",
                            description="The command to execute abacus",)
    
    bohrium_machine: String = Field(default=RECOMMAND_MACHINE,
                            title="Bohrium Machine",
                            description="The bohrium machine type to run abacus",)


def post_process(job_type, job_path, work_path):
    """
    Post process the job
    """
    if job_type in ["scf", "relax", "cell-relax"]:
        pass
    elif job_type == "band":
        pwd = os.getcwd()
        os.chdir(work_path)
        PostBand(job_path).run()
        json.dump({
            "band structure of " + i: {"type": "image", "file": os.path.join(i, "band.png")} for i in job_path
        },
        open("supermetrics.json", "w"),indent=4)
        os.chdir(pwd)
    else:
        raise ValueError(f"Unknown job type: {job_type}")

class AutoABACUSModel(
    NewSetting,
    comm_class_exampleSource.ExampleSourceSet,
                    comm_class.OutputSet,
                    comm_class.ConfigSet,
                    BaseModel):
    ...  

def AutoABACUSRunner(opts:AutoABACUSModel) -> int:
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
        
        job_type = parse_job_type(opts.job_type)
        pppath, orbpath = parse_pp_type(opts.pp_type)
        allparams, job_path = PrepInput(allfiles, opts.stru_type, 
                           pp_path=pppath,
                           orb_path=orbpath, 
                           abacus_command=opts.abacus_command,
                           machine=opts.bohrium_machine,
                           image=opts.abacus_image,
                           jobtype=job_type,
                           lcao=True if opts.basis_type == "lcao" else False,
                           ).run()
        
        allparams["config"] = comm_func.read_config(opts)
        allparams["save_path"] = "."
        allparams["bohrium_group_name"] = "AutoAbacus"
        allparams["post_dft"] = {
            "image": "registry.dp.tech/dptech/abacustest:latest",
	        "metrics": {
                "dft_type": "abacus",
                "metrics_name": [
                    {"label": "''.join([i+str({label}.count(i)) for i in set({label})])"},
                    "normal_end",
                    "converge",
                    "nkstot",
                    "ibzk",
                    "nbands",
                    "nelec",
                    "natom",
                    "scf_steps",
                    "total_time",
                    "scf_time",
                    {
                        "scf_time/step": "{scf_time}/{scf_steps}",
                        "ks_solver": "{INPUT}['ks_solver']"
                    },
                    "stress_time",
                    "energy",
                    "energy_per_atom",
                    "force",
                    "stress",
                    "drho_last",
                    "denergy_last",
                    "version"
                ],
                "save_file": "metrics.json"
            }
        }
        
        os.chdir(pwd)

        # execut
        stdout,stderr = comm_func.exec_abacustest(allparams,work_path)
        post_process(job_type, job_path, work_path)
        os.chdir(pwd)
        comm_report.gen_report(opts,logs,work_path,output_path,allparams)
    except:
        traceback.print_exc()
        return 1

    return 0
