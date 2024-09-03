import traceback,json
from dp.launching.typing.basic import BaseModel,String
from dp.launching.typing import Field
import os

from abacustest.report import gen_html


from . import (comm_class,
               comm_func,
               comm_report,
               comm_class_exampleSource,
               readsetting)

class CommandSet(BaseModel):
    command: String = Field(default="",
                            title="Command(optional)",
                                    description="bash command if needed",)


class ReportModel(
    #CommandSet,
    #comm_class_exampleSource.ExampleSet,
    comm_class_exampleSource.ExampleSourceSet,
    #comm_class_exampleSource.DatasetSet,
                    comm_class.OutputSet,
                    comm_class.ConfigSet,
                    BaseModel):
    ...  

def ReportModelRunner(opts:ReportModel) -> int:
    try:
        logs = comm_class.myLog()

        paths = comm_func.create_path(str(opts.IO_output_path))
        output_path = paths["output_path"]
        work_path = paths["work_path"]
        download_path = paths["download_path"]

        datas = comm_class_exampleSource.read_source(opts,work_path,download_path,logs.iprint)
        #if opts.command and opts.command.strip() != "":
        #    cwd = os.getcwd()
        #    os.chdir(work_path)
        #    try:
        #        return_code, stdout, stderr = comm_func.run_command(opts.command.strip())
        #        logs.iprint(f"{stdout}\n{stderr}\nrun abacustest over!\n")
        #    except:
        #        traceback.print_exc()
        #        return 1
        #    os.chdir(cwd)
        
        if os.path.isfile(os.path.join(work_path,"setting.json")):
            try:
                setting = json.load(open(os.path.join(work_path,"setting.json")))
            except:
                setting = {}
            if "report" in setting:
                report_setting = setting["report"]  
                pwd = os.getcwd()
                os.chdir(work_path)
                gen_html(report_setting,"abacustest.html")
                os.chdir(pwd)
            
        allparams = {"save_path": "",}


        comm_report.gen_report(opts,logs,work_path,output_path,allparams)
        
    except:
        traceback.print_exc()
        return 1

    return 0
