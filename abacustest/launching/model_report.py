from genericpath import isdir
import traceback
from dp.launching.typing.basic import BaseModel,String
from dp.launching.typing import Field
import os

from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement


from . import (comm_class,
               comm_func,
               comm_class_exampleSource,
               readsetting)

class CommandSet(BaseModel):
    command: String = Field(default="",
                                    description="bash command",)


class ReportModel(
    CommandSet,
    comm_class_exampleSource.ExampleSet,
    comm_class_exampleSource.ExampleSourceSet,
    comm_class_exampleSource.DatasetSet,
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
        if opts.command and opts.command.strip() != "":
            cwd = os.getcwd()
            os.chdir(work_path)
            try:
                return_code, stdout, stderr = comm_func.run_command(opts.command.strip())
                logs.iprint(f"{stdout}\n{stderr}\nrun abacustest over!\n")
            except:
                traceback.print_exc()
                return 1
            os.chdir(cwd)
            
        allparams = {"save_path": "",}


        reports = comm_func.produce_metrics_superMetrics_reports(allparams,work_path,output_path)

        logfname = "output.log"
        logs.write(os.path.join(str(opts.IO_output_path),logfname))
        log_section = ReportSection(title="",
                                  elements=[AutoReportElement(title='', path=logfname, description="")])
        reports.append(log_section)

        if reports:
            report = Report(title="abacus test report",
                            sections=reports,
                            description="a report of abacustest")
            report.save(output_path)
            
        #move results to output_path
        comm_func.move_results_to_output(work_path,output_path,allparams.get("save_path","results"))
        
        cwd = os.getcwd()
        os.chdir(os.path.join(output_path,allparams.get("save_path","results")))
        allfiles = os.listdir(".")
        alldirs = []
        for f in allfiles:
            if os.path.isdir(f):
                alldirs.append(f)
        packed_file_name = "results.zip"
        comm_func.pack(alldirs,packed_file_name,"zip")
        os.chdir(cwd)
    except:
        traceback.print_exc()
        return 1

    return 0
