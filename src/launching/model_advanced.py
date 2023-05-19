import traceback
from typing import Literal

from zmq import has

from dp.launching.cli import to_runner,SubParser,run_sp_and_exit
from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.cli import to_runner, default_minimal_exception_handler
from dp.launching.typing import InputFilePath, OutputDirectory
from dp.launching.typing import (
    BaseModel,
    Set,
    Boolean,
    Field,
    DflowAccessToken,
    DflowArgoAPIServer,
    DflowK8sAPIServer,
    DflowStorageEndpoint,
    DflowStorageRepository,
    BohriumMachineType,
    BohriumImage,
    BohriumPlatform,
    BohriumJobType,
    BohriumUsername,
    BohriumPassword,
    BohriumProjectId,
    BenchmarkLabels,
    BenchmarkTags,
)
import os,shutil,glob

from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement


from . import (comm_class,
               comm_func,
               comm_class_exampleSource,
               comm_class_rundft,
               comm_class_postdft,
               comm_class_metrics,
               readsetting)

class AdvancedModel(comm_class.TrackingSet,
                    comm_class_metrics.metricsSaveFileSet,
                    comm_class_metrics.MetricsSet,
                    comm_class_postdft.PostdftImageSet,
                    comm_class.PostdftCommandSet,
                    comm_class_exampleSource.PostdftExtraFileSet,
                    comm_class_rundft.RundftImageSet,
                    comm_class.RundftCommandSet,
                    comm_class.NgroupSet,
                    comm_class_exampleSource.RundftExtraFileSet,
                    comm_class_exampleSource.ExampleSourceSet,
                    comm_class.OutputSet,
                    comm_class.ConfigSet,
                    BaseModel):
    ...  

def AdvancedModelRunner(opts: AdvancedModel) -> int:
    try:
        logs = comm_class.myLog()

        paths = comm_func.create_path(str(opts.IO_output_path))
        output_path = paths["output_path"]
        work_path = paths["work_path"]
        download_path = paths["download_path"]

        allparams = readsetting.ReadSetting(logs,opts,work_path,download_path)
        if allparams == None:
            return 1

        stdout,stderr = comm_func.exec_abacustest(allparams,work_path)
        logs.iprint(f"{stdout}\n{stderr}\nrun abacustest over!\n")
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
    except:
        traceback.print_exc()
        return 1

    return 0
