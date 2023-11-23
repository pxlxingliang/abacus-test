import traceback
from dp.launching.typing.basic import BaseModel
import os

from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement


from . import (comm_class,
               comm_func,
               comm_pmetrics,
               comm_class_exampleSource,
               comm_class_rundft,
               comm_class_postdft,
               comm_class_metrics,
               readsetting)

class PostdftModel(comm_class.TrackingSet,
                    comm_class_metrics.metricsSaveFileSet,
                    comm_class_metrics.MetricsSet,
                    comm_class_postdft.PostdftImageSet,
                    comm_class_postdft.PostdftCommandSet,
                    comm_class_exampleSource.PostdftExtraFileSet,
                    comm_class_exampleSource.PostdftExampleSet,
                    comm_class_exampleSource.PostdftExampleSourceSet,
                    comm_class.OutputSet,
                    comm_class.ConfigSet,
                    BaseModel):
    ...  

class PostdftDatasetsModel(comm_class.TrackingSet,
                    comm_class_metrics.metricsSaveFileSet,
                    comm_class_metrics.MetricsSet,
                    comm_class_postdft.PostdftImageSet,
                    comm_class_postdft.PostdftCommandSet,
                    comm_class_exampleSource.PostdftExtraFileNeededSet,
                    comm_class_exampleSource.PostdftExampleSet,
                    comm_class_exampleSource.DatasetSet,
                    comm_class.OutputSet,
                    comm_class.ConfigSet,
                    BaseModel):
    ... 

def PostdftModelRunner(opts) -> int:
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
        reports = comm_pmetrics.produce_metrics_superMetrics_reports(allparams,work_path,output_path)

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
        comm_func.pack_results(output_path,allparams.get("save_path","results"))
    except:
        traceback.print_exc()
        return 1

    return 0
