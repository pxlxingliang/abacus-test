import traceback,json,sys
from dp.launching.typing.basic import BaseModel,String
from dp.launching.typing import Field,DataSet
import os

from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement


from . import (comm_class,
               comm_func,
               comm_pmetrics,
               comm_class_exampleSource,
               readsetting) 

class RundftImage(BaseModel):
    rundft_image: String = Field(default="",
                          title="Rundft Image",
                          description="If you do not want to use the rundft image defined in setting file, please fill in the image address here",)
    
class ReuseModel(
    comm_class.TrackingSet,
    RundftImage,
    comm_class_exampleSource.ScriptExampleSet,
    comm_class_exampleSource.ScriptSourceSet,
    comm_class_exampleSource.ScriptDatasetSet,
    comm_class_exampleSource.ExampleSet,
    comm_class_exampleSource.ExampleSourceSet,
    comm_class_exampleSource.DatasetSet,
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
            print("\nERROR: Please supply the setting information!!!!")
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

        tracking_set = comm_class.TrackingSet.parse_obj(opts)
        if tracking_set:
            if "aim_access_token" not in allparams["config"]:
                allparams["config"]["aim_access_token"] = str(tracking_set.get("token"))
            if "post_dft" not in allparams:
                allparams["post_dft"] = {}
            if "upload_tracking" not in allparams["post_dft"]:
                allparams["post_dft"]["upload_tracking"] = {}
            if "name" not in allparams["post_dft"]["upload_tracking"]:
                allparams["post_dft"]["upload_tracking"]["name"] = tracking_set.get("name")
            if "experiment" not in allparams["post_dft"]["upload_tracking"]:
                allparams["post_dft"]["upload_tracking"]["experiment"] = tracking_set.get("experiment")
            if "tags" not in allparams["post_dft"]["upload_tracking"]:
                allparams["post_dft"]["upload_tracking"]["tags"] = tracking_set.get("tags")

        # modify the image of rundft
        if opts.rundft_image != None and opts.rundft_image.strip() != "":
            new_image = opts.rundft_image.strip()
            if "run_dft" in allparams:
                if isinstance(allparams["run_dft"],list):
                    for idx in range(len(allparams["run_dft"])):
                        if "image" in allparams["run_dft"][idx]:
                            allparams["run_dft"][idx]["image"] = new_image
                elif isinstance(allparams["run_dft"],dict):
                    if "image" in allparams["run_dft"]:
                        allparams["run_dft"]["image"] = new_image
                
        #execut
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
