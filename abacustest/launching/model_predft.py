import traceback
from dp.launching.typing.basic import BaseModel
import os

from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement


from . import (comm_class,
               comm_func,
               comm_class_exampleSource,
               comm_class_prepare,
               comm_class_predft,
               readsetting)


class PredftModel(
    comm_class_predft.PredftImageSet,
    comm_class_predft.PredftCommandSet,
    comm_class_exampleSource.PredftExtraFileSet,
    #comm_class_predft.PredftGroupSizeSet,
    comm_class_prepare.PrepareSet,
    comm_class_prepare.PrepareOrbLibSet,
    comm_class_prepare.PreparePPLibSet,
    comm_class_prepare.PrepareDPKSDescriptorSet,
    comm_class_prepare.PrepareKptTemplateSet,
    comm_class_prepare.PrepareStruTemplateSet,
    comm_class_prepare.PrepareInputTemplateSet,
    comm_class_exampleSource.PrepareExtraFileSet,
    comm_class_exampleSource.PrepareExampleSourceSet,
    comm_class.OutputSet,
    comm_class.ConfigSet,
        BaseModel):
    ...  

def PredftModelRunner(opts: PredftModel) -> int:
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
        
        #move results to output_path
        comm_func.move_results_to_output(work_path,output_path,allparams.get("save_path","results"))
        cwd = os.getcwd()
        os.chdir(os.path.join(output_path,allparams.get("save_path","results")))
        allfiles = os.listdir(".")
        alldirs = []
        for f in allfiles:
            if os.path.isdir(f):
                alldirs.append(f)
        packed_file_name = "examples.zip"
        comm_func.pack(alldirs,packed_file_name,"zip")
        os.chdir(cwd)
        
    except:
        traceback.print_exc()
        return 1

    return 0
