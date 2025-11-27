import traceback
from dp.launching.typing.basic import BaseModel,String
from dp.launching.typing import Field

from . import comm_report
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE


from . import (comm_class,
               comm_func,
               comm_pmetrics,
               comm_class_exampleSource,
               comm_class_prepare,
               comm_class_predft,
               comm_class_rundft,
               comm_class_postdft,
               comm_class_metrics,
               readsetting)

class NewSetting(BaseModel):
    abacus_image: String = Field(default=RECOMMAND_IMAGE,
                          title="Abacus Image",
                          description="The image to run abaucs.",)
    abacus_command: String = Field(default=RECOMMAND_COMMAND,
                            title="Abacus Command",
                            description="The command to execute abacus",)
    bohrium_machine: String = Field(default=RECOMMAND_MACHINE,
                            title="Bohrium Machine",
                            description="The bohrium machine type to run abacus",)

class ExpertModel(
    #comm_class.TrackingSet,
    #comm_class_metrics.metricsSaveFileSet,
    comm_class_metrics.MetricsSet,
    #comm_class_postdft.PostdftImageSet,
    #comm_class_postdft.PostdftCommandSet,
    #comm_class_exampleSource.PostdftExtraFileSet,
    #comm_class_rundft.RundftImageSet,
    #comm_class_rundft.RundftCommandSet,
    #comm_class_rundft.RundftGroupSizeSet,
    #comm_class_exampleSource.RundftExtraFileSet,
    #comm_class_predft.PredftImageSet,
    #comm_class_predft.PredftCommandSet,
    #comm_class_exampleSource.PredftExtraFileSet,
    NewSetting,
    comm_class_prepare.PrepareSet,
    #comm_class_prepare.PrepareOrbLibSet,
    #comm_class_prepare.PreparePPLibSet,
    #comm_class_prepare.PrepareDPKSDescriptorSet,
    #comm_class_prepare.PrepareKptTemplateSet,
    #comm_class_prepare.PrepareStruTemplateSet,
    #comm_class_prepare.PrepareInputTemplateSet,
    #comm_class_exampleSource.PrepareExtraFileSet,
    #comm_class_exampleSource.PrepareExampleSet,
    comm_class_exampleSource.PrepareExampleSourceSet,
    comm_class.OutputSet,
    comm_class.ConfigSet,
        BaseModel):
    ...  

class ExpertDatasetsModel(
    comm_class.TrackingSet,
    comm_class_metrics.metricsSaveFileSet,
    comm_class_metrics.MetricsSet,
    comm_class_postdft.PostdftImageSet,
    comm_class_postdft.PostdftCommandSet,
    comm_class_exampleSource.PostdftExtraFileNeededSet,
    comm_class_rundft.RundftImageSet,
    comm_class_rundft.RundftCommandSet,
    comm_class_rundft.RundftGroupSizeSet,
    comm_class_exampleSource.RundftExtraFileNeededSet,
    comm_class_predft.PredftImageSet,
    comm_class_predft.PredftCommandSet,
    comm_class_exampleSource.PredftExtraFileNeededSet,
    comm_class_prepare.PrepareSet,
    comm_class_prepare.PrepareOrbLibPathSet,
    comm_class_prepare.PreparePPLibPathSet,
    comm_class_prepare.PrepareDPKSDescriptorPathSet,
    comm_class_prepare.PrepareKptTemplatePathSet,
    comm_class_prepare.PrepareStruTemplatePathSet,
    comm_class_prepare.PrepareInputTemplatePathSet,
    comm_class_exampleSource.PrepareExtraFileNeededSet,
    comm_class_exampleSource.PrepareExampleSet,
    comm_class_exampleSource.DatasetSet,
    comm_class.OutputSet,
    comm_class.ConfigSet,
        BaseModel):
    ...  


def ExpertModelRunner(opts) -> int:
    try:
        logs = comm_class.myLog()

        paths = comm_func.create_path(str(opts.IO_output_path))
        output_path = paths["output_path"]
        work_path = paths["work_path"]
        download_path = paths["download_path"]

        allparams = readsetting.ReadSetting(logs,opts,work_path,download_path)
        if allparams == None:
            return 1
        
        allparams["run_dft"] = {
            "command": opts.abacus_command,
            "image": opts.abacus_image,
            "bohrium": {
                "scass_type": opts.bohrium_machine,
                "job_type": "container",
                "platform": "ali",
            },
        }
        allparams["post_dft"]["image"] = "registry.dp.tech/dptech/abacustest:latest"
        

        stdout,stderr = comm_func.exec_abacustest(allparams,work_path)
        logs.iprint(f"{stdout}\n{stderr}\nrun abacustest over!\n")
        
        comm_report.gen_report(opts,logs,work_path,output_path,allparams)
        
        
    except:
        traceback.print_exc()
        return 1

    return 0
