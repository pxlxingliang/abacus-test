import traceback,json,sys,os
from dp.launching.typing.basic import BaseModel,String,Float,Int
from dp.launching.typing import Field,DataSet
from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement


from . import (comm_class,
               comm_func,
               comm_pmetrics,
               comm_class_exampleSource) 
from abacustest.lib_model import model_007_FDStress as FDStress


class NewSetting(BaseModel):
    fd_step: Float = Field(default=0.0001,
                            title="FD Step (ratio of cell deformation)",
                            description="The finite difference step size. ",)
    
    fd_number: Int = Field(default=5,
                           titile= "FD Number",
                            description="The number of finite difference steps. Will calculate extra 2*Fd_Number structures.",)
    
    abacus_image: String = Field(default="registry.dp.tech/deepmodeling/abacus-intel:latest",
                          title="Abacus Image",
                          description="The image to run abaucs.",)
    abacus_command: String = Field(default="OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log",
                            title="Abacus Command",
                            description="The command to execute abacus",)
    bohrium_machine: String = Field(default="c32_m64_cpu",
                            title="Bohrium Machine",
                            description="The bohrium machine type to run abacus",)


class FDStressModel(
    #comm_class.ConfigSetGithub,
    #comm_class.TrackingSet,
    NewSetting,
    #comm_class_exampleSource.ScriptExampleSet,
    #comm_class_exampleSource.ScriptSourceSet,
    #comm_class_exampleSource.ScriptDatasetSet,
    #comm_class_exampleSource.ExampleSet,
    comm_class_exampleSource.ExampleSourceSet,
    #comm_class_exampleSource.DatasetSet,
    #MyModel,
    #ReuseDataset,
                    comm_class.OutputSet,
                    comm_class.ConfigSet,
                    BaseModel):
    ...  

def FDStressModelRunner(opts:FDStressModel) -> int:
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
        allexamples = [i for i in allfiles if os.path.isdir(i)]
        if len(allexamples) == 0:
            logs.iprint("no example found, exit! Please prepare the inputs of each example in each folder.")
            return 1
        subfolders = FDStress.PrepareFDStress(allexamples,opts.fd_step,opts.fd_number).run()
        if len(subfolders) == 0:
            logs.iprint("no example found, exit! Please prepare the inputs of each example in each folder.")
            return 1

        allparams = {
            "config": comm_func.read_config(opts),
            "save_path": ".",
            "bohrium_group_name": "FDStress",
            "run_dft": {
                "example": subfolders,
                "command": opts.abacus_command,
                "image": opts.abacus_image,
                "bohrium": {
                    "scass_type": opts.bohrium_machine,
                    "job_type": "container",
                    "platform": "ali",
                },
            },
        }
        os.chdir(pwd)

        # execut
        stdout,stderr = comm_func.exec_abacustest(allparams,work_path)
        os.chdir(work_path)
        FDStress.PostProcessFDStress(allexamples,False).run()
        os.chdir(pwd)
        

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

        # move results to output_path
        comm_func.move_results_to_output(work_path,output_path,allparams.get("save_path","results"))
        # comm_func.pack_results(output_path,allparams.get("save_path","results"))
    except:
        traceback.print_exc()
        return 1

    return 0
