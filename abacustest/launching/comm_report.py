import os
from .comm_class import OutputSet,ConfigSet,ConfigSetGithub
from . import comm_pmetrics,comm_func

from dp.launching.report import Report,AutoReportElement,ReportSection,ChartReportElement

def gen_report_param(params,show_param=None):
    """
    generate a report section for the parameters
    """
    # get all member of params
    report_member = "Input Parameters:\n"
    if show_param is None:
        config_members = list(ConfigSet.__fields__.keys()) + list(ConfigSetGithub.__fields__.keys()) + list(OutputSet.__fields__.keys())

        
        for im in params.__fields__.keys():
            if im not in config_members:
                report_member += f"{im}={getattr(params,im)}\n"
    else:
        for im in show_param:
            if im in params.__fields__.keys():
                report_member += f"{im}={getattr(params,im)}\n"
    
    report_member += "\n"
    return report_member

def gen_report(opts,logs,work_path,output_path,abacustest_param):
    param_c = gen_report_param(opts)
    with open(os.path.join(output_path,"ABACUSTEST_INPUT_PARAM.txt"),"w") as f:
        f.write(param_c)
    
    reports = [
        ReportSection(title="INPUT PARAMETER",
                      elements=[AutoReportElement(title='', path="ABACUSTEST_INPUT_PARAM.txt", description="")])
    ]
    if abacustest_param is not None:
        reports += comm_pmetrics.produce_metrics_superMetrics_reports(abacustest_param,work_path,output_path)
    logfname = "output.log"
    logs.write(os.path.join(output_path,logfname))
    log_section = ReportSection(title="",
                              elements=[AutoReportElement(title='', path=logfname, description="")])
    reports.append(log_section)
    if reports:
        report = Report(title="abacus test report",
                        sections=reports,
                        description="a report of abacustest")
        report.save(output_path)
    # move results to output_path
    if abacustest_param is not None:
        comm_func.move_results_to_output(work_path,output_path,abacustest_param.get("save_path","results"))
    else:
        os.system(f"mv {work_path}/* {output_path}/")
    #comm_func.pack_results(output_path,abacustest_param.get("save_path","results"))