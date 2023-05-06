from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.typing import InputFilePath, OutputDirectory
from dp.launching.typing import (
    Set,
    Boolean,
    Field
)
from dp.launching.report import Report

from . import comm_class,comm_func,get_aim_data
import json,traceback,datetime,os,re,pickle
from abacustest.outresult import Table2FeishuInteractive
from dp.launching.report import ChartReportElement,ReportSection,AutoReportElement

NORMAL = {"iGM(SCF_steps)": 100*50,
          "iGM(total_time)": {
              "pw-cg": 100 / 1.555e-03,
              "pw-dav": 100 / 1.555e-03,
              "lcao-genelpa": 100 / 2.605e-03,
              "lcao-scalapack": 100 / 2.605e-03,
              "abacus-gnu-pw-cg": 100 / 1.555e-03,
              "abacus-gnu-pw-david": 100 / 1.555e-03,
              "gnu-lcao-elpa": 100 / 2.605e-03,
              "gnu-lcao-scalapack": 100 / 2.605e-03,
              "OMP-intel-cg": 100 / 1.555e-03,
              "OMP-intel-dav": 100 / 1.555e-03,
              "OMP-intel-elpa": 100 / 2.605e-03,
              "OMP-intel-scalapack": 100 / 2.605e-03,
              "OMP-gnu-cg": 100 / 1.555e-03,
              "OMP-gnu-dav": 100 / 2.605e-03,
              "OMP-gnu-elpa": 100 / 2.605e-03,
              "OMP-gnu-scalapack": 100 / 2.605e-03,
          }
}
TABLEHEAD = {
    "iGM(SCF_steps)": "SCFConverge Score",
    "iGM(total_time)": "Performance Score"
}
PROFILE_NAME = {
    "pw-cg": "intel-cg",
    "pw-dav": "intel-dav",
    "lcao-genelpa": "intel-elpa",
    "lcao-scalapack": "intel-scalapack",
    "abacus-gnu-pw-cg": "gnu-cg",
    "abacus-gnu-pw-david": "gnu-dav",
    "gnu-lcao-elpa": "gnu-elpa",
    "gnu-lcao-scalapack": "gnu-scalapack",
    "OMP-intel-cg": "omp-intel-cg",
    "OMP-intel-dav": "omp-intel-dav",
    "OMP-intel-elpa": "omp-intel-elpa",
    "OMP-intel-scalapack": "omp-intel-scalapack",
    "OMP-gnu-cg": "omp-gnu-cg",
    "OMP-gnu-dav": "omp-gnu-dav",
    "OMP-gnu-elpa": "omp-gnu-elpa",
    "OMP-gnu-scalapack": "omp-gnu-scalapack",
}
DIGIT = {"iGM(SCF_steps)": 0,
         "iGM(total_time)": 0,
         "NormalEnd_ratio": 2,
         'converge_ratio': 2
}

def get_profile_value(profile,allruninfos):
    '''
    allruninfos = {
        "run_name": run["props"]["name"],
        "experiment_name": run["props"]["experiment"]["name"],
        "tags": [i["name"] for i in run["props"]["tags"]],
        "creation_time": run["props"]["creation_time"],
        "end_time": run["props"]["end_time"],
        "metric": metric
    },

    metric = {metric_name: value}

    allvalues = {
        profile_i: {metric_i:[[create_time,metric_values,runname],...],}
    }
    '''
    allvalues = {}
    for i in profile: allvalues[i] = {}
    for run in allruninfos:
        #only get run_name endswith ".summary"
        runname = run["run_name"]
        if not run["run_name"].endswith(".summary"):
            continue
        tags = run["tags"]
        create_time = datetime.datetime.utcfromtimestamp(run["creation_time"]+8*3600)
 
        for i in profile:
            tag = 'benchmark-profile-' + i
            if tag in tags:
                for metric_name,metric_values in run["metric"].items():
                    if metric_name not in allvalues[i]:
                        allvalues[i][metric_name] = []
                    allvalues[i][metric_name].append([create_time,metric_values,runname])
                break
    return allvalues

def produce_outtable(allvalues):
    '''
    outtable = [["profile","date","profile_link",metrics_value]]
    metrics_value = a string for the final display
    This function will collect the last value, which may not be today
    '''

    allmetrics = []
    for k,v in allvalues.items():
        for ik,iv in v.items():
            if not ik.startswith("__system__") and ik not in allmetrics:
                allmetrics.append(ik)
            iv.sort()
    allmetrics.sort()
    
    outtable = [["profile","date","profile_link"]]
    for ii in allmetrics:
        outtable[0].append(TABLEHEAD.get(ii,ii))
    for iprofile,profile_value in allvalues.items():
        outtable.append([PROFILE_NAME.get(iprofile,iprofile),None,None])
        not_set_date = True
        for ik in allmetrics:
            value = "--"
            ivalue = profile_value.get(ik)
            if ivalue:
                if not_set_date:
                    outtable[-1][1] = ivalue[-1][0].strftime('%Y/%m/%d %H:%M')
                    outtable[-1][2] = "https://benchmark.mlops.dp.tech/?request=GET%3A%2Fprojects%2Fabacustest%2Fruns%2F" + ivalue[-1][2].split(".")[-2]
                    not_set_date = False
                
                if ik in NORMAL:
                    if isinstance(NORMAL[ik],dict):
                        coef = NORMAL[ik].get(iprofile,1.0)
                    else:
                        coef = NORMAL[ik]
                else:
                    coef = 1.0

                value0 = ivalue[-1][1][0]
                value_abs = ("%." + "%df"%DIGIT.get(ik,0)) % float(value0*coef)
                value = value_abs + "[--]"
                if len(ivalue) > 1:
                    value1 = ivalue[-2][1][0]
                    if value1 != 0 and isinstance(value0,(int,float)) and isinstance(value1,(int,float)):
                        value_rel = ("%" + ".2f" + "%" + "%") % ((value0 - value1)*100/value1)
                        if (value0 - value1) >= 0:
                            value_rel = "+" + value_rel
                        value = value_abs + f"[{value_rel}]"
                        if (value0 - value1) > 0:
                            value = "<font color='green'>" + value + "</font>"
                        elif (value0 - value1) < 0:
                            value = "<font color='red'>" + value + "</font>"
            outtable[-1].append(value)
    return outtable
        
def send_to_feishu(outtable,webhook,comment):
    today = datetime.datetime.now().strftime('%Y%m%d')
    new_table = [["profile"] + [TABLEHEAD.get(i,i) for i in outtable[0][3:]]]  #table head

    for itable in outtable[1:]:
        iprofile = itable[0]
        date0 = itable[1]
        link = itable[2]
        today_data = False
        if date0 != None and date0.split(":")[0] == today:
            today_data = True
            if link != None:
                iprofile = f"[{iprofile}]({link})"

        new_table.append([iprofile])
        for value in itable[3:]:
            if today_data: 
                new_table[-1].append(value)
            else:
                new_table[-1].append("--")

    Table2FeishuInteractive(new_table,webhook,title="abacus performance test %s" % today,comment=comment)
    return outtable

def produce_html_table(outtable,comment):
    new_table = [["profile"] + [TABLEHEAD.get(i,i) for i in outtable[0][3:]] + ["date"]]  #table head

    for itable in outtable[1:]:
        iprofile = itable[0]
        date0 = itable[1]
        link = itable[2]
        if link != None:
            iprofile = f"<a href=\"{link}\" target=\"detail\">{iprofile}</a>"

        new_table.append([iprofile])
        for value in itable[3:]:
            new_table[-1].append(value)
        new_table[-1].append(date0)

    html = '''<html><head><meta charset="UTF-8"></head><body>'''
    html += comm_func.produce_html_table(new_table) 
    html += "<p>" + "<br>".join(comment.split("\n")) + "</p>"
    html += '</body></html>'

    return html

class SummaryModel(BaseModel):
    feishu_webhook:  String = Field(default = None,
                                    title="FeiShu Webhook")
    AIM_TOKEN: String = Field(title="AIM tracking token")
    #experiment: String = Field(title="abacustest/benchmark")
    #experiment_id: String = Field(title="7ab4e46a-43fb-440a-828d-4fbdef5b4709")
    IO_output_path: OutputDirectory = Field(default="./output")

def echart_report(allvalues):
    #transfer the allvalues format
    
    chart_section = []
    for iprofile,profile_value in allvalues.items():
        metric_name = []
        metric_value = []
        for imetric,imetric_value in profile_value.items():
            ivalue = []
            if imetric in NORMAL:
                if isinstance(NORMAL[imetric],dict):
                    coef = NORMAL[imetric].get(iprofile,1.0)
                else:
                    coef = NORMAL[imetric]
            else:
                coef = 1.0
            for ii in imetric_value:
                ivalue.append([ii[0].strftime('%Y/%m/%d\n%H:%M'),ii[1][0]*coef])
                #ivalue.append([ii[0],ii[1][0]*coef])
            metric_name.append(TABLEHEAD.get(imetric,imetric))
            metric_value.append(ivalue)

        nmetrics = len(metric_name)
        inter = 5
        height = (0.9 - inter/100.0 * (nmetrics-1))/nmetrics*100
        interv = height + inter

        option = {
            "title": {
                "text": PROFILE_NAME.get(iprofile,iprofile),
                "left": 'center',
                "top": 0
            },
            "legend": {
                "data": metric_name,
                "top":"10%"
            },
            "toolbox": {
                "feature": {
                    "dataZoom": {
                        "yAxisIndex": 'none'
                    },
                    "restore": {},
                    "saveAsImage": {}
                }
            },
            "tooltip": {
                "trigger": 'axis',
                "axisPointer": {
                    "type": 'cross',
                    "animation": False,
                    "label": {
                        "backgroundColor": '#505765'
                    }
                }
            },
            "dataZoom": 
                {
                    "show": True,
                    "realtime": True,
                    "start":80,
                    "end":100,
                    "xAxisIndex":list(range(nmetrics)),
                    "left":"5%",
                    "right":"5%"
                }
            ,
            "grid": [{"left": f"{i*interv+5}%","right":f"{(nmetrics-i-1)*interv+5}%"} for i in range(nmetrics)],
            "xAxis": [{"type": "category", "gridIndex": i} for i in range(nmetrics)],
            "yAxis": [{"gridIndex": i,"scale": True,"name":metric_name[i],"nameLocation":"middle","axisLabel":{"inside":True}} for i in range(nmetrics)],
            "series": [
                {
                    "name": metric_name[i],
                    "type": 'line',
                    "xAxisIndex": i,
                    "yAxisIndex": i,
                    "data": metric_value[i],
                    "markLine": "markLineOpt"
                } for i in range(nmetrics)
            ]
        }
        chart_section.append(ChartReportElement(options=option,title=PROFILE_NAME.get(iprofile,iprofile)))

    return ReportSection(title="metrics chart",elements=chart_section,ncols=1)

def SummaryModelRunner(opts:SummaryModel):
    benchmark_profile = ["pw-cg", 
                         "pw-dav", 
                         "lcao-genelpa", 
                         "lcao-scalapack",
                         "abacus-gnu-pw-david",
                         "gnu-lcao-elpa",
                         "OMP-intel-dav",
                         "OMP-intel-elpa"]
    comment = "\n --- \n说明：\n- 括号里面的值为与上次结果相比增加的百分比\n"
    comment += "- cg/dav使用PW测试集，elpa/scalapack使用LCAO测试集\n"
    comment += "- omp指计算命令为OMP_NUM_THREADS=16 mpirun -np 1 abacus\n"
    comment += "- 不含omp指计算命令为OMP_NUM_THREADS=1 mpirun -np 16 abacus\n"
    comment += "- %s以50步作为100分，值越大，表明收敛需要的步数越少\n" % TABLEHEAD.get("iGM(SCF_steps)","iGM(SCF_steps)")
    comment += "- %s以abacus v3.0.0 的结果为基准，值越大表明总耗时越少，计算效率越高：\n" % TABLEHEAD.get("iGM(total_time)","iGM(total_time)")
    comment += "    - cg/dav以v3.0.0 intel-cg的结果为100分\n"
    comment += "    - genelpa/scalapack以v3.0.0 intel-elpa的结果为100分\n"
    comment += "- 可点击profile查看Benchmark run的详细结果"
    
    output_path = opts.IO_output_path
    os.makedirs(output_path,exist_ok=True)
    html_file_name = "result.html"
    html_file = os.path.join(output_path,html_file_name)

    token = opts.AIM_TOKEN
    experiment = "abacustest/benchmark"
    experiment_id = "7ab4e46a-43fb-440a-828d-4fbdef5b4709"

    if 1:
        allruns,allruninfos = get_aim_data.get_runs(token,experiment,experiment_id,collect_metrics=False)
        pickle.dump((allruns,allruninfos),open("tracking.pkl",'wb'))
    else:
        allruns,allruninfos = pickle.load(open("tracking.pkl",'rb'))

    allvalues = get_profile_value(benchmark_profile,allruninfos)

    #send to feishu and produce html file
    outtable = produce_outtable(allvalues)
    if opts.feishu_webhook:
        send_to_feishu(outtable,opts.feishu_webhook,comment)
    html_content = produce_html_table(outtable,comment)
    with open(html_file,'w') as f1: f1.write(html_content)

    #plot
    html_section = ReportSection(title="metrics chart",
                             elements=[AutoReportElement(title='metrics', path=html_file_name, description="")])
    chart_section = echart_report(allvalues)

    report = Report(title="abacus test report",
                        sections=[html_section,chart_section],
                        description="a report of abacustest")
    report.save(output_path)

    #parse inputs  
    return 0



