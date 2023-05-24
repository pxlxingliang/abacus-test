from sympy import intersection
from dp.launching.typing.basic import BaseModel, Int, String
from dp.launching.typing import (
    Set,
    Boolean,
    Field,
    BenchmarkLabels
)
from dp.launching.report import Report

from . import comm_class,comm_func,get_aim_data
import json,traceback,datetime,os,re,pickle
from abacustest.outresult import Table2FeishuInteractive
from dp.launching.report import ChartReportElement,ReportSection,AutoReportElement


def transfer_tag(old_tag):
    """
    transfer:
        "aim_tag":{
        "intel-cg": ["benchmark-profile-pw-cg","benchmark-schedule-intel-cg"],
        "intel-dav": ["benchmark-profile-pw-dav","benchmark-schedule-intel-dav"]}
    to:
        "aim_tag":{
        "benchmark-profile-pw-cg": "intel-cg",
        "benchmark-schedule-intel-cg": "intel-cg",
        "benchmark-profile-pw-dav": "intel-dav",
        "benchmark-schedule-intel-dav": "intel-dav"}
    """
    new_tag = {}
    for k,v in old_tag.items():
        for i in v:
            new_tag[i] = k
    return new_tag

def clean_metrics_coef(profile,metrics,metrics_coef):
    """
    return {(iprofile,imetrics):coef}
    """
    new_metrics_coef = {}
    for i in profile:
        for j in metrics:
            if j not in metrics_coef:
                new_metrics_coef[(i,j)] = 1.0
            elif isinstance(metrics_coef[j],dict):
                new_metrics_coef[(i,j)] = metrics_coef[j].get(i,1.0)
            else:
                new_metrics_coef[(i,j)] = metrics_coef[j]
    return new_metrics_coef

def get_profile_value(allruninfos,profile,aim_tag,metrics,metrics_coef,metrics_name,token):
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
    value is a list

    allvalues = {
        profile_i: {metric_i:[[create_time,metric_values,runname],...],}
    }
    profile_i and metric_i are the names in profile and metrics, which are the same as the name in final table.
    metric_values = value * coef
    '''
    new_tag = transfer_tag(aim_tag)
    new_metrics_coef = clean_metrics_coef(profile,metrics,metrics_coef)
    
    allvalues = {}
    for i in profile: allvalues[i] = {}
    for run in allruninfos:
        #only get run_name endswith ".summary"
        if not run["run_name"].endswith(".summary"):
            continue
        runname = run["run_name"].split(".")[1]
        if runname.startswith("schedule-"):
            url = "https://launching.mlops.dp.tech/?request=GET%3A%2Fapplications%2Fabacustest%2Fjobs%2F" + runname
        else:
            url = "https://benchmark.mlops.dp.tech/?request=GET%3A%2Fprojects%2Fabacustest%2Fruns%2F" + runname
        tags = run["tags"]
        create_time = datetime.datetime.utcfromtimestamp(run["creation_time"]+8*3600)

        for itag in tags:
            if itag in new_tag and new_tag[itag] in allvalues:
                iprofile = new_tag[itag]
                for imetric in metrics:
                    for metric_in_aim in metrics_name.get(imetric,[imetric]):
                        if metric_in_aim in run["metric"]:
                            if imetric not in allvalues[iprofile]:
                                allvalues[iprofile][imetric] = []
                            allvalues[iprofile][imetric].append([create_time,
                                                                 run["metric"][metric_in_aim][0] * new_metrics_coef[(iprofile,imetric)],
                                                                 url, run["run_hash"]]) #the last one should be version, store hash now
                            break
                break
    versions = {}
    for k,v in allvalues.items():
        hash = None
        for ik,iv in v.items():
            iv.sort()
            hash = iv[-1][-1]
            #find the version of last iv
            if hash not in versions:
                versions[hash] = get_aim_data.get_version(hash,token)
            iv[-1][-1] = None if not versions[hash] else versions[hash][0]
            #set version to be None, except for the last iv
            for i in iv[:-1]:
                i[-1] = None
    return allvalues

def produce_outtable(allvalues,profile,metrics,digit):
    '''
    outtable = [["profile","date","profile_link",metrics_value]]
    metrics_value = a string for the final display
    This function will collect the last value, which may not be today
    '''
    
    outtable = [["profile","date","profile_link","version(commit)"]]
    for ii in metrics:
        outtable[0].append(ii)
        
    for iprofile in profile:
        if iprofile not in allvalues:
            continue
        profile_value = allvalues[iprofile]
        
        outtable.append([iprofile,None,None,None])
        not_set_date = True
        for ik in metrics:
            value = "--"
            ivalue = profile_value.get(ik)
            if ivalue:
                if not_set_date:
                    outtable[-1][1] = ivalue[-1][0].strftime('%Y/%m/%d %H:%M')
                    outtable[-1][2] = ivalue[-1][2]
                    outtable[-1][3] = ivalue[-1][3]
                    not_set_date = False

                value0 = ivalue[-1][1]
                value_abs = ("%." + "%df"%digit.get(ik,0)) % float(value0)
                value = value_abs + "[--]"
                if len(ivalue) > 1:
                    value1 = ivalue[-2][1]
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
    today = datetime.datetime.now().strftime('%Y/%m/%d')
    new_table = [["profile"] + [i for i in outtable[0][4:]]]  #table head

    version_list = []
    profile_list = []
    for itable in outtable[1:]:
        iprofile = itable[0]
        date0 = itable[1]
        link = itable[2]
        profile_list.append(iprofile)
        today_data = False
        if date0 != None and date0.split()[0] == today:
            today_data = True
            version_list.append(itable[3])
            if link != None:
                iprofile = f"[{iprofile}]({link})"

        new_table.append([iprofile])
        for value in itable[4:]:
            if today_data: 
                new_table[-1].append(value)
            else:
                new_table[-1].append("--")

    new_comment = ""
    if len(set(version_list)) == 1 and version_list[0] != None:
        new_comment += f"version: {version_list[0]}\n"
    elif len(set(version_list)) > 1:
        version_list_remove_none = set(version_list) - set([None])
        if len(version_list_remove_none) == 1:
            new_comment += f"version: {list(version_list_remove_none)[0]}\n"
        else:
            new_comment += f"version: {version_list}\n"
        
    
    Table2FeishuInteractive(new_table,webhook,title="abacus performance test %s" % today,comment=new_comment + comment)
    return outtable

def produce_html_table(outtable,comment):
    new_table = [["profile"] + [i for i in outtable[0][4:]] + ["test_date","version"]]  #table head

    for itable in outtable[1:]:
        iprofile = itable[0]
        date0 = itable[1]
        link = itable[2]
        version = itable[3]
        if link != None:
            iprofile = f"<a href=\"{link}\" target=\"detail\">{iprofile}</a>"

        new_table.append([iprofile])
        for value in itable[4:]:
            new_table[-1].append(value)
        new_table[-1].append(date0)
        
        if version != None:
            pattern = r"v([\d.]+)\((\w+)\s\((.+)\)\)"
            result = re.match(pattern,str(version))
            if result:
                commit = result.group(2)
                github_link = "https://github.com/deepmodeling/abacus-develop/commit/" + commit
                version = f"<a href=\"{github_link}\" target=\"detail\">{version}</a>"
        new_table[-1].append(str(version))

    html = '''<html><head><meta charset="UTF-8"></head><body>'''
    html += comm_func.produce_html_table(new_table) 
    html += "<p>" + "<br>".join(comment.split("\n")) + "</p>"
    html += '</body></html>'

    return html

class Summary(BaseModel):
    feishu_webhook:  String = Field(default = None,
                                    title="FeiShu Webhook")
    AIM_TOKEN: String = Field(title="AIM tracking token")
    setting: String = Field(default = None,title="setting json string")
    Config_dflow_labels: BenchmarkLabels
    #experiment: String = Field(title="abacustest/benchmark")
    #experiment_id: String = Field(title="7ab4e46a-43fb-440a-828d-4fbdef5b4709")

class SummaryModel(Summary,comm_class.OutputSet,BaseModel):
    ...

def echart_report(allvalues):
    #transfer the allvalues format
    
    chart_section = []
    for iprofile,profile_value in allvalues.items():
        metric_name = []
        metric_value = []
        for imetric,imetric_value in profile_value.items():
            ivalue = []
            for ii in imetric_value:
                ivalue.append([ii[0].strftime('%Y/%m/%d\n%H:%M'),ii[1]])
                #ivalue.append([ii[0],ii[1][0]*coef])
            metric_name.append(imetric)
            metric_value.append(ivalue)

        nmetrics = len(metric_name)
        inter = 5
        height = (0.9 - inter/100.0 * (nmetrics-1))/nmetrics*100
        interv = height + inter

        option = {
            "title": {
                "text": iprofile,
                "left": 'center',
                "top": 0
            },
            "legend": {
                "data": metric_name,
                "top":"10%"
            },
            "toolbox": {
                "feature": {
                    "dataZoom": {},
                    "dataView": {"readOnly": False},
                    "magicType": {"type": ['line', 'bar']},
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
        chart_section.append(ChartReportElement(options=option,title=iprofile))

    return ReportSection(title="metrics chart",elements=chart_section,ncols=1)

def SummaryModelRunner(opts:SummaryModel):
    comment = """
说明：
- 括号里面的值为与上次结果相比增加的百分比
- cg/dav使用PW测试集，elpa/scalapack使用LCAO测试集
- omp指计算命令为OMP_NUM_THREADS=16 mpirun -np 1 abacus
- 不含omp指计算命令为OMP_NUM_THREADS=1 mpirun -np 16 abacus
- SCFConverge Score以50步作为100分，值越大，表明收敛需要的步数越少
- Performance Score以abacus v3.0.0 的结果为基准，值越大表明总耗时越少，计算效率越高：
    - cg/dav以v3.0.0 intel-cg的结果为100分
    - genelpa/scalapack以v3.0.0 intel-elpa的结果为100分
- 可点击profile查看Benchmark run的详细结果    
    """
    setting = {
        "experiment": "abacustest/benchmark",
        "experiment_id": "7ab4e46a-43fb-440a-828d-4fbdef5b4709",
        "profile": [
            "intel-cg",
            "intel-dav",
            "intel-elpa",
            "intel-scalapack",
            #            "gnu-cg",
            "gnu-dav",
            "gnu-elpa",
            #            "gnu-scalapack",
            #            "omp-intel-cg",
            "omp-intel-dav",
            "omp-intel-elpa",
            #           "omp-intel-scalapack",
            #           "omp-gnu-cg",
            #           "omp-gnu-dav",
            #           "omp-gnu-elpa",
            #          "omp-gnu-scalapack",
            "exx-test"
        ],
        "aim_tag": {
            "intel-cg": ["benchmark-profile-pw-cg", "benchmark-schedule-intel-cg"],
            "intel-dav": ["benchmark-profile-pw-dav", "benchmark-schedule-intel-dav"],
            "intel-elpa": ["benchmark-profile-lcao-genelpa", "benchmark-schedule-intel-elpa"],
            "intel-scalapack": ["benchmark-profile-lcao-scalapack", "benchmark-schedule-intel-scalapack"],
            "gnu-cg": ["benchmark-profile-abacus-gnu-pw-cg", "benchmark-schedule-gnu-cg"],
            "gnu-dav": ["benchmark-profile-abacus-gnu-pw-david", "benchmark-schedule-gnu-dav"],
            "gnu-elpa": ["benchmark-profile-gnu-lcao-elpa", "benchmark-schedule-gnu-elpa"],
            "gnu-scalapack": ["benchmark-profile-gnu-lcao-scalapack", "benchmark-schedule-gnu-scalapack"],
            "omp-intel-cg": ["benchmark-profile-OMP-intel-cg", "benchmark-schedule-omp-intel-cg"],
            "omp-intel-dav": ["benchmark-profile-OMP-intel-dav", "benchmark-schedule-omp-intel-dav"],
            "omp-intel-elpa": ["benchmark-profile-OMP-intel-elpa", "benchmark-schedule-omp-intel-elpa"],
            "omp-intel-scalapack": ["benchmark-profile-OMP-intel-scalapack", "benchmark-schedule-omp-intel-scalapack"],
            "omp-gnu-cg": ["benchmark-profile-OMP-gnu-cg", "benchmark-schedule-omp-gnu-cg"],
            "omp-gnu-dav": ["benchmark-profile-OMP-gnu-dav", "benchmark-schedule-omp-gnu-dav"],
            "omp-gnu-elpa": ["benchmark-profile-OMP-gnu-elpa", "benchmark-schedule-omp-gnu-elpa"],
            "omp-gnu-scalapack": ["benchmark-profile-OMP-gnu-scalapack", "benchmark-schedule-omp-gnu-scalapack"],
            "exx-test": ["benchmark-schedule-exx-test"]
        },
        "metrics": ["NormalEnd_ratio", "converge_ratio", "SCFConverge Score", "Performance Score"],
        "metrics_name": {
            "NormalEnd_ratio": ["NormalEnd_ratio","TrueRatio(normal_end)"],
            "converge_ratio": ["converge_ratio","TrueRatio(converge)"],
            "SCFConverge Score": ["iGM(SCF_steps)","iGM(scf_steps)"],
            "Performance Score": ["iGM(total_time)"]
        },
        "metrics_coef": {
            "NormalEnd_ratio": 1,
            "converge_ratio": 1,
            "SCFConverge Score": 100*50,
            "Performance Score": {
                "intel-cg": 100 / 1.555e-03,
                "intel-dav": 100 / 1.555e-03,
                "intel-elpa": 100 / 2.605e-03,
                "intel-scalapack": 100 / 2.605e-03,
                "gnu-cg": 100 / 1.555e-03,
                "gnu-dav": 100 / 1.555e-03,
                "gnu-elpa": 100 / 2.605e-03,
                "gnu-scalapack": 100 / 2.605e-03,
                "omp-intel-cg": 100 / 1.555e-03,
                "omp-intel-dav": 100 / 1.555e-03,
                "omp-intel-elpa": 100 / 2.605e-03,
                "omp-intel-scalapack": 100 / 2.605e-03,
                "omp-gnu-cg": 100 / 1.555e-03,
                "omp-gnu-dav": 100 / 1.555e-03,
                "omp-gnu-elpa": 100 / 2.605e-03,
                "omp-gnu-scalapack": 100 / 2.605e-03
            }
        },
        "digit": {
            "NormalEnd_ratio": 2,
            "converge_ratio": 2,
            "SCFConverge Score": 0,
            "Performance Score": 0
        },
        "comment": comment
    }
    output_path = opts.IO_output_path
    os.makedirs(output_path,exist_ok=True)
    html_file_name = "result.html"
    html_file = os.path.join(output_path,html_file_name)
    
    token = opts.AIM_TOKEN
    if opts.setting != None and opts.setting.strip() != "":
        setting = json.loads(opts.setting)
    experiment = setting.get("experiment","abacustest/benchmark")
    experiment_id = setting.get("experiment_id","7ab4e46a-43fb-440a-828d-4fbdef5b4709")
    profile = setting.get("profile",[])
    aim_tag = setting.get("aim_tag",{})
    metrics = setting.get("metrics",[])
    metrics_name = setting.get("metrics_name",{})
    metrics_coef = setting.get("metrics_coef",{})
    digit = setting.get("digit",{})
    comment = setting.get("comment","")
    
    json.dump(setting,open(os.path.join(output_path,"setting.json"),'w'),indent=4)

    if 1:
        alltags = []
        for i in profile:
            alltags += aim_tag.get(i,[i])
        allruns,allruninfos = get_aim_data.get_runs(token,experiment,experiment_id,alltags,collect_metrics=False,GetVersion=False)
        pickle.dump((allruns,allruninfos),open("tracking.pkl",'wb'))
    else:
        allruns,allruninfos = pickle.load(open("tracking.pkl",'rb'))

    allvalues = get_profile_value(allruninfos,profile,aim_tag,metrics,metrics_coef,metrics_name,token)

    
    #send to feishu and produce html file
    outtable = produce_outtable(allvalues,profile,metrics,digit)
    if opts.feishu_webhook:
        current_job = opts.Config_dflow_labels["benchmark-job"]
        current_url = "https://launching.mlops.dp.tech/?request=GET%3A%2Fapplications%2Fabacustest%2Fjobs%2F" + current_job
        comment_add = f"\n\nClick [here]({current_url}) to view detailed reports and historical trend graphs."
        send_to_feishu(outtable,opts.feishu_webhook,comment+comment_add)
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

'''


'''

