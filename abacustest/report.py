import os,sys,argparse,json,copy
import traceback
from abacustest.lib_report import table as tb
from abacustest.myflow import globV

HTML_HEAD = """
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {
            text-align: left;
        }
        
        .head1{
            font-size: 20px;
            font-weight: bold;
            white-space: pre-wrap;
            line-height: 2; 
            margin-top: 1rem;   
        }
        
        .head2{
            font-size: 18px;
            font-weight: bold;
            white-space: pre-wrap;
            line-height: 2;
        }
        
        .head3{
            font-size: 16px;
            font-weight: bold;
            white-space: pre-wrap;
            line-height: 2;
        }

        .tabletitle{
            font-size: 16px;
            font-weight: bold;
            word-wrap: break-word;
            line-height: 2;
        }
        
        .imagetitle{
            font-size: 16px;
            font-weight: bold;
            word-wrap: break-word;
            line-height: 2;
        }

        .doc {
            font-family: Verdana, sans-serif;
            text-align: left;
            display: inline-block;
            font-size: 16px;
            width: 100%;
            word-wrap: break-word;
            white-space: pre-wrap;
            line-height: 1.5;
            margin-bottom: 0.5rem;
        }
        
        #keys table {
          border-collapse: collapse;
          width: 100%;
        }
        #keys td {
          border: none;
          padding: 5px;
          text-align: left;
          line-height: 0.8;
        }
        
        img {
            max-width: 600px;
            max-height: 600px;
            height: auto;
            cursor: zoom-in;
        }

        .overlay {
            display: none;
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background-color: rgba(0, 0, 0, 0.8);
            z-index: 1000;
            overflow: auto;
        }

        .overlay img {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            max-width: 100%;
            max-height: 100%;
            object-fit: contain;
            cursor: zoom-out;
        }
    </style>
</head>
"""

def _table2html(table,has_head=True):    
    # add title
    html = '''\t<table border="2px">\n'''

    # add table head
    start_row = 0
    if has_head:
        html += '\t\t<thead><tr>'
        for i in range(len(table[0])):
            html += '<th>%s</th>' % table[0][i]
        html += '</tr></thead>\n'
        start_row = 1
    
    # add table body
    html += '\t\t<tbody>'
    for i in range(start_row,len(table)):
        html += '\t\t\t<tr>'
        for j in range(len(table[i])):
            html += '<td>%s</td>' % table[i][j]
        html += '</tr>\n'
    html += '\t\t</tbody>\n'
    html += '\t</table>\n'
    
    return html
    

def metrics2html(metrics_set):
    metric_file = metrics_set.get("content","")
    criteria = metrics_set.get("criteria",{})
    title = metrics_set.get("title","")
    metrics = metrics_set.get("metrics",[])
    sort = metrics_set.get("sort",[])
    center = metrics_set.get("center",True)
    
    if not os.path.exists(metric_file):
        print(f"Error: {metric_file} does not exist!")
        return ""
    
    table = tb.file2table(metric_file)
    if not table:
        print(f"Error: transfer {metric_file} to table failed!")
        return ""
    if table in [[],None]:
        return ""
    if not sort: sort = [table[0][0]]
    table, pass_num = tb.format_table(table, metrics, sort, criteria)
    
    html = ""
    if title:
        html += f'''\t<div class="tabletitle">{title}</div>\n'''
        
    if criteria:
        passnum = pass_num["all"]["pass"]
        totalnum = pass_num["all"]["total"]
        icolor = "green" if passnum == totalnum else "red"
        # if all passed, then color the number to green, else red
        html += f'''\t<div class="head2">Pass/Total: <font color="{icolor}">{passnum}/{totalnum} ({passnum/totalnum*100:.2f}%)</font></div>\n'''
        html += tb.gen_criteria(criteria,pass_num)
    html += _table2html(table,has_head=True)
    
    if center:
        html = "\t<center>\n" + html + "\t</center>\n"
    return html

def supermetrics2html(supermetrics_set):
    metric_file = supermetrics_set.get("content","")
    criteria = supermetrics_set.get("criteria",{})
    title = supermetrics_set.get("title","")
    center = supermetrics_set.get("center",True)
    
    if not os.path.exists(metric_file):
        print(f"Error: {metric_file} does not exist!")
        return ""
    
    sm = json.load(open(metric_file))
    if not sm:
        print("Error: read supermetrics failed!")
        return ""
    
    html = ""
    if title:
        html += f'''\t<div class="tabletitle">{title}</div>\n'''
    html += tb.gen_criteria_sm(criteria,sm)
    if center:
        html = "\t<center>\n" + html + "\t</center>\n"
    return html

def table2html(table_set):
    filename = table_set.get("content","")
    title = table_set.get("title","")
    center = table_set.get("center",True)
    
    if not os.path.exists(filename):
        print(f"Error: {filename} does not exist!")
        return ""
    
    html = ""
    if title:
        html += f'''\t<div class="tabletitle\">{title}</div>\n'''
        
    filetype = os.path.splitext(filename)[1]
    if filetype == ".csv":
        table = tb.csv2table(filename)
    else:
        print(f"Error: file type '{filetype}' of table is not supported!")
        return ""

    html += _table2html(table,has_head=True)
    if center:
        html = "\t<center>\n" + html + "\t</center>\n"
    return html

def image2html(image_set):
    image_file = image_set.get("content","")
    title = image_set.get("title","")
    center = image_set.get("center",True)
    
    if isinstance(image_file,str):
        image_file = [image_file]
        
    html = ""

    for ifile in image_file:
        if not os.path.exists(ifile):
            print(f"Error: image '{ifile}' does not exist!")
        html += f"""\t<img class="thumbnail" src="{ifile}" onclick="openFullscreen(this)">\n"""
    
    if title != "":
        html += f'''\t<div class="imagetitle">{title}</div>\n'''
    if center:
        html = "\t<center>\n" + html + "\t</center>\n"
    return html
    
def text2html(text_set):
    '''
    Transform a text to html format
    '''
    text = text_set.get("content","")
    if isinstance(text, str):
        text = text.split("\n")
    elif isinstance(text, list):
        pass
    else:
        print(f"Error: the content of text should be str or list, but {type(text)} is given!")
        return ""
    html = ""
    for ic in text:    
        html += f'''\t<div class="doc">    {ic}</div>\n'''
    return html   

def get_job_address():
    job = globV.get_value("JOB_ADDRESS","")
    if job != "":
        job = f'''<a href="{job}">link</a>'''
    return job

def get_version():
    if os.path.isfile("version.dat"):
        with open("version.dat") as f:
            version = f.read().strip()
    else:
        version = ""
    return version
    
def get_test_date():
    from datetime import datetime
    run_date = globV.get_value("START_TIME",datetime.now()).strftime("%Y-%m-%d")
    return run_date
         
def keys2html(keys):
    html = ""
    html += """\t<table id="keys">\n"""
    comm_keys = ["test_date","version"]
    #if not keys.get("job_address",""):
    #    keys["job_address"] = get_job_address()
        
    if not keys.get("version",""):
        keys["version"] = get_version()
    if not keys.get("test_date",""):
        keys["test_date"] = get_test_date()
        
    # write the comman keys
    for ikey in comm_keys:
        html += f'''\t\t<tr><td><strong>{ikey.capitalize()}</strong></td><td>:</td><td>{keys.get(ikey,"")}</td></tr>\n'''
    # write the other keys
    for ikey in keys.keys():
        if ikey not in comm_keys:
            html += f'''\t\t<tr><td><strong>{ikey.capitalize()}</strong></td><td>:</td><td>{keys.get(ikey,"")}</td></tr>\n'''
    html += """\t</table>\n"""
    return html

def gen_script(has_image=False):
    if not has_image:
        return ""
    
    html = ""
    
    if has_image:
        # add script for image zoom
        html += '''
    <div class="overlay" id="overlay" onclick="closeFullscreen()">
        <img id="fullscreenImage" src="">
    </div>
    '''
    
    html += '''\t<script>\n'''
    
    if has_image:
        html += '''
        var overlay = document.getElementById("overlay");
        var fullscreenImage = document.getElementById("fullscreenImage");
        
        function openFullscreen(imgElement) {
            fullscreenImage.src = imgElement.src; 
            overlay.style.display = "block";
        }
        
        function closeFullscreen() {
            overlay.style.display = "none";
            fullscreenImage.src = ""
        }
        '''
        
    html += '''
    </script>
    '''
    
    return html

def gen_html(report_setting, output):
    '''
    A report should be like this:
    Test Date/Version/Targets/Datasets/Properties/Criteria/Job Address:
    
    Content Section:
    '''
    
    '''
    Keys shold be a dict:
    {
        "test_date": "2020-01-01",
        "version": "1.0",
        "targets": "target1",
        "datasets": "dataset1",
        "properties": "property1",
        "criteria": "criteria1",
        "job_address": "job_address1"
    }
    
    The content section should has the following format:
    head1:
    head2:
    head3:
    text:
    image:
    table:
    metrics:
    supermetrics:
    
    If the content is a table, it can have an extra criteria to color the value, and if the value pass the criteria, the value will be colored to green, else red.
    If the content is a table or image, it can have an extra title.
    If the content is text, it can be a list of str, and each str will be a paragraph.
    
    Should be a list of dict:
    {
        "type": "head1",
        "content": "This is a head1"
    },
    {
        "type": "text",
        "content": "This is a text" or ["This is a text1","This is a text2"]
    }
    {
        "type": "image",
        "content": "image1.png",
        "title": "This is a image"
    },
    {
        "type": "table",
        "title": "This is a table",
        "content": "table1.csv"
    },
    {
        "type": "metrics",  # the row of a metrics table should be the example name, and the column should be the metric name
        "title": "This is a metrics table",
        "content": "table1.csv"
        "criteria":{
            "key1": "x > 0",
            "key2": "x < 1", # key is the name of the column, and the criteria is a string, should be evaluated by python
        },
        "sort": ["key1","key2"] # sort the table based on the key, default is the first column
        "metrics": ["metric1","metric2",] # only show the metrics in the list, default is all
    },
    {
        "type": "supermetrics",  # the row of a metrics table should be the example name, and the column should be the metric name
        "title": "This is a metrics table",
        "content": "table1.csv"
        "criteria":{
            "key1": "x > 0",
            "key2": "x < 1", # key is the name of the column, and the criteria is a string, should be evaluated by python
        }
    }
    '''
    
    keys = report_setting.get("keys",{})
    content = report_setting.get("content",[])
    
    html = HTML_HEAD + "\n<body>\n"
    html += keys2html(keys) + "\n"

    # write content
    has_image = False
    for item in content:
        itype = item.get("type","text")
        icontent = item.get("content","")
        if itype in ["head1","head2","head3"]:
            html += f'''\t<div class="{itype}">{icontent}</div>\n'''
        elif itype == "text":
            html += text2html(item)
        elif itype == "image":
            html += image2html(item)
            has_image = True
        elif itype == "table":
            html += table2html(item)
        elif itype == "metrics":
            html += metrics2html(item)
        elif itype == "supermetrics":
            html += supermetrics2html(item)
        else:
            continue
        html += "\n"
    
    # add script for image zoom
    html += gen_script(has_image)
         
    html += """\n</body>\n</html>"""
    
    with open(output,"w") as f:
        f.write(html)
    return html

def ReportArgs(parser):  
    parser.description = "Read metrics.json and generate the report"
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type',required=True)
    parser.add_argument('-o', '--output', type=str,  default="abacustest.html",help='The output file name, default is abacustest.html')
    return parser

def Report(param):
    globV._init()
    param_file = param.param
    if not os.path.exists(param_file):
        print(f"Error: {param_file} does not exist!")
        sys.exit(1)
    report_setting = json.load(open(param_file)).get("report",{})

    if report_setting == {}:
        print("Error: report section is empty!")
        sys.exit(1)
    output = param.output
    gen_html(report_setting,output)
    

def main():
    parser = argparse.ArgumentParser()
    param = ReportArgs(parser).parse_args()
    Report(param)
    
if __name__ == "__main__":
    main()