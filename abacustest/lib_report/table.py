import os, json, traceback, copy

def csv2table(csvfile):
    '''
    Transform a csv file to a table
    '''
    if not os.path.exists(csvfile):
        print(f"Error: {csvfile} does not exist!")
        return []
    with open(csvfile) as f:
        lines = f.readlines()
        table = []
        for line in lines:
            table.append(line.strip().split(","))
    return table

def json2table(jsonfile):
    '''
    Transform a json file to a table
    The key should be the first column that is the example name, and the value should be a dict, which are the other columns, and the keys are the metric name.
    
    '''
    if not os.path.exists(jsonfile):
        print(f"Error: {jsonfile} does not exist!")
        return []
    values = json.load(open(jsonfile))
    table = []
    metrics = []
    for k,v in values.items():
        for ik in v.keys():
            if ik not in metrics:
                metrics.append(ik)

    table.append(["example"] + metrics)
    for ikey in values.keys():
        table.append([ikey] + [values[ikey].get(imetric,None) for imetric in metrics])
    
    return table

def json2table_sm(jsonfile):
    '''
    Transform a json file to a table
    The key should be the first column that is supermetrics name, and the value is the value.
    '''
    if not os.path.exists(jsonfile):
        print(f"Error: {jsonfile} does not exist!")
        return []
    values = json.load(open(jsonfile))
    table = []
    for k,v in values.items():
        if isinstance(v,(int,float,str,bool,type(None))):
            table.append([k,v])
        else:
            print(f"Warning: type of the value of {k} is {type(v)}, which is not supported. Ignored.")
    
    return table

def file2table(metricfile):
    # based on the file type to read the file
    if not os.path.exists(metricfile):
        print(f"Error: {metricfile} does not exist!")
        return None
    
    filetype = os.path.splitext(metricfile)[1]
    if filetype == ".csv":
        try:
            table = csv2table(metricfile)
        except:
            traceback.print_exc()
            print(f"Error: transfer {metricfile} to table failed!")
            return None
    elif filetype == ".json":
        try:
            table = json2table(metricfile)
        except:
            traceback.print_exc()
            print(f"Error: transfer {metricfile} to table failed!")
            return None
    else:
        print(f"Error: file '{filetype}' of table is not supported!")
        return None
    return table

def rotate_table(table):
    '''
    Rotate a table
    '''
    new_table = []
    for i in range(len(table[0])):
        new_table.append([table[j][i] for j in range(len(table))])
    return new_table

def isort(itable_input,head_list):
    itable = copy.deepcopy(itable_input)
    heads = copy.deepcopy(itable[0])
    sort_idx = []
    for i in head_list:
        if i not in heads:
            print("isort() ERROR:",i,"is not a head of input table. Table head is:",heads)
            return itable
        else:
            sort_idx.append(heads.index(i))
    jtable = []
    for i in itable[1:]:
        jtable.append([i[j] for j in sort_idx] + [i])
    jtable.sort()
    return [heads] + [i[-1] for i in jtable]

def output_float(f, prec=2):
    '''
    Output a string of float with precision, default 2
    If f is a small number, output in scientific notation
    '''
    if f == None:
        return "---"
    elif isinstance(f, str):
        return f
    elif isinstance(f, int):
        return str(f)
    
    try:
        f = float(f)
    except:
        return f
    
    if abs(f) < pow(10,-1*prec):
        return '%.2e' % f
    else:
        return '%.*f' % (prec, f)
    
def judge_metric(x, criteria):
    '''
    Judge if a metric is good or bad based on criteria
    '''
    try:
        x = float(x)
        sm_pass = eval(criteria)
        return bool(sm_pass)
    except:
        return None
    
def format_table(table,metrics_name=None, sort=None, criteria=None, color={True:"green",False:"red"}):
    '''
    table: a list of list, each list is a row of the table. The first row is the head of the table
    metrics_name: a list of metrics name, which will be output in the table
    sort: a list of metrics name, which will be used to sort the table
    criteria: a dict of criteria, the key is the metric name, and the value is the criteria
    
    Do things:
    1. If the value is a number, then round it to 4 digits
    2. If the value pass the criteria, the value will be colored to green, else red.
    3. If the value is None, then transfer to "---"
    
    criteria is a dict:
    {
        "key": "key_name",
        "criteria": "x > 0"   # this is a string, should be evaluated by python
    }
    
    return a table, and a dict of pass number
    {
        "key_metrics": {"pass": 1, "notpass": 2},
        "all": {"pass": 3, "total": 4}
    }
    '''
    print("criteria:",criteria)
    print("table head:",table[0])
    
    if metrics_name in [[],None]:
        metrics_name = table[0]
    if table[0][0] not in metrics_name:
        metrics_name = [table[0][0]] + metrics_name # the first colume should be the example name
    
    metric_idx = []
    for i in metrics_name:
        if i not in table[0]:
            metric_idx.append(None)
        else:
            metric_idx.append(table[0].index(i))
    
    new_table = [metrics_name]
    for i in range(1,len(table)):
        new_table.append([None if j == None else table[i][j] for j in metric_idx])
    
    if sort:
        new_table = isort(new_table,sort)    
    table = new_table
        
    pass_num = {k:{"pass":0,"notpass":0} for k in criteria.keys()}
    pass_num["all"] = {"pass":0,"total":len(table)-1}
    for i in range(1,len(table)):
        allpass = True
        for j in range(len(table[i])):
            metric_name = table[0][j]
            if metric_name in criteria:
                metric_pass = judge_metric(table[i][j], criteria[metric_name])
                if metric_pass == True:
                    table[i][j] = '<font color="%s">%s</font>' % (color[True], output_float(table[i][j]))
                    pass_num[metric_name]["pass"] += 1
                else:
                    table[i][j] = '<font color="%s">%s</font>' % (color[False], output_float(table[i][j]))
                    if metric_pass == False:
                        pass_num[metric_name]["notpass"] += 1
                    allpass = False
            else:
                table[i][j] = output_float(table[i][j])
        if allpass:
            pass_num["all"]["pass"] += 1
                    
    return table, pass_num

def gen_criteria(criteria,pass_num):
    html =  f'''
    <table border="2px">
        <tbody>
            <tr>
                <td>Key</td>
                <td>Pass/Total</td>
                <td>Criteria</td>
            </tr>
        '''
    for key in criteria.keys():
        if key not in pass_num:
            continue
        icolor = "green" if pass_num[key]["pass"] == pass_num[key]["pass"]+pass_num[key]["notpass"] else "red"
        html += f'''<tr>
                <td>{key}</td>
                <td style="color:{icolor}">{pass_num[key]["pass"]}/{pass_num[key]["pass"]+pass_num[key]["notpass"]}</td> 
                <td>{criteria[key]}</td>
            </tr>
        '''
        
    html += '''</tbody>
    </table>\n\n'''
    return html

def gen_criteria_sm(criteria,sm):
    '''
    sm = {key: value}
    '''
    html =  f'''
    <table border="2px">
        <tbody>
            <tr>
                <td>Metric</td>
                <td>Value</td>
                <td>Criteria</td>
            </tr>
        '''
    for ik,iv in sm.items():
        if ik in criteria:
            icolor = "green" if judge_metric(iv, criteria[ik]) else "red"
            value = output_float(iv)
            html += f'''<tr>
                    <td>{ik}</td>
                    <td style="color:{icolor}">{value}</td> 
                    <td>{criteria[ik]}</td>
                </tr>
            '''
        else:
            html += f'''<tr>
                    <td>{ik}</td>
                    <td>{output_float(iv)}</td> 
                    <td>---</td>
                </tr>
            '''
        
    html += '''</tbody>
    </table>\n\n'''
    return html