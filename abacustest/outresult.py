import os,sys,argparse,glob,json,traceback
import numpy as np

def CheckFile(ifile):
    if os.path.isfile(ifile):
        return True
    else:
        print("WARNING: Can not find %s!!!" % ifile)
        return False

def TableOutput(datalist,maxlen=18,digitmax=3,digit=None,left=None,outframe=True,scintific=None,sep=" "):
    'datalist = [[col1,col2..],[col1,col2..]..]'
    context = ''
    ncolumn = len(datalist[0])
    collen = []
    for i in datalist[0]:
        collen.append(len(str(i)))
    for i in range(1,len(datalist)):
        if len(datalist[i]) != ncolumn:
            print("number of column in each line is not same, line %d" % i)
            print(datalist[i])
            return 0
        for j in range(len(datalist[i])):
            if isinstance(datalist[i][j],(float)):
                if scintific != None and scintific[j]:
                    fe = "e"
                else:
                    fe = "f"
                if digit != None :
                    if digit[j] >= 0:
                        lencol = len(("%."+str(digit[j])+fe)%datalist[i][j])
                    else:
                        lencol = len(str(datalist[i][j]))
                else:
                    lencol = len(("%."+str(digitmax)+fe)%datalist[i][j])
            else:
                lencol = len(str(datalist[i][j]))
            if lencol > collen[j] : collen[j] = lencol
    format1 = ''
    totallength = 0
    for i in range(len(collen)):
        if left != None and left[i]:
            lr = "%-"
        else:
            lr = "%"

        length = str(collen[i]) if collen[i] < maxlen else str(maxlen)
        format1 = format1 + lr + length + '.' + length + 's' + sep
        totallength += int(length)

    totallength += len(collen)-1
    if outframe:
        context += totallength*"-" + "\n"

    for i in range(len(datalist)):
        value = []
        for ij,j in enumerate(datalist[i]):
            if isinstance(j,(float,int)):
                if scintific != None and scintific[ij]:
                    fe = "e"
                else:
                    fe = "f"
                if digit != None:
                    if digit[ij] >= 0:
                        j = ("%."+str(digit[ij])+fe) % j
                    else:
                        j = str(j)
                else:
                    j = ("%."+str(digitmax)+fe) % j
            value.append(str(j) )
        context += format1 % tuple(value)
        context += '\n'

    if outframe:
        context += totallength*"-" + "\n"

    return context

def Table2FeishuInteractive(tablelist,webhook,title=None,comment=None,digitmax=3,digit=None,scintific=None):
    """
    translate a table to json used in feishu interactive
    """
    card = {}
    if title:
        #https://open.feishu.cn/document/ukTMukTMukTM/ukTNwUjL5UDM14SO1ATN
        card["header"] = {
            "title": {
                "tag":"plain_text",
                "content": title
            },
            "template":"blue"  #
        }

    #table
    card["elements"] = []
    
    #context
    #https://open.feishu.cn/document/ukTMukTMukTM/ucTNwUjL3UDM14yN1ATN/column-set
    for i,itable in enumerate(tablelist):
        context = []
        for ii,ivalue in enumerate(itable):
            if i > 0 and isinstance(ivalue,(int,float)):
                if digit:
                    idigit = digit[ii]
                else:
                    if isinstance(ivalue,float):
                        idigit = digitmax
                    else:
                        idigit = 0
                if idigit < 0: 
                    idigit = 0
                if scintific:
                    fe = "e" if scintific[ii] else "f"
                else:
                    fe = "f"
                format1 = "%%.%d%s" % (idigit,fe)
                if isinstance(ivalue,bool):
                    value = str(ivalue)
                else:
                    value = format1 % float(ivalue)
            else:
                value = str(ivalue)

            context.append({
                "tag": "column",
                "width": "weighted",
                "weight": 1,
                "elements": [
                {
                  "tag": "markdown",
                  "text_align": "center",
                  "content": value
                }]
            })
        card["elements"].append({
        "tag": "column_set",
        "flex_mode": "stretch",  #stretch, bisect, flow, trisect
        "background_style": "grey" if i == 0 else "default",
        "horizontal_spacing": "default",
        "columns": context})

    if comment:
        card["elements"].append({
        "tag": "markdown",
        "content": comment
        })

    feishu = {
        "msg_type": "interactive",
        "card" :card
    }

    import requests
    response = requests.post(webhook,
        headers={'Content-Type': 'application/json'},
        data=json.dumps(feishu))
    print(response.url,response.status_code,response.text,response.headers)

def GetParamValue(result,param,examplename):
    #result is a dict, whose key is the name of some parameters, and value is the value of the parameter
    #param is the param name, or a formula containing the param name which is identified by '
    #such as: param = "total_energy" or "'total_energy'/'natom'"

    def getvalue(iparam):
        if iparam in result:
            return result.get(iparam,None)
        else:
            if "/" in iparam:
                param_list = iparam.split("/")
                tmp_result = result
                for iiparam in param_list[:-1]:
                    if iiparam not in tmp_result:
                        return None
                    else:
                        tmp_result = tmp_result[iiparam]
                return tmp_result.get(param_list[-1],None)
            else:
                return None


    if "'" in param:
        formula = ""
        readparam = False
        paramname = ""
        tmp_dic = {}
        for i in param:
            if i == "'":
                if not readparam:
                    readparam = True
                    paramname = ""
                else:
                    value = getvalue(paramname) 
                    if value == None:
                        #print(paramname,type(paramname),result)
                        print("%s: value of '%s' is None, set %s to be None" % (examplename,paramname,param))
                        return None
                    
                    if isinstance(value,(int,float)):
                        formula += str(value)
                    elif isinstance(value,list):
                        key = 'np%d'%len(tmp_dic)
                        tmp_dic[key] = np.array(value)
                        formula += "tmp_dic['%s']" % key
                    else:
                        print("%s: value of %s is " % paramname, value,", whose type is %s, is not int or float, set %s to be None" % (examplename,type(value),param))
                        return None
                    readparam = False
            elif readparam:
                paramname += i
            else:
                formula += i
        try:
            #print(formula)
            #print(eval(formula))
            return eval(formula)
        except:
            traceback.print_exc()
            return None       
    else:
        return getvalue(param)

def ProduceExpandDic(allresult,example_idx):
    #construct the expand dict, the key is the key in "outparams_expand"
    #the value is a dict, whose key is "0", "1", .. "ntype"
    # and whose value is the value defined in outparams.
    expand_dic = {}
    outparams = allresult.get("outparams",[])
    outparams_paramname = [i[0] for i in outparams]
    example_name = allresult.get("example_name",[])
    results = allresult.get("result",[[]])
    type_idx = allresult.get("type_idx",[])
    #print(outparams_paramname)
    for key in allresult.get("outparams_expand",{}):
        if key not in outparams_paramname:
            print("key '%s' defined in 'outparams_expand' is not defined in 'outparams', skip it" % key)
            continue
        expand_dic[key] = {}
        for i in range(len(results)):
            result = GetParamValue(results[i][example_idx],
                                   outparams[outparams_paramname.index(key)][1][type_idx[i]],
                                   example_name[example_idx])
            expand_dic[key]['%d'%i] = result
    #print(expand_dic)
    return expand_dic

def OutParam(allresult,split_example="----"):
    '''
    allresult = {
        "result": [[type1_result1,type1_result2,...],[type2_result1,type2_result2,..],...],
        "type_name": [str],
        "example_name": ["example1","example2",...],
        "outparams": [[param1_name_in_table, [param1_name_in_typeAresult,param_name1_in_typeBresult,...], significant_digit_number, False],
                      [param2_name_in_table, [param2_name_in_typeAresult,param_name2_in_typeBresult,...], significant_digit_number, False],
                      ...],
        "type_idx": [],
        "outparams_expand": {str:[]},
        "outparams_comment": ""
    }
    
    type_name: the name of each job type.
    type_result: is a dictionary, whose key is param name.
    type_idx: ntype elemens. Each element is the index outparams[][1], which indicate the param name in this job type.
    outparams_expand: a dictionary. The key is outparams[][0], and value is a list of each job type. You can define the calculation between different job types.
                        Such as: to calculate the enery difference between job type 0 and 1, and show in line of type0, you can write as {"energy": ["'0' - '1'",""]}
    
    '''
    outtable = [['example',"jobType"]]
    param = allresult.get("outparams")
    for i in param:
        if i[-1] != None: outtable[0].append(i[0])
           
    ncol = len(outtable[0]) 
    outtable.append(ncol * ["----"])
    
    results = allresult.get("result",[[]])
    example_name = allresult.get("example_name",[])
    type_idx = allresult.get("type_idx",[])
    allparam_value = [[] for i in range(len(results))] #[[type1_example1_dict, type1_example2_dict],[type2_...]]
    for i in range(len(results[0])):  #i for example
        expand_dic = ProduceExpandDic(allresult,i)
        for j in range(len(results)):  # j for job type
            iexample = [example_name[i]] if j==0 else [" "]   #add example name
            iexample.append(allresult["type_name"][j])
            allparam_value[j].append({})
            for iparam in param:
                if iparam[0] not in expand_dic:
                    value = GetParamValue(results[j][i],iparam[1][type_idx[j]],example_name[i])
                else:
                    expand_formula = allresult['outparams_expand'][iparam[0]]
                    if j < len(expand_formula) and expand_formula[j].strip() != "":
                        value = GetParamValue(expand_dic[iparam[0]],expand_formula[j],example_name[i]+"/"+iparam[0])
                    else:
                        value = " "
                if iparam[-1] != None:
                    iexample.append(value) 
                #print(allparam_value[j][-1])
                allparam_value[j][-1][iparam[0]] = value if value != " " else None 
            outtable.append(iexample)
        if split_example != None: outtable.append(ncol * [split_example])
    if split_example != None: outtable = outtable[:-1]
    
    
    digit = [-1,-1] + [i[2] for i in param]
    left = [True,True]
    for i in param:
        if len(i) == 4:
            left.append(i[3])
        else:
            left.append(False)
            
    cc = '\nSome key results of ' + ", ".join(allresult.get("type_name")) + "\n"
    comment = allresult.get("outparams_comment",[])
    if len(comment) > 0:
        cc += "\n".join(comment) + "\n"
    cc += TableOutput(outtable,maxlen=18,digit=digit,left=left)
    if allresult.get("webhook"):
        Table2FeishuInteractive(outtable,webhook=allresult.get("webhook"),title="abacustest",digit=digit)
    return cc,allparam_value

class MetricsMethod:
    @classmethod
    def allmethod(cls):
        a = cls()
        method = {}
        for k,v in a.method_list().items():
            method[k] = v[0]
        return method
    
    @classmethod
    def allmethod_str(cls):
        a = cls()
        method = ''
        for k,v in a.method_list().items():
            method += "%10s: %s\n" % (k,v[1])
        return method
    
    def method_list(self):
        return {"GM":[self.GM,"calculate the geometric mean value"],
                "iGM":[self.iGM, "calculate the geometric mean of the inverse of value"],
                "MEAN":[self.MEAN,"mean of all values"],
                "TrueRatio":[self.TrueRatio,"the ratio of True"]
                }
    def MEAN(self,valuelist):
        return np.array(valuelist).mean()
        
    def CalGM(self,valuelist,func=lambda x:x):
        n = len(valuelist)
        gm = 1
        for i in valuelist:
            gm *= func(i)
        return (gm)**(1/n)
    
    def GM(self,valuelist):
        return self.CalGM(valuelist)
    
    def iGM(self,valuelist):
        return self.CalGM(valuelist,lambda x: 1/x)
    
    def TrueRatio(self,valuelist):
        nTrue = 0
        for i in valuelist:
            if i == True:
                nTrue += 1
        
        return float(nTrue)/float(len(valuelist))
        
def OutMetrics(allresults,allparam_value):
    type_name = allresults.get("type_name")
    example_name = allresults.get("example_name",[])
    results = allresults.get("result",[[]])
    
    outtable = [["metrics"] + type_name + ['example_number'],(len(type_name)+2)*["----"]]
    metrics = allresults.get("metrics",[])

    method_list = MetricsMethod.allmethod()
    hasnotsupportmethod = False
    notsupportmethod = []
    comment = []
    allmetric_value = {}
    for metric in metrics:
        name = metric.get("name")
        paramname = metric.get("param_name") 
        method = metric.get("method")
        condition = metric.get("condition","").strip() #if condition is not meet, value will be set to None
        doclean = metric.get("doclean",False)  #if remove the value of None
        normalization = metric.get("normalization",True)
        if "comment" in metric:
            icomment = metric["comment"].strip()
            if icomment != "":
                comment.append("%s:%s\n" % (name,icomment))
        #if paramname not in outparams_paramname:
        #    print("key '%s' defined in 'metrics/param_name' is not defined in 'outparams', skip it" % paramname) 
        #    continue
        if  method not in method_list:
            if method not in notsupportmethod: notsupportmethod.append(method)
            hasnotsupportmethod = True
            continue
        tmp_result = [[] for i in range(len(type_name))]
        for i in range(len(results[0])):
            tmp = []
            hasnone = False
            for j in range(len(results)):
                tmp_name = type_name[j] + "/" + example_name[i]
                if condition == "":
                    value = GetParamValue(allparam_value[j][i],paramname,tmp_name)
                else:
                    conditionvalue = GetParamValue(allparam_value[j][i],condition,tmp_name)
                    if conditionvalue:
                        value = GetParamValue(allparam_value[j][i],paramname,tmp_name)
                    else:
                        value = None
                if value == None:
                    hasnone = True
                tmp.append(value)
            if doclean and hasnone:
                continue
            else:
                for j,itmp in enumerate(tmp): 
                    if itmp!=None: tmp_result[j].append(itmp)
        if doclean and len(tmp_result[0]) == 0:
            print("%s: values of all examples are None for some job type, skipt it!" % name)
            continue
        
        nresult = [str(len(i)) for i in tmp_result]
        nexample = nresult[0] if doclean else "/".join(nresult)
               
        outtable.append([name])
        metric_value = []
        for i in range(len(tmp_result)):
            try:
                metric_value.append(method_list[method](tmp_result[i]))
            except:
                traceback.print_exc()
                metric_value.append(None)

            if normalization:
                try:
                    outtable[-1].append(metric_value[i]/metric_value[0]) #do normalization by divide by the first value
                except:
                    traceback.print_exc()
                    outtable[-1].append(None)
            else:
                outtable[-1].append(metric_value[i])        

        if normalization and metric_value[0] != None:
            comment += "%-20s is normalized, value of %s is %.3e.\n" % (name,type_name[0],metric_value[0])

        outtable[-1] += [nexample]
        allmetric_value[name] = metric_value
        
    if hasnotsupportmethod:
        print("Method %s are not supportted now.\nSupported methods are:\n%s" % (" ".join(notsupportmethod),MetricsMethod.allmethod_str()))

    digit = [-1] + len(type_name)*[3] + [0]
    left = [True] + (len(type_name)+1)*[False]
    scintific = (len(type_name)+2)*[False]
    
    cc = "\nSome key metrics\n"
    cc += TableOutput(outtable,maxlen=50,digit=digit,left=left,scintific=scintific)
    if allresults.get("webhook"):
        Table2FeishuInteractive(outtable[0:1] + outtable[2:],webhook=allresults.get("webhook"),title="abacustest",digit=digit)
    cc += "".join(comment)
    #cc += "Notice: the exmaples with value of None for some job type will be excluded.\n"
    return cc,allmetric_value    

def GetAllResults(result_setting):
    #allresults = {"type_name":[types],"result_file": [result1,result2,result3,....] }
    #result is a two dimension list, 1st dimension is for result of each example in each type,
    #second dimension is types
    #len("example_name") = len("result"[i]) = len("result"[j])
    #len("type_name") = len("result")
    #result123 = [example1_result, example2_result, example3_result]

    resultf = result_setting["result_file"]
    type_name = result_setting.get("type_name",["jobtype"])
    
    #deal with type_name
    if isinstance(type_name,str):
        if len(resultf) > 1:
            type_name = [type_name+str(i) for i in range(len(resultf))]
        else:
            type_name = [type_name]
    elif isinstance(type_name,list) and len(type_name) == 1 and len(resultf) > 1:
            type_name = [type_name[0]+str(i) for i in range(len(resultf))]
    if not isinstance(type_name,list):
        print("type_name should be a list of type names")
        return False
    if (len(resultf) != len(type_name)):     
        print("number of result_fiel is not equal to type_name")
        return False
    
    #deal with example_name_idx
    example_name_idx = result_setting.get("example_name_idx",-1) 
    if isinstance(example_name_idx,int):
        example_name_idx = len(type_name) * [example_name_idx]
    elif isinstance(example_name_idx,list) and len(example_name_idx) == 1 and len(type_name) > 1:
        example_name_idx = len(type_name) * example_name_idx
    if not isinstance(example_name_idx,list):
        print("example_name_idx should be a list of int")
        return False        
    if (len(example_name_idx) != len(type_name)):
        print("number of example_name_idx is not equal to type_name")
        return False
    
    #deal with type_idx
    type_idx = result_setting.get("type_idx",[0 for i in range(len(type_name))])
    if isinstance(type_idx,int):
        type_idx = len(type_name) * [type_idx]
    elif isinstance(type_idx,list) and len(type_idx) == 1 and len(type_name) > 1:
        type_idx = len(type_name) * type_idx
    if not isinstance(type_idx,list):
        print("type_idx should be a list of int")
        return False 
    if (len(type_idx) != len(type_name)):
        print("number of type_idx is not equal to type_name")
        return False
    
    #deal with results
    result = []
    example_name = []
    for i,iresultf in enumerate(resultf):
        results = json.load(open(iresultf))
        result.append([])
        if i==0:
            for k,v in results.items():
                example_name.append(k.split("/")[example_name_idx[i]])
                result[0].append(v)
        else:
            tmp_dict = {}
            for k,v in results.items():
                tmp_dict[k.split("/")[example_name_idx[i]]] = v
            for iexample in example_name:
                result[-1].append(tmp_dict.get(iexample,{}))         
    
    #print("example number: %d" % (len(example_name)))
    #print("type number: %d" % (len(type_name)))
    #print("type idx: %d" % (len(type_idx)))

    hasfalse = False
    if len(result) != len(type_name):
        print("type number is %d, but the type of result is %d" % (len(type_name),len(result)))
        hasfalse = True
        
    if len(type_idx) != len(type_name):
        print("type idx is %d, but the type number is %d" % (len(type_idx),len(type_name)))
        hasfalse = True
        
    for i,ir in enumerate(result):
        if len(ir) != len(example_name):
            print("type %d has %d results, which is not equal to example number %d" % (i+1,len(ir),len(example_name)))
            hasfalse = True

    if hasfalse:
        return False

    outparams = result_setting.get("outparams",[])
    metrics = result_setting.get("metrics",[])
    outparams_expand = result_setting.get("outparams_expand",{})
    outparams_comment = result_setting.get("outparams_comment",{})
    plot = result_setting.get("plot",[])
    webhook = result_setting.get("webhook",None)
    
    return {"example_name": example_name, 
            "type_name": type_name, 
            "result": result,
            "type_idx": type_idx, 
            "outparams": outparams,
            "outparams_expand":outparams_expand,
            "outparams_comment":outparams_comment,
            "metrics": metrics,
            "plot":plot,
            "webhook":webhook}
    
def check_file(ifile):
    # check if ifile exist,
    # if exist, then modify the file name to add a number at the end until the file does not exist
    # return the new file name
    ifile = os.path.abspath(ifile)
    if os.path.isfile(ifile):
        ifile_prefix = os.path.splitext(ifile)[0]
        ifile_suffix = os.path.splitext(ifile)[1]
        i = 1
        while os.path.isfile(ifile):
            ifile = ifile_prefix + "_%d" % i + ifile_suffix
            i += 1
    return ifile

def pandas_out(allresult,savefile = None, report_sample_max = 2,print_result=True, print_list_seperate = False, float_prec=None):
    """
    allresult = {sample1: {key1:value,key2:value},
                 sample2: {key1:value,key2:value}}
    If value is a list, will print out separately.
    """
    import pandas as pd
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 200)
    if float_prec != None:
        s = "{:." + str(float_prec) + "f}"
        pd.set_option('display.float_format', s.format)

    normal_result = {}
    list_result = []
    allsamples = [i for i in allresult.keys()]
    allkeys = []
    for v in allresult.values():
        for iv in v.keys():
            if iv not in allkeys:
                allkeys.append(iv)
        
    allkeys_seperate = []
    savefile_prefix = "" if not savefile else os.path.splitext(savefile)[0]
    savefile_names = []
    for ikey in allkeys:
        allkeys_seperate.append(False)
        if not print_list_seperate: # if not print list separately, then print out the list together
            continue
        for i in allsamples:
            if isinstance(allresult[i].get(ikey,None),dict):
                allkeys_seperate[-1] = True
                break
            elif isinstance(allresult[i].get(ikey,None),list):
                #if len(allresult[i][ikey]) > 0 and isinstance(allresult[i][ikey][0],(list,dict)):
                    allkeys_seperate[-1] = True
                    break
    nsample = 0    
    for isample in allsamples:
        nsample += 1
        normal_result[isample] = {}
        for i,ikey in enumerate(allkeys):
            sfname = savefile_prefix +"_" + isample.replace("/","_")+"_"+ikey.replace("/","_")
            if allkeys_seperate[i]:
                if isinstance(allresult[isample].get(ikey,None),(list,dict)):
                    try:
                        list_result.append("%s\t%s:\n" % (isample,ikey) + str(pd.DataFrame.from_dict(allresult[isample][ikey])))
                        if savefile:
                            fname_final = check_file(sfname+".csv")
                            pd.DataFrame.from_dict(allresult[isample][ikey]).to_csv(fname_final)
                            if nsample <= report_sample_max:
                                savefile_names.append(fname_final)
                            
                    except:
                        list_result.append("%s\t%s:\n" % (isample,ikey) + str(allresult[isample][ikey]))
                        if savefile:
                            fname_final = check_file(sfname+".txt")
                            with open(fname_final,'w') as f:
                                f.write(str(allresult[isample][ikey]))
                                if nsample <= report_sample_max:
                                    savefile_names.append(fname_final)
                            
                else:
                    if allresult[isample].get(ikey,None) != None:
                        list_result.append("%s\t%s:\n" % (isample,ikey) + str(allresult[isample].get(ikey,None)))
                        if savefile:
                            fname_final = check_file(sfname+".txt")
                            with open(fname_final,'w') as f:
                                f.write(str(allresult[isample][ikey]))
                                if nsample <= report_sample_max:
                                    savefile_names.append(fname_final)
            else:
                normal_result[isample][ikey] = allresult[isample].get(ikey,None)
                
    if False not in allkeys_seperate:
        normal_result = ""
    else:
        pddata = pd.DataFrame.from_dict(normal_result,orient='index')
        if savefile:
            savefile_name_final = check_file(savefile)
            pddata.to_csv(savefile_name_final, index_label="example")
            savefile_names.insert(0,savefile_name_final)
        normal_result = str(pddata)
        
    if True not in allkeys_seperate:
        list_result = ""
    else:
        list_result = "\n\n".join(list_result)
    
    if print_result:
        print("%s\n\n%s" % (normal_result,list_result))      
    return normal_result,list_result,savefile_names

def OutResultArgs(parser):  
    parser.description = "This script is used to output the summary of results"
    parser.add_argument('-r', '--result', type=str, help='the result file from collectdata, should be .json type',action="extend",nargs="*")
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type')
    parser.add_argument('-m', '--metrics', default=["ALL"], help='The metrics that needed to be shown. Use ALL to show all metrics, default is ALL.', action="extend",nargs="*")
    parser.add_argument('-o', '--output', type=str, help='output the selected metrics to a json file')
    parser.add_argument('--csv', type=str, help='if output the selected metrics to csv file')
    parser.add_argument('--prec', type=int, default=None, help='the precision of float number in pandas output, default is None')
    parser.add_argument("--list_sep", default=0, const=1, nargs='?', type=int, help="if print the list separately, default is 0")
    
    return parser

def outresult(param):
    if param.result != None:
        allresult_files = param.result #if len(param.result) == 1 else param.result[1:]
        metric_name = param.metrics if len(param.metrics) == 1 else param.metrics[1:]
        allresult = {}
        allfiles = []
        for ifile in allresult_files:
            for iif in glob.glob(ifile):
                allfiles.append(iif) 
        for ifile in allfiles:
            result = json.load(open(ifile))
            for k,v in result.items():
                key = "%s:%s" % (ifile,k) if len(allfiles) > 1 else k
                allresult[key] = v
        if "ALL" not in metric_name:
            need_metric = []
            for i in metric_name:
                if i not in need_metric:
                    need_metric.append(i)
            new_result = {}
            for ik,iv in allresult.items():
                new_result[ik] = {}
                for i in need_metric:
                    new_result[ik][i] = iv.get(i,None)
            allresult = new_result
        if param.output != None:
            json.dump(allresult,open(param.output,'w'),indent=4)
        pandas_out(allresult,param.csv,print_list_seperate=param.list_sep,float_prec=param.prec)
    
    if param.param!= None:
        if not CheckFile(param.param):
            sys.exit(1)
        
        report_part = json.load(open(param.param))
        if "allresults" in report_part:
            report_part = report_part["allresults"]
        #elif "report" in report_part:
        #    report_part = report_part["report"]
        allresults = GetAllResults(report_part)
        if len(allresults['type_name']) == 1: split_example = None
        else: split_example = "----"
        cc_outparam,allparam_value = OutParam(allresults,split_example=split_example)
        cc_outmetrics,allmetric_value = OutMetrics(allresults,allparam_value)
        
        output_value = {}
        output_value['type_name'] = allresults['type_name']
        output_value['example_name'] =  allresults['example_name']
        output_value['param_value'] = allparam_value
        output_value['metrics_value'] = allmetric_value
        
        json.dump(allmetric_value,open("superMetrics.json","w"),indent=4)
        if param.output != None:
            json.dump(output_value,open(param.output,'w'),indent=4)
        if allmetric_value:
            print(cc_outmetrics)
        print(cc_outparam)

def main():
    parser = argparse.ArgumentParser()
    param = OutResultArgs(parser).parse_args()
    outresult(param)
    
if __name__ == "__main__":
    main()