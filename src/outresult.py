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

def GetParamValue(result,param,examplename):
    #result is a dict, whose key is the name of some parameters, and value is the value of the parameter
    #param is the param name, or a formula containing the param name which is identified by '
    #such as: param = "total_energy" or "'total_energy'/'natom'"
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
                    value = result.get(paramname,None)
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
        return result.get(param,None)

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
        "outparams_command": ""
    }
    
    type_name: the name of each job type.
    type_result: is a dictionary, whose key is param name.
    type_idx: ntype elemens. Each element is the index outparams[][1], which indicate the param name in this job type.
    outparams_expand: a dictionary. The key is outparams[][0], and value is a list of each job type. You can define the calculation between different job types.
                        Such as: to calculate the enery difference between job type 0 and 1, and show in line of type0, you can write as {"energy": ["'0' - '1'",""]}
    
    '''
    outtable = [['example']]
    param = allresult.get("outparams")
    for i in param:
        if i[-1] != None: outtable[0].append(i[0])
           
    ncol = len(outtable[0])
    if split_example != None: outtable.append(ncol * [split_example])
    
    results = allresult.get("result",[[]])
    example_name = allresult.get("example_name",[])
    type_idx = allresult.get("type_idx",[])
    allparam_value = [[] for i in range(len(results))] #[[type1_example1_dict, type1_example2_dict],[type2_...]]
    for i in range(len(results[0])):  #i for example
        expand_dic = ProduceExpandDic(allresult,i)
        for j in range(len(results)):  # j for job type
            iexample = [example_name[i]] if j==0 else [" "]   #add example name
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
    
    digit = [-1] + [i[2] for i in param]
    left = [True]
    for i in param:
        if len(i) == 4:
            left.append(i[3])
        else:
            left.append(False)
            
    cc = '\nSome key results of ' + ", ".join(allresult.get("type_name")) + "\n"
    command = allresult.get("outparams_command",[])
    if len(command) > 0:
        cc += "\n".join(command) + "\n"
    cc += TableOutput(outtable[:-1],maxlen=18,digit=digit,left=left)
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
                "MEAN":[self.MEAN,"mean of all values"]
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
        
def OutMetrics(allresults,allparam_value):
    type_name = allresults.get("type_name")
    example_name = allresults.get("example_name",[])
    results = allresults.get("result",[[]])
    
    outtable = [["metrics"] + type_name + ['example_number'] + [type_name[0]],(len(type_name)+3)*["----"]]
    metrics = allresults.get("metrics")

    method_list = MetricsMethod.allmethod()
    hasnotsupportmethod = False
    notsupportmethod = []
    command = []
    allmetric_value = {}
    for metric in metrics:
        name = metric.get("name")
        paramname = metric.get("param_name") 
        method = metric.get("method")
        condition = metric.get("condition","").strip() 
        doclean = metric.get("doclean",True)
        if "command" in metric:
            icommand = metric["command"].strip()
            if icommand != "":
                command.append("%s:%s\n" % (name,icommand))
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
            print("%s: All examples with values of None for some job type, skipt it!" % name)
            continue
        
        nresult = [str(len(i)) for i in tmp_result]
        nexample = nresult[0] if doclean else "/".join(nresult)
               
        outtable.append([name])
        metric_value = []
        for i in range(len(tmp_result)):
            metric_value.append(method_list[method](tmp_result[i]))
            if doclean:
                outtable[-1].append(metric_value[i]/metric_value[0]) #do normalization by divide by the first value
            else:
                outtable[-1].append(metric_value[i])
        outtable[-1] += [nexample] + [metric_value[0]]
        allmetric_value[name] = metric_value
        
    if hasnotsupportmethod:
        print("Method %s are not supportted now.\nSupported methods are:\n%s" % (" ".join(notsupportmethod),MetricsMethod.allmethod_str()))

    digit = [-1] + len(type_name)*[3] + [0] + [3]
    left = [True] + (len(type_name)+2)*[False]
    scintific = (len(type_name)+2)*[False] + [True]
    
    cc = "\nSome key metrics\nThe middle %d columns are relative value devided by %s\n" % (len(type_name),type_name[0])
    cc += "The last column is the calculated value of %s\n" % type_name[0]
    cc += TableOutput(outtable,maxlen=50,digit=digit,left=left,scintific=scintific)
    cc += "".join(command)
    #cc += "Notice: the exmaples with value of None for some job type will be excluded.\n"
    return cc,allmetric_value    

def GetAllResults(jsonf):
    #allresults = {"example_name":[examplenames], "type_name":[types],"result": [result1,result2,result3,....] }
    #result is a two dimension list, 1st dimension is for result of each example in each type,
    #second dimension is types
    #len("example_name") = len("result"[i]) = len("result"[j])
    #len("type_name") = len("result")
    #result123 = [example1_result, example2_result, example3_result]

    allresults = json.load(open(jsonf)).get('allresults')
    example_name = []
    example_name_layer = allresults.get("example_name_layer",-1)
    for ie in allresults.get("example_name",[]):
        example_name += [i.split('/')[example_name_layer] for i in glob.glob(ie)]
    
    type_name = allresults.get("type_name",[])
    type_idx = allresults.get("type_idx",[0 for i in range(len(type_name))])
    
    result = []
    for ir in allresults.get("result",[]):
        result.append([])
        for iir in ir:
            alliir = glob.glob(iir)
            alliir.sort()
            for iiir in alliir:
                iresult = json.load(open(iiir))
                if len(iresult.keys()) > 0:
                    iresult = list(iresult.values())[0]
                result[-1].append(iresult)
    
    print("example number: %d" % (len(example_name)))
    print("type number: %d" % (len(type_name)))
    print("type idx: %d" % (len(type_idx)))

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
        sys.exit(1)

    outparams = allresults.get("outparams",[])
    metrics = allresults.get("metrics",[])
    outparams_expand = allresults.get("outparams_expand",{})
    outparams_command = allresults.get("outparams_command",{})
    plot = allresults.get("plot",[])
    
    return {"example_name": example_name, 
            "type_name": type_name, 
            "result": result,
            "type_idx": type_idx, 
            "outparams": outparams,
            "outparams_expand":outparams_expand,
            "outparams_command":outparams_command,
            "metrics": metrics,
            "plot":plot}

def OutResultArgs(parser):  
    parser.description = "This script is used to output the summary of results"
    parser.add_argument('-p', '--param', type=str, default=None, help='the parameter file, should be .json type')
    return parser

def outresult(param):
    if not CheckFile(param.param):
        sys.exit(1)
        
    allresults = GetAllResults(param.param)
    cc_outparam,allparam_value = OutParam(allresults)
    cc_outmetrics,allmetric_value = OutMetrics(allresults,allparam_value)
    print(cc_outmetrics,cc_outparam)

def main():
    parser = argparse.ArgumentParser()
    outresult(OutResultArgs(parser).parse_args())
    
if __name__ == "__main__":
    main()


