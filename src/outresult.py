#!/usr/bin/env python
import os,sys,argparse,glob,json,traceback

def CheckFile(ifile):
    if os.path.isfile(ifile):
        return True
    else:
        print("WARNING: Can not find %s!!!" % ifile)
        return False

def TableOutput(datalist,maxlen=18,digitmax=3,digit=None,left=None,outframe=True,sep=" "):
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
            if type(datalist[i][j]) == float:
                if digit != None :
                    if digit[j] >= 0:
                        lencol = len(("%."+str(digit[j])+"f")%datalist[i][j])
                    else:
                        lencol = len(str(datalist[i][j]))
                else:
                    lencol = len(("%."+str(digitmax)+"f")%datalist[i][j])
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
            if type(j) in [float,int]:
                if digit != None:
                    if digit[ij] >= 0:
                        j = ("%."+str(digit[ij])+"f") % j
                    else:
                        j = str(j)
                else:
                    j = ("%."+str(digitmax)+"f") % j
            value.append(str(j) )
        context += format1 % tuple(value)
        context += '\n'

    if outframe:
        context += totallength*"-" + "\n"

    return context

def GetParamValue(result,param):
    if "'" in param:
        formula = ""
        readparam = False
        paramname = ""
        for i in param:
            if i == "'":
                if not readparam:
                    readparam = True
                    paramname = ""
                else:
                    value = result.get(paramname,None)
                    if value == None:
                        print("value of %s is None, set %s to be None" % (paramname,param))
                        return None
                    if not isinstance(value,(int,float)):
                        print("value of %s is " % paramname, value,", whose type is %s, is not int or float, set %s to be None" % (type(value),param))
                        return None
                    formula += str(value)
                    readparam = False
            elif readparam:
                paramname += i
            else:
                formula += i
        try:
            return eval(formula)
        except:
            traceback.print_exc()
            return None       
    else:
        return result.get(param,None)

def OutParam(allresult,split_example="----"):
    param = allresult.get("params")
    ncol = len(param) + 1
    outtable = [['example'] + [i[0] for i in param]]  # table head
    if split_example != None: outtable.append(ncol * [split_example])
    
    results = allresult.get("result",[[]])
    example_name = allresult.get("example_name",[])
    type_idx = allresult.get("type_idx",[])
    for i in range(len(results[0])):
        for j in range(len(results)):
            iexample = [example_name[i]] if j==0 else [" "]   #add example name
            for iparam in param:
                iexample.append(GetParamValue(results[j][i],iparam[1][type_idx[j]])) #add the value
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
    cc += TableOutput(outtable[:-1],maxlen=18,digit=digit,left=left)
    return cc

def CalGM(valuelist,func=None):
    def a(i):
        return i
    if func == None:
        func = a

    n = len(valuelist)
    gm = 1
    for i in valuelist:
        gm *= func(i)
    return (gm)**(1/n)

def GetValueList(allresult,idx,keyname):
    value = []
    skip_idx = []
    for i,ir in enumerate(allresult):
        v = ir[idx].get(keyname,None)
        if v == None:
            print("WARNING: %s of %s is None, the calculation of related index will exclude this sample" % (keyname,ir[0]))
            skip_idx.append(i)
        value.append(v)
    return value,skip_idx

def GetInverseKeyValue(allresult,keyname,outname):
    totaltime1,skip_idx1 = GetValueList(allresult,1,keyname)
    totaltime2,skip_idx2 = GetValueList(allresult,2,keyname)
    totaltime1_clean = []
    totaltime2_clean = []
    skip_idx = skip_idx1 + skip_idx2
    for i in range(len(totaltime1)):
        if i not in skip_idx:
            totaltime1_clean.append(totaltime1[i])
            totaltime2_clean.append(totaltime2[i])
    def inverse(i):
        return 1.0/i
    return [outname,CalGM(totaltime1_clean,inverse)/CalGM(totaltime2_clean,inverse)]

def KeyIndex(allresult):
    outtable = [['Index','GM(ABACUS)/GM(QE)']]
    cc = 'GM: Geometric Mean\n'
    # GM
    for keyname,outname in [['total_time','inverse(total time)'],\
                            ['scf_steps','inverse(SCF steps)'],\
                            ['force_time','inverse(force time)'],\
                            ['stress_time','inverse(stress time)'],\
                            ['ibzk','inverse(irreducible K points)']]:
        outtable.append(GetInverseKeyValue(allresult,keyname,outname))
    
    cc += TableOutput(outtable,maxlen=30,digitmax=3,left=[True,False])
    return cc

def GetAllResults():
    #allresults = {"example_name":[examplenames], "type_name":[types],"result": [result1,result2,result3,....] }
    #result is a two dimension list, 1st dimension is for result of each example in each type,
    #second dimension is types
    #len("example_name") = len("result"[i]) = len("result"[j])
    #len("type_name") = len("result")
    #result123 = [example1_result, example2_result, example3_result]

    allresults = json.load(open('plusu.json')).get('allresults')
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

    params = allresults.get("outparams",[])
    metrics = allresults.get("metrics",[])
    
    return {"example_name": example_name, 
            "type_name": type_name, 
            "result": result,
            "type_idx": type_idx, 
            "params": params,
            "metrics": metrics}

def CleanData(list2):
    length = len(list2[0])
    for i in range(len(list2[0])):
        for j in range(len(list2)):
            if list2[j][i] == None:
                for k in range(len(list2)):
                    del list2[k][i-length]
                break
    return list2

def main():
    allresults = GetAllResults()
    print(OutParam(allresults))
    
    total_time = []
    for ir in allresults['result']:
        total_time.append([i.get('total_time',None) for i in ir])
    total_time = CleanData(total_time)
    print(len(total_time[0]))
    inverse = lambda x : 1.0/x
    print("GM of inverse of total_time:")
    gm = []
    for i,it in enumerate(total_time):
        gm.append(CalGM(it,inverse))
    print(allresults['type_name'])
    print([i/gm[0] for i in gm])
if __name__ == "__main__":
    main()


