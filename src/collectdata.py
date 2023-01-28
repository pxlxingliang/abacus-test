#!/usr/bin/env python3
import os,sys,json
sys.path.append(os.path.split(__file__)[0])
from lib_collectdata.collectdata import RESULT
import argparse

def parse_param(paramf):
    if os.path.isfile(paramf):
        print("Read the parameters from %s" % paramf)
        allparams = json.load(open(paramf))
        if 'PARAM' not in allparams:
            print("PARAM is not defined in %s" % paramf)
            return []
        else:
            return allparams["PARAM"]
    else:
        print("ERROR: can not find file %s" % paramf)
        return []

def parse_value(abacus_result,allparams):
    allresult = {}
    for param in allparams:
        if isinstance(param,str):
            allresult[param] = abacus_result[param]
        elif isinstance(param,dict):
            for k,v in param.items():
                if not isinstance(k,str):
                    print("Error: %s should be a string, skip!" % str(k))
                    continue

                value = abacus_result[k]
                if not isinstance(value,dict):
                    print("the value of %s is not a dictionary" % k)
                    allresult[k] = value
                else:
                    if isinstance(v,str):
                       allresult["%s/%s"% (k,v)] = value.get(v,None)
                    elif isinstance(v,(list,tuple)):
                        for iv in v:
                            allresult["%s/%s"% (k,iv)] = value.get(iv,None)
                    else:
                        print("%s should be str, list or tuple" % str(v))
        else:
            print("%s should be str or dict" % str(param))
    return allresult                        

def Parser():
    parser = argparse.ArgumentParser(description="This script is used to collect some key values from the output of ABACUS/QE/VASP jobs")
    parser.add_argument('-j', '--jobs', default=["."], help='the path of jobs', action="extend",nargs="*")
    parser.add_argument('-t', '--type', type=int, default=0, help='0:abacus, 1:qe, 2:vasp. Default: 0',choices=[0,1,2])
    parser.add_argument('-p', '--param', type=str, default=None, help='the parameter file, should be .json type')
    parser.add_argument('-o', '--output', type=str, default="result.json",help='the file name to store the output results, default is "result.json"')
    parser.add_argument('--outparam', nargs='?',type=int, const=1, default=0,help='output the registed parameters, you can set the type by -t or --type to choose abacus/qe/vasp. 0: No, 1: yes')
    return parser.parse_args()

def main():
    param = Parser()
    outputf = param.output
    paramf = param.param
    alljobs = param.jobs if len(param.jobs) == 1 else param.jobs[1:]

    alltype = {0:"abacus",1:"qe",2:"vasp"}
    jobtype = alltype.get(param.type)
    
    if param.outparam:
        RESULT(fmt=jobtype,outparam=True)
        return

    if paramf == None:
        print("ERROR: you have not define the parameter file. You can specify by -p or --param")
        return
    if not os.path.isfile(paramf):
        print("ERROR: can not find parameter file %s!!!" % paramf)
        return

    allparams = parse_param(paramf)
    allresult = {}
    for ipath in alljobs:
        if not os.path.isdir(ipath):
            print("ERROR: %s is not a directory, skip it!" % ipath)
            continue

        print("Handle %s" % ipath)
    
        result = RESULT(fmt=jobtype, path=ipath)
        allresult[ipath] = parse_value(result,allparams)

    print("Write the results to %s" % outputf)
    json.dump(allresult,open(outputf,"w"))


if __name__ == "__main__":
    main()



