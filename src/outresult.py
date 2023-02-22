import os,sys,argparse,glob,json,traceback
import numpy as np

def pandas_out(allresult):
    """
    allresult = {sample1: {key1:value,key2:value},
                 sample2: {key1:value,key2:value}}
    If value is a list, will print out separately.
    """
    import pandas as pd
    normal_result = {}
    list_result = []
    allsamples = [i for i in allresult.keys()]
    allkeys = [i for i in allresult[allsamples[0]].keys()]
    allkeys_seperate = []
    for ikey in allkeys:
        allkeys_seperate.append(False)
        for i in allsamples:
            if isinstance(allresult[i][ikey],dict):
                allkeys_seperate[-1] = True
                break
            elif isinstance(allresult[i][ikey],list):
                if isinstance(allresult[i][ikey][0],(list,dict)):
                    allkeys_seperate[-1] = True
                    break
        
    for isample in allsamples:
        normal_result[isample] = {}
        for i,ikey in enumerate(allkeys):
            if allkeys_seperate[i]:
                if isinstance(allresult[isample][ikey],(list,dict)):
                    try:
                        list_result.append("%s\t%s:\n" % (isample,ikey) + str(pd.DataFrame.from_dict(allresult[isample][ikey])))
                    except:
                        list_result.append("%s\t%s:\n" % (isample,ikey) + str(allresult[isample][ikey]))
                else:
                    list_result.append("%s\t%s:\n" % (isample,ikey) + str(allresult[isample][ikey]))
            else:
                normal_result[isample][ikey] = allresult[isample][ikey]
                
    if False not in allkeys_seperate:
        normal_result = ""
    else:
        normal_result = str(pd.DataFrame.from_dict(normal_result,orient='index'))
        
    if True not in allkeys_seperate:
        list_result = ""
    else:
        list_result = "\n\n".join(list_result)
        
    print("%s\n\n%s" % (normal_result,list_result))
           
    return normal_result,list_result

def OutResultArgs(parser):  
    parser.description = "This script is used to output the summary of results"
    parser.add_argument('-r', '--result', type=str, default=["result.json"], help='the result file from collectdata, should be .json type',action="extend",nargs="*")
    return parser

def outresult(param):
    allresult_files = param.result if len(param.result) == 1 else param.result[1:]
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
    pandas_out(allresult)

def main():
    param = argparse.ArgumentParser()
    outresult(param)
    
if __name__ == "__main__":
    main()


