import os,sys,json
sys.path.append(os.path.split(__file__)[0])
from .lib_collectdata.collectdata import RESULT
from .lib_collectdata.comm import get_metric_from_str
import argparse
import traceback
from abacustest.lib_prepare.abacus import ReadInput

def parse_param(paramf):
    if os.path.isfile(paramf):
        print("Read the parameters from %s" % paramf)
        allparams = json.load(open(paramf))
        if 'PARAM' not in allparams:
            print("PARAM is not defined in %s" % paramf)
            if "post_dft" in allparams and isinstance(allparams["post_dft"],dict) and "metrics" in allparams["post_dft"] and \
                isinstance(allparams["post_dft"]["metrics"],dict) and "metrics_name" in allparams["post_dft"]["metrics"]:
                    print("Find post_dft/metrics/metrics_name, and will collectdata by these metrics")
                    return allparams["post_dft"]["metrics"]["metrics_name"]
            else:
                return []
        else:
            return allparams["PARAM"]
    else:
        print("ERROR: can not find file %s" % paramf)
        return []

def get_metric_value(result,metric_str):
    if "{" not in metric_str:
        return result[metric_str]
    else:
        metrics = get_metric_from_str(metric_str)
        values = {}
        for i in metrics:
            values[i] = result[i]
        metrics_real_str = metric_str.format(**values)
        try:
            value = eval(metrics_real_str)
        except:
            print("ERROR: can not evaluate %s" % metric_str)
            print("     The real str is: %s" % metrics_real_str)
            traceback.print_exc()
            value = None
        return value

def parse_value(abacus_result,allparams):
    allresult = {}
    for param in allparams:
        if isinstance(param,str):
            allresult[param] = get_metric_value(abacus_result,param)
        elif isinstance(param,dict):
            # if is a dict, then the key is a rename of param, and the value is an eval string, and the key of abacus_result can used by {key}
            for k,v in param.items():
                if not isinstance(k,str):
                    print("Error: %s should be a string, skip!" % str(k))
                    continue
                
                # remain when key is INPUT, indicate want to get the sub value of INPUT that is a dict,
                # such as INPUT["ecutwfc"], INPUT["basis"]
                # this can be achieved by type "{INPUT}['ecutwfc']"
                
                if k in ["INPUT"]:
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
                    if isinstance(v,str):
                        allresult[k] = get_metric_value(abacus_result,v)
                    else:
                        print("%s should be str" % str(v))
        else:
            print("%s should be str or dict" % str(param))
    return allresult                        

def CollectDataArgs(parser):
    parser.description = "This script is used to collect some key values from the output of ABACUS/QE/VASP jobs"
    parser.add_argument('-j', '--jobs', default=["."], help='the path of jobs', action="extend",nargs="*")
    parser.add_argument('-t', '--type', type=int, default=0, help='0:abacus, 1:qe, 2:vasp. Default: 0',choices=[0,1,2])
    parser.add_argument('-p', '--param', type=str, default=None,nargs="*", help='the parameter file, or parameter name.')
    parser.add_argument('-o', '--output', type=str, default="metrics.json",help='the file name to store the output results, default is "metrics.json"')
    parser.add_argument('-m', '--modules',help='add extra modules. Default only module \'job-type\' will be loaded, such as: \'abacus\' for abacus type. You can check all modules by --outparam', action="extend",nargs="*")
    parser.add_argument('--newmethods', help='the self-defined python modules, and shuold be format of import, such as "abc"(the file name is abc.py), "a.b.c" (teh file is a/b/c.py)', action="extend",nargs="*")
    parser.add_argument('--outparam', nargs='?',type=int, const=1, default=0,help='output the registed parameters, you can set the type by -t or --type to choose abacus/qe/vasp. 0: No, 1: yes')
    parser.add_argument('--ref', type=str, nargs='?',default=None,const="resultREF.json",help='A json file includes the reference value of some keys. Generally, get values of keys start with \"delta_\" require this file. Default is resultREF.json')
    return parser

NO_PARAM_WARNING = """
WARNING: you have not defined the parameters to collect.
         You can specify by -p or --param. The context can be like: 
         {
             "PARAM":["natom","total_time"]
         }
         Only collect some common parameters, such as natom, total_time, converge, normal_end, etc.
         
         You can execute below commnd to get all the parameters that can be collected:
             abacustest collectdata --outparam
"""

def collectdata(param):    
    outputf = param.output
    paramf = param.param
    alljobs = param.jobs if len(param.jobs) == 1 else param.jobs[1:]

    alltype = {0:"abacus",1:"qe",2:"vasp"}
    jobtype = alltype.get(param.type)
    
    if param.outparam:
        RESULT(fmt=jobtype,outparam=True,newmethods=param.newmethods,modules=param.modules)
        return   
    if paramf == None:
        allparams = []
    else:
        allparams = []
        for iparam in paramf:
            if os.path.isfile(iparam) and iparam.endswith(".json"):
                allparams += parse_param(iparam)
            else:
                allparams.append(iparam)
    
    if len(allparams) == 0:
        print(NO_PARAM_WARNING)
        
    allresult = {}
    for ipath in alljobs:
        if not os.path.isdir(ipath):
            print("ERROR: %s is not a directory, skip it!" % ipath)
            continue

        print("Handle %s" % ipath)
    
        result = RESULT(fmt=jobtype, path=ipath,newmethods=param.newmethods,modules=param.modules,resultREF=param.ref)

        if os.path.isfile(os.path.join(ipath, "INPUT")):
            input_param = ReadInput(os.path.join(ipath, "INPUT"))
        else:
            input_param = {}
        
        job_type = input_param.get("calculation", "scf")
        if len(allparams) == 0:
            allparams = ["normal_end","converge","nkstot","ibzk",
                "nbands","nelec","natom","scf_steps","total_time",
                "scf_time",
                {
                    "scf_time/step": "{scf_time}/{scf_steps}",
                    "ks_solver": "{INPUT}['ks_solver']"
                },
                "stress_time",
                "energy",
                "energy_per_atom",
                "force",
                "stress",
                "drho_last",
                "denergy_last",
                "version"]
            #allparams = list(result.AllMethod().keys())
            if job_type in ["relax","cell-relax"]:
                allparams.append("relax_steps")
                allparams.append("relax_converge")
                allparams.append("largest_gradient")

                if job_type == "cell-relax":
                    allparams.append("largest_gradient_stress")
                    allparams.append("lattice_constants")
            
            if input_param.get("nspin",1) in [2,4]:
                allparams.append("total_mag")
                allparams.append("absolute_mag")
                allparams.append("atom_mag")

        allresult[ipath] = parse_value(result,allparams)

    print("Write the results to %s" % outputf)
    json.dump(allresult,open(outputf,"w"),indent=4)
    
    from .outresult import pandas_out
    pandas_out(allresult)

def main():
    parser = argparse.ArgumentParser()
    param = CollectDataArgs(parser).parse_args()
    collectdata(param)
    
if __name__ == "__main__":
    main()



