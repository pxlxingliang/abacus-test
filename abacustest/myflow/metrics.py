import os,glob,json,traceback,sys
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.collectdata import parse_value

class Metrics:
    def __init__(self,dft_type="abacus",metrics_name=[],newmethods=[],path=["."],modules=[]):
        self.dft_type = dft_type
        self.metrics_name = metrics_name
        self.newmethods = newmethods
        self.path = path
        self.modules = modules
        pass
    
    def get_metrics(self,save_file = None):
        allvalue = {}
        print("metrics setting (getcwd(), path=, newmethods=, dft_type=, modules=):",os.getcwd(),self.path,self.newmethods,self.dft_type,self.modules)
        print("os.listdir:",os.listdir("."))
        for ipath in self.path:
            for iipath in glob.glob(ipath):
                if not os.path.isdir(iipath):
                    continue
                try:
                    print("get metrics from:",iipath)
                    result = RESULT(fmt=self.dft_type,newmethods=self.newmethods,path=iipath,modules=self.modules)
                    if len(self.metrics_name) == 0:
                        self.metrics_name = result.AllMethod().keys()
                    allvalue[iipath] = parse_value(result,self.metrics_name) 
                    
                    #write version.dat, this file is used by report section
                    if not os.path.isfile("version.dat"):
                        version = result["version"]
                        version = "unknow" if version == None else version
                        with open("version.dat",'w') as f:
                            f.write(version)
                except:
                    traceback.print_exc()
        if save_file != None:
            json.dump(allvalue,open(save_file,'w'),sort_keys=True,indent=4) 
        return allvalue
    
    @staticmethod
    def TransferMetricsOPIO(metrics_io):
        if isinstance(metrics_io,dict):
            poin_metrics = [metrics_io]
        elif isinstance(metrics_io,list):
            poin_metrics = metrics_io
        else:
            poin_metrics = []
        return poin_metrics 

    @staticmethod
    def ParseMetricsOPIO(metrics_io,example_path):
        all_dfttype = {0:"abacus",1:"qw",2:"vasp"}
        dft_type = metrics_io.get("dft_type",None)
        dftt = None
        if dft_type == None:
            return None
        elif isinstance(dft_type,int):
            if dft_type not in all_dfttype:
                print("Unknow dft type '%d'. Supportted:" % dft_type,str(all_dfttype))
            else:
                dftt = all_dfttype[dft_type]
        elif isinstance(dft_type,str):
            if dft_type.lower() == "abacus":
                dftt = "abacus"
            elif dft_type.lower() == "qe":
                dftt = "qe"
            elif dft_type.lower() == "vasp":
                dftt = "vasp"
            else:
                print("Unknow dft type '%s'. Supportted:" % dft_type,str(all_dfttype))
        else:
            print("Unknow dft type '",dft_type,"'. Supportted:",str(all_dfttype))

        if "path" not in metrics_io:
            epath = example_path
        else:
            epath = metrics_io["path"]
            
        if dftt == None:
            return None
        else:
            print(metrics_io)
            return Metrics(dft_type=dftt,
                           path= epath,
                       metrics_name= metrics_io.get("metrics_name",[]),
                       newmethods=metrics_io.get("newmethods",[]),
                       modules = metrics_io.get("modules",[]))       

    @staticmethod
    def SuperMetricsResult(super_metrics_setting):
        if not super_metrics_setting.get("result_file"):
            return None,None,None
        from abacustest import outresult
        allresults = outresult.GetAllResults(super_metrics_setting)
        if not allresults:
            return None,None,None
        split_example=None if len(allresults["type_name"]) == 1 else "----"
        cc_outparam,allparam_value = outresult.OutParam(allresults,split_example=split_example)
        cc_outmetrics,allmetric_value = outresult.OutMetrics(allresults,allparam_value)
        report = ""
        if len(allresults["metrics"]) > 0:
            report += cc_outmetrics
            save_file = super_metrics_setting.get("save_file","superMetric.json")
            json.dump(allmetric_value,open(save_file,'w'),indent=4)
        report += cc_outparam

        type_name = allresults["type_name"]
        if allmetric_value and len(type_name) > 1:
            #if there has more than 1 type, transfer allmetric_value to one type one value
            allmetric_value_final = {}
            for k,v in allmetric_value.items():
                for i,itype in enumerate(type_name):
                    allmetric_value_final[k+"_" +itype] = v[i]
        else:
            allmetric_value_final = allmetric_value     

        return allparam_value,allmetric_value_final,report

    
    @staticmethod
    def rotate_metrics(dict1):
        '''
        dict1 = {example1: {metric1:,metric2:,...}, example2: {}}
        will be rotated to:
        {metric1:[example1_value,example2_value],metric2:[example1_value,example2_value,...]}
        '''
        allkey = [i for i in dict1.keys()]
        allmetrics = dict1[allkey[0]].keys()
        newdict = {}
        for imetric in allmetrics:
            newdict[imetric] = []
            for iexample in allkey:
                newdict[imetric].append(dict1.get(iexample,{}).get(imetric,None))
        return allkey,newdict
    
    @staticmethod
    def Transfer2Table(dict1):
        '''
        dict1 = {example1: {metric1:,metric2:,...}, example2: {}}
        will be transfered to:
        [{"sample":example1,"metric1":,"metric2":,...},{}]
        '''
        from dp.tracking import Table
        allkey = [i for i in dict1.keys()]
        allmetrics = dict1[allkey[0]].keys()
        newlist = []
        for ikey in allkey:
            newlist.append({"sample":ikey})
            for imetric in allmetrics:
                newlist[-1][imetric] = dict1[ikey].get(imetric,None)
        return Table(newlist)
    
    
def ReadMetrics(poin_metrics,do_upload_tracking,example_path):
    """
    poin_metrics is a list of dict
    poin_metrics = [
        {
            "value_from_file": str,     #a json file that has the metrics and value, which will also be uploaded to tracking, should be {"example":{"key":value,...},...} 
            "group_name": str,  #self-defined group name, only used in context for tracking
            "path": ["*"],      #the path of dft jobs
            "dft_type": "abacus",   #dft type, can be abacus/qe/vasp
            "metrics_name": [],     #the name of metrics that will collected from dft job
            "newmethods": [str],    # self-defined methods used in abacustest collectdata
            "modules": [str],       # the modules that abacustest collectdata will load
            "save_file": "metrics.json",  #file to store the value of metrics
        },
        {"path": ["*"], ...}, ...
    ]

    return tracking_values
    tracking_values = [(metrics_value,context),...]
    metrics_value = {metric_name:value}
    
    """
    tracking_values = []
    for im,imetric in enumerate(poin_metrics):
        metrics_value = {}

        metric_form_calc = Metrics.ParseMetricsOPIO(imetric,example_path)
        if metric_form_calc:
            value = metric_form_calc.get_metrics(save_file=imetric.get("save_file","metrics.json"))
        else:
            value = {}

        if do_upload_tracking:
            metric_from_file = imetric.get("value_from_file",None)
            if metric_from_file: 
                if os.path.isfile(metric_from_file):
                    try:
                        for k,v in json.load(open(metric_from_file)).items():
                            if k not in value:
                                value[k] = {}
                            for ik,iv in v.items():
                                value[k][ik] = iv
                    except:
                        traceback.print_exc()
                else:
                    print("Can not find file %s" % metric_from_file,file=sys.stderr)
                    print("Current path: %s, listdir:" % os.path.abspath("."),os.path.listdir("."),file=sys.stderr)
                    
            if value: 
                try:
                    examples,new_dict = Metrics.rotate_metrics(value)
                    for k,v in new_dict.items():
                        metrics_value[k] = v
                    metrics_value["sample_name_%d" % im] = examples
                    metrics_value["metrics_%d"%im] = Metrics.Transfer2Table(value) 
                except:
                    traceback.print_exc()
                    
            context = {"subset":"metrics%d"%im,"datatype":"metrics"}
            if "group_name" in imetric:
                context["group_name"] = imetric.get("group_name")

            if metrics_value:
                tracking_values.append((metrics_value,context))

    return tracking_values

def ReadSuperMetrics(poin_super_metrics,do_upload_tracking):
    """
    poin_super_metrics = [{
        "group_name": str,  # a json file that has the metrics and value, which will also be uploaded to tracking
        "value_from_file": str, # self-defined super metrics, can also define some files
                                # {key:{"type":"image","file":"file.png"}}
        "save_file": "superMetric.json",    #same as that defined in outresult report
        "result_file": ["metrics.json"],
        "example_name_idx":0,
        "outparams":[
            ["ScfSteps", ["scf_steps"], 0],...
        ],
        "metrics":[{}]
    },...]

    return tracking_summary
    tracking_summary = [(metrics_value,context),...]
    metrics_value = {metric_name:value}

    """
    tracking_summary = []  
    log = ""             
    for isuper,super_metrics_setting in enumerate(poin_super_metrics):
        super_metric_value = {}
        allparam_value,allmetric_value,report = Metrics.SuperMetricsResult(super_metrics_setting)

        if do_upload_tracking:
            if allmetric_value:
                for k,v in allmetric_value.items():
                    super_metric_value[k] = v
                from dp.tracking import Table
                super_metric_value["super_metrics_%d"%isuper] = Table([allmetric_value])
            if report:
                super_metric_value["report"] = report
                log += report

            super_metric_from_file = super_metrics_setting.get("value_from_file",None)
            if super_metric_from_file:
                if os.path.isfile(super_metric_from_file):
                    try:
                        for k,v in json.load(open(super_metric_from_file)).items():
                            super_metric_value[k] = v
                    except:
                        traceback.print_exc()
                else:
                    print("Can not find file %s" % super_metrics_setting.get("value_from_file"),file=sys.stderr)
                    print("Current path: %s, listdir:" % os.path.abspath("."),os.path.listdir("."),file=sys.stderr)
            if super_metric_value:
                context = {"subset":"super_metrics%d"%isuper,"datatype":"super_metrics"}
                if "group_name" in super_metrics_setting:
                    context["group_name"] = super_metrics_setting.get("group_name")
                tracking_summary.append((super_metric_value,context))
    
    return tracking_summary,log