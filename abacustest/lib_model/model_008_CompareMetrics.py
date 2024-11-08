from ..model import Model
import os, glob, json
from . import comm
import numpy as np
from abacustest.outresult import pandas_out


class CompareMetricsModel(Model):
    @staticmethod
    def model_name(): # type: ignore
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "comparem"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "Compare two metrics.json files"
    
    @staticmethod
    def add_args(parser):
        parser.add_argument("-f","--file",type=str,help="two metrics json files",action="extend",nargs=2,)
        parser.add_argument("-m","--metric",default=[],type=str,help="the metrics to be compared",action="extend",nargs="*",)
        parser.add_argument("-c","--const",default = [], type=str,help="the metrics of first metrics file will be shown. Usally be some comm metrics.",action="extend",nargs="*",)
        parser.add_argument("--show", default=0,const=1,help="show the example and metrics",nargs="?",)

        
    def run(self,params):
        m1f,m2f = params.file
        compare_metrics = params.metric
        const_metrics = params.const
        cm = CompareMetrics(m1f,m2f, compare_metrics, const_metrics,params.show)
        cm.run()
        
                    
class CompareMetrics:
    def __init__(self, m1f,m2f, compare_metrics, const_metrics,show_m):
        self.m1f = m1f
        self.m2f = m2f
        self.compare_metrics = compare_metrics
        self.const_metrics = const_metrics
        self.show_metrics = show_m
    
    def run(self):
        if self.show_metrics:
            return self.show()
        
        if not os.path.exists(self.m1f):
            print(f"{self.m1f} does not exist.")
            return
        if not os.path.exists(self.m2f):
            print(f"{self.m2f} does not exist.")
            return

        m1 = json.load(open(self.m1f))
        m2 = json.load(open(self.m2f))
        
        examples = list(set(m1.keys()) & set(m2.keys()))
        examples.sort()
        
        report = {}
        
        for ie in examples:
            report[ie] = {}
            for ic in self.const_metrics:
                v1 = m1.get(ie,{}).get(ic,None)
                v2 = m2.get(ie,{}).get(ic,None)
                if v1 == v2:
                    report[ie][ic] = v1
                else:
                    report[ie][ic] = [v1,v2]
            for ic in self.compare_metrics:
                value1 = m1.get(ie,{}).get(ic,None)
                value2 = m2.get(ie,{}).get(ic,None)
                itype,ivalue,ivalue_ratio = self.max_dev(value1,value2)
                report[ie][f"{itype}({ic})"] = ivalue
                if ivalue_ratio is not None:
                    report[ie][f"{itype}Ratio({ic})"] = ivalue_ratio
        json.dump(report,open("report.json","w"),indent=4)    
        print("The report is saved in report.json")
        print(f"The results of {self.m1f} and {self.m2f} are compared, and the deviation between them.") 
        print("  Dev means the difference between two values.")
        print("  MaxDev means the maximum absolue difference between two values.")
        print(f"  DevRatio means the difference ratio between two values from {self.m2f}")
        pandas_out(report)    
    
    def show(self):
        for mf in [self.m1f,self.m2f]:
            print(f"File: {mf}")
            m = json.load(open(mf))
            examples = list(m.keys())
            examples.sort()
            examples = " ".join(examples)
            print(f"Examples: {examples}")
            metrics = []
            for ie in m.keys():
                metrics += m[ie].keys()
            metrics = list(set(metrics))
            metrics.sort()
            metrics = "\" \"".join(metrics)
            print(f"Metrics: \"{metrics}\"")

    def max_dev(self,v1,v2):
        if None in [v1,v2]:
            return "Dev",None,None
        try:
            if isinstance(v1,str) or isinstance(v2,str):
                return "", f"{v1}|{v2}",None
            elif isinstance(v1,(int,float)) and isinstance(v2,(int,float)):
                if v2 == 0:
                    return "Dev",v1-v2, None
                else:
                    return "Dev",v1-v2, (v1-v2)/v2
            else:
                return "MaxDev",abs(np.array(v1)-np.array(v2)).max().tolist(), None
        except:
            return "Dev",None  ,None          
    
            
    
    