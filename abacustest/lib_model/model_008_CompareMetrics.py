from ..model import Model
import os, glob, json
from . import comm
import numpy as np
from abacustest.outresult import pandas_out
from abacustest.report import gen_html
from abacustest.lib_report.table import output_float
from .comm_plot import auto_set_yaxis


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
        parser.add_argument("--name",type=str,help="The name of these two object, which is used in plot",default=None,action="extend",nargs=2,)
        parser.add_argument("--sort",default=1,const=1,help="If sort the example",nargs="?",)

        
    def run(self,params):
        m1f,m2f = params.file
        compare_metrics = params.metric
        const_metrics = params.const
        cm = CompareMetrics(m1f,m2f, compare_metrics, const_metrics,params.show, params.name,params.sort)
        cm.run()
        
                    
class CompareMetrics:
    def __init__(self, m1f,m2f, compare_metrics, const_metrics,show_m,name, sort):
        self.m1f = m1f
        self.m2f = m2f
        self.compare_metrics = compare_metrics
        self.const_metrics = const_metrics
        self.show_metrics = show_m
        self.name = name
        self.sort = sort
        print(f"Compare {m1f} and {m2f}")
        print(f"Metrics to compare: {compare_metrics}")
        print(f"Metrics to show: {const_metrics}")
    
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
        self.m1 = m1
        self.m2 = m2
        
        report, sm = self.cal_dev(m1,m2)
        self.plot_metrics(report,"compare_metrics.png")
        json.dump({"fig_results":{"type": "image", "file": "compare_metrics.png"}}, open("supermetrics.json", "w"), indent=4)
        
        json.dump(report,open("report.json","w"),indent=4)    
        json.dump(sm,open("summary.json","w"),indent=4)
        print( "The report is saved in report.json")
        print(f"The results of {self.m1f} and {self.m2f} are compared, and the deviation between them.") 
        print( "  Dev means the difference between '{self.m1f}' and '{self.m2f}'.")
        print( "  MaxDev means the maximum absolue difference between two values.")
        print(f"  Ratio means '{self.m1f}' devided by '{self.m2f}'")
        print("")
        pandas_out(report) 
        print("")
        print("Summary:")
        print(sm)
        
        # gen html
        gen_html({"content": [{"type": "supermetrics", "content": "summary.json"},
                              {"type": "metrics","content": "report.json"}]
                  },
                 "abacustest.html")
        print("The html file is saved in abacustest.html")   
    
    def is_float(self, v):
        try:
            float(v)
            return True
        except:
            return False
    
    def plot_metrics(self,report,fname):
        nmetrics = len(self.compare_metrics)
        if nmetrics == 0:
            print("No metrics to plot")
            return

        examples = list(report.keys())
        if self.sort:
            examples.sort()
        
        if self.name is None:
            label1 = self.m1f
            label2 = self.m2f
        else:
            label1 = self.name[0]
            label2 = self.name[1]
            
        # now, we need to calculate how many figs we need
        plot_values = []
        for im in self.compare_metrics:
            x1 = range(len(examples))
            x2 = range(len(examples))
            
            y1 = [self.m1[ie].get(im,None) for ie in examples]
            y2 = [self.m2[ie].get(im,None) for ie in examples]
            x1, y1 = comm.clean_none_list(x1,y1)
            x2, y2 = comm.clean_none_list(x2,y2)
            if len(x1) == 0 or len(x2) == 0:
                continue
            # check if y is str
            #if any([isinstance(iv,str) for iv in y1+y2]):
            #    continue
            
            # only plot the scalar
            if any([not self.is_float(iv) for iv in y1+y2]):
                continue
            
            plot_values.append({"name":im,"x1": x1, "y1": y1, "x2": x2, "y2": y2})
            if f"Dev({im})" in report[examples[0]]:
                dev_x = range(len(examples))
                dev_y = [report[ie][f"Dev({im})"] for ie in examples]
                dev_x, dev_y = comm.clean_none_list(dev_x,dev_y)
                plot_values[-1]["sub_x"] = dev_x
                plot_values[-1]["sub_y"] = dev_y
                plot_values[-1]["sub_name"] = f"Dev({im})"
            if f"Ratio({im})" in report[examples[0]]:
                dev_ratio_x = range(len(examples))
                dev_ratio_y = [report[ie][f"Ratio({im})"] for ie in examples]
                dev_ratio_x, dev_ratio_y = comm.clean_none_list(dev_ratio_x,dev_ratio_y)
                plot_values.append({"name":im, "x1": x1, "y1": y1, "x2": x2, "y2": y2, 
                                        "sub_x": dev_ratio_x, "sub_y": dev_ratio_y, "sub_name": f"Ratio({im})"})
            
        ncol = min(2, len(plot_values))
        nrow = len(plot_values)//ncol
        while nrow*ncol < len(plot_values): nrow += 1   
        
        width = max(6, 0.4*len(examples))
        if width > 12:
            ncol = 1
            nrow = len(plot_values)
        
        import matplotlib.pyplot as plt
        
        fig, axs = plt.subplots(nrow,ncol,figsize=(width*ncol,6*nrow))
        fontsize = 18
        for i,plot_value in enumerate(plot_values):
            metric = plot_values[i]["name"]
            if nrow == 1 and ncol == 1:
                ax = axs
            elif nrow == 1 or ncol == 1:
                ax = axs[i]
            else:
                ax = axs[i//ncol,i%ncol]
            ax.plot(plot_value["x1"],plot_value["y1"],"r-o",label=label1)
            ax.plot(plot_value["x2"],plot_value["y2"],"b-^",label=label2)
            ax.legend(loc="upper left",fontsize=fontsize)
            title = metric
            if "sub_name" in plot_value:
                title += f", {plot_value['sub_name']}"
                if plot_value['sub_name'].startswith("Dev"):
                    title += f" (MeanDev: {output_float(self.cal_mean(plot_value['sub_y']))})"
                elif plot_value['sub_name'].startswith("Ratio"):
                    title += f" (GeoMeanRatio: {output_float(self.cal_mean(plot_value['sub_y'],geom=True))})"
            ax.set_title(title,fontsize=fontsize+2)
            ax.set_xlabel("Example",fontsize=fontsize+2)
            ax.set_ylabel(metric,fontsize=fontsize+2)
            ax.set_xticks(range(len(examples)))
            ax.set_xticklabels(examples,rotation=30,fontsize=fontsize)
            ax.tick_params(axis='y',labelsize=fontsize)
            auto_set_yaxis(ax, plot_value["y1"]+plot_value["y2"], margin1=0.1, margin2=0.3)
            
            if "sub_x" in plot_value:
                ax2 = ax.twinx()
                ax2.plot(plot_value["sub_x"],plot_value["sub_y"],"g--",label=plot_value["sub_name"],marker='o',markerfacecolor='none')
                ax2.set_ylabel(plot_value["sub_name"],fontsize=fontsize+2)
                ax2.legend(loc="upper right",fontsize=fontsize)
                ax2.tick_params(axis='y',labelsize=fontsize)
                auto_set_yaxis(ax2, plot_value["sub_y"], margin1=0.1, margin2=0.3)
                
            
        plt.tight_layout()
        plt.savefig(fname,dpi=400)
        plt.close()

    def cal_mean(self, values, geom=False):
        values = comm.clean_none_list(values)
        if len(values) == 0:
            return None
        if geom:
            return np.exp(np.mean(np.log(values)))
        else:
            return np.mean(values)

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

    def cal_dev(self, m1, m2):
        examples = list(set(m1.keys()) & set(m2.keys()))
        if self.sort:
            examples.sort()
        
        e1 = [[m1.get(im,{}).get(ik, None) for im in examples] for ik in self.compare_metrics]
        e2 = [[m2.get(im,{}).get(ik, None) for im in examples] for ik in self.compare_metrics]
        dev = [[] for i in range(len(self.compare_metrics))]
        ratio = [[] for i in range(len(self.compare_metrics))]
        for i in range(len(self.compare_metrics)):
            for j in range(len(examples)):
                v1 = e1[i][j]
                v2 = e2[i][j]
                idev, idevratio = self.cal_one_dev(v1,v2)
                dev[i].append(idev)
                ratio[i].append(idevratio)
        
        # check if the dev or ratio is None
        has_dev = [any([iv is not None for iv in idev]) for idev in dev]
        has_dev_ratio = [any([iv is not None for iv in idev]) for idev in ratio]
        
        report = {} 
        for idx_e, ie in enumerate(examples):
            report[ie] = {}
            for im in self.const_metrics:
                report[ie][im] = m1.get(ie,{}).get(im,None)
            for idx_m, im in enumerate(self.compare_metrics):
                report[ie][f"{im}_1"] = m1.get(ie,{}).get(im,None)
                report[ie][f"{im}_2"] = m2.get(ie,{}).get(im,None)
                if has_dev[idx_m]:
                    report[ie][f"Dev({im})"] = dev[idx_m][idx_e]
                if has_dev_ratio[idx_m]:
                    report[ie][f"Ratio({im})"] = ratio[idx_m][idx_e]
        
        sm = {}
        for i, im in enumerate(self.compare_metrics):
            if has_dev[i]:
                dev_nonone = [iv for iv in dev[i] if iv is not None]
                if len(dev_nonone) > 0:
                    sm[f"MeanDev({im})"] = self.cal_mean(dev_nonone)
                    sm[f"MaxDev({im})"] = max(dev_nonone)
                    sm[f"MinDev({im})"] = min(dev_nonone)
            if has_dev_ratio[i]:
                ratio_nonone = [iv for iv in ratio[i] if iv is not None]
                if len(ratio_nonone) > 0:
                    sm[f"GeoMeanRatio({im})"] = self.cal_mean(ratio_nonone,geom=True)
                    sm[f"MaxRatio({im})"] = max(ratio_nonone)
                    sm[f"MinRatio({im})"] = min(ratio_nonone)
                                          
        return report, sm             

    def cal_one_dev(self,v1,v2):
        if None in [v1,v2]:
            return None,None
        try:
            if isinstance(v1,str) or isinstance(v2,str):
                return f"{v1}|{v2}",None
            elif isinstance(v1,(int,float)) and isinstance(v2,(int,float)):
                if v2 == 0:
                    return v1-v2, None
                else:
                    return v1-v2, v1/v2
            else:
                return abs(np.array(v1)-np.array(v2)).max().tolist(), None
        except:
            return None, None          
    

    