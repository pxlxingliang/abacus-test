from ..model import Model
import os, glob, json, traceback
from . import comm
import numpy as np
from abacustest.outresult import pandas_out
from abacustest.report import gen_html
from abacustest.lib_report.table import output_float
from .comm_plot import auto_set_yaxis
from abacustest.lib_collectdata.collectdata import RESULT


class SPTestModel(Model):
    @staticmethod
    def model_name(): # type: ignore
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "sptest"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "Show the results of FP test"
    
    @staticmethod
    def add_args(parser):
        parser.add_argument("-j","--job",default=[], action="extend",nargs="*" ,help='the path of calculated jobs, or a json file containing results. If the jobpath is in json file, the results in job path is prefered.')
        parser.add_argument("--jobtype", default="abacus",type=str, help="the job type, should be abacus/qe/vasp. Default is abacus")
        parser.add_argument("--ref",type=str,default=None, help="the reference results in json file",)
        parser.add_argument("--testn",type=str,default="TEST", help="the label of test jobs",)
        parser.add_argument("--refn",type=str,default="REF", help="the label of reference results",)

        
    def run(self,params):
        cm = SPTest(params.job, params.jobtype, params.ref, params.testn, params.refn)
        cm.run()
        
                    
class SPTest:
    def __init__(self, jobs, jobtype, ref, testn, refn,):
        self.test_r = self.read_test(jobs, jobtype)
        self.ref_r = self.read_reference(ref)
        self.testn = testn
        self.refn = refn
    
    def run(self):
        if len(self.test_r) == 0:
            print("No results found")
            return
        # save test_r to metrics.json
        json.dump(self.test_r, open("metrics.json", "w"), indent=4)
        
        # 1. gen report.json
        report, cretria_set = self.gen_report()
        json.dump(report, open("report.json", "w"), indent=4)
        
        # 2. plot the statistical results
        self.plot_report(report, cretria_set, "sptest.png")
        json.dump({"fig_results":{"type": "image", "file": "sptest.png"}}, open("supermetrics.json", "w"), indent=4)
        
        # 3. generate the abacus.html
        gen_html({"content": [{"type": "metrics",
                               "content": "report.json",
                               "title": "Table of results.",
                               "criteria": cretria_set},
                              {"type": "image","content": "sptest.png"}]
                  },
                 "abacustest.html")
        
    
    def clean_example(self, r):
        return {k.rstrip("/"): r[k] for k in r.keys()}  
    
    def write_testr(self, fname="metrics.json"):
        json.dump(self.test_r, open(fname, "w"), indent=4)
    
    def read_reference(self, ref):
        if ref is None:
            return None
        else:
            try:
                return self.clean_example(json.load(open(ref)))
            except:
                traceback.print_exc()
                return None
    
    def read_test(self, jobs, jobtype):
        job_path = []
        result_file = []
        for job in jobs:
            if os.path.isfile(job) and job.endswith(".json"):
                result_file.append(job)
            elif os.path.isdir(job):
                job_path.append(job)
            else:
                continue
    
        r = {}
        for i in result_file:
            try:
                r.update(json.load(open(i)))  
            except:
                traceback.print_exc()
                print(f"Error in reading {i}, skip it")
            
        for i in job_path:
            try:
                result = RESULT(path = i, fmt = jobtype)
                ir = {key: result[key] for key in ["natom","label","nelec","normal_end", "converge", "total_time", "scf_steps", 
                                                   "energy", "force", "stress", "virial"]}
                r.update({i: ir})
            except:
                traceback.print_exc()
                print(f"Error in reading {i}, skip it")
        
        return self.clean_example(r)
    
    
    def gen_label(self, ilabel):
        if ilabel is None:
            return None
        elif isinstance(ilabel, str):
            return ilabel
        elif isinstance(ilabel, list):
            label_uniq = []
            for i in ilabel:
                if i not in label_uniq:
                    label_uniq.append(i)
            return "".join([i+str(ilabel.count(i)) for i in label_uniq])
        else:
            return None
    
    def cal_dev(self, r1, r2):
        if None in [r1, r2]:
            return None
        if isinstance(r1, (int, float)) and isinstance(r2, (int, float)):
            return r1 - r2
        if isinstance(r1, list) and isinstance(r2, list):
            r1 = np.array(r1).flatten()
            r2 = np.array(r2).flatten()
            if r1.shape != r2.shape:
                return None
            return np.max(np.abs(r1 - r2)).tolist()
            
    def remove_none(self, r):
        new_r = {}
        value_keys = list(r.values())[0].keys()
        values = {key: all([r[i][key] is None for i in r]) for key in value_keys}
        
        for example, value in r.items():
            new_r[example] = {key: value[key] for key in value_keys if not values[key]}
        return new_r
      
    def sort_report(self, r, key_list):
        new_r = {}
        for example, value in r.items():
            new_r[example] = {key: value[key] for key in key_list if key in value}
        return new_r              
        
    def gen_report(self):
        report = {}
        for example, result in self.test_r.items():
            report[example] = {
                "label": self.gen_label(result.get("label")),
                "natom": result.get("natom"),
                "nelec": result.get("nelec"),
                "NormalEnd": result.get("normal_end"),
                "Converge": result.get("converge"),
                "TotTime(s)": result.get("total_time"),
                "SCFSteps": result.get("scf_steps"),
            }
            
            if self.ref_r is not None:
                for key,key_result in zip(["NormalEnd", "Converge", "TotTime(s)", "SCFSteps"],
                                          ["normal_end", "converge", "total_time", "scf_steps"]):
                    report[example][f"{key}({self.testn})"] = report[example].pop(key)
                    report[example][f"{key}({self.refn})"] = self.ref_r.get(example, {}).get(key_result)
                
                for key, name in zip(["energy", "force", "stress"],
                                     ["DevE(eV)", "MaxDevF(eV/A)", "MaxDevS(kbar)"]):
                    report[example][name] = self.cal_dev(result.get(key), self.ref_r.get(example, {}).get(key))
                
                if None in [report[example]["natom"], report[example]["DevE(eV)"]]:
                    report[example]["DevE(meV/atom)"] = None
                else:
                    report[example]["DevE(meV/atom)"] = report[example]["DevE(eV)"] * 1000 / report[example]["natom"]
            
            for k1, k2, k3 in [["TotTime(s)", "SCFSteps", "Time/step(s)"],
                                  [f"TotTime(s)({self.testn})", f"SCFSteps({self.testn})", f"Time/step(s)({self.testn})"],
                                  [f"TotTime(s)({self.refn})", f"SCFSteps({self.refn})", f"Time/step(s)({self.refn})"]]:
                if None not in [report[example].get(k1), report[example].get(k2)]:
                    report[example][k3] = report[example][k1]/report[example][k2]
                else:
                    report[example][k3] = None

        # remove the key with all None
        report = self.remove_none(report)
        # now sort the keys, if key is not in report, then it will be ignored
        report = self.sort_report(report, ["label", "natom", "nelec", 
                                           "NormalEnd", "Converge", "TotTime(s)","SCFSteps", "Time/step(s)",
                                           f"NormalEnd({self.testn})", f"Converge({self.testn})", f"TotTime(s)({self.testn})", f"SCFSteps({self.testn})", f"Time/step(s)({self.testn})",
                                           f"NormalEnd({self.refn})", f"Converge({self.refn})", f"TotTime(s)({self.refn})", f"SCFSteps({self.refn})", f"Time/step(s)({self.refn})",
                                           "DevE(eV)", "DevE(meV/atom)", "MaxDevF(eV/A)", "MaxDevS(kbar)"])
        metric_keys = list(report[list(report.keys())[0]].keys())
        cretria_set = {
            "NormalEnd": "bool(x)",
            "Converge": "bool(x)",
            f"NormalEnd({self.testn})": "bool(x)",
            f"Converge({self.testn})": "bool(x)",
            "DevE(meV/atom)": "x < 1",
            "MaxDevF(eV/A)": "x < 0.01",
            "MaxDevS(kbar)": "x < 1",
        }
        cretria_set = {k: v for k, v in cretria_set.items() if k in metric_keys}
        return report, cretria_set
    
    def plot_report(self, report, cretria_set, fname):
        metric_keys = list(report[list(report.keys())[0]].keys())
        example_name = list(report.keys())
        example_idx = np.arange(len(example_name))
        
        # accuracy results
        acc_keys = ["DevE(meV/atom)", "MaxDevF(eV/A)", "MaxDevS(kbar)"]
        acc_data = [None if i not in metric_keys else [report[k][i] for k in example_name] for i in acc_keys]
        acc_ifplot = [True if i is not None else False for i in acc_data]
        
        # performance results
        per_keys = ["TotTime(s)", "SCFSteps", "Time/step(s)"]
        per_data1 = [None if f"{i}({self.testn})" not in metric_keys else [report[k][f"{i}({self.testn})"] for k in example_name] for i in per_keys]
        per_data2 = [None if f"{i}({self.refn})" not in metric_keys else [report[k][f"{i}({self.refn})"] for k in example_name] for i in per_keys]
        per_ifplot = [True if per_data1[i] is not None else False for i in range(len(per_keys))]
        
        nfig = sum(acc_ifplot) + sum(per_ifplot)
        if nfig == 0:
            print("No data to plot")
            return
        
        ncol = min(2, nfig)
        nrow = nfig//ncol
        while nrow*ncol < nfig: nrow += 1   
        
        width = max(6, 0.4*len(example_name))
        if width > 12:
            ncol = 1
            nrow = nfig
        fontsize = 18
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(nrow,ncol,figsize=(width*ncol,6*nrow))
        axs = np.array(axs).flatten()
        for i in range(len(acc_keys)):
            if acc_ifplot[i]:
                criteria_line = float(cretria_set[acc_keys[i]].split("<")[-1])
                x, y = comm.clean_none_list(example_idx, acc_data[i])
                axs[i].bar(x, y, label=acc_keys[i])
                axs[i].set_xticks(example_idx)
                axs[i].set_xticklabels(example_name, rotation=30, fontsize=fontsize-2, ha="right")
                axs[i].set_title(acc_keys[i], fontsize=fontsize)
                auto_set_yaxis(axs[i], y, log_scale_threshold=100)
                axs[i].axhline(criteria_line, color="red", linestyle="--", label=f"{criteria_line}")
                
                axs[i].set_xlabel("Example", fontsize=fontsize)
                axs[i].set_ylabel(acc_keys[i], fontsize=fontsize)
                axs[i].legend(fontsize=fontsize-2)
                
        for i in range(len(per_keys)):
            if per_ifplot[i]:
                idx = i + sum(acc_ifplot)
                if per_data2[i] is None:
                    x, tot_y = comm.clean_none_list(example_idx, per_data1[i])
                    axs[idx].bar(x, y, label=f"{self.testn}")
                else:
                    tot_y = []
                    x, y = comm.clean_none_list(example_idx, per_data1[i])
                    tot_y += y
                    axs[idx].bar(np.array(x)-0.2, y, label=f"{self.testn}", width=0.4)
                    x, y = comm.clean_none_list(example_idx, per_data2[i])
                    tot_y += y
                    axs[idx].bar(np.array(x)+0.2, y, label=f"{self.refn}", width=0.4)
                    
                axs[idx].set_xticks(example_idx)
                axs[idx].set_xticklabels(example_name, rotation=30, fontsize=fontsize-2, ha="right")
                axs[idx].set_title(per_keys[i], fontsize=fontsize)
                axs[idx].set_xlabel("Example", fontsize=fontsize)
                axs[idx].set_ylabel(per_keys[i], fontsize=fontsize)
                axs[idx].legend(fontsize=fontsize-2, loc="upper left")
                auto_set_yaxis(axs[idx], tot_y,log_scale_threshold=100)
                
                if per_data2[i] is not None:
                    ratio = [None if per_data1[i][j] is None or per_data2[i][j] is None else per_data1[i][j]/per_data2[i][j] for j in range(len(example_name))]
                    x, y = comm.clean_none_list(example_idx, ratio)
                    if len(y) == 0: continue
                    # calculate the geometric mean of ratio
                    ratio = np.prod(y)**(1/len(y))
                    axs1 = axs[idx].twinx()
                    axs1.plot(x, y, color="red", linestyle="--", label=f"{self.testn}/{self.refn}={output_float(ratio)}")
                    axs1.legend(fontsize=fontsize-2, loc="upper right")
                    auto_set_yaxis(axs1, y,log_scale_threshold=100)
        plt.tight_layout()
        plt.savefig(fname)
        plt.close()
        
        
        
        
    
        
        
                    
                
                
        
        
        

    