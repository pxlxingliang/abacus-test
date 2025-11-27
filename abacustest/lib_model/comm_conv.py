import glob,os,json
import numpy as np
from abacustest.lib_collectdata.collectdata import RESULT
from . import comm_plot,comm

def shift_data(data_list, base_index=0):
    # shift the data to the base_index
    # data_list should a list of float or list of list of float
    if not isinstance(data_list,list):
        return None
    if base_index < 0 or base_index >= len(data_list):
        return None
    
    if isinstance(data_list[base_index],list):
        # it is a list of list of float
        new_list = []
        if None in data_list[base_index]:
            return None
        for i in range(len(data_list)):
            if len(data_list[i]) != len(data_list[base_index]):
                new_list.append(None)
            else:
                new_list.append([data_list[i][j] - data_list[base_index][j] for j in range(len(data_list[i]))])
    else:
        if data_list[base_index] is None:
            return None
        new_list = [data_list[i] - data_list[base_index] for i in range(len(data_list))]
    return new_list

class PostConv:
    """
    Used to do the postprocess of convergence test of ecutwfc and kspacing.
    
    It will firstly collect the metrics of the jobs, and save the results to metrics.json file.
    If the metric_file is defined, it will also read the results in file, and merge the results. This file should have the same format as metrics.json.
    
    The it will rearrnge the results to plotdata, and plot the results.
    
    For default, energy/force/stress/atom_mag/mag_force/band_gap will be plotted.
    
    Args:
        test_key: str, the key of the test, ecutwfc or kspacing. It should be a parameter name of INPUT. Special case is kpt, which will be used to plot the k-point mesh.
        jobs: list of str, the list of the job folders.
        metric_file: str, the file of the metrics. or a liat of metrics files.
        extra_y: dict, the extra y to plot. The key is the metric name, and the value is the y title.
        job_type: str, abacus or qe or vasp, the type of the job.
        shift_data: max/min/False, whether to shift the data. Default is None.
            None: will shift the data to VALUE with max ecutwfc or min kspacing.
            max: will shift the data to the max test_key.
            min: will shift the data to the min test_key
            False: will not shift the data.
    """
    def __init__(self, test_key="ecutwfc",jobs=None, metric_file=None, extra_y=None, job_type="abacus",shift_data=None,x_name=None):
        self.test_key = test_key
        self.jobs = jobs
        self.metric_file = None
        if metric_file is not None:
            if isinstance(metric_file,str):
                self.metric_file = [metric_file]
            elif isinstance(metric_file,list):
                self.metric_file = metric_file
                
        self.extra_y = [extra_y] if isinstance(extra_y,str) else extra_y
        self.job_type = job_type

        self.plot_keys = {"energy_per_atom":"Energy (meV/atom)", 
                          "force":"Force (eV/$\mathrm{\AA}$)", 
                          "stress":"Stress (kbar)", 
                          "atom_mag":"Atomic maganetic moment ($\mathrm{\mu}$B)", 
                          "ds_mag_force": "Maganetic force ($\mathrm{\mu}$B/$\mathrm{\AA}$)",
                          "band_gap":"Band gap (eV)"}
        # key is the metric name, value is the y title
        
        if self.extra_y is not None:
            for iey,iey_unit in self.extra_y.items():
                self.plot_keys[iey] = iey_unit
        
        self.shift_idx = None  # shift the energy to value of shift_idx (0 means the first value, -1 means the last value)
        if shift_data == "max":
            self.shift_idx = -1
        elif shift_data == "min":
            self.shift_idx = 0
        elif shift_data == None:
            if self.test_key in ["ecutwfc", "kpt"]:
                self.shift_idx = -1
            elif self.test_key == "kspacing":
                self.shift_idx = 0
        
        if x_name is not None:
            self.x_name = x_name.capitalize()
        else:
            if test_key == "ecutwfc":
                self.x_name = "Ecutwfc (Ry)"
            elif test_key == "kspacing":
                self.x_name = "Kspacing (1/bohr)"
            elif test_key == "kpt":
                self.x_name = "K-point mesh"
            else:
                self.x_name = test_key.capitalize()

    def run(self):
        metrics_jobs = None
        metrics_file = []
        if self.jobs is not None:
            metrics_jobs = self._collect_metrics()
            json.dump(metrics_jobs, open("metrics.json", "w"), indent=4)
        if self.metric_file is not None:
            metrics_file = [json.load(open(i)) for i in self.metric_file if os.path.isfile(i)]
        
        if metrics_jobs is None and len(metrics_file)==0:
            print("No metrics are found.")
            return
        
        plotdata = self.rearrage_data(metrics_jobs, metrics_file)
        plotdata = self.sort_plotdata(plotdata)
        fnames = self.plot_data(plotdata)
        
        # write a supermetrics.json
        json.dump(plotdata, open("result.json", "w"), indent=4)
        if fnames:
            json.dump({k:{"type": "image", "file": k}  for k in fnames}, open("supermetrics.json","w"),indent=4)
    
    def _collect_metrics(self):
        jobs = self.jobs
        test_key = self.test_key
        job_type = self.job_type
        
        allmetrics = {}
        for job in jobs:
                for ijob in glob.glob(os.path.join(job,"*")):
                    if os.path.isdir(ijob):
                        results = RESULT(path=ijob,fmt=job_type)
                        if test_key == "kpt":
                            test_key_v = results["kpt"]
                        else:
                            if job_type == "abacus":
                                test_key_v = results["INPUT"].get(test_key,None)
                            else:
                                test_key_v = results[test_key]
                                
                        allmetrics[ijob] = {
                            test_key: test_key_v,
                            "energy_per_atom": results["energy_per_atom"],
                            "energy": results["energy"],
                            "natom": results["natom"],
                            "normal_end": results["normal_end"],
                            "denergy_last": results["denergy_last"],
                            "drho_last": results["drho_last"],
                            "converge": results["converge"],
                            "band_gap": results["band_gap"],
                            "scf_steps": results["scf_steps"],
                            "total_time": results["total_time"],
                            "scf_time": results["scf_time"],
                            "force": results["force"],
                            "stress": results["stress"],
                            "ds_mag_force": results["ds_mag_force"],
                            "version": results["version"],
                            "ks_solver": None if not isinstance(results["INPUT"],dict) else results["INPUT"].get("ks_solver",None),
                            "scf_time/step": None if not results["scf_time"] or not results["scf_steps"] else results["scf_time"] / results["scf_steps"]
                        }
                        for ikey in self.plot_keys:
                            if ikey not in allmetrics[ijob]:
                                allmetrics[ijob][ikey] = results[ikey]         
        return allmetrics

    def get_xkey_v(self, iv):
        # we need deal with the x_key value to transfer it to a float or int or string
        x_key_v = iv.get(self.test_key,None)
        if x_key_v is None:
            return None
        
        if self.test_key == "kspacing":
            kspacing = x_key_v
            if isinstance(kspacing,str):
                try:
                    ksplit = [float(ik) for ik in kspacing.split()]
                    if len(ksplit) == 1:
                        kspacing = ksplit[0]
                    elif len(ksplit) == 3:
                        if ksplit[0] == ksplit[1] == ksplit[2]:
                            kspacing = ksplit[0]
                except:
                    return kspacing
            return kspacing
        elif isinstance(x_key_v,str):
            x_key_v_split = x_key_v.split()
            if len(x_key_v_split) == 1:
                try:
                    return float(x_key_v)
                except:
                    return x_key_v
            elif len(x_key_v_split) > 1:
                try:
                    return tuple([float(i) for i in x_key_v_split])
                except:
                    return x_key_v

        return x_key_v
                
    def get_ykey_v(self, iv, y_key):
        y_key_v = iv.get(y_key,None)
        if y_key_v is None:
            return None
        
        if y_key in ["energy","energy_per_atom"]:
            # here will actually use energy_per_atom in unit meV/atom
            if "energy_per_atom" in iv:
                y_key_v = iv["energy_per_atom"]
            else:
                e = iv.get("energy",None)
                na = iv.get("natom",None)
                if e is not None and na is not None:
                    y_key_v = e / na
                elif e is not None and na is None:
                    y_key_v = e
                    self.plot_keys[y_key] = "meV"
                else:
                    y_key_v = None
            if y_key_v is not None:
                y_key_v *= 1000 # transfer to meV/atom
        elif y_key == "mag_force":
            return iv.get("ds_mag_force",None)
        
        # if y_key_v is a list of list, we should transfer it to one list
        if isinstance(y_key_v,list):
            y_key_v = np.array(y_key_v).flatten().tolist()
                
        return y_key_v   
    
    def rearrage_data(self, metrics_jobs, metrics_file):
        # transfer the metrics to a dict with key is example name and value is a dict of subfolder, x_key, and y_keys
        '''
        {
            "example_name": {
                "subfolder": [],
                "x_key": [],  # x_key should be a list of float or int or string
                "y_key1": [], # y_key should be a list of float or int or list, element of y_key is matched with x_key
                "y_key2": [],
                ...
            }
        }
        '''
        plotdata = {}
        x_key = self.test_key
        for allmetrics in [metrics_jobs] + metrics_file:
            if allmetrics is None:
                continue
            
            for i, iv in allmetrics.items():
                x_value = self.get_xkey_v(iv)
                if x_value is None: # do not save the data if the x_key is None
                    continue
                example_name = os.path.dirname(i.rstrip("/"))
                if example_name not in plotdata:
                    plotdata[example_name] = {
                        "subfolder": [],
                        x_key: [],
                    }
                    for ik in self.plot_keys:
                        plotdata[example_name][ik] = []

                if x_value in plotdata[example_name][x_key]: # if one x_value is already in the plotdata, skip it. For case value from jobs and file are the same.
                    continue
                plotdata[example_name]["subfolder"].append(os.path.basename(i.rstrip("/")))
                plotdata[example_name][x_key].append(x_value)
                for ik in self.plot_keys:
                    plotdata[example_name][ik].append(self.get_ykey_v(iv,ik))   
        return plotdata
    
    def sort_kpt(self, kpts):
        # kpts is a list of k-point meshs, such as [[2,2,2],[3,3,3],[4,4,4]]
        # we need sort the kpts by the number of k-points, and then by the k-point mesh
        # the number of k-points is the product of the three elements
        kpoints = [ np.prod(kpt) for kpt in kpts ]
        sort_idx = sorted(range(len(kpoints)), key=lambda k: (kpoints[k], kpts[k]))
        
        # if 3 k-points are the same, then return the first one
        if all([len(set(k)) == 1 for k in kpts]):
            kpts = [i[0] for i in kpts]
        else:
            kpts = ["_".join([str(j) for j in i]) for i in kpts]
        kpts = [kpts[i] for i in sort_idx]
        return sort_idx, kpts
        
    
    def sort_plotdata(self, plotdata):
        # we need sort the plotdata by x_key
        # and if shift_idx is not False, we need shift the data to the shift_idx
        new_plotdata = {}
        for ik,iv in plotdata.items():
            
            if self.test_key == "kpt":
                # if the test_key is kpt, we need sort the kpts
                sort_idx, x_values = self.sort_kpt(iv[self.test_key])
            else:
                x_values = iv[self.test_key]
                sort_idx = sorted(range(len(x_values)), key=lambda k: x_values[k])
                x_values = [x_values[i] for i in sort_idx]
                
            new_plotdata[ik] = {
                "subfolder": [iv["subfolder"][i] for i in sort_idx],
                self.test_key: x_values
            }
            for ikey in self.plot_keys:
                if self.shift_idx is None or ikey in ["band_gap"]:
                    new_plotdata[ik][ikey] = [iv[ikey][i] for i in sort_idx]
                else:
                    iv_shifted = shift_data(iv[ikey],base_index=sort_idx[self.shift_idx])
                    new_plotdata[ik][ikey] = [iv_shifted[i] for i in sort_idx] if iv_shifted is not None else None
        return new_plotdata
    
    def plot_data(self, plotdata):
        allexamples = list(plotdata.keys())
        allexamples.sort()
        fnames = []
        ykeys = list(self.plot_keys.keys())
        ykeys.sort()
        
        for ik in allexamples:
            
            x = plotdata[ik][self.test_key]
            
            for ykey in ykeys:
                fname = ik.replace("/","_") + f"-{ykey}.png"
                y = plotdata[ik][ykey]
                if y is None:
                    continue
                
                ix, iy = comm.clean_none_list(x,y)
                if len(ix) == 0:
                    continue
                if isinstance(ix[0],(tuple,list)):
                    ix = [str(i) for i in ix]
                
                # if iy is a list of list, we should rotate it
                if isinstance(iy[0],list):
                    iys = np.array(iy).T
                else:
                    iys = [iy]

                title = f"{ik}({ykey})"
                ytitle = self.plot_keys[ykey]
                if self.shift_idx is None or ykey in ["band_gap"]:
                    ytitle = ytitle.capitalize()
                else:
                    ytitle = "Relative " + ytitle
                        
                comm_plot.plot_line_point(ix, iys,xtitle=self.x_name,ytitle=ytitle,title=title, 
                                      fname=fname,figsize=(8,6),fontsize=16,grid=True)
                fnames.append(fname)
        
        return fnames        
                
                
            