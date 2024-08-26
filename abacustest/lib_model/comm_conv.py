'''
some comm functions for convegence test of kspacing and ecutwfc
'''
import glob,os,json
from abacustest.lib_collectdata.collectdata import RESULT
import numpy as np

def collect_metrics(jobs):
    allmetrics = {}
    for job in jobs:
            for ijob in glob.glob(os.path.join(job,"*")):
                if os.path.isdir(ijob):
                    results = RESULT(path=ijob,fmt="abacus")
                    if results["INPUT"].get("kspacing",None) is None:
                        kspacing = None
                    else:
                        kspacing = results["INPUT"]["kspacing"]
                        if isinstance(kspacing,str):
                            ksplit = [float(ik) for ik in kspacing.split()]
                            if len(ksplit) == 1:
                                kspacing = ksplit[0]
                            elif len(ksplit) == 3:
                                if ksplit[0] == ksplit[1] == ksplit[2]:
                                    kspacing = ksplit[0]
                    ecutwfc = results["INPUT"].get("ecutwfc",None)
                    
                    allmetrics[ijob] = {
                        "ecutwfc": ecutwfc,
                        "kspacing": kspacing,
                        "energy_per_atom": results["energy_per_atom"],
                        "energy": results["energy"],
                        "natom": results["natom"],
                        "normal_end": results["normal_end"],
                        "converge": results["converge"],
                        "scf_steps": results["scf_steps"],
                        "total_time": results["total_time"],
                        "scf_time": results["scf_time"],
                        "force": results["force"],
                        "stress": results["stress"],
                        "atom_mag": None if results["ds_mag"] is None else np.array(results["ds_mag"]).flatten().tolist(),
                        "mag_force": None if results["ds_mag_force"] is None else np.array(results["ds_mag_force"]).flatten().tolist(),
                        "version": results["version"],
                        "ks_solver": None if not isinstance(results["INPUT"],dict) else results["INPUT"].get("ks_solver",None),
                        "scf_time/step": None if not results["scf_time"] or not results["scf_steps"] else results["scf_time"] / results["scf_steps"]
                    }
    return allmetrics

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
    def __init__(self, jobs, metric=None):
        self.jobs = jobs
        self.metric = metric
    
    def run(self):
        if self.metric:
            allmetrics = json.load(open(self.metric))
        else:
            allmetrics = self._collect_metrics()
            json.dump(allmetrics, open("metrics.json", "w"), indent=4)
        
        plotdata = self._metrics2plotdata(allmetrics)
        fnames = self._plot_ecutwfc(plotdata)
        # write a supermetrics.json
        json.dump(plotdata, open("result.json", "w"), indent=4)
        json.dump({k:{"type": "image", "file": v}  for k,v in fnames.items()}, open("supermetrics.json","w"),indent=4)
    
    def _collect_metrics(self):
        jobs = self.jobs
        if len(jobs) == 0:
            print(
                "Warning: have not set the jobs -j. Try to find the results in current folder."
            )
            jobs = ["."]
        return comm_conv.collect_metrics(jobs)

    def _metrics2plotdata(self, allmetrics):
        # transfer the allmetrics to plotdata, which is a dict where key is example name 
        # and value is a dict of ecutwfc, energy_per_atom, and subfolder.
        plotdata = {}
        for i, iv in allmetrics.items():
            example_name = os.path.dirname(i.rstrip("/"))
            if example_name not in plotdata:
                plotdata[example_name] = {
                    "subfolder": [],
                    "ecutwfc": [],
                    "energy_per_atom": [],
                    "force": [],
                    "stress": [],
                    "atom_mag": [],
                    "mag_force": [],
                }
            ecutwfc = iv.get("ecutwfc", None)
            if "energy_per_atom" in iv:
                energy_per_atom = iv.get("energy_per_atom", None)
            elif "energy" in iv and "natom" in iv:
                try:
                    energy_per_atom = iv["energy"] / iv["natom"]
                except:
                    print(f"Warning: {i} has no energy_per_atom, and energy ({iv['energy']}) can not devided by natom ({iv['natom']}).")
                    energy_per_atom = None
            else:
                energy_per_atom = None
                print(f"Warning: {i} has no energy_per_atom or energy and natom.")
            
            if ecutwfc ==None or energy_per_atom == None:
                print(f"Warning: {i} has no ecutwfc or energy_per_atom.")
                continue
            plotdata[example_name]["subfolder"].append(os.path.basename(i.rstrip("/")))
            plotdata[example_name]["ecutwfc"].append(iv.get("ecutwfc"))
            plotdata[example_name]["energy_per_atom"].append(energy_per_atom)
            plotdata[example_name]["force"].append(iv.get("force"))
            plotdata[example_name]["stress"].append(iv.get("stress"))
            plotdata[example_name]["atom_mag"].append(iv.get("atom_mag"))
            plotdata[example_name]["mag_force"].append(iv.get("mag_force"))
            

        # sort the subfolder/energy_per_atom by ecutwfc
        allkeys = list(plotdata.keys())
        allkeys.sort()
        sorteddata = {}
        for ik,iv in plotdata.items():
            if len(iv["ecutwfc"]) == 0:
                print(f"Warning: {ik} has no ecutwfc.")
                continue
            idx = sorted(range(len(iv["ecutwfc"])), key=lambda k: iv["ecutwfc"][k])
            sorteddata[ik] = {
                "subfolder": [iv["subfolder"][i] for i in idx],
                "ecutwfc": [iv["ecutwfc"][i] for i in idx],
                "energy_per_atom": [(iv["energy_per_atom"][i] - iv["energy_per_atom"][idx[-1]])*1000  for i in idx], # shift energy to the last one
            }
            for ikk in ["force","stress","atom_mag","mag_force"]:
                if len(iv[ikk]) == 0:
                    continue
                idata = comm_conv.shift_data(iv[ikk],base_index=idx[-1])
                if idata is None:
                    continue
                # rotate the data
                sorteddata[ik][ikk] = [[idata[j][i] for j in idx] for i in range(len(idata[0]))]
        return sorteddata
    
    def _plot_ecutwfc(self, plotdata):
        allkeys = list(plotdata.keys())
        allkeys.sort()
        
        fnames = {}
        # plot each example in a subplot
        for i, ik in enumerate(allkeys):
            fname = ik.replace("/","_")+".png"
            comm_plot.plot_line_point(plotdata[ik]["ecutwfc"], [plotdata[ik]["energy_per_atom"]],xtitle="Kspacing (1/bohr)",ytitle="Energy (meV/atom)",title=ik, 
                                      fname=fname,figsize=(8,6),fontsize=16,grid=True)
            fnames[ik] = fname
            
            for ikey,iunit in zip(["force", "stress","atom_mag", "mag_force"],
                                  ["Force (eV/Ang)", "Stress (kbar)","Atom magnetic moment (uB)","Magnetic force (uB/A)"]):
                if ikey not in plotdata[ik] or len(plotdata[ik][ikey]) == 0:
                    continue
                fname = ik.replace("/","_") + "-" + ikey + ".png"
                comm_plot.plot_line_point(plotdata[ik]["ecutwfc"], plotdata[ik][ikey],xtitle="Kspacing (1/bohr)",ytitle=iunit,title=ik+"-"+ikey, 
                                      fname=fname,figsize=(8,6),fontsize=16,grid=True)
                fnames[ik + "-" + ikey] = fname
        return fnames   