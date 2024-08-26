from ..model import Model
import numpy as np
import os, glob, json
from . import comm,comm_plot,comm_conv
from abacustest.lib_collectdata.collectdata import RESULT

SETTING_TMP = {
    "save_path": "results",
    "bohrium_group_name": "convergence-kspacing",
    "prepare": {"example_template": [], "mix_input": {"kspacing": []}},
    "run_dft": {
        "command": "OMP_NUM_THREADS=1 mpirun -n 16 abacus | tee out.log",
        "image": "registry.dp.tech/deepmodeling/abacus-intel:latest",
        "bohrium": {
            "scass_type": "c32_m64_cpu",
            "job_type": "container",
            "platform": "ali",
        },
    },
    "post_dft": {
        "image": "registry.dp.tech/dptech/abacustest:latest",
        "metrics": {
            "dft_type": "abacus",
            "metrics_name": [
                "normal_end",
                "converge",
                "natom",
                "scf_steps",
                "total_time",
                "scf_time",
                {
                    "scf_time/step": "{scf_time}/{scf_steps}",
                    "ks_solver": "{INPUT}['ks_solver']",
                    "kspacing": "{INPUT}['kspacing']",
                },
                "energy",
                "energy_per_atom",
                "force",
                "stress",
                "version",
            ],
            "save_file": "metrics.json",
            "path": ["*/*"],
        },
    },
}


class ConvEcutwfc(Model):
    @staticmethod
    def description():  # type: ignore
        return "Do a convergence test of the kspacing"

    @staticmethod
    def model_name():  # type: ignore
        return "convkspacing"

    @staticmethod
    def prepare_args(parser):
        parser.add_argument(
            "-j",
            "--jobs",
            default=[],
            help="the path of jobs",
            action="extend",
            nargs="*",
        )
        parser.add_argument(
            "-v",
            "--value",
            default=[],
            type=float,
            help="the value of kspacing. Default is 0.08, 0.09, 0.1, 0.12, 0.15, 0.2",
            action="extend",
            nargs="*",
        )
        parser.add_argument(
            "-r",
            "--run",
            default=0,
            type=int,
            help="run the calculation after prepare. Default is 0, not run. 1 is run."
        )

    def run_prepare(self, params):
        jobs = params.jobs
        if len(params.jobs) == 0:
            print(
                "Warning: have not set the jobs -j. Will find all folders in the current directory."
            )
            jobs = ["*"]

        value = params.value if len(params.value) > 0 else [0.08, 0.09, 0.1, 0.12, 0.15, 0.2]
        real_jobs = comm.get_job_list(jobs)
        print("The jobs are:", ", ".join(real_jobs))

        setting = SETTING_TMP
        setting["prepare"]["example_template"] = real_jobs
        setting["prepare"]["mix_input"]["kspacing"] = value
        setting["post_dft"]["metrics"]["path"] = [
            os.path.join(i, "*") for i in real_jobs
        ]
        comm.dump_setting(setting)
        
        comm.doc_after_prepare("convergence test of kspacing", real_jobs, ["setting.json"])
        print("After finish the calculation, you can run below command to do the postprocess:")
        print(f"    abacustest model {self.model_name()} post -j {' '.join(real_jobs)}\n")
        
        if params.run:
            bash_script = "kspacing_convergence.sh"
            with open(bash_script,"w") as f:
                f.write("abacustest submit -p setting.json\n")
                f.write("cd results\n")
                f.write(f"abacustest model {self.model_name()} post -m metrics.json\n")
            os.system(f"bash {bash_script} &")

    
    @staticmethod
    def postprocess_args(parser):
        group = parser.add_mutually_exclusive_group()
        
        group.add_argument(
            "-j",
            "--jobs",
            default=[],
            help="the path of jobs. There should have several subfolders with different kspacing for each job.",
            action="extend",
            nargs="*",
        )
        group.add_argument(
            "-m",
            "--metric",
            default=None,
            help="the metrics.json file. Key is the examplename/subfolder and value is a dict of metrics value containing kspacing, energy_per_atom (or energy and natom).",
        )
        
    def run_postprocess(self, params):
        PostConvKspacing(params.jobs, params.metric).run()     


class PostConvKspacing:
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
        fnames = self._plot_kspacing(plotdata)
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
        # and value is a dict of kspacing, energy_per_atom, and subfolder.
        plotdata = {}
        for i, iv in allmetrics.items():
            example_name = os.path.dirname(i.rstrip("/"))
            if example_name not in plotdata:
                plotdata[example_name] = {
                    "subfolder": [],
                    "kspacing": [],
                    "energy_per_atom": [],
                    "force": [],
                    "stress": [],
                    "atom_mag": [],
                    "mag_force": [],
                }
            kspacing = iv.get("kspacing", None)
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
            
            if kspacing ==None or energy_per_atom == None:
                print(f"Warning: {i} has no kspacing or energy_per_atom.")
                continue
            plotdata[example_name]["subfolder"].append(os.path.basename(i.rstrip("/")))
            plotdata[example_name]["kspacing"].append(iv.get("kspacing"))
            plotdata[example_name]["energy_per_atom"].append(energy_per_atom)
            plotdata[example_name]["force"].append(iv.get("force"))
            plotdata[example_name]["stress"].append(iv.get("stress"))
            plotdata[example_name]["atom_mag"].append(iv.get("atom_mag"))
            plotdata[example_name]["mag_force"].append(iv.get("mag_force"))
            

        # sort the subfolder/energy_per_atom by kspacing
        allkeys = list(plotdata.keys())
        allkeys.sort()
        sorteddata = {}
        for ik,iv in plotdata.items():
            if len(iv["kspacing"]) == 0:
                print(f"Warning: {ik} has no kspacing.")
                continue
            idx = sorted(range(len(iv["kspacing"])), key=lambda k: iv["kspacing"][k])
            sorteddata[ik] = {
                "subfolder": [iv["subfolder"][i] for i in idx],
                "kspacing": [iv["kspacing"][i] for i in idx],
                "energy_per_atom": [(iv["energy_per_atom"][i] - iv["energy_per_atom"][idx[0]])*1000  for i in idx], # shift energy to the last one
            }
            for ikk in ["force","stress","atom_mag","mag_force"]:
                if len(iv[ikk]) == 0:
                    continue
                idata = comm_conv.shift_data(iv[ikk],base_index=idx[0])
                if idata is None:
                    continue
                # rotate the data
                sorteddata[ik][ikk] = [[idata[j][i] for j in idx] for i in range(len(idata[0]))]
        return sorteddata
    
    def _plot_kspacing(self, plotdata):
        allkeys = list(plotdata.keys())
        allkeys.sort()
        
        fnames = {}
        # plot each example in a subplot
        for i, ik in enumerate(allkeys):
            fname = ik.replace("/","_")+".png"
            comm_plot.plot_line_point(plotdata[ik]["kspacing"], [plotdata[ik]["energy_per_atom"]],xtitle="Kspacing (1/bohr)",ytitle="Energy (meV/atom)",title=ik, 
                                      fname=fname,figsize=(8,6),fontsize=16,grid=True)
            fnames[ik] = fname
            
            for ikey,iunit in zip(["force", "stress","atom_mag", "mag_force"],
                                  ["Force (eV/Ang)", "Stress (kbar)","Atom magnetic moment (uB)","Magnetic force (uB/A)"]):
                if ikey not in plotdata[ik] or len(plotdata[ik][ikey]) == 0:
                    continue
                fname = ik.replace("/","_") + "-" + ikey + ".png"
                comm_plot.plot_line_point(plotdata[ik]["kspacing"], plotdata[ik][ikey],xtitle="Kspacing (1/bohr)",ytitle=iunit,title=ik+"-"+ikey, 
                                      fname=fname,figsize=(8,6),fontsize=16,grid=True)
                fnames[ik + "-" + ikey] = fname
        return fnames