from ..model import Model
import os, glob, json
from . import comm
from abacustest.lib_collectdata.collectdata import RESULT

SETTING_TMP = {
    "save_path": "results",
    "bohrium_group_name": "convergence-ecutwfc",
    "prepare": {"example_template": [], "mix_input": {"ecutwfc": []}},
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
                    "ecutwfc": "{INPUT}['ecutwfc']",
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
        return "Do a convergence test of the ecutwfc"

    @staticmethod
    def model_name():  # type: ignore
        return "convecutwfc"

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
            help="the value of ecutwfc. Default is 50 to 100 in step 10",
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

        value = params.value if len(params.value) > 0 else list(range(50, 101, 10))
        real_jobs = comm.get_job_list(jobs)
        print("The jobs are:", ", ".join(real_jobs))

        setting = SETTING_TMP
        setting["prepare"]["example_template"] = real_jobs
        setting["prepare"]["mix_input"]["ecutwfc"] = value
        setting["post_dft"]["metrics"]["path"] = [
            os.path.join(i, "*") for i in real_jobs
        ]
        comm.dump_setting(setting)
        
        comm.doc_after_prepare("convergence test of ecutwfc", real_jobs, ["setting.json"])
        print("After finish the calculation, you can run below command to do the postprocess:")
        print(f"    abacustest model {self.model_name()} post -j {' '.join(real_jobs)} -r result.json\n")
        
        if params.run:
            bash_script = "ecutwfc_convergence.sh"
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
            help="the path of jobs. There should have several subfolders with different ecutwfc for each job.",
            action="extend",
            nargs="*",
        )
        group.add_argument(
            "-m",
            "--metric",
            default=None,
            help="the metrics.json file. Key is the examplename/subfolder and value is a dict of metrics value containing ecutwfc, energy_per_atom (or energy and natom).",
        )
        parser.add_argument(
            "-o",
            "--output",
            default="ecutwfc.png",
            help="the output picture name, default is ecutwfc.png",
        )
        parser.add_argument(
            "-s",
            "--save",
            default="metrics.json",
            help="Save the metrics to a file when defining the jobs. Default is metrics.json",
        )
        parser.add_argument(
            "-r",
            "--result",
            default="result.json",
            help="Save the data of the plot to a json file. Default is result.json.",
        )
    
    def _collect_metrics(self, params):
        jobs = params.jobs
        if len(jobs) == 0:
            print(
                "Warning: have not set the jobs -j. Try to find the results in current folder."
            )
            jobs = ["."]
        allmetrics = {}
        for job in jobs:
            for ijob in glob.glob(os.path.join(job,"*")):
                if os.path.isdir(ijob):
                    results = RESULT(path=ijob,fmt="abacus")
                    
                    allmetrics[ijob] = {
                        "ecutwfc": None if not isinstance(results["INPUT"],dict) else results["INPUT"].get("ecutwfc",None),
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
                        "version": results["version"],
                        "ks_solver": None if not isinstance(results["INPUT"],dict) else results["INPUT"].get("ks_solver",None),
                        "scf_time/step": None if not results["scf_time"] or not results["scf_steps"] else results["scf_time"] / results["scf_steps"]
                    }
        return allmetrics
    
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
                "energy_per_atom": [iv["energy_per_atom"][i] - iv["energy_per_atom"][idx[-1]]  for i in idx], # shift energy to the last one
            }
        return sorteddata
    
    def _plot_ecutwfc(self, plotdata, savef):
        import matplotlib.pyplot as plt
        allkeys = list(plotdata.keys())
        allkeys.sort()
        
        # plot each example in a subplot
        fig, ax = plt.subplots(len(allkeys), 1, figsize=(8, 6*len(allkeys)))
        for i, ik in enumerate(allkeys):
            ax_tmp = ax[i] if len(allkeys) > 1 else ax
            ax_tmp.plot(plotdata[ik]["ecutwfc"], plotdata[ik]["energy_per_atom"], "o-")
            ax_tmp.set_title(ik)
            ax_tmp.set_xlabel("ecutwfc")
            ax_tmp.set_ylabel("energy_per_atom (eV)")
            ax_tmp.grid()
        plt.tight_layout()
        plt.savefig(savef)    
        
    def run_postprocess(self, params):
        if params.metric:
            allmetrics = json.load(open(params.metric))
        else:
            allmetrics = self._collect_metrics(params)
            json.dump(allmetrics, open("metrics.json", "w"), indent=4)
        
        plotdata = self._metrics2plotdata(allmetrics)
        self._plot_ecutwfc(plotdata, params.output)
        # write a supermetrics.json
        json.dump({"ecutwfc_converge":{"type": "image", "file": params.output}}, open("supermetrics.json","w"),indent=4)      
