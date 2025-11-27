from ..model import Model
import os, glob, json
from . import comm,comm_conv,comm_plot
from abacustest.lib_collectdata.collectdata import RESULT
from abacustest.constant import RECOMMAND_IMAGE, RECOMMAND_COMMAND, RECOMMAND_MACHINE

SETTING_TMP = {
    "save_path": "results",
    "bohrium_group_name": "convergence-test",
    "prepare": {"example_template": [], "mix_input": {}},
    "run_dft": {
        "command": "OMP_NUM_THREADS=1 mpirun -n 16 abacus | tee out.log",
        "image": RECOMMAND_IMAGE,
        "bohrium": {
            "scass_type": RECOMMAND_MACHINE,
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
                "energy",
                "energy_per_atom",
                "force",
                "stress",
                "version",
                {
                    "scf_time/step": "{scf_time}/{scf_steps}",
                    "ks_solver": "{INPUT}['ks_solver']",
                    "ecutwfc": "{INPUT}['ecutwfc']",
                }
            ],
            "save_file": "metrics.json"
        },
    },
}


class ConvEcutwfc(Model):
    @staticmethod
    def description():  # type: ignore
        return "Do a convergence test"

    @staticmethod
    def model_name():  # type: ignore
        return "conv"

    @staticmethod
    def prepare_args(parser):
        parser.add_argument(
            "-j",
            "--jobs",
            help="the path of jobs",
            action="extend",
            nargs="*",
        )
        parser.add_argument(
            "-k",
            "--key",
            type=str,
            help="the key name to be tested, should be a INPUT para. Such as ecutwfc, kspacing, etc.",
        )
        parser.add_argument(
            "-v",
            "--value",
            help="the value to be tested",
            action="extend",
            nargs="*",
        )

    def run_prepare(self, params):
        real_jobs = comm.get_job_list(params.jobs)
        if len(real_jobs) == 0:
            print("No jobs found.")
            return
        print("The jobs are:", ", ".join(real_jobs))
        setting = SETTING_TMP
        setting["prepare"]["example_template"] = real_jobs
        setting["prepare"]["mix_input"][params.key] = params.value
        setting["post_dft"]["metrics"]["metrics_name"][-1][params.key] = f"{{INPUT}}['{params.key}']"
        
        comm.dump_setting(setting)
        comm.doc_after_prepare(f"convergence test of {params.key}", real_jobs, ["setting.json"])

    
    @staticmethod
    def postprocess_args(parser):
        #group = parser.add_mutually_exclusive_group()
        parser.add_argument(
            "-k",
            "--key",
            help="the key name to be tested, should be a INPUT para. Such as ecutwfc, kspacing, etc.",
        )
        parser.add_argument(
            "-j",
            "--jobs",
            help="the path of jobs. There should have several subfolders with different ecutwfc for each job.",
            action="extend",
            nargs="*",
        )
        parser.add_argument(
            "-m",
            "--metric",
            default=None,
            help="the metrics.json file. Key is the examplename/subfolder and value is a dict of metrics value containing kspacing, energy_per_atom (or energy and natom).",
        )
        parser.add_argument(
            "-e",
            "--extra",
            nargs="*",
            action="extend",
            help="the extra metrics that need to be ploted. The format is \"metrics_name:title of metrics\"",
        )
        
    def run_postprocess(self, params):
        extra_y = None
        if params.extra:
            extra_y = {}
            for e in params.extra:
                e = e.split(":")
                extra_y[e[0]] = e[1].strip()
        comm_conv.PostConv(test_key=params.key,jobs=params.jobs,metric_file=params.metric,extra_y=extra_y).run()