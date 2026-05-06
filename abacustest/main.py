import argparse
import sys

def print_head():
    from importlib.metadata import version
    __version__ = version("abacustest")
    print("\n")
    print("--"*30)
    print("+        ABACUSTEST")
    print(f"+ version: {__version__}")
    print(f"+ GITHUB: https://github.com/pxlxingliang/abacus-test/tree/develop")
    print(f"+ BohriumApp: https://bohrium.dp.tech/apps/abacustest")
    print("--"*30+"\n")

COMMANDS = {
    "submit": {
        "help": "Submit a workflow",
        "args_func": "SubmitArgs",
        "module": "abacustest",
        "handler": "submitflow",
    },
    "mlops-submit": {
        "help": "Just submit a workflow",
        "args_func": "SubmitArgs",
        "module": "abacustest",
        "handler": "submitflow",
    },
    "status": {
        "help": "Check the status of the dflow job",
        "args_func": "CheckStatusArgs",
        "module": "abacustest",
        "handler": "checkstatus",
    },
    "download": {
        "help": "Download the results of the dflow job",
        "args_func": "DownloadArgs",
        "module": "abacustest",
        "handler": "downloadflow",
    },
    "collectdata": {
        "help": "Collect the specified key from the result files",
        "args_func": "CollectDataArgs",
        "module": "collectdata",
        "handler": "collectdata.collectdata",
    },
    "outresult": {
        "help": "Screen the result files",
        "args_func": "OutResultArgs",
        "module": "outresult",
        "handler": "outresult.outresult",
    },
    "prepare": {
        "help": "Prepare the inputs",
        "args_func": "PrepareArgs",
        "module": "prepare",
        "handler": "prepare.PrepareInput",
    },
    "report": {
        "help": "Generate a html report for the results",
        "args_func": "ReportArgs",
        "module": "report",
        "handler": "report.Report",
    },
    "remote": {
        "help": "Push the results to the remote server",
        "args_func": "RemoteArgs",
        "module": "remote",
        "handler": "remote.Remote",
    },
    "model": {
        "help": "Prepare and post-process the specified model",
        "module": "model",
        "handler": "model.RunModel",
    },
    "skills": {
        "help": "Show the path of abacustest skills",
        "module": "skills",
        "handler": "show_skills",
    },
}

def main():
    print_head()
    
    my_parser = argparse.ArgumentParser(description="abacustest")
    subparser = my_parser.add_subparsers(dest="command")
    
    if len(sys.argv) > 1 and sys.argv[1] == "model":
        from . import model as model_mod
        
        model_parser = subparser.add_parser("model", help=COMMANDS["model"]["help"])
        
        if len(sys.argv) == 2 or (len(sys.argv) > 2 and sys.argv[2] in ["-h", "--help"]) or \
           (len(sys.argv) > 2 and sys.argv[2] not in model_mod.MODEL_ARGS):
            model_mod.ModelArgs(model_parser, list_all_models=True)
        else:
            model_subcommand = sys.argv[2] if len(sys.argv) > 2 else None
            model_mod.ModelArgs(model_parser, list_all_models=False, model_subcommand=model_subcommand)
    else:
        from . import arguments as args_mod
        for cmd, cfg in COMMANDS.items():
            if cmd in ["model", "skills"]:
                subparser.add_parser(cmd, help=cfg["help"])
            else:
                args_func = getattr(args_mod, cfg["args_func"])
                args_func(subparser.add_parser(cmd, help=cfg["help"]))
    
    param = my_parser.parse_args()
    
    if not param.command:
        my_parser.print_help()
        return
    
    cfg = COMMANDS[param.command]
    
    if param.command == "model":
        if not hasattr(param, 'model') or param.model is None:
            model_parser.print_help()
            return
        from . import model as model_mod
        model_mod.RunModel(param)
    else:
        mod = __import__(f"abacustest.{cfg['module']}", fromlist=[cfg['handler'].split('.')[-1]])
        handler = getattr(mod, cfg['handler'].split('.')[-1])
        handler(param)

if __name__ == '__main__':
    main()
