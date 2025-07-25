import argparse,os
from . import abacustest,collectdata,outresult,prepare,report,remote,model
from importlib.metadata import version
__version__ = version("abacustest")

def parser():
    my_parser = argparse.ArgumentParser(description="abacustest")
    subparser = my_parser.add_subparsers(dest="command")
    
    abacustest.AbacusTestArgs(subparser.add_parser("submit",help="Submit a workflow"))
    abacustest.AbacusTestCheckStatusArgs(subparser.add_parser("status", help="Check the status of the dflow job"))
    abacustest.AbacusTestDownloadArgs(subparser.add_parser("download", help="Download the results of the dflow job"))
    collectdata.CollectDataArgs(subparser.add_parser("collectdata", help="Collect the specified key from the result files"))
    outresult.OutResultArgs(subparser.add_parser("outresult", help="Screen the result files"))
    prepare.PrepareArgs(subparser.add_parser("prepare", help="Prepare the inputs"))
    report.ReportArgs(subparser.add_parser("report", help="Generate a html report for the results"))
    report.ReportArgs(subparser.add_parser("remote", help = "Push the results to the remote server"))
    model.ModelArgs(subparser.add_parser("model", help = "Prepare and post-process the specified model"))
    return my_parser

def print_head():
    print("\n")
    print("--"*30)
    print("+        ABACUSTEST")
    print(f"+ version: {__version__}")
    print(f"+ GITHUB: https://github.com/pxlxingliang/abacus-test/tree/develop")
    print(f"+ BohriumApp: https://bohrium.dp.tech/apps/abacustest")
    print("--"*30+"\n")
    
    
def main():
    print_head()
    my_parser = parser()
    param = my_parser.parse_args()
    if param.command in ['submit','mlops-submit']:
        abacustest.abacustest(param)
    elif param.command == 'status':
        abacustest.checkstatus(param)
    elif param.command == 'download':
        abacustest.downloadflow(param)
    elif param.command == 'collectdata':
        collectdata.collectdata(param)
    elif param.command == 'outresult':
        outresult.outresult(param)
    elif param.command == 'prepare':
        prepare.PrepareInput(param)
    elif param.command == "report":
        report.Report(param)
    elif param.command == "remote":
        remote.Remote(param)
    elif param.command == "model":
        model.RunModel(param)
    else:
        print(my_parser.parse_args(['-h']))
    
if __name__ == '__main__':
    main()    