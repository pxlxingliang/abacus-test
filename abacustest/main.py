import argparse
from . import abacustest,collectdata,outresult,prepare,report,remote

def parser():
    my_parser = argparse.ArgumentParser(description="abacustest")
    subparser = my_parser.add_subparsers(dest="command")
    
    abacustest.AbacusTestArgs(subparser.add_parser("submit"))
    abacustest.AbacusTestArgs(subparser.add_parser("mlops-submit"))
    abacustest.AbacusTestCheckStatusArgs(subparser.add_parser("status"))
    collectdata.CollectDataArgs(subparser.add_parser("collectdata"))
    outresult.OutResultArgs(subparser.add_parser("outresult"))
    prepare.PrepareArgs(subparser.add_parser("prepare"))
    report.ReportArgs(subparser.add_parser("report"))
    report.ReportArgs(subparser.add_parser("remote"))
    
    
    return my_parser
    
def main():
    my_parser = parser()
    param = my_parser.parse_args()
    if param.command in ['submit','mlops-submit']:
        abacustest.abacustest(param)
    elif param.command == 'status':
        abacustest.checkstatus(param)
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
    else:
        print(my_parser.parse_args(['-h']))
    
if __name__ == '__main__':
    main()    