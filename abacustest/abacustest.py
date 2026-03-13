#!/usr/bin/env python3
import os,sys,argparse
sys.path.append(os.path.split(__file__)[0])
from .myflow import flow,globV,comm

def submitflow(param):
    from dflow.python import upload_packages
    upload_packages.append(os.path.split(__file__)[0])
    globV._init()
    flow.RunJobs(param)
    comm.printinfo("\nAll jobs are finished!!!")

def checkstatus(param):
    globV._init()
    status = flow.CheckStatus(param)
    print(status)

def downloadflow(param):
    globV._init()
    flow.DownloadFlow(param)

def main():
    from abacustest.arguments import SubmitArgs,CheckStatusArgs

    parser = argparse.ArgumentParser(description="abacustest")
    subparser = parser.add_subparsers(dest="command")
    SubmitArgs(subparser.add_parser("submit",help="Submit a workflow"))
    SubmitArgs(subparser.add_parser("mlops-submit",help="Just submit a workflow and will not wait for the result"))
    CheckStatusArgs(subparser.add_parser("status", help="Check the status of the dflow job"))
    parser = parser.parse_args()
    if parser.command == "submit":
        submitflow(parser)
    elif parser.command == "status":
        checkstatus(parser)
    elif parser.command == "download":
        downloadflow(parser)

if __name__ == '__main__':
    main()
