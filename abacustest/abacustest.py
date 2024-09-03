#!/usr/bin/env python3
import os,sys,argparse
sys.path.append(os.path.split(__file__)[0])
from .myflow import flow,globV,comm
from dflow.python import upload_packages
#import numpy,pandas
#upload_packages.append(numpy.__path__[0])
#upload_packages.append(pandas.__path__[0])
#upload_packages.append(pymatgen.__path__[0])
#upload_packages.append(os.path.join(os.path.split(__file__)[0],'myflow'))
#upload_packages.append(os.path.join(os.path.split(__file__)[0],'lib_collectdata'))
upload_packages.append(os.path.split(__file__)[0])

def AbacusTestArgs(parser):
    parser.description = "This script is used to run a testing"
    parser.add_argument('-p', '--param', type=str, default="job.json",help='the job setting file, default is job.json')
    parser.add_argument('-s', '--save', type=str, default=None,help='the folder where the results will be put in, default: result')
    parser.add_argument('--override', type=int, default=1,help="when the save folder exists, if override it. 0: no, 1: yes. ")
    parser.add_argument('--debug',nargs='?',type=int, const=1, default=0,help="debug mode for dflow")
    parser.add_argument('--download',default=1,type=int, help="if wait the finish of flow and download the results, 0: no, 1: yes. Default is 1")
    return parser

def AbacusTestCheckStatusArgs(parser):
    parser.description = "Check the status of the dflow job"
    parser.add_argument('-p', '--param', type=str, default="job.json",help='the file for bohrium account information, default is "job.json"')
    parser.add_argument("job_id", help="the job id of dflow")
    return parser

def AbacusTestDownloadArgs(parser):
    parser.description = "Download the results of the dflow job"
    parser.add_argument("job_id", help="the job id of your workflow")
    parser.add_argument('-p', '--param', type=str, default="job.json",help='the file for bohrium account information, default is "job.json"')
    parser.add_argument('-s', '--save', type=str, default=None,help='the folder where the results will be put in, default: result')
    return parser

def abacustest(param):
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
    parser = argparse.ArgumentParser(description="abacustest")
    subparser = parser.add_subparsers(dest="command")
    AbacusTestArgs(subparser.add_parser("submit",help="Submit a workflow"))
    AbacusTestArgs(subparser.add_parser("mlops-submit",help="Just submit a workflow and will not wait for the result"))
    AbacusTestCheckStatusArgs(subparser.add_parser("status", help="Check the status of the dflow job"))
    parser = parser.parse_args()
    if parser.command == "submit":
        abacustest(parser)
    elif parser.command == "status":
        checkstatus(parser)
    elif parser.command == "download":
        downloadflow(parser)

if __name__ == '__main__':
    main()
