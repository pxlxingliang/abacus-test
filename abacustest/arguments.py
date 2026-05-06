import argparse

def _job_type(string):
    string = string.lower()
    type_map = {'0': 'abacus', '1': 'qe', '2': 'vasp'}
    if string in type_map:
        return type_map[string]
    valid_choices = list(type_map.values())
    if string not in valid_choices:
        raise argparse.ArgumentTypeError(f"Invalid job type: {string}. Valid options are: {valid_choices} or 0/1/2")
    return string

def SubmitArgs(parser):
    parser.description = "This script is used to run a testing"
    parser.add_argument('-p', '--param', type=str, default="job.json",help='the job setting file, default is job.json')
    parser.add_argument('-s', '--save', type=str, default=None,help='the folder where the results will be put in, default: result')
    parser.add_argument('--override', type=int, default=1,help="when the save folder exists, if override it. 0: no, 1: yes. ")
    parser.add_argument('--debug',nargs='?',type=int, const=1, default=0,help="debug mode for dflow")
    parser.add_argument('--download',default=1,type=int, help="if wait the finish of flow and download the results, 0: no, 1: yes. Default is 1")
    return parser

def CheckStatusArgs(parser):
    parser.description = "Check the status of the dflow job"
    parser.add_argument('-p', '--param', type=str, default="job.json",help='the file for bohrium account information, default is "job.json"')
    parser.add_argument("job_id", help="the job id of dflow")
    return parser

def DownloadArgs(parser):
    parser.description = "Download the results of the dflow job"
    parser.add_argument("job_id", default=None, nargs="?",help="the job id of your workflow")
    parser.add_argument('-p', '--param', type=str, default="job.json",help='the file for bohrium account information, default is "job.json"')
    parser.add_argument('-s', '--save', type=str, default=None,help='the folder where the results will be put in, default: result')
    return parser

def CollectDataArgs(parser):
    parser.description = "This script is used to collect some key values from the output of ABACUS/QE/VASP jobs"
    parser.add_argument('-j', '--jobs', default=["."], help='the path of jobs', action="extend",nargs="*")
    parser.add_argument('-t', '--type', type=_job_type, default='abacus', help='abacus/qe/vasp or 0/1/2. Default: abacus', choices=['abacus', 'qe', 'vasp'])
    parser.add_argument('-p', '--param', type=str, default=None,nargs="*", help='the parameter file, or parameter name.')
    parser.add_argument('-o', '--output', type=str, default="metrics.json",help='the file name to store the output results, default is "metrics.json"')
    parser.add_argument('-m', '--modules',help='add extra modules. Default only module \'job-type\' will be loaded, such as: \'abacus\' for abacus type. You can check all modules by --outparam', action="extend",nargs="*")
    parser.add_argument('--newmethods', help='the self-defined python modules, and shuold be format of import, such as "abc"(the file name is abc.py), "a.b.c" (teh file is a/b/c.py)', action="extend",nargs="*")
    parser.add_argument('--outparam', nargs='?',type=int, const=1, default=0,help='output the registed parameters, you can set the type by -t or --type to choose abacus/qe/vasp. 0: No, 1: yes')
    parser.add_argument('--ref', type=str, nargs='?',default=None,const="resultREF.json",help='A json file includes the reference value of some keys. Generally, get values of keys start with \"delta_\" require this file. Default is resultREF.json')
    return parser

def OutResultArgs(parser):  
    parser.description = "This script is used to output the summary of results"
    parser.add_argument('-r', '--result', type=str, help='the result file from collectdata, should be .json type',action="extend",nargs="*")
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type')
    parser.add_argument('-m', '--metrics', default=["ALL"], help='The metrics that needed to be shown. Use ALL to show all metrics, default is ALL.', action="extend",nargs="*")
    parser.add_argument('-o', '--output', type=str, help='output the selected metrics to a json file')
    parser.add_argument('--csv', type=str, help='if output the selected metrics to csv file')
    parser.add_argument('--prec', type=int, default=None, help='the precision of float number in pandas output, default is None')
    parser.add_argument("--list_sep", default=0, const=1, nargs='?', type=int, help="if print the list separately, default is 0")
    return parser

def PrepareArgs(parser):
    parser.description = "This script is used to prepare the INPUTS OF ABACUS JOB"
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type',required=True)
    parser.add_argument('-s', '--save', type=str,  default="abacustest",help='where to store the inputs, default is abacustest ')
    parser.add_argument('--nolink',  nargs='?',type=int, const=1, default=0,help='if link the files in the example folder, default is 0')
    return parser

def ReportArgs(parser):  
    parser.description = "Read metrics.json and generate the report"
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type',required=True)
    parser.add_argument('-o', '--output', type=str,  default="abacustest.html",help='The output file name, default is abacustest.html')
    return parser


def RemoteArgs(parser):  
    parser.description = '''push the local files to remote, like github
    Need to specify the parameter file, which should be .json type, and the content should be like:
    "config":{
        "github_username": "",
        "github_email": "",
        "github_token": ""
    },
    "add_date": "bda/dates",
    "remote":{
        "files":{
            "file1_local": "file1_remote",
            "file2_local": "file2_remote",
        },
        "github":{
            "repo": "",
            "branch": "",
            "commit_msg": ""
        }
    }
    the key of "files" is the local file/folder, and the value is the remote file/folder.
    '''
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type',required=True)
    return parser

def SkillsArgs(parser):
    parser.description = "Show the path of abacustest skills"
    return parser
