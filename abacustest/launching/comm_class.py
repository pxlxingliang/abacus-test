from dp.launching.typing.basic import BaseModel, Int, String, Float, List, Optional, Union, Dict
from dp.launching.typing import InputFilePath, OutputDirectory
from dp.launching.typing import (
    BaseModel,
    Set,
    Boolean,
    Field,
    DflowAccessToken,
    DflowArgoAPIServer,
    DflowK8sAPIServer,
    DflowStorageEndpoint,
    DflowStorageRepository,
    BohriumMachineType,
    BohriumImage,
    BohriumPlatform,
    BohriumJobType,
    BohriumUsername,
    BohriumPassword,
    BohriumProjectId,
    BohriumTicket,
    BenchmarkLabels,
    BenchmarkTags,
    DflowLabels 
)
from typing import Literal

import re,os

class OutputSet(BaseModel):
    IO_output_path: OutputDirectory = Field(default="./outputs")

class ConfigSet(BaseModel):
    #Bohrium config
    Config_bohrium_username:   BohriumUsername
#    Config_bohrium_password:   BohriumPassword
    Config_bohrium_project_id:     BohriumProjectId
    Config_bohrium_ticket: BohriumTicket

    #dflow set
    Config_dflow_host: DflowArgoAPIServer
    Config_dflow_s3_config_endpoint: DflowStorageEndpoint
    Config_dflow_k8s_api_server: DflowK8sAPIServer
    Config_dflow_token: DflowAccessToken
    Config_dflow_labels: DflowLabels

class ConfigSetGithub(BaseModel):
    Config_github_username: String = Field(default="",title="github_username",description="github username")
    Config_github_email: String = Field(default="",title="github_email",description="github email")
    Config_github_token: String = Field(default="",title="github_token",description="github token")

class ImageDispatcher(BaseModel):
    type: Literal["Use dp-dispatcher"]

    '''
    "machine_dict": {
                    "remote_root": "xxxx",
                    "remote_profile": {
                        "hostname": "xxxx",
                        "username": "xxx",
                        "password": "xxx",
                        "port": 22
                    }
    },
    "resources_dict": {
        "number_node": 1,
        "cpu_per_node": 8,
        "gpu_per_node": 1,
        "queue_name": "normal"
    }
    '''
    host: String = Field(default="",title="host",description="host of remote server")
    username: String = Field(default="",title="username",description="username of remote server")
    password: String = Field(default="",title="password",description="password of remote server")
    port: Int = Field(default=22,title="port",description="port of remote server")
    remote_root: String = Field(default="",title="remote_root",description="an absolute path on remote server to run the calculation")
    number_node: Int = Field(default=1,title="number_node",description="number of nodes to run the calculation")
    cpu_per_node: Int = Field(default=8,title="cpu_per_node",description="number of cpus per node to run the calculation")
    gpu_per_node: Int = Field(default=0,title="gpu_per_node",description="number of gpus per node to run the calculation")
    queue_name: String = Field(default="",title="queue_name",description="queue name to run the calculation")

    @classmethod
    def construct_dispatcher(cls,opts):
        param = {
            "dispatcher":{
                "machine_dict":{
                    "remote_root":opts.remote_root,
                    "remote_profile":{
                        "hostname":opts.host,
                        "username":opts.username,
                        "password":opts.password,
                        "port":opts.port
                    }
                },
                "resources_dict":{
                    "number_node":opts.number_node,
                    "cpu_per_node":opts.cpu_per_node,
                    "gpu_per_node":opts.gpu_per_node,
                    "queue_name":opts.queue_name
                }
            }
        }
        if opts.gpu_per_node < 1:
            del param["dispatcher"]["resources_dict"]["gpu_per_node"]
        if opts.queue_name.strip() == "":
            del param["dispatcher"]["resources_dict"]["queue_name"]
        return param
    
class NgroupSet(BaseModel):
    ngroup: Int = Field(
        default=0, description="Number of groups to run in parallel. If not set, all examples will be run in parallel.", ge=0)

class TrackingSet(BaseModel):
    Tracking_metrics: Boolean = Field(
        default=False, description="If tracking, will display historical metrics values based on test and experience")
    Tracking_token: String = Field(
        default=None, description="If want to track metrics, please enter your token to access AIM")
    Tracking_tags: BenchmarkTags
    #tags = [f"benchmark-application-{application.name}", f"benchmark-version-{job.version}", f"benchmark-schedule-{job.properties.get('source_name', 'none')}",
    #        f"benchmark-job-{job.name}"]

    @classmethod
    def parse_obj(cls, opts):
        if opts.Tracking_metrics:
            default_tags = opts.Tracking_tags
            schedule = "_".join(
                re.split("-", default_tags[2], 2)[-1].strip().split())
            application_name = re.split("-", default_tags[0], 2)[-1]
            job_name = re.split("-", default_tags[3], 2)[-1]
            if opts.Tracking_token != None and opts.Tracking_token.strip() != "":
                token = opts.Tracking_token
            else:
                token = os.environ.get("AIM_ACCESS_TOKEN")

            return {
                "name": schedule + "." + job_name,
                "experiment": application_name + "/benchmark",
                "tags": default_tags,
                "token": token
            }
        else:
            return None


class myLog:
    def __init__(self):
        self.logs = ""
    
    def __add__(self, other):
        if isinstance(other, myLog):
            self.logs += other.logs
        else:
            self.logs += str(other)
        return self

    def iprint(self, mess, *args):
        allmess = " ".join([str(mess)]+[str(i) for i in args])
        print(allmess)
        self.logs += allmess + "\n"

    def write(self, filename):
        with open(filename, 'w') as f1:
            f1.write(self.logs)
