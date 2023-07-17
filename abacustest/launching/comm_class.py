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
    BenchmarkTags
)

import re,os

class OutputSet(BaseModel):
    IO_output_path: OutputDirectory = Field(default="./output")

class ConfigSet(BaseModel):
    #Bohrium config
    Config_lbg_username:   BohriumUsername
#    Config_lbg_password:   BohriumPassword
    Config_project_id:     BohriumProjectId
    Config_bohrium_ticket: BohriumTicket

    #dflow set
    Config_config_host: DflowArgoAPIServer
    Config_s3_config_endpoint: DflowStorageEndpoint
    Config_config_k8s_api_server: DflowK8sAPIServer
    Config_config_token: DflowAccessToken

    Config_dflow_labels: BenchmarkLabels


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

    def iprint(self, mess, *args):
        allmess = " ".join([str(mess)]+[str(i) for i in args])
        print(allmess)
        self.logs += allmess + "\n"

    def write(self, filename):
        with open(filename, 'w') as f1:
            f1.write(self.logs)
