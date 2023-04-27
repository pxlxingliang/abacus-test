from dp.launching.cli import to_runner,SubParser,run_sp_and_exit
from dp.launching.typing.basic import BaseModel, Int, String, Float,List,Optional,Union,Dict
from dp.launching.cli import to_runner, default_minimal_exception_handler
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
    BohriumProjectId
)
from enum import Enum

class ConfigSet(BaseModel):
    #Bohrium config
    Config_lbg_username:   BohriumUsername
    Config_lbg_password:   BohriumPassword
    Config_project_id:     BohriumProjectId

    #dflow set
    Config_config_host: DflowArgoAPIServer
    Config_s3_config_endpoint: DflowStorageEndpoint
    Config_config_k8s_api_server: DflowK8sAPIServer
    Config_config_token: DflowAccessToken 

class AbacusMetricEnum(String,Enum):
    AbacusMetric_metric1 = 'version'
    AbacusMetric_metric2 = 'ncore'
    AbacusMetric_metric3 = 'normal_end'
#    AbacusMetric_metric4 = 'INPUT'
#    AbacusMetric_metric5 = 'kpt'
    AbacusMetric_metric6 = 'nbands'
    AbacusMetric_metric7 = 'converge'
    AbacusMetric_metric8 = 'total_mag'
    AbacusMetric_metric9 = 'absolute_mag'
    AbacusMetric_metric10 = 'nkstot'
    AbacusMetric_metric11 = 'ibzk'
    AbacusMetric_metric12 = 'natom'
    AbacusMetric_metric13 = 'nelec'
    AbacusMetric_metric14 = 'energy'
    AbacusMetric_metric15 = 'volume'
#    AbacusMetric_metric16 = 'fft_grid'
    AbacusMetric_metric17 = 'efermi'
    AbacusMetric_metric18 = 'energy_per_atom'
#    AbacusMetric_metric19 = 'stress'
#    AbacusMetric_metric20 = 'force'
    AbacusMetric_metric21 = 'band_gap'
    AbacusMetric_metric22 = 'total_time'
    AbacusMetric_metric23 = 'stress_time'
    AbacusMetric_metric24 = 'force_time'
    AbacusMetric_metric25 = 'scf_time'
    AbacusMetric_metric26 = 'scf_time_each_step'
    AbacusMetric_metric27 = 'step1_time'
    AbacusMetric_metric28 = 'scf_steps'
    AbacusMetric_metric29 = 'atom_mag'
#    AbacusMetric_metric30 = 'drho'
    AbacusMetric_metric31 = 'lattice_constant'
    AbacusMetric_metric32 = 'cell'
#    AbacusMetric_metric33 = 'coordinate'
    AbacusMetric_metric34 = 'element_list'
    AbacusMetric_metric35 = 'atomlabel_list'
#    AbacusMetric_metric36 = 'delta_energy'
#    AbacusMetric_metric37 = 'delta_energyPerAtom'
    AbacusMetric_metric38 = 'relax_converge'
    AbacusMetric_metric39 = 'relax_steps'


class SuperMetricMethodEnum(String,Enum):
    SuperMetricMethod_method1 = "iGM"
    SuperMetricMethod_method2 = "GM"
    SuperMetricMethod_method3 = "TrueRatio"
    SuperMetricMethod_method4 = "MEAN"

'''
class RunSetUpload(BaseModel):
    uploaddatahub_project: String =  Field(default="",description="Please set the project where dataset will be uploaded to. If not set, will not upload")
    uploaddatahub_datasetname: String =  Field(default="",description="Please set the path in project where dataset will be uploaded to. If not set, will not upload")
    uploaddatahub_properties: Dict[String,String] =  Field(default={},description="")
    uploaddatahub_tags: List[String] =  Field(default=[],description="")
    uploadtracking_test_name: String =  Field(default="",description="Please set the name of test. If not set, will not upload metrics to tracking")
    uploadtracking_experiment_name: String =  Field(default="",description="Please set the name of experiment. If not set, will not upload metrics to tracking")
    uploadtracking_tags: List[String] =  Field(default=[],description="")
'''