from dp.launching.cli import to_runner, SubParser, run_sp_and_exit
from dp.launching.typing.basic import BaseModel, Int, String, Float, List, Optional, Union, Dict
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
    BohriumProjectId,
    BenchmarkLabels,
    BenchmarkTags
)
from enum import Enum
from typing import Literal
import re,copy

class AbacusMetricEnum(String, Enum):
    AbacusMetric_metric1 = 'version'
    AbacusMetric_metric2 = 'ncore'
    AbacusMetric_metric3 = 'normal_end'
    AbacusMetric_metric4 = 'INPUT:ks_solver'
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

class SuperMetricMethodEnum(String, Enum):
    SuperMetricMethod_method1 = "iGM"
    SuperMetricMethod_method2 = "GM"
    SuperMetricMethod_method3 = "TrueRatio"
    SuperMetricMethod_method4 = "MEAN"

class SuperMetrics(BaseModel):
    #name : String
    param_name: AbacusMetricEnum 
    method: SuperMetricMethodEnum
    #normalization: Boolean

class MetricsSet(BaseModel):
    metrics: Set[AbacusMetricEnum]

class SuperMetricsSet(BaseModel):
    super_metrics: List[SuperMetrics] = Field(default=[])

def convert_metrics(metrics_list):
    #the metrics in launching may have the format of KEY1:KEY2
    #this should be transfer to KEY1:{[KEY2]} in abacustest metrics
    #and transfer to KEY1/KEY2 in abacustest supermetrics
    new_metrics = []
    dict_tmp = {}
    for imetric in metrics_list:
        if ":" in imetric:
            allkeys = imetric.split(":")
            if allkeys[0] not in dict_tmp:
                dict_tmp[allkeys[0]] = []
            dict_tmp[allkeys[0]].append(allkeys[1])
        else:
            new_metrics.append(imetric)
    if dict_tmp:
        new_metrics.append(copy.deepcopy(dict_tmp))
    return new_metrics

def convert_supermetrics_metrics_name(imetric):
    if ":" in imetric:
        allkeys = imetric.split(":")
        return "%s/%s" % (allkeys[0].strip(),allkeys[1].strip())
    else:
        return imetric.strip()

def parse_metrics_set(metrics_set: MetricsSet, super_metrics_set: SuperMetricsSet):
    out_dict = {}
    metrics = list(metrics_set.metrics)
    has_super_metrics = False
    if len(super_metrics_set.super_metrics) > 0:
        has_super_metrics = True
        #complete metrics
        for i in super_metrics_set.super_metrics:
            if i.param_name not in metrics:
                metrics.append(i.param_name)

    #convert some special metrics (such as: KEY1:KEY2,..)
    metrics = convert_metrics(metrics)

    if len(metrics) == 0:
        return out_dict
    
    #set metrics
    out_dict["metrics"] = {
        "dft_type": "abacus",
        "metrics_name": metrics,
        "save_file": "metrics.json"
    }
    
    if has_super_metrics:
        out_dict["super_metrics"] = [{
            "save_file": "superMetrics.json",
            "result_file": ["metrics.json"],
            "metrics":[],
            "outparams":[]
        }]  
        for i in super_metrics_set.super_metrics:
            metric_name = convert_supermetrics_metrics_name(i.param_name)
            out_dict["super_metrics"][-1]["metrics"].append({
                "name": f"{i.method}({metric_name})" ,
                "param_name": metric_name,
                "method": i.method,
                "normalization": False})
            out_dict["super_metrics"][-1]["outparams"].append([metric_name, [metric_name], -1])
    
    return out_dict