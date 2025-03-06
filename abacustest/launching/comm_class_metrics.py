from dp.launching.typing.basic import BaseModel, Int, String, Float, List, Optional, Union, Dict
from dp.launching.typing import (
    BaseModel,
    Set,
    Boolean,
    Field,
    )
from enum import Enum
from typing import Literal
import re,copy

class AbacusMetricEnum(String, Enum):
    AbacusMetric_metric1 = 'version'
    AbacusMetric_metric2 = 'ncore'
    AbacusMetric_metric3 = 'normal_end'
    AbacusMetric_metric4_1 = 'INPUT:ks_solver'
    AbacusMetric_metric4_2 = 'INPUT:ecutwfc'
    AbacusMetric_metric4_3 = 'INPUT:kspacing'
    AbacusMetric_metric4_4 = 'INPUT:lcao_ecut'
    AbacusMetric_metric5 = 'kpt'
    AbacusMetric_metric6_0 = 'nbase'
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
    AbacusMetric_metric16 = 'fft_grid'
    AbacusMetric_metric17 = 'efermi'
    AbacusMetric_metric18 = 'energy_per_atom'
    AbacusMetric_metric19 = 'stress'
    AbacusMetric_metric19_1 = 'virial'
    AbacusMetric_metric19_2 = 'pressure'
    AbacusMetric_metric20 = 'force'
    AbacusMetric_metric20_1 = 'energies'
    AbacusMetric_metric20_2 = 'forces'
    AbacusMetric_metric20_3 = 'stresses'
    AbacusMetric_metric20_4 = 'virials'
    AbacusMetric_metric20_5 = 'pressures'
    AbacusMetric_metric21 = 'band_gap'
    AbacusMetric_metric22 = 'total_time'
    AbacusMetric_metric23 = 'stress_time'
    AbacusMetric_metric24 = 'force_time'
    AbacusMetric_metric25 = 'scf_time'
    AbacusMetric_metric26 = 'scf_time_each_step'
    AbacusMetric_metric27 = 'step1_time'
    AbacusMetric_metric28 = 'scf_steps'
    AbacusMetric_metric29 = 'atom_mag'
    AbacusMetric_metric29_1 = 'atom_mag_u'
    AbacusMetric_metric29_2 = 'atom_mags'
    AbacusMetric_metric29_3 = 'atom_elec'
    AbacusMetric_metric29_4 = 'atom_orb_elec'
    AbacusMetric_metric29_5 = 'atom_elec_u'
    AbacusMetric_metric30 = 'drho'
    AbacusMetric_metric30_1 = 'drho_last'
    AbacusMetric_metric30_2 = 'denergy'
    AbacusMetric_metric30_3 = 'denergy_last'
    AbacusMetric_metric30_4 = 'denergy_womix'
    AbacusMetric_metric30_5 = 'denergy_womix_last'
    AbacusMetric_metric31 = 'lattice_constant'
    AbacusMetric_metric31_1 = 'lattice_constants'
    AbacusMetric_metric32 = 'cell'
    AbacusMetric_metric32_1 = 'cells'
    AbacusMetric_metric32_2 = 'cell_init'
    AbacusMetric_metric33 = 'coordinate'
    AbacusMetric_metric33_1 = 'coordinate_init'
    AbacusMetric_metric34 = 'element_list'
    AbacusMetric_metric35 = 'atomlabel_list'
    AbacusMetric_metric36 = 'delta_energy'
    AbacusMetric_metric37 = 'delta_energyPerAtom'
    AbacusMetric_metric38 = 'relax_converge'
    AbacusMetric_metric39 = 'relax_steps'
    AbacusMetric_metric40 = 'bda_mag_moment'
    AbacusMetric_metric41 = 'bda_bond_length'
    AbacusMetric_metric42 = 'point_group'
    AbacusMetric_metric43 = 'point_group_in_space_group'
    AbacusMetric_metric44 = 'largest_gradient'
    AbacusMetric_metric45 = 'ds_lambda_step'
    AbacusMetric_metric46 = 'ds_lambda_rms'
    AbacusMetric_metric47 = 'ds_mag'
    AbacusMetric_metric48 = 'ds_mag_force'
    
    AbacusMetric_metric_mem1 = 'mem_vkb'
    AbacusMetric_metric_mem2 = 'mem_psipw'

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
    metrics: Set[AbacusMetricEnum] = Field(default=("version",'INPUT:ks_solver',"normal_end",'converge','energy','total_time','scf_steps'))
    #super_metrics: List[SuperMetrics] = Field(default=None)

class metricsSaveFileSet(BaseModel):    
    metrics_savefile: String = Field(default=None,
                                             description = "If you need to display or upload the metrics generated by executing the postdft command, please enter the file name here. \
Should be a json file, and the key is the name of the metric, and the value is the value of the metric.")
    
    super_metrics_savefile: String = Field(default=None,
                                             description = "If you need to display or upload the super metrics generated by executing the postdft command, please enter the file name here. \
Should be a json file, and the key is the name of the metric, and the value is the value of the metric.")

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

def parse_metrics_set(metrics_set:MetricsSet):
    out_dict = {}
    if hasattr(metrics_set,"metrics"):
        if metrics_set.metrics == None:
            metrics = []
        else:
            metrics = list(metrics_set.metrics)
        has_super_metrics = False
        if hasattr(metrics_set,"super_metrics") and metrics_set.super_metrics != None and len(metrics_set.super_metrics) > 0:
            has_super_metrics = True
            #complete metrics
            for i in metrics_set.super_metrics:
                if i.param_name not in metrics:
                    metrics.append(i.param_name)

        #convert some special metrics (such as: KEY1:KEY2,..)
        metrics = convert_metrics(metrics)

        if len(metrics) > 0:
            #set metrics
            out_dict["metrics"] = {
                "dft_type": "abacus",
                "metrics_name": metrics,
                "save_file": "metrics_default.json",
                "modules":["bda"]
            }

            if has_super_metrics:
                out_dict["super_metrics"] = [{
                    "save_file": "superMetrics.json",
                    "result_file": ["metrics_default.json"],
                    "metrics":[],
                    "outparams":[]
                }]  
                for i in metrics_set.super_metrics:
                    metric_name = convert_supermetrics_metrics_name(i.param_name)
                    out_dict["super_metrics"][-1]["metrics"].append({
                        "name": f"{i.method}({metric_name})" ,
                        "param_name": metric_name,
                        "method": i.method,
                        "normalization": False})
                    out_dict["super_metrics"][-1]["outparams"].append([metric_name, [metric_name], -1])
                    
    if hasattr(metrics_set,"metrics_savefile"):
        if metrics_set.metrics_savefile != None and metrics_set.metrics_savefile.strip() != "":
            if "metrics" not in out_dict:
                out_dict["metrics"] = {}
            out_dict["metrics"]["value_from_file"] = metrics_set.metrics_savefile.strip()
        if metrics_set.super_metrics_savefile != None and metrics_set.super_metrics_savefile.strip() != "":
            if "super_metrics" not in out_dict:
                out_dict["super_metrics"] = [{}]
            out_dict["super_metrics"][-1]["value_from_file"] = metrics_set.super_metrics_savefile.strip()
    
    return out_dict

