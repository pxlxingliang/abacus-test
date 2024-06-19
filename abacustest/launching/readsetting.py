

from . import (comm_class,
               comm_class_exampleSource,
               comm_class_predft,
               comm_class_rundft,
               comm_class_postdft,
               comm_class_metrics,
               comm_class_prepare,
               comm_func)

def ReadSetting(logs:comm_class.myLog,opts,work_path,download_path):
    """
    {
        "config":{},
        "pre_dft":{
            "image":
            "bohrium":{"scass_type","job_type","platform"},
            "example":[],
            "extra_files":[],
            "command":
            "work_directories_filename":
        }
        "run_dft":[{
            "image":
            "bohrium":{"scass_type","job_type","platform"},
            "example":[],
            "group_size":,
            "extra_files":[],
            "command":
        }],
        "post_dft":{
            "command":,
            "extra_files":,
            "image":,
            "bohrium":,
            "metrics":{
                "path":,
                "metric_name":
                "value_from_file":,},
            "super_metric":{
                "value_from_file"},
            "upload_datahub":,
            "upload_tracking":
        }
    }
    """
    logs.iprint("read config setting ...")
    config = comm_func.read_config(opts)

    # download examples/predft_examples/rundft_examples/postdft_examples/
    # predft_extrafiles/rundft_extrafiles/postdft_extrafiles
    logs.iprint("read source setting ...")
    datas = comm_class_exampleSource.read_source(opts,work_path,download_path,logs.iprint)
    if datas == None:
        logs.iprint("Error: download examples or rundft_extrafiles or postdft_extrafiles failed!")
        return None  
    
    need_prepare, prepare = comm_class_prepare.construct_input(datas,opts,work_path,download_path,logs)
    #prepare is a dict
    
    need_predft, pre_dft = comm_class_predft.construct_input(datas,opts,logs)
    #pre_dft is a dict

    need_rundft, run_dft = comm_class_rundft.construct_input(datas,opts,logs)
    #run_dft is a list of dict

    need_postdft, post_dft = comm_class_postdft.construct_input(datas,opts,logs)
    #post_dft is a dict
    
    #read metrics setting
    logs.iprint("read metrics setting ...")
    metrics_set = comm_class_metrics.parse_metrics_set(opts)
    if "metrics" in metrics_set:
        post_dft["metrics"] = metrics_set["metrics"]
        if need_prepare:
            # read INPUT parameters setted in prepare, and set it to metrics
            if "mix_input" in prepare and len(prepare["mix_input"]) > 0:
                if "metrics_name" not in post_dft["metrics"]:
                    post_dft["metrics"]["metrics_name"] = []
                post_dft["metrics"]["metrics_name"].append({"INPUT":[]})
                for ikey in prepare["mix_input"]:
                    if "|" in ikey:
                        post_dft["metrics"]["metrics_name"][-1]["INPUT"] += ikey.split("|")
                    else:
                        post_dft["metrics"]["metrics_name"][-1]["INPUT"].append(ikey)
            # do not set path for metrics 
        #elif run_dft[-1].get("example"):
        #    post_dft["metrics"]["path"] = run_dft[-1]["example"]
        #elif post_dft.get("example"):
        #    post_dft["metrics"]["path"] = post_dft.get("example")
        #else:
        #    post_dft["metrics"]["path"] = ["*"]
        
        
 
        need_postdft = True
    if "super_metrics" in metrics_set:
        post_dft["super_metrics"] = metrics_set["super_metrics"]
        need_postdft = True

    #read tracking setting
    if need_postdft and hasattr(opts,"Tracking_metrics"):
        tracking_set = comm_class.TrackingSet.parse_obj(opts)
        if tracking_set:
            config["aim_access_token"] = tracking_set.get("token")
            post_dft["upload_tracking"] = {
                "tags": tracking_set.get("tags"),
                "name": tracking_set.get("name"),
                "experiment": tracking_set.get("experiment")
            }

    allparams = {"config": config,
            "save_path": "results"}
    if need_prepare:
        allparams["prepare"] = prepare
    if need_predft:
        allparams["pre_dft"] = pre_dft
    if need_rundft:
        allparams["run_dft"] = run_dft
    if need_postdft:
        allparams["post_dft"] = post_dft
    logs.iprint("read setting over!\n")
    return allparams
