

from . import (comm_class,
               comm_class_exampleSource,
               comm_class_rundft,
               comm_class_postdft,
               comm_class_metrics,
               comm_func)

def ReadSetting(logs:comm_class.myLog,opts,work_path,download_path):
    """
    {
        "config":{},
        "run_dft":[{
            "image":
            "bohrium":{"scass_type","job_type","platform"},
            "example":[],
            "extra_files":[],
            "command":
        }],
        "post_dft":{
            "command":,
            "extra_files":,
            "iamge":,
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

    #download examples/rundft_extrafiles/postdft_extrafiles
    logs.iprint("read source setting ...")
    datas = comm_class_exampleSource.read_source(opts,work_path,download_path,logs.iprint)
    if datas == None or not datas.get("example"):
        logs.iprint("Error: download examples or rundft_extrafiles or postdft_extrafiles failed!")
        return None  
    
    #read rundft
    logs.iprint("read run dft setting ...")
    run_dft = [{}]
    
    #read rundft example
    if datas.get("example"):
        run_dft[-1]["example"] = datas.get("example")    
        logs.iprint("\texample:",run_dft[-1]["example"])
    
    #read rundft extra files
    if datas.get("rundft_extrafiles"):
        run_dft[-1]["extra_files"] = datas.get("rundft_extrafiles")
        logs.iprint("\trundft_extrafiles:",run_dft[-1]["extra_files"])
        
    #read rundft command
    if hasattr(opts,"rundft_command"):
        logs.iprint("\tcommand:",opts.rundft_command)
        run_dft[-1]["command"] = opts.rundft_command
    
    #read rundft ngroup        
    if hasattr(opts,"ngroup") and opts.ngroup > 0:
        logs.iprint("\tngroup:",opts.ngroup)
        run_dft[-1]["ngroup"] = opts.ngroup
    
    #read rundft image
    if hasattr(opts,"rundft_image_set"):
        logs.iprint("\timage:",opts.rundft_image_set.image)
        for k,v in comm_class_rundft.parse_image_set(opts.rundft_image_set).items():
            run_dft[-1][k] = v
        if "bohrium" in run_dft[-1]:   
            logs.iprint("\tbohrium:",run_dft[-1]["bohrium"])

    #read postdft
    logs.iprint("read post dft setting ...")
    need_post_dft = False
    post_dft = {}
    
    #read postdft image
    if hasattr(opts,"postdft_image_set"):
        logs.iprint("\timage:",opts.postdft_image_set.image)
        post_dft = comm_class_postdft.parse_image_set(opts.postdft_image_set)
        if "bohrium" in post_dft:
            logs.iprint("\tbohrium:",post_dft["bohrium"])
    
    #read postdft command
    if hasattr(opts,"postdft_command"): 
        if opts.postdft_command != None and opts.postdft_command.strip() != "":
            post_dft["command"] = opts.postdft_command.strip()
            need_post_dft = True
            logs.iprint("\tcommand:",post_dft["command"])
        
    #read postdft extra files
    if datas.get("postdft_extrafiles"):
        post_dft["extra_files"] = datas.get("postdft_extrafiles")
        logs.iprint("\tpostdft_extrafiles:",post_dft["extra_files"])
        
    #read metrics setting
    logs.iprint("read metrics setting ...")
    metrics_set = comm_class_metrics.parse_metrics_set(opts)
    if "metrics" in metrics_set:
        post_dft["metrics"] = metrics_set["metrics"]
        post_dft["metrics"]["path"] = run_dft[-1]["example"]
        need_post_dft = True
    if "super_metrics" in metrics_set:
        post_dft["super_metrics"] = metrics_set["super_metrics"]
        need_post_dft = True

    #read tracking setting
    if need_post_dft and hasattr(opts,"Tracking_metrics"):
        tracking_set = comm_class.TrackingSet.parse_obj(opts)
        if tracking_set:
            config["AIM_ACCESS_TOKEN"] = tracking_set.get("token")
            post_dft["upload_tracking"] = {
                "tags": tracking_set.get("tags"),
                "name": tracking_set.get("name"),
                "experiment": tracking_set.get("experiment")
            }

    allparams = {"config": config,
            "run_dft": run_dft,
            "save_path": "results"}
    if need_post_dft:
        allparams["post_dft"] = post_dft
    logs.iprint("read setting over!\n")
    return allparams
