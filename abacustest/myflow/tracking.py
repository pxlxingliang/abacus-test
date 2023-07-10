
import os,traceback
            
def upload_to_tracking(tracking_setting,tracking_values,tracking_summary,AIM_ACCESS_TOKEN=None):
    """
    tracking_values = [(tracking_value1, context1),
                        (tracking_value2,context2),...]
    tracking_value = {name:value}
    tracking_setting = {"name",
                        "experiment",
                        "tags":[]}
    """
    all_tracking_values = []
    if tracking_values:
        all_tracking_values += tracking_values
    if tracking_summary:
        all_tracking_values += tracking_summary
    if all_tracking_values:
        name = tracking_setting.get("name","") + ".summary"
        experiment = tracking_setting.get("experiment")
        tags = tracking_setting.get("tags")
        
        if AIM_ACCESS_TOKEN:
            os.environ["AIM_ACCESS_TOKEN"] = AIM_ACCESS_TOKEN
            
        if "AIM_ACCESS_TOKEN" not in os.environ:
            print("Upload tracking error. Please set 'AIM_ACCESS_TOKEN' information.")
            return None

        from dp.tracking import Run, Text, Table, Image, HTML
        tracking_run = Run(repo='aim://tracking-api.dp.tech:443')
        run_hash = tracking_run.hash
        tracking_run.name = name
        tracking_run.experiment = experiment
        for tag in tags: tracking_run.add_tag(tag)

        def my_track(value,name,context):
            try:
                if isinstance(value,(int,float,Table)):
                    tracking_run.track(value,name=name,context=context)
                elif isinstance(value,str):
                    tracking_run.track(Text(value),name=name,context=context)  
                elif isinstance(value,dict):
                    # upload a file
                    file_type = value.get("type")
                    file_name = value.get("file")
                    if file_type and file_name and os.path.isfile(file_name):
                        if file_type == "image":
                            tracking_run.track(Image(file_name,format="jpeg",optimize=True,quality=50),
                                               name=name,
                                               context=context)
                        elif file_type == "html":
                            with open(file_name) as f: lines = f.read()
                            tracking_run.track(HTML(lines),name=name,context=context)
                else:
                    print(type(value))  
                    tracking_run.track(value,name=name,context=context)
            except:
                traceback.print_exc()
                print("upload tracking failed, name:",name,"value:",value)

        for tracking_value,context in all_tracking_values:
            for k,v in tracking_value.items():
                k = k.replace("/",".")
                print("upload to tracking: %s" % k)
                if isinstance(v,list):
                    for iv in v: 
                        my_track(iv,k,context)
                else:
                    my_track(v,k,context)

        tracking_run.close()
    return True

