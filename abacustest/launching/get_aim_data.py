import traceback
import requests,json
import struct
from pprint import pprint
from typing import Iterator
from tqdm  import tqdm
#from loguru import logger
from aim.storage.treeutils import decode_tree

def decode_encoded_tree_stream(stream: Iterator[bytes], concat_chunks=False) -> bytes:
    prev_chunk_tail = b''
    if concat_chunks:
        data = b''
        for chunk in stream:
            data += chunk
        while data:
            (key_size,), data_tail = struct.unpack('I', data[:4]), data[4:]
            key, data_tail = data_tail[:key_size], data_tail[key_size:]
            (value_size,), data_tail = struct.unpack('I', data_tail[:4]), data_tail[4:]
            value, data_tail = data_tail[:value_size], data_tail[value_size:]
            data = data_tail
            yield key, value
    else:
        for chunk in stream:
            data = prev_chunk_tail + chunk
            prev_chunk_tail = b''
            while data:
                try:
                    (key_size,), data_tail = struct.unpack('I', data[:4]), data[4:]
                    key, data_tail = data_tail[:key_size], data_tail[key_size:]
                    (value_size,), data_tail = struct.unpack('I', data_tail[:4]), data_tail[4:]
                    value, data_tail = data_tail[:value_size], data_tail[value_size:]
                    data = data_tail
                except Exception:
                    prev_chunk_tail = data
                    break
                yield key, value
        assert prev_chunk_tail == b''
 
def get_version(run_hash,token):
    params = {
        "record_density": 50,
        "index_density": 5,
    }
    response = requests.post(f'https://tracking.mlops.dp.tech/api/runs/{run_hash}/texts/get-batch',
                        params=params,
                        json=[{'name': "version", 'context': {'datatype': "metrics",'subset': "metrics0"}}],
                        headers={
                            "Content-Type": "application/json",
                            "Authorization": f"Bearer {token}",
                        })
    if response.status_code != 200:
        return None
    try:
        response.raise_for_status()
        decoded_response = decode_tree(decode_encoded_tree_stream(response.iter_content(chunk_size=512*1024), concat_chunks=True))
        version = [i[0]["data"] for i in decoded_response["values"]]
        index = decoded_response["iters"]
        return [x[0] for x in sorted(zip(version,index),key=lambda x:x[1])]
    except:
        #print("hash:",run_hash)
        #traceback.print_exc()
        return None
        
def get_experiment_info(experiment_id, token):
    response = requests.get(f'https://tracking.mlops.dp.tech/api/experiments/{experiment_id}/',
                        headers={
                            "Content-Type": "application/json",
                            "Authorization": f"Bearer {token}",
                        })
    response.raise_for_status()
    return response.json()
def get_run_info_batch(experiment, token, offset=None, batch_size=10):
    query = f"run.experiment == '{experiment}' and not run.archived"
    params = {'q': query, 'report_progress': False, 'limit': batch_size, 'exclude_params': True}
    if offset:
        params['offset'] = offset
    response = requests.get('https://tracking.mlops.dp.tech/api/runs/search/run',
                        params=params,
                        headers={
                            "Content-Type": "application/json",
                            "Authorization": f"Bearer {token}",
                        })
    response.raise_for_status()
    decoded_response = decode_tree(decode_encoded_tree_stream(response.iter_content(chunk_size=512*1024), concat_chunks=True))
    return list(sorted(decoded_response.items(), key=lambda i: i[1]['props']['creation_time'], reverse=True))

def collect_from_runinfo(run,token,hash,GetVersion=False):
        metric = {}
        for imetric in run["traces"]["metric"]:
            if imetric["name"].startswith("__"):
                continue
            #only collect super_metrics
            #if imetric["context"].get("subset",None) == "super_metrics0" or \
            #   imetric["context"].get("datatype",None) == "super_metrics":
            metric[imetric["name"]] = [imetric["last_value"]["last"]]
        if GetVersion:
            version = get_version(hash,token)
            if version:
                metric["version"] = version
        return {
            "run_hash": hash,
            "run_name": run["props"]["name"],
            "experiment_name": run["props"]["experiment"]["name"],
            "tags": [i["name"] for i in run["props"]["tags"]],
            "creation_time": run["props"]["creation_time"],
            "end_time": run["props"]["end_time"],
            "metric": metric
            }

def do_collect_info(all_run,all_run_infos, run_infos,token,needed_tags=None,GetVersion=False):
    for hash, run in run_infos:
        if needed_tags:
            tags_in_run = set([i["name"] for i in  run["props"]["tags"]] )
            if not needed_tags.intersection(tags_in_run):
                continue
        all_run_infos.append(collect_from_runinfo(run,token,hash,GetVersion)) #collect selected infos
        all_run.append(run)  #store all infos

def do_collect_metrics(all_run,all_run_infos, run_info, token,needed_tags=None,GetVersion=False):
    if needed_tags:
        tags_in_run = set([i["name"] for i in  run_info[1]["props"]["tags"]] )
        if not needed_tags.intersection(tags_in_run):
            return
    resp = requests.post(f'https://tracking.mlops.dp.tech/api/runs/{run_info[0]}/metric/get-batch',
        data=json.dumps(run_info[1]['traces']['metric']), 
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {token}",
            })
    run_info_tmp = run_info[1]
    run_info_tmp["metric"] = resp.json()
    all_run.append(run_info_tmp)  #store extra metric values

    selected_run_info = collect_from_runinfo(run_info[1],token,run_info[0],GetVersion)
    for imetric in run_info_tmp["metric"]:
        if imetric["name"].startswith("__"):
            continue
        selected_run_info["metric"][imetric["name"]] = imetric["values"]
    all_run_infos.append(selected_run_info)

def get_runs(token,experiment,experiment_id,needed_tags=None,collect_metrics=False,GetVersion=False):
    import time
    start = time.time()
    batch_size = 45
    offset = None
    all_run_infos = []
    all_run = []
    experiment_info = get_experiment_info(experiment_id, token)
    #logger.info(experiment_info)
    total = experiment_info['run_count']
    needed_tags = set(needed_tags) if needed_tags else None
    with tqdm(total=total) as pbar:
        pbar.set_description(f"fetch next batch offset={offset}, batch_size={batch_size}")
        run_infos = get_run_info_batch(experiment, token, batch_size=batch_size)
        while run_infos:
            if collect_metrics:
                for run_info in run_infos:
                    do_collect_metrics(all_run,all_run_infos, run_info, token,needed_tags,GetVersion)
                    pbar.update(1)
            else:
                do_collect_info(all_run,all_run_infos, run_infos, token,needed_tags,GetVersion)
                pbar.update(batch_size)
            try:
                if len(run_infos) < batch_size:
                    break
                offset = run_infos[-1][0]
                pbar.set_description(f"fetch next batch offset={offset}, batch_size={batch_size}")
                run_infos = get_run_info_batch(experiment, token, offset=offset, batch_size=batch_size)
            except KeyError: # already empty
                break
    elapsed = time.time() - start
    #logger.info(f'All done! elapsed: {elapsed}')
    return all_run,all_run_infos

'''
all_run is a list of all runs (dict), is the original format
one of the run is like:
{
    "props": {
        "name": "lcao-scalapack.lcao-scalapack-86073e.summary",
        "description": null,
        "experiment": {
            "id": "7ab4e46a-43fb-440a-828d-4fbdef5b4709",
            "name": "abacustest/benchmark",
            "description": null
        },
        "tags": [
            {
                "id": "0c23152e-ad9a-4612-aef7-4f94f4a914e4",
                "name": "benchmark-project-abacustest",
                "color": null,
                "description": null
            },
            ...
        ],
        "archived": false,
        "creation_time": 1682010699.0,
        "end_time": 1682010699.462929,
        "active": false
    },
    "traces": {
        "metric": [
            {
                "context": {},
                "name": "__system__cpu",
                "last_value": {
                    "dtype": "float",
                    "first_step": 0,
                    "last_step": 0,
                    "last": 0.0,
                    "version": 2
                }
            },
            ...
        ]
    }
}

The run_infos is a list of the dict of critical informations of each run
one of the dict is like:
{
            "run_name": run["props"]["name"],
            "experiment_name": run["props"]["experiment"]["name"],
            "tags": [i["name"] for i in run["props"]["tags"]],
            "creation_time": run["props"]["creation_time"],
            "end_time": run["props"]["end_time"],
            "metric": metric
}
metric: {metric_nam1:[],
        metric_name2:[],
        ...}

Be noticed, if the collect_metrics is False, the the metric of run_infos is only
the last value. And only set collect_metrics to be True, we can get the all elements of 
each metric.
'''