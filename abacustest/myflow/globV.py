#!/usr/bin/env python
from datetime import datetime
def _init():
    global _global_dict
    _global_dict = {"START_TIME":datetime.now()}

def set_value(key,value):
    # check if _global_dict has been initialized
    if "_global_dict" not in globals():
        _init()
    
    if key in _global_dict:
        print("WARNING: key %s has been used, the value will be modified from '" % key,_global_dict[key],"' to '",value,"'")
    _global_dict[key] = value

def get_value(key,defValue=None):
    if "_global_dict" not in globals():
        _init()
    return _global_dict.get(key,defValue)




