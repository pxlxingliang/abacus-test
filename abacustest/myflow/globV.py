#!/usr/bin/env python

def _init():
    global _global_dict
    _global_dict = {}

def set_value(key,value):
    if key in _global_dict:
        print("WARNING: key %s has been used, the value will be modified from '" % key,_global_dict[key],"' to '",value,"'")
    _global_dict[key] = value

def get_value(key,defValue=None):
    return _global_dict.get(key,defValue)




