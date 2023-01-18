#!/usr/bin/env python
import sys,os
#import importlib

def printAllMethod(allmethod):
    for k,v in allmethod.items():
        print("%20s:\t%30s\t%s"% (k,v[0],v[1]))

def RESULT(fmt="abacus",outparam=False,**kwargs):
    if fmt == "abacus":
        from .resultAbacus import ResultAbacus
        from .abacus import abacus

        if outparam:
            printAllMethod(ResultAbacus.AllMethod())
        else:
            return ResultAbacus(**kwargs)

    elif fmt == "qe":
        from .resultQe import ResultQe
        from .qe import qe

        if outparam:
            printAllMethod(ResultQe.AllMethod())
        else:
            return ResultQe(**kwargs)

    elif fmt == 'vasp':
        from .resultVasp import ResultVasp
        from .vasp import vasp

        if outparam:
            printAllMethod(ResultVasp.AllMethod())
        else:
            return ResultVasp(**kwargs)
    else:
        print("ERROR: unkonw format %s" % fmt)
        sys.exit(1)
