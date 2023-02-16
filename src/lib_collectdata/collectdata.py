#!/usr/bin/env python
import sys,os
import traceback
#import importlib

def printAllMethod(allmethod):
    for k,v in allmethod.items():
        print("%20s:\t%30s\t%s"% (k,v[0],v[1]))

def import_new_method(newmethods=None):
    if newmethods != None:
        import importlib
        for imethod in newmethods:
            try:
                importlib.import_module(imethod)
                print("Import module: %s" % imethod)
            except:
                traceback.print_exc()
                print("Import module %s failed, skip it!" % imethod)
        importlib.invalidate_caches()


def RESULT(fmt="abacus",outparam=False,newmethods=None,**kwargs):

    if fmt == "abacus":
        from .resultAbacus import ResultAbacus
        from .abacus import abacus
        import_new_method(newmethods)
            
        if outparam:
            printAllMethod(ResultAbacus.AllMethod())
        else:
            return ResultAbacus(**kwargs)

    elif fmt == "qe":
        from .resultQe import ResultQe
        from .qe import qe
        import_new_method(newmethods)

        if outparam:
            printAllMethod(ResultQe.AllMethod())
        else:
            return ResultQe(**kwargs)

    elif fmt == 'vasp':
        from .resultVasp import ResultVasp
        from .vasp import vasp
        import_new_method(newmethods)

        if outparam:
            printAllMethod(ResultVasp.AllMethod())
        else:
            return ResultVasp(**kwargs)
    else:
        print("ERROR: unkonw format %s" % fmt)
        sys.exit(1)
