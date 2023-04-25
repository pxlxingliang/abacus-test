import sys,os,glob
import traceback
import importlib

def printAllMethod(allmethod,fmt):
    ipath = os.path.split(os.path.abspath(__file__))[0]
    print("Job type: %s, all modules: " % fmt,[os.path.split(i)[1][:-3] for i in glob.glob(os.path.join(ipath,"%s/*.py"%fmt))])
    print("\nmetric_name,\tmodule:function,\tdescription")
    for k,v in allmethod.items():
        filename = v[1] 
        funcname = v[0]
        out = filename.split("/")[-1][:-3] + ":" + funcname 
        print("%20s:\t%-30s\t%s"% (k,out,v[2]))
    print("\nNOTICE1: keys start with \"delta_\" require a json file includes the reference values. You can specify the file by \"--ref\"")
    print("NOTICE2: by default, only keys in module '%s' are available. You can use other modules by '--modules'" % fmt)

def import_new_method(newmethods=[]):
    '''import the self-defined methods'''
    if newmethods == None:
        return
    sys.path.append(os.getcwd())
    for imethod in newmethods:
        try:
            importlib.import_module(imethod)
            print("Import module: %s" % imethod)
        except:
            traceback.print_exc()
            print("Import module %s failed, skip it!" % imethod)
    #importlib.invalidate_caches()

def import_modules(fmt,modules):
    '''
    Can only import the modules in lib_collectdata.{fmt}
    '''
    if modules == None:
        return
    for module in modules:
        try:
            importlib.import_module("abacustest.lib_collectdata.%s.%s"%(fmt,module))
            print("import module abacustest.lib_collectdata.%s.%s" % (fmt,module))
        except:
            traceback.print_exc()
            print("Import module %s failed, skip it!" % module)

def RESULT(fmt="abacus",outparam=False,newmethods=[],modules=[],**kwargs):

    if fmt == "abacus":
        from .resultAbacus import ResultAbacus
        from .abacus import abacus
        import_modules(fmt,modules)
        import_new_method(newmethods)
            
        if outparam:
            printAllMethod(ResultAbacus.AllMethod(),fmt)
        else:
            return ResultAbacus(**kwargs)

    elif fmt == "qe":
        from .resultQe import ResultQe
        from .qe import qe
        import_modules(fmt,modules)
        import_new_method(newmethods)

        if outparam:
            printAllMethod(ResultQe.AllMethod(),fmt)
        else:
            return ResultQe(**kwargs)

    elif fmt == 'vasp':
        from .resultVasp import ResultVasp
        from .vasp import vasp
        import_modules(fmt,modules)
        import_new_method(newmethods)

        if outparam:
            printAllMethod(ResultVasp.AllMethod(),fmt)
        else:
            return ResultVasp(**kwargs)
    else:
        print("ERROR: unkonw format %s" % fmt)
        sys.exit(1)
