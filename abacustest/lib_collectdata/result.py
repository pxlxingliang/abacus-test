import traceback,inspect

class Result:
    def __init__(self):
        self._PARAM_VALUE={}
    
    @classmethod
    def register(cls,**key):
        def aa(method):
            
            for ikey,descript in key.items():
                #print("register",ikey)
                if ikey in cls._PARAM_DIC:
                    method_org = cls._PARAM_DIC[ikey][0].__name__
                    print("WARNING: '%s' has been defined in %s(), modify to %s()" % 
                          (ikey,method_org,method.__name__))
                cls._PARAM_DIC[ikey] = (method,descript)
                #print(cls)

                #for line in inspect.getsourcelines(method):
                #    if ("self.[\'" in line and "\']" in line) or \
                #        ("self.[\"" in line and "\"]" in line):
                            
            return method
        return aa

    def __setitem__(self,key,value):
        if key not in self._PARAM_DIC:
            print("WARNING: the method to catch '%s' is not registered, but is doing the assignment" % key)
            print(inspect.currentframe().f_back)

        try:
            value_org = self._PARAM_VALUE[key]
            if value_org != value:
                print("WARNING: value of '%s' has been catched, but value is different:" % key,value_org,value)
        except:
            self._PARAM_VALUE[key] = value

    def __getitem__(self,key):
        try:
            return self._PARAM_VALUE[key]
        except:
            pass

        try:
            func = self._PARAM_DIC[key][0]
        except:
            print("ERROR: no method to catch value of '%s', or the method is not registered!" % key)
            return None
        
        func_source_file = inspect.getsource(func)
        fback_source_file = inspect.getsource(inspect.currentframe().f_back)
        if func_source_file == fback_source_file:
            framinfo = inspect.getframeinfo(inspect.currentframe().f_back)
            print("WARNING: try to catch the value of '%s' by %s()" % (key,func.__qualname__))
            print("  This is a recursion call, and should be avoid!")
            print("  File:'%s', line %d, %s()" % framinfo[0:3])
            print("%s" % framinfo[3][0])
            return None

        try:
            func(self)
        except:
            print("ERROR: excecute function %s() failed, skip it!" % func.__qualname__)
            traceback.print_exc()
        
        try:
            return self._PARAM_VALUE[key]
        except:
            print("WARNING: try to catch the value of '%s' by %s(), but failed" % 
                    (key,self._PARAM_DIC[key][0].__qualname__))
            return None

    @classmethod
    def AllMethod(cls):
        #return a dict, whose key is the parameter name, 
        #and value is a tuple of functionname, filename, and the description
        paramdic = {}
        for key,value in cls._PARAM_DIC.items():
            filename = value[0].__code__.co_filename
            value = ("%s()"% (value[0].__qualname__),filename,value[1])
            paramdic[key] = value
        return paramdic
    
    def AllParamValue(self):
        return self._PARAM_VALUE

