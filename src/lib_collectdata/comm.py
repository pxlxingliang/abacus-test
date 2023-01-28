import os,sys,traceback
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element

HARTREE2EV = 27.211396132
EV2RY = 2.0 / HARTREE2EV
RY2EV = HARTREE2EV / 2.0

def ReadFile(ifile,warn=True):
    if ifile == None:
        return []

    if os.path.isfile(ifile):
        with open(ifile) as f1:
            context = f1.readlines()
    else:
        if warn: print("WARNING: can not find file %s" % ifile)
        context = []

    return context

def FindOutput(path,keyinfo):
    'find the file that has keyinfo'
    cwd = os.getcwd()
    os.chdir(path)
    allfiles = os.listdir(".")
    for ifile in allfiles:
        if os.path.isfile(ifile):
            with open(ifile) as f1: lines = f1.read()
            if keyinfo in lines:
                os.chdir(cwd)
                return os.path.join(path,ifile)
    os.chdir(cwd)
    return None

def ReadXmlFile(ifile,warn=True):
    if ifile == None:
        return None

    if os.path.isfile(ifile):
        tree = ET.parse(ifile)
        return tree.getroot()
    else:
        if warn: print("WARNING: can not find file %s" % ifile)
        return None

def XmlFindMultiLayer(root,layerlist):
    tmp = root
    for i in layerlist:
        if tmp == None: return None
        tmp = tmp.find(i)
    return tmp

def XmlFindMultiLayerText(root,layerlist):
    tmp = XmlFindMultiLayer(root,layerlist)
    if isinstance(tmp,Element):
        return tmp.text
    return None

def XmlGetText(tmp1,func=lambda x:x,idx: int=None):
    funcname = sys._getframe().f_code.co_name
    if tmp1 == None:
        return None
    elif isinstance(tmp1,list):
        if idx == None:
            return [func(i.text) for i in tmp1]
        else:
            try:
                return func(tmp1[idx].text)
            except:
                print("ERROR: %s, try to get the text of %d-th element, the length of list is %d" 
                        % (funcname,idx,len(tmp1)))
                traceback.print_exc()
                return None
    elif isinstance(tmp1,Element):
        return func(tmp1.text)
    else:
        print("Not support type '%s' in function '%s' now." % 
              (type(tmp1),funcname))
        return None

def istr(x,n=None):
    if type(x) == float and n != None:
            return "%%.%df" % n % x
    return str(x)

def iint(x):
    try:
        return int(x)
    except:
        return None

def ifloat(x):
    try:
        return float(x)
    except:
        return None

def ibool(x):
    try:
        return bool(x)
    except:
        return None

def imath(a,b,symbol):
    if a == None or b == None:
        return None
    elif symbol == '*':
        return a*b
    elif symbol == '+':
        return a+b
    elif symbol == '-':
        return a-b
    elif symbol == '/':
        return a/b
    else:
        print("Not support '%s' now" % symbol)
        return None

