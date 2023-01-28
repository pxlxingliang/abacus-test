from . import globV

def printinfo(istr):
    if globV.get_value("OUTINFO"):
        print(istr)
