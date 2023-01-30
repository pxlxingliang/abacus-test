from . import globV

def printinfo(istr):
    LOGFILE = "abacustest.log"
    with open(LOGFILE,'a+') as f1:
        f1.write(istr + "\n")
    if globV.get_value("OUTINFO"):
        print(istr,flush=True)
