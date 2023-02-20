from . import globV
import os,shutil

def printinfo(istr):
    LOGFILE = "abacustest.log"
    with open(LOGFILE,'a+') as f1:
        f1.write(str(istr) + "\n")
    if globV.get_value("OUTINFO"):
        print(istr,flush=True)
        
def GetBakFile(sfile):
    n = 1
    bk = sfile + ".bak%d" % n
    while os.path.exists(bk):
        n += 1
        bk = sfile + ".bak%d" % n
    return bk

def CopyFiles(path1,path2,move = False):
    '''copy the files in path1 to path2'''
    abspath1 = os.path.abspath(path1)
    abspath2 = os.path.abspath(path2)
    
    def CopyFile(old_path,new_path):
        for ifile in os.listdir(old_path):
            iold_path = os.path.join(old_path,ifile)
            if os.path.isfile(iold_path):
                shutil.copy(iold_path,new_path)
            else:
                if not os.path.isdir(os.path.join(new_path,ifile)):
                    os.makedirs(os.path.join(new_path,ifile))
                CopyFile(iold_path,os.path.join(new_path,ifile))

    if abspath2.startswith(abspath1):
        '''
        If path2 is a son path of path1,
        we will firstly create a tmp path, and copy/move
        files in path1 to the tmp path, and then replace tmp path
        to path2.
        '''
        tmp_path = GetBakFile(os.path.split(abspath1)[0])
        os.makedirs(tmp_path)
        if move:
            for i in os.listdir(abspath1):
                shutil.move(os.path.join(abspath1,i),tmp_path)
        else:
            #shutil.copytree(abspath1,tmp_path,dirs_exist_ok=True)
            CopyFile(abspath1,tmp_path)
        
        if not os.path.isdir(abspath2):
                os.makedirs(abspath2)                
        #shutil.copytree(tmp_path,abspath2,dirs_exist_ok=True)
        CopyFile(tmp_path,abspath2)
        shutil.rmtree(tmp_path)
        
    else:
        if not os.path.isdir(abspath2):
                os.makedirs(abspath2)
        if move:
            for i in os.listdir(abspath1):
                shutil.move(os.path.join(abspath1,i),abspath2)
        else:
            #shutil.copytree(abspath1,abspath2,dirs_exist_ok=True)
            CopyFile(abspath1,abspath2)
            
