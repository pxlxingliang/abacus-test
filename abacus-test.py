#!/usr/bin/env python3
import os,sys
sys.path.append('lib')
from lib.dflow import dflow,globV
from dflow.python import upload_packages
import numpy
upload_packages.append(os.path.split(numpy.__file__)[0])
upload_packages.append(os.path.join(os.path.split(__file__)[0],'lib'))

def main():
    globV._init()
    dflow.RunJobs()
    print("\nAll jobs are finished!!!")

if __name__ == '__main__':
    main()
