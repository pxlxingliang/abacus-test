#!/usr/bin/env python3
from lib.dflow import dflowOP,globV
from dflow.python import upload_packages
import os,sys
import numpy
upload_packages.append(os.path.split(numpy.__file__)[0])
upload_packages.append('lib')

def main():
    globV._init()
    dflowOP.RunJobs()
    print("\nAll jobs are finished!!!")

if __name__ == '__main__':
    main()
