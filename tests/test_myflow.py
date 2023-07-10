import unittest,os,shutil,re
from pathlib import Path

from abacustest.myflow import flow,dflowOP,globV
from dflow import (
    Workflow)
import defaultset
from abacustest.myflow import PredftOP,RundftOP,PostdftOP

class TestPredft(unittest.TestCase):
    def setUp(self):
        globV._init()
        globV.set_value("ABBREVIATION",{})
        globV.set_value("PARAM_FNAME","param.json")
        globV.set_value("HOST","test.test")
        self.cwd = os.getcwd()
        self.work_path = Path("testflow")
        os.makedirs("testflow",exist_ok=True)
        os.chdir(self.work_path)
        self.examples = ["a","b","c"]
        self.exatrafiles = ["a.txt","b.txt","c.txt"]
        self.exatrafiles1 = ["a1.txt","b1.txt","c1.txt"]
        self.exatrafiles2 = ["a2.txt","b2.txt","c2.txt"]
        for i in self.examples: 
            os.makedirs(i,exist_ok=True) 
            Path(os.path.join(i,f"{i}.inp")).touch()
        for i in self.exatrafiles: Path(i).touch()
        for i in self.exatrafiles1: Path(i).touch()
        for i in self.exatrafiles2: Path(i).touch()

    def tearDown(self):
        os.chdir(self.cwd)
        if os.path.isdir("testflow"):
            shutil.rmtree("testflow")

    
    def test_predft(self):
        param = {
            "save_path": "result",
            "pre_dft": {
                "example": self.examples,
                "command": "echo 1 > log",
                "image": "python:3.8"
            }   
        }
        globV.set_value("PARAM_CONTEXT",param)
        allstep,stepname,allsave_path = dflowOP.ProduceAllStep(param)
        wf = Workflow(name="abacustest")
        wf.add(allstep)
        wf.submit()
        flow.waitrun(wf,stepname,allsave_path)
        
        expect_dir = ["result", "result/a", "result/b", "result/c"]
        expect_file = ["result/a/log", "result/b/log", "result/c/log",
                       "result/a/a.inp", "result/b/b.inp", "result/c/c.inp"
                       ]
        
        for i in expect_dir:
            self.assertTrue(os.path.isdir(i),f"{i} is not exist!")
        for i in expect_file:
            self.assertTrue(os.path.isfile(i),f"{i} is not exist!")     
    
    def test_rundft(self):
        param = {
            "save_path": "result",
            "run_dft": {
                "example": [["a","b"],"c"],
                "command": "echo run > run.log",
                "image": "python:3.8",
                "extra_files": self.exatrafiles1
            },  
        }
        globV.set_value("PARAM_CONTEXT",param)
        allstep,stepname,allsave_path = dflowOP.ProduceAllStep(param)
        wf = Workflow(name="abacustest")
        wf.add(allstep)
        wf.submit()
        flow.waitrun(wf,stepname,allsave_path)
        
        expect_dir = ["result", "result/a", "result/b", "result/c"]
        expect_file = [
                       "result/a/a.inp", "result/b/b.inp", "result/c/c.inp",
                       "result/a/run.log", "result/b/run.log", "result/c/run.log",
                       "result/a/a1.txt", "result/a/b1.txt", "result/a/c1.txt",
                       "result/b/a1.txt", "result/b/b1.txt", "result/b/c1.txt",
                       "result/c/a1.txt", "result/c/b1.txt", "result/c/c1.txt"
                       ]
        
        for i in expect_dir:
            self.assertTrue(os.path.isdir(i),f"{i} is not exist!")
        for i in expect_file:
            self.assertTrue(os.path.isfile(i),f"{i} is not exist!")

    def test_postdft(self):
        param = {
            "save_path": "result",
            "post_dft": {
                "example": ["a","b","c"],
                "command": "echo post > post.log",
                "image": "python:3.8",
                "extra_files": self.exatrafiles2
            },  
        }
        globV.set_value("PARAM_CONTEXT",param)
        allstep,stepname,allsave_path = dflowOP.ProduceAllStep(param)
        wf = Workflow(name="abacustest")
        wf.add(allstep)
        wf.submit()
        flow.waitrun(wf,stepname,allsave_path)
        
        expect_dir = ["result", "result/a", "result/b", "result/c"]
        expect_file = [
                       "result/a/a.inp", "result/b/b.inp", "result/c/c.inp",
                       "result/post.log", 
                       "result/a2.txt", "result/b2.txt", "result/c2.txt"
                       ]
        
        for i in expect_dir:
            self.assertTrue(os.path.isdir(i),f"{i} is not exist!")
        for i in expect_file:
            self.assertTrue(os.path.isfile(i),f"{i} is not exist!")
        
    def test_predft_rundft_postdft(self):
        param = {
            "save_path": "result",
            "pre_dft": {
                "example": self.examples,
                "command": "mkdir -p d1 d2; echo d1 > workdirs;echo d2 >> workdirs; echo predft > d1/log; cp *.txt d1;cp *.inp d2",
                "image": "python:3.8",
                "extra_files": self.exatrafiles,
                "work_directories_filename": "workdirs"
            } ,
            "run_dft": {
                "command": "echo run > run.log",
                "image": "python:3.8",
                "extra_files": self.exatrafiles1,
                "group_size":2
            },  
            "post_dft": {
                "command": "echo post > post.log",
                "image": "python:3.8",
                "extra_files": self.exatrafiles2,
            },
        }
        globV.set_value("PARAM_CONTEXT",param)
        allstep,stepname,allsave_path = dflowOP.ProduceAllStep(param)
        wf = Workflow(name="abacustest")
        wf.add(allstep)
        wf.submit()
        flow.waitrun(wf,stepname,allsave_path)
        
        expect_dir = ["result", "result/a", "result/b", "result/c",
                      "result/a/d1", "result/b/d1", "result/c/d1",
                      "result/a/d2", "result/b/d2", "result/c/d2",]
        expect_file = [
            "result/a/d1/log", "result/b/d1/log", "result/c/d1/log",
            "result/a/d1/a.txt", "result/a/d1/b.txt", "result/a/d1/c.txt",
            "result/b/d1/a.txt", "result/b/d1/b.txt", "result/b/d1/c.txt",
            "result/c/d1/a.txt", "result/c/d1/b.txt", "result/c/d1/c.txt",
            "result/a/d2/a.inp", "result/b/d2/b.inp", "result/c/d2/c.inp",
            "result/a/d1/run.log", "result/b/d1/run.log", "result/c/d1/run.log",
            "result/a/d2/run.log", "result/b/d2/run.log", "result/c/d2/run.log",
            "result/a/d1/a1.txt", "result/a/d1/b1.txt", "result/a/d1/c1.txt",
            "result/b/d1/a1.txt", "result/b/d1/b1.txt", "result/b/d1/c1.txt",
            "result/c/d1/a1.txt", "result/c/d1/b1.txt", "result/c/d1/c1.txt",
            "result/a/d2/a1.txt", "result/a/d2/b1.txt", "result/a/d2/c1.txt",
            "result/b/d2/a1.txt", "result/b/d2/b1.txt", "result/b/d2/c1.txt",
            "result/c/d2/a1.txt", "result/c/d2/b1.txt", "result/c/d2/c1.txt",
            "result/post.log", "result/a2.txt", "result/b2.txt", "result/c2.txt"
                       ]
        
        for i in expect_dir:
            self.assertTrue(os.path.isdir(i),f"{i} is not exist!")
        for i in expect_file:
            self.assertTrue(os.path.isfile(i),f"{i} is not exist!") 
    
    
    def test_prepare_predft(self):
        os.makedirs("prepare",exist_ok=True)
        for iexample in ["a","b","c"]:
            os.makedirs(os.path.join("prepare",iexample),exist_ok=True)
            Path(os.path.join("prepare",iexample,"STRU")).write_text(defaultset.STRU1)
            Path(os.path.join("prepare",iexample,"INPUT")).write_text(defaultset.INPUT1)
            Path(os.path.join("prepare",iexample,"KPT")).write_text(defaultset.KPT1)
            Path(os.path.join("prepare",iexample,"Si.upf")).write_text("Si.upf")
            Path(os.path.join("prepare",iexample,"Si.orb")).write_text("Si.orb")
        # the new prepared a/b/c should be in current work path, 
        # and has 3*2=6 subdirs named 00000,00001,00002,00003,00004,00005
        param = {
            "save_path": "result",
            "prepare":{
                "example_template":["prepare/a","prepare/b","prepare/c"],
                "mix_input":{
                    "kspacing": [0.1,0.2,0.3],
                    "ks_solver":["cg","dav"]
                }
            },
            "pre_dft": {
                "command": "echo 1 > log",
                "image": "python:3.8"
            }   
        }
        globV.set_value("PARAM_CONTEXT",param)
        allstep,stepname,allsave_path = dflowOP.ProduceAllStep(param)
        wf = Workflow(name="abacustest")
        wf.add(allstep)
        wf.submit()
        flow.waitrun(wf,stepname,allsave_path)
        
        examples = ["a","b","c"]
        sub_dir = ["00000", "00001", "00002", "00003", "00004", "00005"]
        expect_file = ["log", "INPUT", "STRU", "KPT", "Si.upf", "Si.orb"]
        input_set_ref = [[0.1,"cg"],[0.1,"dav"],[0.2,"cg"],[0.2,"dav"],[0.3,"cg"],[0.3,"dav"]]
        input_set_ref.sort()
        
        self.assertTrue(os.path.isdir("result"))
        for iexample in examples:
            input_set = []
            for isubdir in sub_dir:
                self.assertTrue(os.path.isdir(os.path.join("result",iexample,isubdir)))
                for ifile in expect_file:
                    self.assertTrue(os.path.isfile(os.path.join("result",iexample,isubdir,ifile)))
                    
                input_text = Path(os.path.join("result",iexample,isubdir,"INPUT")).read_text()
                kspacing = re.search(r"kspacing\s+0.[1,2,3]",input_text)
                solver = re.search(r"ks_solver\s+(cg|dav)",input_text)
                self.assertTrue(kspacing)
                self.assertTrue(solver)
                input_set.append([float(kspacing.group().split()[1]),solver.group().split()[1]])
            input_set.sort()
            self.assertEqual(input_set,input_set_ref)
                
        
        
        
        