import unittest,os,shutil,re
from pathlib import Path

from abacustest import prepare
import defaultset

class TestPredft(unittest.TestCase):
    def setUp(self):
        self.cwd = os.getcwd()
        self.work_path = Path("testprepare")
        os.makedirs(self.work_path,exist_ok=True)
        os.chdir(self.work_path)
        
        # create example template
        examples = ["a","b","c"]
        for i in examples:
            os.makedirs(i,exist_ok=True)
            Path(os.path.join(i,"STRU")).write_text(defaultset.STRU1)
            Path(os.path.join(i,"INPUT")).write_text(defaultset.INPUT1)
            Path(os.path.join(i,"Si.upf")).write_text("Si.upf")
            Path(os.path.join(i,"Si.orb")).write_text("Si.orb")
            Path(os.path.join(i,"KPT")).write_text(defaultset.KPT1)
        
        # create pp lib
        os.makedirs("pplib",exist_ok=True)
        Path("pplib/Si.upf1").write_text("Si.upf1")
        
        # create orb lib
        os.makedirs("orblib",exist_ok=True)
        Path("orblib/Si.orb1").write_text("Si.orb1")
        
        # create INPUT/KPT/STRU template
        Path("INPUT").write_text(defaultset.INPUT2)
        Path("KPT").write_text(defaultset.KPT2)
        Path("STRU").write_text(defaultset.STRU2)

    def tearDown(self):
        os.chdir(self.cwd)
        if os.path.isdir(self.work_path):
            shutil.rmtree(self.work_path)
    
    
    def test_prepare_example_template(self):
        #case1: use example template, and mix input
        param_setting = {
            "example_template": ["a","b","c"],
            "mix_input": {
                "ecutwfc": [50, 60, 70],
                "kspacing": [0.1, 0.12, 0.13]
            }
        }
        all_path_setting = prepare.DoPrepare(param_setting,"abacustest")
        
        # check result folder
        self.assertTrue(os.path.isdir("abacustest"))
        
        # check new create examples, should create 9 examples for each template
        self.assertEqual(len(all_path_setting), 3)
        self.assertEqual(len(all_path_setting[0]), 9)
        self.assertEqual(len(all_path_setting[1]), 9)
        self.assertEqual(len(all_path_setting[2]), 9)
        
        # check files in each example, and ecutwfc/kspacing value in INPUT
        ecut_kspacing_ref = [[50,0.1],[50,0.12],[50,0.13],[60,0.1],[60,0.12],[60,0.13],[70,0.1],[70,0.12],[70,0.13]]
        ecut_kspacing_ref.sort()
        for i in all_path_setting:
            ecut_kspacing = []
            for j in i:
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.upf")))
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.orb")))
                self.assertTrue(os.path.isfile(os.path.join(j,"STRU")))
                self.assertTrue(os.path.isfile(os.path.join(j,"INPUT")))
                self.assertTrue(os.path.isfile(os.path.join(j,"KPT")))
                input_text = Path(os.path.join(j,"INPUT")).read_text()
                ecutwfc = re.search(r"ecutwfc\s+[5,6,7]0",input_text)
                kspacing = re.search(r"kspacing\s+0\.1[0,2,3]?",input_text)
                self.assertTrue(ecutwfc)
                self.assertTrue(kspacing)
                ecut_kspacing.append([float(ecutwfc.group().split()[1]),float(kspacing.group().split()[1])])
            ecut_kspacing.sort()
            self.assertEqual(ecut_kspacing,ecut_kspacing_ref)
    
    def test_prepare_input_template(self):
        #case2: use INPUT template, and mix input/KPT 
        # check mix_kpt by using list or int
        param_setting = {
            "example_template": ["a","b","c"],
            "input_template": "INPUT",
            "mix_input": {
                "ecutwfc": [50, 60, 70]
            },
            "mix_kpt": [[1,2,3],2]
        }
        all_path_setting = prepare.DoPrepare(param_setting,"abacustest")  
        
        # check result folder
        self.assertTrue(os.path.isdir("abacustest"))
        
        # check new create examples, should create 2*3 = 6 examples for each template
        self.assertEqual(len(all_path_setting), 3)
        self.assertEqual(len(all_path_setting[0]), 6)
        self.assertEqual(len(all_path_setting[1]), 6)
        self.assertEqual(len(all_path_setting[2]), 6)
        
        # check files in each example, and ecutwfc/kspacing value in INPUT
        ecut_kpt_ref = [[50,[1,2,3]],[50,[2,2,2]],
                        [60,[1,2,3]],[60,[2,2,2]],
                        [70,[1,2,3]],[70,[2,2,2]]]
        
        for i in all_path_setting:
            ecut_kpt = []
            for j in i:
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.upf")))
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.orb")))
                self.assertTrue(os.path.isfile(os.path.join(j,"STRU")))
                self.assertTrue(os.path.isfile(os.path.join(j,"INPUT")))
                self.assertTrue(os.path.isfile(os.path.join(j,"KPT")))
                
                input_text = Path(os.path.join(j,"INPUT")).read_text()
                self.assertTrue(re.search(r"basis_type\s+pw",input_text))
                ecutwfc = re.search(r"ecutwfc\s+[5,6,7]0",input_text)
                self.assertTrue(ecutwfc)
                
                kpt_text = Path(os.path.join(j,"KPT")).read_text()
                kpt = [ int(i) for i in kpt_text.strip().split("\n")[-1].split()[:3]]
                ecut_kpt.append([float(ecutwfc.group().split()[1]),kpt])
            ecut_kpt.sort()
            self.assertEqual(ecut_kpt,ecut_kpt_ref) 

    def test_prepare_kpt_template(self):
        # case3: use kpt template, and mix STRU 
        param_setting = {
            "example_template": ["a","b","c"],
            "kpt_template": "KPT",
            "mix_stru": ["STRU","a/STRU"]
        }
        all_path_setting = prepare.DoPrepare(param_setting,"abacustest")  
        
        # check result folder
        self.assertTrue(os.path.isdir("abacustest"))
        
        # check new create examples, should create 2*3 = 6 examples for each template
        print(all_path_setting)
        self.assertEqual(len(all_path_setting), 3)
        self.assertEqual(len(all_path_setting[0]), 2)
        self.assertEqual(len(all_path_setting[1]), 2)
        self.assertEqual(len(all_path_setting[2]), 2)

        lc_ref = [5.1,10.2]
        for i in all_path_setting:
            lcs = []
            for j in i:
                self.assertTrue((os.path.isfile(os.path.join(j,"Si.upf")) or os.path.isfile(os.path.join(j,"Si.upf1"))))
                self.assertTrue((os.path.isfile(os.path.join(j,"Si.orb")) or os.path.isfile(os.path.join(j,"Si.orb1"))))
                self.assertTrue(os.path.isfile(os.path.join(j,"STRU")))
                self.assertTrue(os.path.isfile(os.path.join(j,"INPUT")))
                self.assertTrue(os.path.isfile(os.path.join(j,"KPT")))

                kpt_text = Path(os.path.join(j,"KPT")).read_text()
                kpt = [ int(i) for i in kpt_text.strip().split("\n")[-1].split()[:3]]
                self.assertEqual(kpt,[2,2,2]) 
                
                stru_text = Path(os.path.join(j,"STRU")).read_text().split("\n")
                #read LATTICE_CONSTANT, the value is in the next line of key word "LATTICE_CONSTANT"
                lc = 0
                for i,line in enumerate(stru_text):
                    if re.search(r"LATTICE_CONSTANT",line):
                        lc = float(stru_text[i+1].strip().split()[0])
                        break
                self.assertTrue(lc in lc_ref)
                lcs.append(lc)
            lcs.sort()
            self.assertEqual(lcs,lc_ref)
    
    def test_prepare_stru_template(self):
        # case4: use STRU template, and USE PP/ORB LIB
        # when specify pp_path or orb_path, prepare will firstly search pporb in lib
        # and if there is no pporb, then use pporb defined STRU
        param_setting = {
            "example_template": ["a","b","c"],
            "stru_template": "a/STRU",
            "pp_path": "pplib",
            "orb_path": "orblib"
        }
        all_path_setting = prepare.DoPrepare(param_setting,"abacustest")  
        
        # check result folder
        self.assertTrue(os.path.isdir("abacustest"))
        
        # check new create examples, should create 3 examples for each template
        print(all_path_setting)
        self.assertEqual(len(all_path_setting), 3)
        self.assertEqual(len(all_path_setting[0]), 1)
        self.assertEqual(len(all_path_setting[1]), 1)
        self.assertEqual(len(all_path_setting[2]), 1)

        for i in all_path_setting:
            for j in i:
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.upf1")))
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.orb1")))
                self.assertTrue(os.path.isfile(os.path.join(j,"STRU")))
                self.assertTrue(os.path.isfile(os.path.join(j,"INPUT")))
                self.assertTrue(os.path.isfile(os.path.join(j,"KPT")))
                
                self.assertEqual(Path(os.path.join(j,"Si.upf1")).read_text(),"Si.upf1")
                self.assertEqual(Path(os.path.join(j,"Si.orb1")).read_text(),"Si.orb1")
                stru_text = Path(os.path.join(j,"STRU")).read_text().split("\n")
                #read LATTICE_CONSTANT, the value is in the next line of key word "LATTICE_CONSTANT"
                lc = 0
                pp = orb = ""
                for i,line in enumerate(stru_text):
                    if re.search(r"LATTICE_CONSTANT",line):
                        lc = float(stru_text[i+1].strip().split()[0])
                    elif re.search(r"ATOMIC_SPECIES",line):
                        pp = stru_text[i+1].strip().split()[2]
                    elif re.search(r"NUMERICAL_ORBITAL",line):
                        orb = stru_text[i+1].strip().split()[0]
            
                self.assertEqual(lc, 10.2)
                self.assertTrue(pp in ["Si.upf1","./Si.upf1"])
                self.assertTrue(orb in ["Si.orb1","./Si.orb1"])
    
    def test_prepare_specify_pporb(self):
        # case4: use kpt template, and USE PP/ORB LIB
        param_setting = {
            "example_template": ["a","b","c"],
            "pp_dict": {"Si": "pplib/Si.upf1"},
            "orb_dict": {"Si": "a/Si.orb"}
        }
        all_path_setting = prepare.DoPrepare(param_setting,"abacustest")  
        
        # check result folder
        self.assertTrue(os.path.isdir("abacustest"))
        
        # check new create examples, should create 3 examples for each template
        print(all_path_setting)
        self.assertEqual(len(all_path_setting), 3)
        self.assertEqual(len(all_path_setting[0]), 1)
        self.assertEqual(len(all_path_setting[1]), 1)
        self.assertEqual(len(all_path_setting[2]), 1)

        for i in all_path_setting:
            for j in i:
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.upf1")))
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.orb")))
                self.assertTrue(os.path.isfile(os.path.join(j,"STRU")))
                self.assertTrue(os.path.isfile(os.path.join(j,"INPUT")))
                self.assertTrue(os.path.isfile(os.path.join(j,"KPT")))
                
                self.assertEqual(Path(os.path.join(j,"Si.upf1")).read_text(),"Si.upf1")
                self.assertEqual(Path(os.path.join(j,"Si.orb")).read_text(),"Si.orb")
                
                stru_text = Path(os.path.join(j,"STRU")).read_text().split("\n")
                #read LATTICE_CONSTANT, the value is in the next line of key word "LATTICE_CONSTANT"
                lc = 0
                pp = orb = ""
                for i,line in enumerate(stru_text):
                    if re.search(r"LATTICE_CONSTANT",line):
                        lc = float(stru_text[i+1].strip().split()[0])
                    elif re.search(r"ATOMIC_SPECIES",line):
                        pp = stru_text[i+1].strip().split()[2]
                    elif re.search(r"NUMERICAL_ORBITAL",line):
                        orb = stru_text[i+1].strip().split()[0]
            
                self.assertEqual(lc, 10.2)
                self.assertTrue(pp in ["Si.upf1","./Si.upf1"])
                self.assertTrue(orb in ["Si.orb","./Si.orb"])
    
    def test_prepare_dpks(self):
        # case4: use kpt template, and USE PP/ORB LIB
        Path("INPUT").write_text("INPUT_PARAMETERS\nbasis_type lcao\ndeepks_scf 1\ndeepks_model dpks.ptg")
        Path("dpks.ptg").write_text("dpks.ptg")
        Path("jle.orb").write_text("jle.orb")
        param_setting = {
            "example_template": ["a","b","c"],
            "input_template": "INPUT",
            "dpks_descriptor": "jle.orb",
            "extra_files": ["dpks.ptg"]
        }
        all_path_setting = prepare.DoPrepare(param_setting,"abacustest")  
        
        # check result folder
        self.assertTrue(os.path.isdir("abacustest"))
        
        # check new create examples, should create 3 examples for each template
        print(all_path_setting)
        self.assertEqual(len(all_path_setting), 3)
        self.assertEqual(len(all_path_setting[0]), 1)
        self.assertEqual(len(all_path_setting[1]), 1)
        self.assertEqual(len(all_path_setting[2]), 1)

        for i in all_path_setting:
            for j in i:
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.upf")))
                self.assertTrue(os.path.isfile(os.path.join(j,"Si.orb")))
                self.assertTrue(os.path.isfile(os.path.join(j,"STRU")))
                self.assertTrue(os.path.isfile(os.path.join(j,"INPUT")))
                self.assertTrue(os.path.isfile(os.path.join(j,"KPT")))
                self.assertTrue(os.path.isfile(os.path.join(j,"dpks.ptg")))
                self.assertTrue(os.path.isfile(os.path.join(j,"jle.orb")))
                
                stru_text = Path(os.path.join(j,"STRU")).read_text().split("\n")
                #read LATTICE_CONSTANT, the value is in the next line of key word "LATTICE_CONSTANT"
                dpks_descriptor = ""
                for i,line in enumerate(stru_text):
                    if re.search(r"NUMERICAL_DESCRIPTOR",line):
                        dpks_descriptor = stru_text[i+1].strip().split()[0]
                        break
                self.assertTrue(dpks_descriptor in ["jle.orb","./jle.orb"])       
    