import os,sys,argparse,glob,json,traceback,re,shutil
import numpy as np
from typing import Dict, List, Optional, Union
import copy
from pathlib import Path
from abacustest import constant
from abacustest.myflow import comm
from abacustest.lib_prepare import abacus as MyAbacus
from abacustest.lib_prepare import abacus2qe as Aba2Qe
from abacustest.lib_prepare import abacus2vasp as Aba2Vasp
from abacustest.lib_prepare import abacus2cp2k as Aba2Cp2k
from abacustest.lib_prepare.comm import translate_strus, collect_pp
 
class PrepareAbacus:
    def __init__(self,
                 save_path: str = Path("."),
                 example_template: str = None,
                 input_template: str = None,
                 kpt_template: str = None,
                 stru_template: str = None,
                 mix_input: Dict[str,any] = {},
                 mix_kpt: List[Union[int,List[int]]] = [],
                 mix_stru: List[str]=[],
                 pert_stru: Dict[str,any] = {},
                 pp_dict: Dict[str,str]= {},
                 orb_dict: Dict[str,str]= {},
                 paw_dict: Dict[str,str]= {},
                 pp_path:str = None,
                 orb_path:str = None,
                 dpks_descriptor:str=None,
                 extra_files: List[str] = [],
                 bak_file = True,
                 no_link = False,
                 link_example_template_extra_files=True,):
        """To prepare the inputs of abacus

        Parameters
        ----------
        save_path : str, optional
            the path to store inputs, by default Path(".")

        example_template : str, optional
            folder name of an abacus inputs template, by default None.
            If specify input/kpt/stru_template, the related file in exmaple_template will
            be ingored. The pp/orb defined in example_template/STRU will also added to 
            pp/orb_dict.

        input_template : str, optional
            file name of INPUT template, by default None

        kpt_template : str, optional
            file name of KPT template, by default None

        stru_template : str, optional
            file name of STRU template, by default None

        mix_input : Dict[str,any], optional
            mix setting of input parameters, by default {}. 
            Such as:{"ecutwfc":[50, 60, 70], 
                     "mixing_beta":[0.3, 0.4, 0.5]}, 
            then will prepare 3*3 = 9 inputs.

        mix_kpt : List, optional
            mixing setting of k point list, by default []. If this parameter
            is defined, then kpt_template will be invalide.
            Such as: [[1,1,1,0,0,0],
                      [2,2,2,0,0,0]], 
            then will prepare 2 KPT files.
            The kpoint list can be a list of 6 values (3 int number for K point
            and 3 float values for the shift value) like: [1,1,1,0,0,0], 
            or 3 values like: [1,1,1] which means [1,1,1,0,0,0], or one int number
            like: 3 which means [3,3,3,0,0,0].

        mix_stru : List, optional
            a list of STRU files, by default []. If this parameter is defined
            then stru_template will be invalide.
            Such as: ["a/STRU","b/STRU"]
            
        pert_stru : Dict, optional
            perturbation setting of structure, by default {}.
            
        pp_dict : Dict, optional
            a dictionary specify the pseudopotential files, by default {}

        orb_dict : Dict, optional
            a dictionary specify the orbital files, by default {}
        
        paw_dict: Dict, optional
            a dictionary specify the paw files, by default {}

        pp_path : str, optional
            the path of your pseudopotential lib, by default None. 
            The file name shuold be started with the element name and followed by "_".
            Such as: Ag_ONCV_PBE-1.0.upf.

        orb_path : str, optional
            the path of your orbital lib, by default None. 
            The file name shuold be started with the element name and followed by "_".
            Such as: Ag_gga_7au_100Ry_4s2p2d1f.orb.

        dpks_descriptor : str, optional
            file name of deepks descriptor, by default None

        extra_files : List, optional
            a list of somes extra files that will be putted in each input, by default []   
            
        bak_file: bool, optional
            when save path is exist, if bak_file is True, then bak the old files,
            
        no_link: bool, optional
            if True, then will not link the files in example_template to save_path,
            
        link_example_template_extra_files: bool, optional
            if True, then will link the extra files in example_template to save_path,
        """        
        self.save_path = save_path
        self.example_template = example_template
        self.input_template = input_template
        self.kpt_template = kpt_template
        self.stru_template = stru_template
        self.mix_input = mix_input
        self.pp_dict = pp_dict
        self.orb_dict = orb_dict
        self.paw_dict = paw_dict
        self.extra_files = self.CheckExtrafile(extra_files)
        self.mix_kpt = mix_kpt
        self.mix_stru = mix_stru
        self.pert_stru = pert_stru
        self.dpks_descriptor = dpks_descriptor if dpks_descriptor else None
        
        self.pp_path = None if not pp_path else pp_path.strip()
        self.orb_path = None if not orb_path else orb_path.strip()
        self.CollectPP()
        self.CollectOrb()
        
        self.input_list,self.input_mix_param = self.Construct_input_list()
        self.kpt_list = self.Construct_kpt_list()
        self.stru_list = self.Construct_stru_list()
        
        self.bak_file = bak_file
        self.no_link = no_link
        
        # when read the template file, will check if the template file folder is same with save_path
        # if same and final structures only one, then will not bak the template folder
        self.template_is_save_path = self.CheckIfTemplateIsSavePath()
        self.example_template_extra_files = self.CheckExampleTemplateExtraFiles(link_example_template_extra_files)
    
    def CheckExampleTemplateExtraFiles(self,link_example_template_extra_files):
        # check the extra files in example_template except for INPUT, KPT, STRU
        if (not link_example_template_extra_files) or (not self.example_template):
            return []
        if not os.path.isdir(self.example_template):
            return []
        allfiles = os.listdir(self.example_template)
        extrafiles = []
        for ifile in allfiles:
            if ifile in ["INPUT","STRU","KPT"]:
                continue
            if os.path.isfile(os.path.join(self.example_template,ifile)):
                extrafiles.append(os.path.abspath(os.path.join(self.example_template,ifile)))
        return extrafiles
        
    def CheckIfTemplateIsSavePath(self):
        if not os.path.exists(self.save_path):
            return False
        CheckedPath = []
        if self.input_template:
            CheckedPath.append(Path(os.path.split(self.input_template)[0]))
        if self.kpt_template:
            CheckedPath.append(Path(os.path.split(self.kpt_template)[0]))
        if self.stru_template:
            CheckedPath.append(Path(os.path.split(self.stru_template)[0]))
        if self.example_template:
            CheckedPath.append(Path(self.example_template))
        for ipath in CheckedPath:
            if Path(self.save_path).samefile(ipath):
                return True
        return False

    def CheckExtrafile(self,extrafiles):
        if not extrafiles:
            return []
        filelist = []
        for ifile in extrafiles:
            if os.path.isfile(ifile):
                filelist.append(ifile)
            else:
                print(f"ERROR: Can not find {ifile}!")
        return filelist

    def CollectPP(self):
        for k,v in collect_pp(self.pp_path).items():
            if k not in self.pp_dict:
                self.pp_dict[k] = v

    def CollectOrb(self):
        for k,v in collect_pp(self.orb_path).items():
            if k not in self.orb_dict:
                self.orb_dict[k] = v
             
    def CollectPaw(self):
        for k,v in collect_pp(self.paw_path).items():
            if k not in self.paw_dict:
                self.paw_dict[k] = v
            
    def parse_combine_input(self,input_dict, combine_str="|"):
        #parse the combine input, such as "ecutwfc|kspacing": "50|0.2"
        new_input_dict = {}
        for key,value in input_dict.items():
            if combine_str in key and isinstance(value,str):
                keys = key.split(combine_str)
                vs = value.split(combine_str)
                if len(keys) != len(vs):
                    print(f"ERROR: number of keys and values are not equal in {key}={value}")
                    new_input_dict[key] = value
                else:
                    for i in range(len(keys)):
                        new_input_dict[keys[i]] = vs[i]
            else:
                new_input_dict[key] = value
        return new_input_dict
            
    def Construct_input_list(self):
        all_inputs = []
        inputf = None
        if self.example_template != None:
            if os.path.isfile(os.path.join(self.example_template,"INPUT")):
                inputf = os.path.join(self.example_template,"INPUT")
        if self.input_template != None:
            if os.path.isfile(self.input_template):
                inputf = self.input_template 
            else:
                print("WARNING: File %s is not found" % self.input_template)

        print("INPUT: %s" % str(inputf))
        input_constant = {}
        if inputf != None:
            input_constant = PrepareAbacus.ReadInput(inputf)
        
        list_param = {}
        input_constant_common = {}
        for k,v in self.mix_input.items():
            #if value is list type, then we need prepare INPUT for each value
            #input_constant stores the parameters whose value is constant
            #list_param stores the parameters that has several values.
            if isinstance(v,(int,float,str)):
                input_constant[k] = v
                input_constant_common[k] = v
            elif isinstance(v,list):
                if len(v) == 1:
                    input_constant[k] = v[0]
                    input_constant_common[k] = v[0]
                else:    
                    list_param[k] = v
                    if k in input_constant:
                        del input_constant[k]
            else:
                print("WARNING: type of '%s' is" % str(v),type(v),"will not add to INPUT")
                input_constant[k] = v
                
        input_constant = self.parse_combine_input(input_constant) # parse the combine input
        print("Invariant INPUT setting:",str(input_constant))

        all_inputs.append(input_constant)
        for k,v in list_param.items():
            for i,iv in enumerate(v):
                if i == 0:
                    #if v only one value, add in all_inputs directly
                    for iinput in all_inputs:
                        iinput[k] = iv   
                else:
                    #if v has more than two values, deepcopy all_inputs and modify v, and then add to all_inputs
                    if i == 1:  
                        tmp_inputs = copy.deepcopy(all_inputs)
                    for iinput in tmp_inputs:
                        iinput[k] = iv
                    all_inputs += copy.deepcopy(tmp_inputs)
        list_param_new = []
        for k in list_param:
            if "|" in k:
                list_param_new += k.split("|")
            else:
                list_param_new.append(k)
        return [self.parse_combine_input(i) for i in all_inputs],list_param_new

    def Construct_kpt_list(self):
        """
        all_kpt is a list of list (6 element) or filename of KPT.
        If not specify mix_kpt and has specify the KPT template, will retrun a list
        of KPT filename, or will return a null list if no KPT template.
        """
        all_kpt = []

        kptf = None
        if self.example_template != None:
            if os.path.isfile(os.path.join(self.example_template,"KPT")):
                kptf = os.path.join(self.example_template,"KPT")
        if self.kpt_template != None:
            if os.path.isfile(self.kpt_template):
                kptf = self.kpt_template 
            else:
                print("WARNING: File %s is not found" % self.kpt_template)

        if len(self.mix_kpt) == 0:
            if kptf != None:
                all_kpt.append(kptf)
            template_file = kptf
        else:
            for ikpt in self.mix_kpt:
                if isinstance(ikpt,int):
                    all_kpt.append([ikpt,ikpt,ikpt,0,0,0])
                elif isinstance(ikpt,list):
                    if len(ikpt) == 3:
                        all_kpt.append(ikpt+[0,0,0])
                    elif len(ikpt) == 6:
                        all_kpt.append(ikpt)
                    else:
                        print("mix_kpt should be a int or list of 3/6 elements, but not ",ikpt)
                elif isinstance(ikpt,str):
                    allkpt = glob.glob(ikpt)
                    allkpt.sort()
                    for iikpt in allkpt:
                        all_kpt.append(iikpt)
                else:
                    print("mix_kpt should be int or list of 3/6 elements, but not ",ikpt)
            template_file = self.mix_kpt
        print(f"KPT: {template_file}")
        return all_kpt
    
    def Construct_stru_list(self):
        all_stru = []
        struf = None
        if self.example_template != None:
            if os.path.isfile(os.path.join(self.example_template,"STRU")):
                struf = os.path.join(self.example_template,"STRU")
        if self.stru_template != None:
            if os.path.isfile(self.stru_template):
                struf = self.stru_template 
            else:
                print("WARNING: File %s is not found" % self.stru_template)

        if len(self.mix_stru) == 0:
            if struf != None:
                all_stru.append(struf)
            else:
                print("Please set the stru_template!")
            template_file = struf
        else:
            for stru in self.mix_stru:
                allstrus = glob.glob(stru)
                allstrus.sort()
                for istru in allstrus:
                    all_stru.append(os.path.abspath(istru))
                if len(allstrus) == 0:
                    print("Structure file '%s' is not exist" % stru)
            template_file = self.mix_stru
        
        print(f"STRU: {template_file}")
        return all_stru
    
    def prepare(self):
        "Return a dict, and the key is the path of inputs,"
        "and the value is list of structure, kpt, and input setting information."
        "The input setting information is a dict of param name and value"

        if not self.stru_list:
            print("No stru files, skip!!!")
            return None
        
        if not self.kpt_list:
            print("WARNING: not set KPT")
            kpt_list = [None]
        else:
            kpt_list = self.kpt_list
        
        if not self.input_list:
            input_list = [None]
        else:
            input_list = self.input_list 
                   
        ipath = -1
        param_setting = {}
        cwd = os.getcwd()
        stru_num = len(self.stru_list) * len(kpt_list) * len(input_list)
        for istru in self.stru_list:  #iteration of STRU
            stru_data = MyAbacus.AbacusStru.ReadStru(istru)
            if not stru_data:
                print(f"ERROR: read STRU failed in {istru}")
                continue
            stru_path = os.path.split(istru)[0]
            if stru_path == "": stru_path = os.getcwd()
            labels = []
            for ilabel in stru_data.get_label():
                if ilabel not in labels:
                    labels.append(ilabel)
            linkstru = True
            skipstru = False
            allfiles = self.example_template_extra_files + self.extra_files  #files that will be linked
            pp_list = []   #pp file name
            orb_list = []  #orb file name
            paw_list = [] #paw file name
            dpks = None    #dpks file name
            input_in_stru_path = {}
            if os.path.isfile(os.path.join(stru_path,"INPUT")):
                input_in_stru_path = PrepareAbacus.ReadInput(os.path.join(stru_path,"INPUT"))
            print(input_in_stru_path)
            for i,ilabel in enumerate(labels):
                #check pp file 
                if ilabel not in self.pp_dict:
                    #print("label '%s' is found in '%s', but not defined in pp_dict." % (ilabel,istru))
                    if stru_data.get_pp():
                        os.chdir(stru_path)
                        pp_in_stru = os.path.join(input_in_stru_path.get("pseudo_dir",""),stru_data.get_pp()[i])
                        
                        if os.path.isfile(pp_in_stru):
                            print("label '%s': link the pseudopotential file '%s' defined in %s" % (ilabel,pp_in_stru,istru))
                            pp_list.append(os.path.basename(pp_in_stru))
                            allfiles.append(os.path.abspath(pp_in_stru))
                            if pp_list[-1] != stru_data.get_pp()[i]:
                                linkstru = False
                        else:
                            pp_list.append(stru_data.get_pp()[i])
                            print("label '%s': the pseudopotential file '%s' defined in %s is not found" % (ilabel,pp_in_stru,istru))
                            #skipstru = True
                            #os.chdir(cwd)
                            #break
                        os.chdir(cwd)
                else:
                    pp_list.append(os.path.basename(self.pp_dict[ilabel]))  #only store the file name to pp_list
                    allfiles.append(self.pp_dict[ilabel]) #store the whole pp file to allfiles                    
                    if not stru_data.get_pp() or pp_list[-1] != stru_data.get_pp()[i]:
                        linkstru = False

                #check orbital file
                if ilabel not in self.orb_dict:
                    if stru_data.get_orb():
                        #print("label '%s' is found in '%s', but not defined in orb_dict." % (ilabel,istru))
                        os.chdir(stru_path)
                        orb_in_stru = os.path.join(input_in_stru_path.get("orbital_dir",""),stru_data.get_orb()[i])
                        if os.path.isfile(orb_in_stru):
                            print("label '%s': link the orb file '%s' defined in %s" % (ilabel,orb_in_stru,istru))
                            orb_list.append(os.path.basename(orb_in_stru))
                            allfiles.append(os.path.abspath(orb_in_stru))
                            if orb_list[-1] != stru_data.get_orb()[i]:  
                                linkstru = False
                        else:
                            orb_list.append(stru_data.get_orb()[i])
                            print("label '%s': the orbital file '%s' defined in %s is not found." % (ilabel,orb_in_stru,istru))
                        os.chdir(cwd)
                else:
                    orb_list.append(os.path.basename(self.orb_dict[ilabel])) 
                    allfiles.append(self.orb_dict[ilabel])
                    if not stru_data.get_orb() or orb_list[-1] != stru_data.get_orb()[i]:  
                        linkstru = False
                
                #check paw file
                if ilabel not in self.paw_dict:
                    if stru_data.get_paw():
                        os.chdir(stru_path)
                        paw_in_stru = stru_data.get_paw()[i]
                        if os.path.isfile(paw_in_stru):
                            print("label '%s': link the paw file '%s' defined in %s" % (ilabel,paw_in_stru,istru))
                            paw_list.append(os.path.basename(paw_in_stru))
                            allfiles.append(os.path.abspath(paw_in_stru))
                            if paw_list[-1] != stru_data.get_paw()[i]:  
                                linkstru = False
                        else:
                            print("label '%s': the paw file '%s' defined in %s is not found." % (ilabel,paw_in_stru,istru))
                        os.chdir(cwd)
                else:
                    paw_list.append(os.path.basename(self.paw_dict[ilabel])) 
                    allfiles.append(self.paw_dict[ilabel])
                    if not stru_data.get_paw() or paw_list[-1] != stru_data.get_paw()[i]:  
                        linkstru = False

            if skipstru:
                continue
            
            #check dpks 
            if self.dpks_descriptor:
                if os.path.isfile(self.dpks_descriptor):
                    dpks = os.path.basename(self.dpks_descriptor)
                    allfiles.append(self.dpks_descriptor)
                else:
                    print("Error: Can not find file %s, skip the prepare of dpks_descriptor" % self.dpks_descriptor)
                    dpks = None
            else:
                os.chdir(stru_path)
                dpks = stru_data.get_dpks()
                if dpks:
                    if os.path.isfile(dpks):
                        allfiles.append(os.path.abspath(dpks))
                        dpks = os.path.basename(dpks)
                    else:
                        print("Error: deepks descriptor is defined in %s/STRU, but can not find the file, skip the prepare" % stru_path)
                        dpks = None
                else:
                    dpks = None
                os.chdir(cwd)
            if dpks != stru_data.get_dpks():
                linkstru = False

            #stru_data will be writen in new folder
            stru_data.set_pp(pp_list)
            stru_data.set_orb(orb_list)
            stru_data.set_paw(paw_list)
            stru_data.set_dpks(dpks)
            
            for ikpt in kpt_list:  #iteration of KPT 
                for iinput in input_list: #iteration of INPUT
                    if self.pert_stru and self.pert_stru.get("pert_number",0) > 0:
                        stru_data_pert = stru_data.perturb_stru(self.pert_stru.get("pert_number",0),
                                                                cell_pert_frac=self.pert_stru.get("cell_pert_frac",None),
                                                                atom_pert_dist=self.pert_stru.get("atom_pert_dist",None),
                                                                mag_rotate_angle=self.pert_stru.get("mag_rotate_angle",None),
                                                                mag_tilt_angle=self.pert_stru.get("mag_tilt_angle",None),
                                                                mag_norm_dist=self.pert_stru.get("mag_norm_dist",None))
                        linkstru = False
                        create_subpath = True
                    else:
                        stru_data_pert = [stru_data]
                        create_subpath = bool(stru_num > 1)
                    
                    for stru_datai in stru_data_pert:
                        ipath += 1
                        self.write_one_example(create_subpath,ipath,cwd,param_setting,istru,ikpt,iinput,linkstru,stru_datai,allfiles)        
        
        return param_setting

    def write_one_example(self,create_subpath,ipath,cwd,param_setting,stru_path,ikpt,iinput,linkstru,stru_data,allfiles):
        """
        create_subpath: bool, if True, then create subpath for each example. This is True when there are more than one example or more than one inputs/kpt/stru setting.
        ipath: int, the index for example_template
        cwd: str, the current working directory, used to check if the save_path is current path and get the relative path of stru file
        param_setting: dict, to store the param setting of each new abacus job
        stru_path: str, the path of STRU file
        ikpt: KPT setting
        iinput: INPUT setting
        linkstru: bool, if True, then link the STRU file, else copy the STRU file
        stru_data: AbacusStru, the structure data
        allfiles: List, the list of all files that will be linked to the save
        """
        if create_subpath:
                save_path = os.path.join(self.save_path,str(ipath).zfill(5))
        else:
            save_path = self.save_path
        
        if os.path.isdir(save_path) and \
            (not Path(save_path).samefile(cwd)) and \
            (self.bak_file) and \
            (not self.template_is_save_path):
                bk = comm.GetBakFile(save_path)
                shutil.move(save_path,bk)
                
        if not os.path.isdir(save_path):
            os.makedirs(save_path)   
        #store the param setting and path
        param_setting[save_path] = [os.path.relpath(stru_path, cwd),ikpt,{}]
        for input_param in self.input_mix_param:
            param_setting[save_path][-1][input_param] = iinput.get(input_param)
        
        #create INPUT   
        if iinput != None: 
            # if pseudo_dir and orbital_dir is defined in INPUT, then remove them
            iinput.pop("pseudo_dir","")
            iinput.pop("orbital_dir","")  
            PrepareAbacus.WriteInput(iinput,os.path.join(save_path,"INPUT"))
        
        #create KPT
        if ikpt != None:
            kptf = os.path.join(save_path,"KPT")
            if isinstance(ikpt,str):
                if os.path.isfile(kptf) and (not Path(ikpt).samefile(Path(kptf))):
                    os.unlink(kptf)
                
                if not os.path.isfile(kptf):
                    if self.no_link:
                        shutil.copy(os.path.abspath(ikpt),kptf)
                    else:
                        os.symlink(os.path.abspath(ikpt),kptf)
            elif isinstance(ikpt,list):
                PrepareAbacus.WriteKpt(ikpt,kptf)
        
        #create STRU
        struf = os.path.join(os.path.join(save_path,"STRU"))
        if linkstru:
            if os.path.isfile(struf) and (not Path(stru_path).samefile(Path(struf))):
                os.unlink(struf)
                    
            if not os.path.isfile(struf):
                if self.no_link:
                    shutil.copy(os.path.abspath(stru_path),struf)
                else:
                    os.symlink(os.path.abspath(stru_path),struf)
        else:
            stru_data.write(struf)
        
        #link other files
        for ifile in allfiles:
            ifile = os.path.abspath(ifile)
            filename = os.path.basename(ifile)
            target_file = os.path.join(save_path,filename)
            if os.path.isfile(target_file) and not Path(ifile).samefile(Path(target_file)):
                os.unlink(target_file)
            if not os.path.isfile(target_file):
                if self.no_link:
                    shutil.copy(ifile,target_file)
                else:
                    os.symlink(ifile,target_file)  
    
    @staticmethod
    def WriteKpt(kpoint_list:List = [1,1,1,0,0,0],file_name:str = "KPT", model:str = "gamma"):
        MyAbacus.WriteKpt(kpoint_list,file_name,model)
    
    @staticmethod
    def ReadInput(INPUTf: str = None, input_lines: str = None) -> Dict[str,any]:
        return MyAbacus.ReadInput(INPUTf,input_lines)

    @staticmethod
    def WriteInput(input_context:Dict[str,any],
                   INPUTf: str = "INPUT"):
        MyAbacus.WriteInput(input_context,INPUTf)           

def CheckExample(example_path,description:str=""):
    print(f"Check ABACUS inputs: {example_path} {description}")
    if not os.path.isdir(example_path):
        print("Can not find example path %s" % example_path)
        return False
    if not os.path.isfile(os.path.join(example_path,"INPUT")):
        print("Can not find INPUT file in %s" % example_path)
        return False
    if not os.path.isfile(os.path.join(example_path,"STRU")):
        print("Can not find STRU file in %s" % example_path)
        return False
    iinput = PrepareAbacus.ReadInput(os.path.join(example_path,"INPUT"))
    pp_path = iinput.get("pseudo_dir","")
    orb_path = iinput.get("orbital_dir","")
    paw_path = iinput.get("paw_dir","")
    use_paw = iinput.get("use_paw",False)
    if use_paw and use_paw.isdigit() and int(use_paw):
        use_paw = True
    else:
        use_paw = False
    istru = MyAbacus.AbacusStru.ReadStru(os.path.join(example_path,"STRU"))
    if not istru:
        print("Read STRU file failed in %s" % example_path)
        return False
    
    cwd = os.getcwd()
    allpass = True
    os.chdir(example_path)
    #check pp
    if not use_paw:
        ppfiles = istru.get_pp()
        if not ppfiles:
            print("Can not find PP in %s" % example_path)
            allpass = False
        else:
            for ppfile in ppfiles:
                real_ppfile = os.path.join(pp_path,ppfile)
                if not os.path.isfile(real_ppfile):
                    print("Can not find PP file %s in %s" % (real_ppfile,example_path))
                    allpass = False
    else:
        pawfiles = istru.get_paw()
        if not pawfiles:
            print("Can not find PAW in %s" % example_path)
            allpass = False
        else:
            for pawfile in pawfiles:
                real_pawfile = os.path.join(paw_path,pawfile)
                if not os.path.isfile(real_pawfile):
                    print("Can not find PAW file %s in %s" % (real_pawfile,example_path))
                    allpass = False
    
    #check orb
    basis = iinput.get("basis_type","pw").lower()
    if basis != "pw" and basis.startswith("lcao"):
        orbfiles = istru.get_orb()
        if not orbfiles:
            print("Can not find ORB in %s/STRU" % example_path)
            allpass = False
        else:
            for orbfile in orbfiles:
                real_orbfile = os.path.join(orb_path,orbfile)
                if not os.path.isfile(real_orbfile):
                    print("Can not find ORB file %s in %s" % (real_orbfile,example_path))
                    allpass = False
    
    #check dpks
    deepks_scf = iinput.get("deepks_scf",False)
    if deepks_scf:
        dpks_model = iinput.get("deepks_model",None)
        if not dpks_model:
            print("deepks_model is not defined in %s/INPUT" % example_path)
            allpass = False
        else:
            if not os.path.isfile(dpks_model):
                print("Can not find dpks model file %s in %s" % (dpks_model,example_path))
                allpass = False
    deepks_out_labels = iinput.get("deepks_out_labels",None)
    if deepks_scf or deepks_out_labels:
        dpks_descriptor = istru.get_dpks()
        if not dpks_descriptor:
            print("Read dpks descriptor failed in %s/STRU" % example_path)
            allpass = False
        elif not os.path.isfile(dpks_descriptor):
            print("Can not find dpks descriptor file %s in %s" % (dpks_descriptor,example_path))
            allpass = False
    
    # check KPT
    # only basis is lcao and use gamma_only, or has define kspacing, KPT file is not needed
    gamma_only = iinput.get("gamma_only",False)
    kspacing = iinput.get("kspacing",None)
    if not ((basis.startswith("lcao") and gamma_only) or kspacing): 
        if not os.path.isfile("KPT"):
            print("Can not find KPT file in %s" % example_path)
            allpass = False 
    os.chdir(cwd)
    return allpass 

def CommPath(pathlist):
    "return (commpath pathlist_without_commpath)"
    if len(pathlist) == 0:
        return None,pathlist
    commpath = os.path.commonpath(pathlist)
    if commpath != "":
        length = len(commpath) + 1
    else:
        length = 0
    return commpath,[i[length:] for i in pathlist]

def DoPrepare(param_setting: Dict[str, any], save_folder: str, no_link: bool = False) -> List[Dict[str, dict]]:  
    """
    param_setting is a dictionary like:
    {
        "example_template":["example_path"]
        "input_template":"INPUT",
        "kpt_template":"KPT",
        "stru_template":"STRU",
        "strus": str/list,  # the stru file of other format, like cif, vasp, etc.
        "stru_format": "cif",  # the format of stru file, should be cif or dpdata supportted format, if strus is setted, then example_template is not needed.
        "mix_input":{
            "ecutwfc":[50,60,70],
            "kspacing":[0.1,0.12,0.13]
        },
        "mix_kpt":[],
        "mix_stru":[],
        "pert_stru":{
            "pert_number": 0,
            "cell_pert_frac": null,
            "atom_pert_dist": null,
            "mag_rotate_angle": null,
            "mag_tilt_angle": null,
            "mag_norm_dist": null
            },
        "pp_dict":{},
        "orb_dict":{},
        "pp_path": str,
        "orb_path": str,
        "dpks_descriptor":"",
        "extra_files":[],
        "link_example_template_extra_files":true,
        "bak_file":true,
        "abacus2qe": bool,
        "qe_setting":{},
        "abacus2vasp": bool,
        "potcar": str,
        "vasp_setting":{},
        "abacus2cp2k": bool,
        "cp2k_setting":{},
        "mix_input_comment":"Do mixing for several input parameters. The diffrent values of one parameter should put in a list.",
        "mix_kpt_comment":"If need the mixing of several kpt setting. The element should be an int (eg 2 means [2,2,2,0,0,0]), or a list of 3 or 6 elements (eg [2,2,2] or [2,2,2,0,0,0]).",
        "mix_stru_comment":"If need the mixing of several stru files. Like: ["a/stru","b/stru"], 
        "bak_file_comment":"If need to bak the old folder if the new folder is exist.",      
    },

    save_folder specifies the destination of the new created examples.
    
    if strus and stru_format is setted, then will firstly convert it to ABACUS format, and save to new folder named by %06d.

    Return a list of dict, which is related to each example_template element. The key of the dict is the newcreated example path, and the value
    is the STRU/KPT/INPUT settings. 
    """
    strus = param_setting.get("strus",None)
    stru_format = param_setting.get("stru_format",None)
    if strus is not None:
        assert stru_format is not None, "If strus is setted, then stru_format should be setted."
        example_template = translate_strus(strus,stru_format)
        if example_template is None:
            print("Convert strus to ABACUS format failed.")
            return []
        if param_setting.get("example_template",None) is not None:
            print("WARNING: strus is setted, example_template is not used.")
    else:
        example_template = param_setting.get("example_template",None)
        
    if isinstance(example_template,str):
        example_template = glob.glob(example_template)
        # if __MACOSX folder is exist, then remove it
        if "__MACOSX" in example_template:
            example_template.remove("__MACOSX")
        if "__MACOSX/" in example_template:
            example_template.remove("__MACOSX/")
        example_template.sort()
    elif isinstance(example_template,list):
        tmp = copy.deepcopy(example_template)
        example_template = []
        for i in tmp:
            alli = glob.glob(i)
            # if __MACOSX folder is exist, then remove it
            if "__MACOSX" in alli:
                alli.remove("__MACOSX")
            if "__MACOSX/" in alli:
                alli.remove("__MACOSX/")
            alli.sort()
            example_template += alli
    else:
        example_template = [example_template]
    
    all_path_setting = []
    commpath,example_template_nocomm = CommPath(example_template)
    print(commpath,example_template_nocomm)
    for idx,iexample in enumerate(example_template):
        if not os.path.isdir(iexample):
            print("Can not find folder %s" % iexample)
            continue
        #if len(example_template) > 1:
        #    save_path = os.path.join(save_folder,example_template_nocomm[idx])
        #    print("\n%s" % iexample)
        #else:
        #    if Path(save_folder) == Path("."):
        #        save_path = commpath
        #    else:
        #        save_path = save_folder
        save_path = os.path.join(save_folder,iexample)
        prepareabacus = PrepareAbacus(save_path=save_path,
                                  example_template=iexample,
                                  input_template=param_setting.get("input_template",None),
                                  kpt_template=param_setting.get("kpt_template",None),
                                  stru_template=param_setting.get("stru_template",None),
                                  mix_input=param_setting.get("mix_input",{}),
                                  mix_kpt=param_setting.get("mix_kpt",[]),
                                  mix_stru=param_setting.get("mix_stru",[]),
                                  pert_stru=param_setting.get("pert_stru",{}),
                                  pp_dict=param_setting.get("pp_dict",{}),
                                  orb_dict=param_setting.get("orb_dict",{}),
                                  pp_path=param_setting.get("pp_path",None),
                                  orb_path=param_setting.get("orb_path",None),
                                  dpks_descriptor=param_setting.get("dpks_descriptor",None),
                                  extra_files=param_setting.get("extra_files",[]),
                                  bak_file = param_setting.get("bak_file",True),
                                  no_link=no_link,
                                  link_example_template_extra_files=param_setting.get("link_example_template_extra_files",True)
                                  )
        all_path_setting.append(prepareabacus.prepare())
            
    if param_setting.get("abacus2qe",False):
        print("\nConvert ABACUS inputs to QE inputs")
        for isetting in all_path_setting:
            if isetting:
                for ipath in list(isetting.keys()):
                    print(ipath)
                    try:
                        Aba2Qe.Abacus2Qe(ipath,save_path=os.path.join(ipath,"input"),qe_param=param_setting.get("qe_setting",{}))
                    except:
                        traceback.print_exc()
    
    if param_setting.get("abacus2vasp",False):
        print("\nConvert ABACUS inputs to VASP inputs")
        potcar = param_setting.get("potcar",None)
        vasp_setting = param_setting.get("vasp_setting",{})
        for isetting in all_path_setting:
            if isetting:
                for ipath in list(isetting.keys()):
                    print(ipath)
                    try:
                        Aba2Vasp.Abacus2Vasp(ipath,save_path=ipath,potcar=potcar,vasp_setting=vasp_setting)
                    except:
                        traceback.print_exc()
    
    if param_setting.get("abacus2cp2k",False):
        print("\nConvert ABACUS inputs to CP2K inputs")
        cp2k_setting = param_setting.get("cp2k_setting",{})
        for isetting in all_path_setting:
            if isetting:
                for ipath in list(isetting.keys()):
                    print(ipath)
                    try:
                        Aba2Cp2k.Abacus2Cp2k(ipath,save_path=os.path.join(ipath,"cp2k.inp"),cp2k_setting=cp2k_setting)
                    except:
                        traceback.print_exc()
                      
    return all_path_setting

def PrepareInput(param):
    param_setting_file = param.param
    save_folder = param.save
    if not os.path.isfile(param_setting_file):
        print("Can not find file %s!!!" % param_setting_file)
        sys.exit(1)

    param_setting = json.load(open(param_setting_file))
    if "prepare" in param_setting:
        param_setting = param_setting["prepare"]
    
    all_path_setting = DoPrepare(param_setting,save_folder,param.nolink)
    for path_setting in all_path_setting:
        if path_setting:
            for k,v in path_setting.items():
                print("%s:%s" % (k,str(v)))

def PrepareArgs(parser):  
    parser.description = "This script is used to prepare the INPUTS OF ABACUS JOB"
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type',required=True)
    parser.add_argument('-s', '--save', type=str,  default="abacustest",help='where to store the inputs, default is abacustest ')
    parser.add_argument('--nolink',  nargs='?',type=int, const=1, default=0,help='if link the files in the example folder, default is 0')
    return parser

def main():
    parser = argparse.ArgumentParser()
    param = PrepareArgs(parser).parse_args()
    PrepareInput(param)
    
if __name__ == "__main__":
    main()
