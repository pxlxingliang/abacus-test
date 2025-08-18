import os,json,glob,shutil,traceback
import subprocess,copy
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,ReadKpt
import select
from typing import List, Dict, Union, Optional, Tuple, Literal
from abacustest.lib_prepare.comm import IsTrue

BOHRIUM_DES = '''If you use Bohrium to accelerate the calculation, you need to set below environment variables:
    export BOHRIUM_USERNAME=<your username> BOHRIUM_PASSWORD=<your password> BOHRIUM_PROJECT_ID=<your project id> 
Or you can add below block in your setting.json:
    "config": {
        "bohrium_username": "<your username>",
        "bohrium_password": "<your password>",
        "bohrium_project_id": "<your project id>"
    } 
Or you can use dispatcher to submit the calculation to the remote server,
and please refer to the dispatcher section for more details in https://github.com/pxlxingliang/abacus-test/tree/develop.'''

def doc_after_prepare(model_name, jobs, files,has_prepare=True):
    print(f"\nPrepare the {model_name} successfully!!!")
    print("\nYou can choose one of below methods to do the calculation:")
    print("1. run with abacustest locally")
    print(BOHRIUM_DES,)
    print("After finishing the setting, you can run below command to submit the calculation:")
    print("        abacustest submit -p setting.json &")
    print("\n2. run with bohrium app abacustest:")
    print("Compress the inputs and submit to https://app.bohrium.dp.tech/abacustest/?request=GET%3A%2Fapplications%2Fabacustest with Reuse Submodel.")
    print("Compress command:")
    print("    tar -zcvf %s.tar.gz %s" % (model_name, " ".join(jobs+files)))
    if has_prepare:
        print("\n3. run with your own script:")
        print("You can execute below command and abacustest will generate all inputs in folder abacustest/, and then you can submit the calculation by yourself.")
        outfolder = "" if len(jobs) > 1 else f"-s {jobs[0]}"
        print("    abacustest prepare -p setting.json " + outfolder + "\n\n")
    
def dump_setting(setting):
    """
    Dump the given setting dictionary to a JSON file named 'setting.json'.
    If 'setting.json' already exists, it will be backed up with a numbered suffix.

    Args:
        setting (dict): The setting dictionary to be dumped.

    Returns:
        None
    """
    if os.path.isfile("setting.json"):
        i = 1
        settingbakf = f"setting.json.bak{i}"
        while os.path.isfile(settingbakf):
            i += 1
            settingbakf = f"setting.json.bak{i}"
        print("Warning: setting.json exists. Will backup it to", settingbakf)
        os.system(f"mv setting.json {settingbakf}")
    json.dump(setting, open("setting.json", "w"), indent=4)

def get_physical_cores():
    """
    Get the number of physical cores in the system.

    Returns:
        int: The number of physical cores.

    Raises:
        CalledProcessError: If the command execution fails.
    """
    cmd = "lscpu | grep 'Core(s) per socket:' | awk '{print $4}'"
    cores = subprocess.check_output(cmd, shell=True).decode().strip()
    return int(cores)

def get_job_list(jobs):
    '''
    Get a list of job directories from the given list of job patterns.

    Args:
        jobs (list): A list of job patterns.

    Returns:
        list: A list of job directories.
    '''
    job_list = []
    for job in jobs:
        for ijob in glob.glob(job):
            if os.path.isdir(ijob):                        
                job_list.append(ijob)
    return job_list
    
def gen_supermetrics(file_name, sm_name=None, image=True):
    # file_type should be image or html
    if sm_name == None:
        basename = os.path.basename(file_name.rstrip("/"))
        if "." in basename:
            sm_name = os.path.splitext(basename)[0]
        else:
            sm_name = basename
    return {
        sm_name: {"type": {True:"image", False: "html"}.get(bool(image)),
                  "file": file_name}
    } 

def bak_file(filename, bak_org=False):
    '''
    Check if the file exists, and if so:
    if bak_org is True, backup the original file to a non-exist filename.bak{n}, n = 1,2,3,...
    else, return the filename.bak{n} that does not exist.
    '''
    if not os.path.exists(filename):
        return filename
    
    i = 1
    bakf = f"{filename}.bak{i}"
    while os.path.exists(bakf):
        i += 1
        bakf = f"{filename}.bak{i}"
        
    if bak_org:
        os.system(f"mv {filename} {bakf}")
        return filename 
    else:
        return bakf

def clean_files(ipath, f_list=[], folder_list=[]):
    '''
    Clean the files in the given folder.
    '''
    for ifile in f_list:
        for jfile in glob.glob(os.path.join(ipath,ifile)):
            if os.path.isfile(jfile):
                os.remove(jfile)
    for ifolder in folder_list:
        for jfolder in glob.glob(os.path.join(ipath,ifolder)):
            if os.path.isdir(jfolder):
                shutil.rmtree(jfolder)       

def get_abacus_inputfiles(ipath):
    """
    Find the input files in the given path and return a dict:
    {
        "input": input_param, # a dict of the input parameters, kpt_file, stru_file, pseudo_dir, orb_dir will be popped out or modified
        "stru": stru, # an AbacusStru object or None. The pp/orb now is only the base name.
        "kptf": kptf, # the absolute path of the kpt file or None if not found
        "pp": pp_abs_path, # a list of the absolute path of the pseudopotential files or [] if defined in STRU
        "orb": orb_abs_path, # a list of the absolute path of the orbital files or [] if defined in STRU
        "extra_files": extra_files # a list of the absolute path of the extra files
    }
    
    1. the kptf/pp/orb files are the absolute path, but not check if the files exist
    2. if the INPUT or STRU is not read, the corresponding value will be None
    3. if INPUT defines the stru_file, pseudo_dir or orb_dir, the returned dict will delete the corresponding key in the input_param
    4. extra_files are the files in the folder but not in the list of INPUT, STRU, KPT, PP, ORB
    """
    if not os.path.isdir(ipath):
        print(f"ERROR: {ipath} is not a directory")
        return {
            "input": None,
            "stru": None,
            "kpt": None,
            "pp": [],
            "orb": [],
            "extra_files": []
        }
    
    pwd = os.getcwd()
    os.chdir(ipath)
    
    if os.path.isfile("INPUT"):
        input_param = ReadInput("INPUT")
        struf = input_param.pop("stru_file","STRU")
        kptf = "KPT"
        if "kpt_file" in input_param:
            kptf = input_param["kpt_file"]
            input_param["kpt_file"] = os.path.basename(kptf)
    else:
        print(f"ERROR: INPUT file not found in {ipath}")
        input_param = {}
        kptf = "KPT"
        struf = "STRU"
        
    if not os.path.isfile(struf):
        print(f"ERROR: STRU file '{struf}' not found in {ipath}")
        stru = None
    else:
        stru = AbacusStru.ReadStru(struf)
        if not stru:
            print(f"ERROR: read STRU failed in {ipath}")
            stru = None
            
    pp_path = input_param.pop("pseudo_dir","")
    orb_path = input_param.pop("orb_dir","")
    
    pp_name = []
    orb_name = []
    pp_abs_path = None
    orb_abs_path = None
    if stru:
        pp = stru.get_pp()
        orb = stru.get_orb()
        pp = pp if pp else []
        orb = orb if orb else []
        pp_abs_path = [os.path.abspath(os.path.join(pp_path,i)) for i in pp]
        orb_abs_path = [os.path.abspath(os.path.join(orb_path,i)) for i in orb]
        pp_name = [os.path.basename(i) for i in pp]
        orb_name = [os.path.basename(i) for i in orb]
        if len(pp_name) > 0:
            stru.set_pp(pp_name)
        if len(orb_name) > 0:
            stru.set_orb(orb_name)    
    
    extra_files = []
    for ifile in os.listdir("."):
        if ifile not in ["INPUT",struf,kptf]+pp_name+orb_name and os.path.isfile(os.path.join(ifile)):
            extra_files.append(os.path.abspath(ifile))
    
    kptf = os.path.abspath(kptf) if os.path.isfile(kptf) else None
    os.chdir(pwd)        
    return {
        "input": input_param if input_param else None,
        "stru": stru,
        "kpt": kptf,
        "pp": pp_abs_path,
        "orb": orb_abs_path,
        "extra_files": extra_files
    }

def clean_none_list(*args):
    '''
    For given same length lists, remove the elements that are None in the same position.
    If the length of the lists are different, return the original lists.
    '''
    if len(args) == 0:
        return []
    n = len(args[0])
    for i in args:
        if len(i) != n:
            return args
    new_args = [[] for i in range(len(args))]
    for i in range(n):
        if all([not isinstance(j[i],type(None)) for j in args]):
            for j in range(len(args)):
                new_args[j].append(args[j][i])
    return tuple(new_args)

def run_command(
        cmd,
        shell=True
):
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=shell,
        executable='/bin/bash'
    )
    out = ""
    err = ""
    while True:
        readable, _, _ = select.select(
            [process.stdout, process.stderr], [], [])

        # 读取已经准备好的输出
        for fd in readable:
            if fd == process.stdout:
                line = process.stdout.readline()
                print(line.decode()[:-1])
                out += line.decode()
            elif fd == process.stderr:
                line = process.stderr.readline()
                print("STDERR:", line.decode()[:-1])
                err += line.decode()

        # 如果子进程已经结束，则退出循环
        return_code = process.poll()
        if return_code is not None:
            break
    return return_code, out, err


def cal_cellparam(cell):
    """
    Calculate the parameters of a cell.
    
    Parameters:
        cell (np.ndarray): The cell matrix.
        
    Returns:
        tuple: A tuple containing the lengths and angles of the cell.
    """
    import numpy as np
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    alpha = np.arccos(np.dot(cell[1], cell[2]) / (b * c)) * 180 / np.pi
    beta = np.arccos(np.dot(cell[0], cell[2]) / (a * c)) * 180 / np.pi
    gamma = np.arccos(np.dot(cell[0], cell[1]) / (a * b)) * 180 / np.pi
    
    return a, b, c, alpha, beta, gamma

def inter_coord(cell1, coord1, 
                cell2, coord2,
                n):
    """
    Interpolate coordinates from one cell to another.
    
    Parameters:
        cell1 (np.ndarray): The first cell matrix.
        coord1 (np.ndarray): Coordinates in the first cell.
        cell2 (np.ndarray): The second cell matrix.
        coord2 (np.ndarray): Coordinates in the second cell.
        n (int): Number of interpolation points.
        
    Returns:
        np.ndarray: Interpolated coordinates in the second cell.
    """
    import numpy as np
    cell1 = np.array(cell1)
    cell2 = np.array(cell2)
    coord1 = np.array(coord1)
    coord2 = np.array(coord2)
    coord1_direct = np.dot(np.linalg.inv(cell1), coord1.T).T
    coord2_direct = np.dot(np.linalg.inv(cell2), coord2.T).T

    coord_diff = coord2_direct - coord1_direct
    # if diff larger than 0.5, then add 1 to the diff
    while np.any(coord_diff > 0.5) or np.any(coord_diff < -0.5):
        coord_diff[coord_diff > 0.5] -= 1
        coord_diff[coord_diff < -0.5] += 1
    coord2_direct = coord1_direct + coord_diff

    maxdiff_idx = np.argmax(np.abs(coord_diff))
    print("Max difference index:", maxdiff_idx)
    print("Max difference value:", np.max(np.abs(coord_diff)))
    
    a1, b1, c1, alpha1, beta1, gamma1 = cal_cellparam(cell1)
    a2, b2, c2, alpha2, beta2, gamma2 = cal_cellparam(cell2)
    print(f"Cell1: a={a1:.2f}, b={b1:.2f}, c={c1:.2f}, alpha={alpha1:.2f}, beta={beta1:.2f}, gamma={gamma1:.2f}")
    print(f"Cell2: a={a2:.2f}, b={b2:.2f}, c={c2:.2f}, alpha={alpha2:.2f}, beta={beta2:.2f}, gamma={gamma2:.2f}")
    
    cells = []
    coords = []
    for i in range(n+1):
        alpha = alpha1 + (alpha2 - alpha1) * i / n
        beta = beta1 + (beta2 - beta1) * i / n
        gamma = gamma1 + (gamma2 - gamma1) * i / n
        
        a = a1 + (a2 - a1) * i / (n - 1)
        b = b1 + (b2 - b1) * i / (n - 1)
        c = c1 + (c2 - c1) * i / (n - 1)
        
        cell = np.array([[a, 0, 0], 
                         [b * np.cos(np.radians(gamma)), b * np.sin(np.radians(gamma)), 0], 
                         [c * np.cos(np.radians(beta)), 
                          c * (np.cos(np.radians(alpha)) - np.cos(np.radians(beta)) * np.cos(np.radians(gamma))) / np.sin(np.radians(gamma)), 
                          c * np.sqrt(1 - np.cos(np.radians(alpha))**2 - np.cos(np.radians(beta))**2 + 2 * np.cos(np.radians(alpha)) * np.cos(np.radians(beta)) * np.cos(np.radians(gamma)))]])

        coord_direct = coord1_direct + coord_diff * i / n
        coord = np.dot(cell, coord_direct.T).T
        
        cells.append(cell)
        coords.append(coord)
    return cells, coords

def read_gaussian_cube(fcube: str):
    """read the Gaussian cube format volumetric data.
    The detailed format information can be found here:
    https://paulbourke.net/dataformats/cube/

    In brief, 

    ```
    [comment line1]
    [comment line2]
    [natom] [x-origin] [y-origin] [z-origin]
    [nx] [e11 in a.u.] [e12 in a.u.] [e13 in a.u.]
    [ny] [e21 in a.u.] [e22 in a.u.] [e23 in a.u.]
    [nz] [e31 in a.u.] [e32 in a.u.] [e33 in a.u.]
    [Z1] [chg1] [x1 in a.u.] [y1 in a.u.] [z1 in a.u.]
    [Z2] [chg2] [x2 in a.u.] [y2 in a.u.] [z2 in a.u.]
    ...
    [Znatom] [xnatom in a.u.] [ynatom in a.u.] [znatom in a.u.]
    [data1] [data2] [data3] ... [data6]
    [data7] [data8] [data9] ... [data12]
    ...
    ```
    
    Args:
        fcube (str): the file name of the cube file.
    
    Returns:
        data (dict): a dictionary containing the following
    """
    import numpy as np
    with open(fcube, 'r') as f:
        lines = f.readlines()
    
    out = {}
    out["comment 1"] = lines[0].strip()
    out["comment 2"] = lines[1].strip()

    natom, x_origin, y_origin, z_origin = map(float, lines[2].split())
    out["natom"] = natom
    out["origin"] = (x_origin, y_origin, z_origin)

    out["nx"] = int(lines[3].split()[0])
    out["ny"] = int(lines[4].split()[0])
    out["nz"] = int(lines[5].split()[0])
    e11, e12, e13 = map(float, lines[3].split()[1:4])
    e21, e22, e23 = map(float, lines[4].split()[1:4])
    e31, e32, e33 = map(float, lines[5].split()[1:4])
    out["R"] = np.array([[e11, e12, e13], [e21, e22, e23], [e31, e32, e33]])

    atomz, chg, coords = [], [], []
    for line in lines[6:6+int(natom)]:
        atomz.append(int(line.split()[0]))
        chg.append(float(line.split()[1]))
        coords.append(list(map(float, line.split()[2:])))
    out["atomz"] = atomz
    out["chg"] = chg
    out["coords"] = coords

    data = []
    for line in lines[6+int(natom):]:
        data.extend(list(map(float, line.split())))
    #out["data"] = np.array(data).reshape(int(nx), int(ny), int(nz)) # is it necessary to reshape?
    out["data"] = np.array(data) # not reshaped
    return out

def write_gaussian_cube(data: dict, fcube: str, ndigits: int = 6):
    """write the Gaussian cube format volumetric data.
    Data should be organized as what is returned by read_cube().
    
    (if present) the ndigits will be used to control the number of digits of volumetric data,
    while the coordination data in header of cube is fixed to 6 digits.
    
    format of each part is defined in the following with format string:
    ```
    %s
     %s
     %4d %11.6f %11.6f %11.6f
     %4d %11.6f %11.6f %11.6f
     %4d %11.6f %11.6f %11.6f
     %4d %11.6f %11.6f %11.6f
     %4d %11.6f %11.6f %11.6f %11.6f
    ...
     %[width].[ndigits]f %[width].[ndigits]f %[width].[ndigits]f %[width].[ndigits]f %[width].[ndigits]f %[width].[ndigits]f
    ...
    ```

    Args:
        data (dict): the dictionary containing the cube data.
        fcube (str): the file name of the cube file.
        ndigits (int): the number of digits to be used for the volumetric data.
                        Default is 6, which is the same as the original cube format.
    
    Returns:
        None
    """
    
    width = ndigits + 7
    # temporarily there is no more format controlling options.

    with open(fcube, 'w') as f:
        f.write(data["comment 1"] + "\n")
        f.write(data["comment 2"] + "\n")
        f.write(" %4d %11.6f %11.6f %11.6f\n" % (data["natom"], data["origin"][0], data["origin"][1], data["origin"][2]))
        f.write(" %4d %11.6f %11.6f %11.6f\n" % (data["nx"], data["R"][0][0], data["R"][0][1], data["R"][0][2]))
        f.write(" %4d %11.6f %11.6f %11.6f\n" % (data["ny"], data["R"][1][0], data["R"][1][1], data["R"][1][2]))
        f.write(" %4d %11.6f %11.6f %11.6f\n" % (data["nz"], data["R"][2][0], data["R"][2][1], data["R"][2][2]))
        for i in range(int(data["natom"])):
            f.write(" %4d %11.6f %11.6f %11.6f %11.6f\n" % (data["atomz"][i], data["chg"][i], data["coords"][i][0], data["coords"][i][1], data["coords"][i][2]))
        for i in range(0, len(data["data"]), 6):
            f.write(" ".join([f"%{width}.{ndigits}e" % x for x in data["data"][i:i+6]]))
            f.write("\n")
    return

def check_abacus_inputs(job: str) -> Tuple[bool, str]:
    """
    Check if the ABACUS input files in the given job directory are valid.
    
    Parameters:
        job (str): The job directory.
        
    Returns:
        Tuple[bool, str]: A tuple containing a boolean indicating if the inputs are valid,
                          and a string message.
    """
    if not os.path.isdir(job):
        return False, f"{job} is not a directory"
    
    # check INPUT file
    if not os.path.isfile(os.path.join(job, "INPUT")):
        return False, f"INPUT file not found in {job}"
    
    input_param = ReadInput(os.path.join(job, "INPUT"))
    struf = input_param.get("stru_file", "STRU")
    # check STRU file
    if not os.path.isfile(os.path.join(job, struf)):
        return False, f"STRU file '{struf}' not found in {job}"
    
    stru = AbacusStru.ReadStru(os.path.join(job, struf))
    if not stru:
        return False, f"Read STRU failed in {job}"
    
    # check K Point file
    kspacing = input_param.get("kspacing", None)
    basis = input_param.get("basis_type", "pw")
    gamma_only = input_param.get("gamma_only", False)
    kptf = input_param.get("kpoint_file", "KPT")
    if kspacing is None and (basis == "pw" or (basis == "lcao" and not gamma_only)) and not os.path.isfile(os.path.join(job, kptf)):
        return False, f"KPT file '{kptf}' not found in {job} and kspacing is not set."
    
    # check pp/orbital files
    pp_path = input_param.get("pseudo_dir", "")
    orb_path = input_param.get("orbital_dir", "")
    pp = stru.get_pp()
    orb = stru.get_orb()
    if pp is None:
        return False, f"Pseudopotential files not defined in STRU file '{struf}' in {job}"
    
    # check pp_file
    pp = [os.path.join(job, pp_path, i) for i in pp]
    for ipp in pp:
        if not os.path.isfile(ipp):
            return False, f"Pseudopotential file '{ipp}' not found in {job}"
    
    # check orbital files
    if basis == "lcao" or input_param.get("onsite_radius", None):
        if orb is None:
            return False, f"Orbital files not defined in STRU file '{struf}' in {job}"
        orb = [os.path.join(job, orb_path, i) for i in orb]
        for iorb in orb:
            if not os.path.isfile(iorb):
                return False, f"Orbital file '{iorb}' not found in {job}"
    
    # check orther INPUT parameters
    ks_solver = input_param.get("ks_solver", None)
    device = input_param.get("device", "cpu")
    if basis == "pw" and ks_solver is not None and ks_solver not in ["cg", "dav", "dav_subspace", "bpcg"]:
        return False, f"Invalid ks_solver '{ks_solver}' in INPUT file in {job}. It should be 'cg', 'dav', 'bpcg', or 'dav_subspace' for PW basis."
    if basis == "lcao" and ks_solver is not None:
        if device == "cpu" and ks_solver not in ["genelpa", "scalapack_gvx","lapack", "elpa"]:
            return False, f"Invalid ks_solver '{ks_solver}' in INPUT file in {job}. It should be 'genelpa', 'scalapack_gvx', 'elpa', or 'lapack' for LCAO basis with CPU device."
        if device == "gpu" and ks_solver not in ["cusolver", "cusolvermp", "elpa"]:
            return False, f"Invalid ks_solver '{ks_solver}' in INPUT file in {job}. It should be 'cusolver', 'cusolvermp', or 'elpa' for LCAO basis with GPU device."
    
    # check nspin and soc
    nspin = input_param.get("nspin", None)
    soc = IsTrue(input_param.get("lspinorb", None))
    noncollinear = IsTrue(input_param.get("noncolin", None))
    if nspin is not None:
        if nspin not in [1, 2, 4, "1", "2", "4"]:
            return False, f"Invalid 'nspin' value: {nspin}. It should be 1 (no spin), 2 (spin polarized), or 4 (non-collinear spin)."
        nspin = int(nspin)
        if nspin in [1, 2] and (soc or noncollinear):
            return False, "Spin-orbit coupling ('lspinorb') or non-collinear spin ('noncolin') is set, but nspin is not 4. Please set nspin to 4 for non-collinear spin calculations."
        if nspin == 4 and noncollinear == False:
            return False, "Non-collinear spin ('noncolin') is set to False, but nspin is 4. Please set noncollinear to True for non-collinear spin calculations."
    elif soc and noncollinear == False:
        return False, "Spin-orbit coupling ('lspinorb') is set, but noncollinear ('noncolin') is set to False. Please set noncollinear to True for non-collinear spin calculations."
    
    return True, "All input files are valid"