import os,json,glob,shutil,traceback
import subprocess,copy
from abacustest.lib_prepare.abacus import AbacusStru,ReadInput,ReadKpt
import select

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
