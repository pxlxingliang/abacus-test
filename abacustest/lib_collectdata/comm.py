import os,sys,traceback
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element

HARTREE2EV = 27.211396132
EV2RY = 2.0 / HARTREE2EV
RY2EV = HARTREE2EV / 2.0
BOHR2A = 0.52917721092
KBAR2HARTREEPERBOHR3 = 3.398927420868445E-6
KBAR2EVPERANGSTROM3 = KBAR2HARTREEPERBOHR3 * HARTREE2EV / BOHR2A**3

def ReadFile(ifile,warn=True):
    if ifile == None:
        return []

    if os.path.isfile(ifile):
        with open(ifile) as f1:
            context = f1.readlines()
    else:
        if warn: print("WARNING: can not find file %s" % ifile)
        context = []

    return context

def FindOutput(path,keyinfo):
    'find the file that has keyinfo'
    cwd = os.getcwd()
    os.chdir(path)
    allfiles = os.listdir(".")
    for ifile in allfiles:
        if os.path.isfile(ifile):
            find = False
            try:
                with open(ifile) as f1: lines = f1.read()
                if keyinfo in lines:
                    find = True
            except:
                pass
                
            if find:
                os.chdir(cwd)
                return os.path.join(path,ifile)
    os.chdir(cwd)
    return None

def ReadXmlFile(ifile,warn=True):
    if ifile == None:
        return None

    if os.path.isfile(ifile):
        try:
            tree = ET.parse(ifile)
            return tree.getroot()
        except:
            traceback.print_exc()
            if warn: print("WARNING: can not parse file %s" % ifile)
            return None
    else:
        if warn: print("WARNING: can not find file %s" % ifile)
        return None

def XmlFindMultiLayer(root,layerlist):
    tmp = root
    for i in layerlist:
        if tmp == None: return None
        tmp = tmp.find(i)
    return tmp

def XmlFindMultiLayerText(root,layerlist):
    tmp = XmlFindMultiLayer(root,layerlist)
    if isinstance(tmp,Element):
        return tmp.text
    return None

def XmlGetText(tmp1,func=lambda x:x,idx: int=None):
    funcname = sys._getframe().f_code.co_name
    if tmp1 == None:
        return None
    elif isinstance(tmp1,list):
        if idx == None:
            return [func(i.text) for i in tmp1]
        else:
            try:
                return func(tmp1[idx].text)
            except:
                print("ERROR: %s, try to get the text of %d-th element, the length of list is %d" 
                        % (funcname,idx,len(tmp1)))
                traceback.print_exc()
                return None
    elif isinstance(tmp1,Element):
        return func(tmp1.text)
    else:
        print("Not support type '%s' in function '%s' now." % 
              (type(tmp1),funcname))
        return None

def istr(x,n=None):
    if type(x) == float and n != None:
            return "%%.%df" % n % x
    return str(x)

def iint(x):
    try:
        return int(x)
    except:
        return None

def ifloat(x):
    try:
        return float(x)
    except:
        return None

def ibool(x):
    try:
        if x.lower() == "true":
            return True
        elif x.lower() == "false":
            return False
        else:
            return bool(int(x))
    except:
        return None

def imath(a,b,symbol):
    if a == None or b == None:
        return None
    elif symbol == '*':
        return a*b
    elif symbol == '+':
        return a+b
    elif symbol == '-':
        return a-b
    elif symbol == '/':
        return a/b
    else:
        print("Not support '%s' now" % symbol)
        return None

def strtime2sec(stim_str):
    # time_str: 1d2h3m4s
    
    day = 0
    hour = 0
    minute = 0
    second = 0
    time_str = stim_str.lower()
    if 'd' in time_str:
        day = int(time_str.split('d')[0])
        time_str = time_str.split('d')[1]
    if 'h' in time_str:
        hour = int(time_str.split('h')[0])
        time_str = time_str.split('h')[1]
    if 'm' in time_str:
        minute = int(time_str.split('m')[0])
        time_str = time_str.split('m')[1]
    if 's' in time_str:
        second = float(time_str.split('s')[0])
    return day*24*3600 + hour*3600 + minute*60 + second

def cal_band_gap(band,efermi):
    # band: [[spin1],[spin2],...]
    # efermi: [efermi_spin1,efermi_spin2]
    if isinstance(efermi,float):
        efermi = [efermi,efermi]
            
    cb = None
    vb = None
    
    for ispin in range(len(band)):
        nband_below_fermi = None
        for ik in range(len(band[ispin])):
            for ib in range(len(band[ispin][ik])):
                if band[ispin][ik][ib] > efermi[ispin]:
                    if nband_below_fermi == None:
                        nband_below_fermi = ib
                    elif nband_below_fermi != ib:
                        return 0
                    icb = band[ispin][ik][ib-1] - efermi[ispin]
                    ivb = band[ispin][ik][ib] - efermi[ispin]
                    if cb == None or icb > cb:
                        cb = icb
                    if vb == None or ivb < vb:
                        vb = ivb
                    break
    
    if cb == None or vb == None:
        band_gap = None
    else:
        band_gap = vb - cb
        if band_gap < 0: band_gap = 0 
    return band_gap

def plot_band(band,fname,efermi=None):
    """
    band: NSPIN X NKPTS X NBAND
    fname: the file name of the band plot
    efermi: float or [float, float]
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    def plot_one(ax, iband, efermi=None, _range=None):
        for j in range(len(iband)):
            ax.plot(iband[j], "--o", linewidth=0.3, markersize=0.3)
        if efermi is not None:
            ax.axhline(y=efermi, color='r', linestyle='--', linewidth=0.3)
            ax.text(int(len(iband[0])/3), efermi, f'Efermi={efermi:.3f} eV', color='r', fontsize=8)
        if _range is not None:
            ax.set_ylim(_range)
        else:
            ax.set_ylim(np.min(iband.flatten()) - 1, np.max(iband.flatten()) + 1)
        ax.set_ylabel('Energy (eV)')
        ax.set_xlabel('k point')
        ax.set_xlim(-1, len(iband[0]))
        
    band = np.array(band)
    nband = len(band)
    if isinstance(efermi, (float, type(None),int)):
        efermi = [efermi] * nband

    fig, axs = plt.subplots(nband, 2, figsize=(8, 4*nband))
    axs = axs.flatten()
    for i in range(nband):
        plot_one(axs[2*i], np.array(band[i]).T, efermi=efermi[i])
        _range = [-5, 5] if efermi is None else [efermi[i]-5, efermi[i]+5]
        plot_one(axs[2*i+1], np.array(band[i]).T, efermi=efermi[i], _range=_range)
        axs[2*i].set_title(f'Spin {i+1} (Full)')
        axs[2*i+1].set_title(f'Spin {i+1} (Zoomed)')
        
    fig.tight_layout()
    plt.savefig(fname, dpi=300)
    plt.close()

def abacus_orb_label(l,m):
    if l == 0:
        return "s"
    elif l == 1:
        if m == 0:
            return "p_z"
        elif m == 1:
            return "p_x"
        elif m == 2:
            return "p_y"
    elif l == 2:
        if m == 0:
            return "d_z^2"
        elif m == 1:
            return "d_xz"
        elif m == 2:
            return "d_yz"
        elif m == 3:
            return "d_x^2-y^2"
        elif m == 4:
            return "d_xy"

    return f"{l}_{m}"

def get_abacus_json(json_dict,key_list):
    # need first find if key_list in json_dict
    v = None
    has_keys = True
    for k in key_list:
        if k in json_dict:
            json_dict = json_dict[k]
        else:
            has_keys = False
            break
    if has_keys:
        v = json_dict
    return v,has_keys

def get_metric_from_str(istr):
    # istr: "{energy}/{natom}"
    # the metric name may in the format of "{metric_name}"
    # return a list of metric name
    
    if "{" not in istr:
        return [istr]

    metric_list = []
    for i in istr.split("{"):
        if "}" in i:
            metric_list.append(i.split("}")[0])
    return metric_list

def get_mulliken(mullikenf):
    """Read mulliken population analysis file
    
    Return:
        atom_mag: list of list, atom_mag[step][atom] = mag or [mx,my,mz]
        atom_elec: list of list, atom_elec[step][atom] = [elec_spin1,(elec_spin2)]
    """
    atom_mag = []
    atom_elec = []
    atom_labels = []
    with open(mullikenf) as f1: lines = f1.readlines()
    for idx,line in enumerate(lines):
        if line[:5] == "STEP:":
            atom_mag.append([])
            atom_elec.append([])
            atom_labels.append([])
        elif "Total Magnetism on atom" in line:
            if "(" in line and ")" in line:
                # for nspin = 4
                sline = line.split("(")[1].split(")")[0].split(",")
                if len(sline) == 3:
                    atom_mag[-1].append([float(sline[0]),float(sline[1]),float(sline[2])])
            else:
                atom_mag[-1].append(float(line.split()[-1]))
        elif "Zeta of" in line:
            atom_labels[-1].append(line.split()[3])
            two_spin = True if "Spin 2" in line else False
            atom_elec[-1].append([])
            j = idx + 1
            while j < len(lines) :
                if "Zeta of" in lines[j]:
                    break
                if "sum over m+zeta" in lines[j] or "SUM OVER M+Zeta" in lines[j] and "M+Zeta+L" not in lines[j]:
                    atom_elec[-1][-1].append([float(lines[j].split()[3])])
                    if two_spin:
                        atom_elec[-1][-1][-1].append(float(lines[j].split()[4]))
                j += 1
    return atom_mag,atom_elec, atom_labels