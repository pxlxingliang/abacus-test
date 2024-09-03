import numpy as np
from regex import F
from . import comm

def eos_fit(volume,energy):
    # refer to https://github.com/molmod/DeltaCodesDFT
    # fit a common third-order Birch-Murnaghan equation of state to the data
    # volume and energy are lists of volumes and energies
    # return volume, bulk_modulus, bulk_deriv, residuals
    if len(volume) != len(energy):
        print('Error: The number of volumes and energies do not match')
        return None
    v, e = comm.clean_none_list(volume,energy)
    if len(v) < 4:
        print('Error: Not enough data points to fit the equation of state')
        return None
    v,e = zip(*sorted(zip(v,e)))
    
    ve = np.array([v,e]).T
    fitdata = np.polyfit(ve[:, 0] ** (-2. / 3.),ve[:, 1], 3, full=True)
    residuals0 = (fitdata[1] / len(ve[:, 0]))**0.5
    deriv0 = np.poly1d(fitdata[0])
    deriv1 = np.polyder(deriv0, 1)
    deriv2 = np.polyder(deriv1, 1)
    deriv3 = np.polyder(deriv2, 1)
    
    volume0 = 0
    x = 0
    for x in np.roots(deriv1):
        if x > 0 and deriv2(x) > 0:
            volume0 = x**(-3./2.)
            break

    if volume0 == 0:
        print('Error: No minimum could be found')
    
    derivV2 = 4./9. * x**5. * deriv2(x)
    derivV3 = (-20./9. * x**(13./2.) * deriv2(x) -
        8./27. * x**(15./2.) * deriv3(x))
    bulk_modulus0 = derivV2 / x**(3./2.)
    bulk_deriv0 = -1 - x**(-3./2.) * derivV3 / derivV2
    b0 = bulk_modulus0 * 160.21766208 # GPa
    
    fit_x = np.linspace(min(min(ve[:,0]),volume0), max(max(ve[:,0]),volume0), 100)
    fit_y = deriv0(fit_x**(-2./3.))
    e0 = deriv0(volume0**(-2./3.))
    return float(volume0), float(e0), fit_x.tolist(), fit_y.tolist(), float(b0), float(bulk_deriv0), float(residuals0)

def cal_delta(v0, b0, bp, v0_ref, b0_ref, bp_ref):
    # v0, b0, bp: volume, bulk modulus, bulk modulus derivative of the material
    # unit of v0: A^3, unit of b0: GPa, unit of bp: unitless
    
    v01 = v0
    v02 = v0_ref
    b01 = b0 * 10.**9. / 1.602176565e-19 / 10.**30.
    b02 = b0_ref * 10.**9. / 1.602176565e-19 / 10.**30.
    b11 = bp
    b12 = bp_ref
    
    Vi = 0.94 * v01
    Vf = 1.06 * v02
    
    a31 = 9. * v01**3. * b01 / 16. * (b11 - 4.)
    a21 = 9. * v01**(7./3.) * b01 / 16. * (14. - 3. * b11)
    a11 = 9. * v01**(5./3.) * b01 / 16. * (3. * b11 - 16.)
    a01 = 9. * v01 * b01 / 16. * (6. - b11)
    
    a32 = 9. * v02**3. * b02 / 16. * (b12 - 4.)
    a22 = 9. * v02**(7./3.) * b02 / 16. * (14. - 3. * b12)
    a12 = 9. * v02**(5./3.) * b02 / 16. * (3. * b12 - 16.)
    a02 = 9. * v02 * b02 / 16. * (6. - b12)
    
    x = [0, 0, 0, 0, 0, 0, 0]
    
    x[0] = (a01 - a02)**2
    x[1] = 6. * (a11 - a12) * (a01 - a02)
    x[2] = -3. * (2. * (a21 - a22) * (a01 - a02) + (a11 - a12)**2.)
    x[3] = -2. * (a31 - a32) * (a01 - a02) - 2. * (a21 - a22) * (a11 - a12)
    x[4] = -3./5. * (2. * (a31 - a32) * (a11 - a12) + (a21 - a22)**2.)
    x[5] = -6./7. * (a31 - a32) * (a21 - a22)
    x[6] = -1./3. * (a31 - a32)**2.
    
    Fi = np.zeros_like(Vi)
    Ff = np.zeros_like(Vf)
    
    for n in range(7):
        Fi = Fi + x[n] * Vi**(-(2.*n-3.)/3.)
        Ff = Ff + x[n] * Vf**(-(2.*n-3.)/3.)
    
    # delta unit in eV
    delta = np.sqrt((Ff - Fi) / (0.12 * v02))
    return float(delta)

def plot_eos_one(ax,results, results_ref=None, label_size=16, legend_size=14,x_label="Volume ($A^3$)",y_label="Energy ($eV/atom$)"):
    """
    plot the EOS
    ax: the axis to plot
    results and results_ref is a dict:
    {
        "volume": [v1,v2,...],
        "energy_per_atom": [e1,e2,...],
    }
    
    will plot on ax with:
    1. volume vs energy with blue circle
    2. if fit the EOS, plot the EOS-fit with blue dashed line, and the v0 with blue dashed line
    3. if has results_ref: 
        3.1 plot the EOS-ref with red star
        3.2 if fit the EOS-ref, plot the EOS-fit-ref with red dashed line, and the v0_ref with red dashed line
    
    return a dict:
    {
        "fit": None or [v0,e0,fitx,fity,b0,bp,res],
        "fit_ref": None or [v0,e0,fitx,fity,b0,bp,res], # optional
        "delta": None or delta, # optional
    }
    """
    if "volume" not in results or "energy_per_atom" not in results:
        print("No volume or energy_per_atom in results.")
        return {}
    
    x = results["volume"]
    y = results["energy_per_atom"]
    if len(x) != len(y):
        print("The length of volume and energy_per_atom is not equal.")
        return {}
    
    fit1 = eos_fit(x,y)
    if fit1:
        v0, e0, fitx, fity, b0, bp, res = fit1
    else:
        print("Failed to fit the EOS.")
        v0 = fitx = fity = b0 = bp = res = None
        e0 = min(y)
    
    y = np.array(y) - e0
    
    ax.scatter(x,y,label="EOS",marker="o",color="blue")
    if not isinstance(fitx,type(None)) and not isinstance(fity,type(None)) and len(fitx) == len(fity):
        ilabel = "EOS-fit"
        if not isinstance(res,type(None)):
            ilabel += f" (RMSD={res:.2e})"
        ax.plot(fitx,np.array(fity)-e0,label=ilabel,color="blue",linestyle="--")
    
    if not isinstance(v0,type(None)):
        ax.axvline(v0,ls="--",color="blue")
        y_offset = max(y)*0.1
        ax.text(v0,y_offset,f"{v0:.2f}",color="blue",fontsize=legend_size)

    maxy = max(y)

    ax.set_xlabel(x_label,fontsize=label_size)
    ax.set_ylabel(y_label,fontsize=label_size)
    ax.tick_params(axis='x', labelsize=label_size)
    ax.tick_params(axis='y', labelsize=label_size)
    ax.legend(loc="upper right",fontsize=legend_size)
    ax.set_ylim(-maxy*0.1,maxy*1.2)
    ax.grid(True)
    
    delta = None
    ref_fit1 = None
    
    if results_ref:
        ref_x = results_ref.get("volume",None)
        ref_y = results_ref.get("energy_per_atom",None)
        if ref_x == None or ref_y == None or len(ref_x) != len(ref_y):
            print("The length of ref volume and ref energy_per_atom is not equal.")
            return {"fit":fit1,"fit_ref":None,"delta":None}
        else:
            ref_fit1 = eos_fit(ref_x,ref_y)
            if ref_fit1:
                ref_v0, ref_e0, ref_fitx, ref_fity, ref_b0, ref_bp, ref_res = ref_fit1
            else:
                print("Failed to fit the EOS-ref.")
                ref_v0 = ref_fitx = ref_fity = ref_b0 = ref_bp = ref_res = None
                ref_e0 = min(ref_y)

            ref_y = np.array(ref_y) - ref_e0
            maxy = max(maxy,max(ref_y))

            ax.scatter(ref_x,ref_y,label="EOS-ref",marker="*",color="red")
            if not isinstance(ref_fitx,type(None)) and not isinstance(ref_fity,type(None)) and len(ref_fitx) == len(ref_fity):
                ilabel = "EOS-fit-ref"
                if ref_res != None:
                    ilabel += f" (RMSD={ref_res:.2e})"
                ax.plot(ref_fitx,np.array(ref_fity)-ref_e0,label=ilabel,color="red",linestyle="--")

            if not isinstance(ref_v0,type(None)):
                ax.axvline(ref_v0,ls="--",color="red")
                y_offset = maxy*0.1
                ax.text(ref_v0,y_offset,f"{ref_v0:.2f}",color="red",fontsize=legend_size)

            if None not in [v0,b0,bp,ref_v0,ref_b0,ref_bp]:
                delta = cal_delta(v0,b0,bp,ref_v0,ref_b0,ref_bp)            
            
            ax.set_ylim(-maxy*0.1,maxy*1.2)
            ax.legend(loc="upper right",fontsize=legend_size)
            
    return {"fit":fit1,"fit_ref":ref_fit1,"delta":delta}
        
    