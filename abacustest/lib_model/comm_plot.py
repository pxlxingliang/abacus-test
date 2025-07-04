from . import comm
import matplotlib.pyplot as plt
from typing import List
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MaxNLocator

#try:
#    from matplotlib import rc
#    rc('font',**{'family':'sans-serif'})
#    rc('text', usetex=True)
#except:
#    print("No latex support in matplotlib.")

def check_font_exists(font_name):
    import matplotlib.font_manager as fm
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    return font_name in available_fonts


def set_font(font_name:str="Times New Roman", 
             font_size:int=None):
    """Set the font for matplotlib.

    Args:
        font_name (str, optional): font name. Defaults to "Times New Roman".
        font_size (int, optional): font size. Defaults to 18.
    """
    try:
        if not check_font_exists(font_name):
            print(f"Font '{font_name}' not found. Using default font.")
            return

        plt.rcParams["font.family"] = font_name  
        
        plt.rcParams["mathtext.fontset"] = "custom"
        plt.rcParams["mathtext.rm"] = font_name
        if check_font_exists(f"{font_name}:italic"):
            plt.rcParams["mathtext.it"] = f"{font_name}:italic"
        if check_font_exists(f"{font_name}:bold"):
            plt.rcParams["mathtext.bf"] = f"{font_name}:bold"
    except Exception as e:
        pass

    if font_size is not None:
        plt.rcParams["font.size"] = font_size


def plot_line_point(x,
                    ys,
                    legends:List[str]=None, 
                    title: str=None, 
                    xtitle: str=None, 
                    ytitle: str=None,
                    xlabels:List[str]=None,
                    fname:str=None,
                    figsize:tuple=(8,6),
                    fontsize:int=18,
                    colors:List[str] =None,
                    markers:List[str] =None,
                    grid:bool=True,
                    legend_outbox:bool=False):
    """Plot x vs ys[i] with line and point.

    Args:
        x (List[str|float|int]): x values.
        ys (List[List[float|int]]): y values.
        legends (List[str], optional): legends. Defaults to None.
        title (str, optional): title. Defaults to None.
        xtitle (str, optional): x title. Defaults to None.
        ytitle (str, optional): y title. Defaults to None.
        xlabels (List[str], optional): x labels. Defaults to None.
        fname (str, optional): file name. Defaults to None.
        figsize (tuple, optional): figure size. Defaults to (8,6).
        fontsize (int, optional): font size. Defaults to 18.
        colors (List[str], optional): colors. Defaults to None.
        markers (List[str], optional): markers. Defaults to None.
        grid (bool, optional): grid. Defaults to True.
        legend_outbox (bool, optional): legend outbox. Defaults to False
    """
    if colors == None:
        colors = ["b","g","r","c","m","y","k"]
    if markers == None:
        markers = ["o","*","+","x","s","p","d"]

    # if x is a list of string, then convert it to a list of int, and set xlabels to x
    if isinstance(x[0],str):
        if xlabels is None:
            xlabels = x
        x = list(range(len(x)))
    
    plt.figure(figsize=figsize)
    for i in range(len(ys)):
        if len(x) != len(ys[i]):
            print("Lengths of x and ys[{}] are different.".format(i))
            continue    

        ix,iy = comm.clean_none_list(x,ys[i])
        if len(ix) == 0:
            print("All elements in x and ys[{}] are None.".format(i))
            continue
        
        if legends is not None:
            ilegend = legends[i]
        else:
            ilegend = None
        
        plt.plot(ix,iy, label=ilegend, color=colors[i%len(colors)], marker=markers[(i//len(colors))%len(markers)])

    if title is not None:
        plt.title(title, fontsize=fontsize)
    if xtitle is not None:
        plt.xlabel(xtitle, fontsize=fontsize)
    if ytitle is not None:
        plt.ylabel(ytitle, fontsize=fontsize)
    if legends is not None:
        if legend_outbox:
            plt.legend(fontsize=fontsize-2, loc="upper left", bbox_to_anchor=(1, 1))
        else:
            plt.legend(fontsize=fontsize-2)
    if xlabels is not None:
        plt.xticks(x, xlabels, fontsize=fontsize-2,rotation=30)
    if grid:
        plt.grid()
    plt.tick_params(axis='x', labelsize=fontsize-2)
    plt.tick_params(axis='y', labelsize=fontsize-2) 
    plt.tight_layout()
    if fname is not None:
        plt.savefig(fname)  
    else:
        plt.show()
    plt.close()
    
def plot_bar(x,
             ys,
             legends:List[str]=None, 
             title: str=None, 
             xtitle: str=None, 
             ytitle: str=None,
             xlabels:List[str]=None,
             fname:str=None,
             figsize:tuple=(8,6),
             fontsize:int=18,
             colors:List[str] =None,
             ref_idx:int=None,
             grid:bool=True,
             legend_outbox:bool=False):
    """Plot x vs ys[i] with bar.

    Args:
        x (List[str|float|int]): x values.
        ys (List[List[float|int]]): y values.
        legends (List[str], optional): legends. Defaults to None.
        title (str, optional): title. Defaults to None.
        xtitle (str, optional): x title. Defaults to None.
        ytitle (str, optional): y title. Defaults to None.
        xlabels (List[str], optional): x labels. Defaults to None.
        fname (str, optional): file name. Defaults to None.
        figsize (tuple, optional): figure size. Defaults to (8,6).
        fontsize (int, optional): font size. Defaults to 18.
        colors (List[str], optional): colors. Defaults to None.
        ref_idx (int, optional): reference index. Defaults to None. If set, then will plot the deviation of ys[i] from ys[ref_idx].
        grid (bool, optional): grid. Defaults to True.
        legend_outbox (bool, optional): legend outbox. Defaults to False
    """
    if colors == None:
        colors = ["b","g","r","c","m","y","k"]
    
    markers = ["o","*","+","x","s","p","d"]

    plt.figure(figsize=figsize)
    for i in range(len(ys)):
        if len(x) != len(ys[i]):
            print("Lengths of x and ys[{}] are different.".format(i))
            continue    

        ix,iy = comm.clean_none_list(x,ys[i])
        if len(ix) == 0:
            print("All elements in x and ys[{}] are None.".format(i))
            continue
        
        if legends is not None:
            ilegend = legends[i]
        else:
            ilegend = None
        
        plt.bar(ix,iy, label=ilegend, color=colors[i%len(colors)])

    if title is not None:
        plt.title(title, fontsize=fontsize)
    if xtitle is not None:
        plt.xlabel(xtitle, fontsize=fontsize)


def plot_boxplot(x, 
                 ys,
                 title: str=None, 
                 xtitle: str=None, 
                 ytitle: str=None,
                 fname:str=None,
                 figsize:tuple=(8,6),
                 fontsize:int=18,):
    """Plot boxplot.

    Args:
        x (List[str|float|int]): x values.
        ys (List[List[float|int]]): y values, the length of ys should be the same as x.
        title (str, optional): title. Defaults to None.
        xtitle (str, optional): x title. Defaults to None.
        ytitle (str, optional): y title. Defaults to None.
        fname (str, optional): file name. Defaults to None.
        figsize (tuple, optional): figure size. Defaults to (8,6).
        fontsize (int, optional): font size. Defaults to 18.
    """
    plt.figure(figsize=figsize)
    plt.boxplot(ys, labels=x)
    if title is not None:
        plt.title(title, fontsize=fontsize)
    if xtitle is not None:
        plt.xlabel(xtitle, fontsize=fontsize)
    if ytitle is not None:
        plt.ylabel(ytitle, fontsize=fontsize)
    
    if fname is not None:
        plt.savefig(fname)
    else:
        plt.show()
    plt.close()


def gen_force_fd_1(ax1, ax2, x, e, ana, fd,
                   x_title = "Position ($\mathrm{\AA}$)",
                   e_title = "Energy (eV)",
                   f_title = "Force (eV/$\mathrm{\AA}$)",
                   e_label = "Energy (eV)",
                   ana_label = "Analytic (eV/$\mathrm{\AA}$)",
                   fd_label = "Finite Difference (eV/$\mathrm{\AA}$)",
                   title = None,
                   font_size = None,
                   grid = True
                   ):
    """
    Generate the force vs distance plot for the first case.
    
    """
    d,e,f,fd = comm.clean_none_list(x, e, ana, fd)
    
    if len(d) == 0 or len(x) == 0:
        return None

    if font_size is None:
        font_size = 18
    fz_label = font_size
    fz_xy_title = font_size + 4
    fz_title = None if font_size is None else font_size + 4
    fz_tick = None if font_size is None else font_size
    fz_marker = None if font_size is None else font_size - 6

    ax1.plot(d,e,"b+-",label=e_label,markersize=fz_marker, linewidth=3)
    ax2.plot(d,f,"ro-",label=ana_label,markersize=fz_marker, linewidth=3)
    ax2.plot(d,fd,"m*--",label=fd_label,markersize=fz_marker, linewidth=3)
    
    ymax1 = max(e)
    ymin1 = min(e)
    ymax2 = max(max(f),max(fd))
    ymin2 = min(min(f),min(fd))
    ax1.set_xlabel(x_title,fontsize=fz_xy_title)
    ax1.set_ylabel(e_title,color="b",fontsize=fz_xy_title)
    ax2.set_ylabel(f_title,color="r",fontsize=fz_xy_title)
    ax1.spines['left'].set_color('blue')
    ax1.tick_params(axis='y', colors='blue')
    ax2.spines['right'].set_color('red')
    ax2.tick_params(axis='y', colors='red')
    ax1.set_ylim(ymin1-(ymax1-ymin1)*0.1,ymax1+(ymax1-ymin1)*0.2)
    ax2.set_ylim(ymin2-(ymax2-ymin2)*0.1,ymax2+(ymax2-ymin2)*0.2)
    ax1.legend(loc="upper left",fontsize=fz_label)
    ax2.legend(loc="upper right",fontsize=fz_label)
    if title is not None:
        ax1.set_title(title,fontsize=fz_title)
    if grid:
        ax2.grid(True)
    if fz_tick is not None:
        ax1.tick_params(axis='x', labelsize=fz_tick, direction='in')
        ax1.tick_params(axis='y', labelsize=fz_tick, direction='in')
        ax2.tick_params(axis='y', labelsize=fz_tick, direction='in')
        
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    ax1.yaxis.get_offset_text().set_visible(False)
    ax1.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    ax2.yaxis.get_offset_text().set_visible(False)
    ax2.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    

def auto_set_yaxis(ax, ys, margin1=0.1, margin2=0.3, log_scale_threshold=1000):
    """Auto set y limit.

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): axes.
        ys (List[float]): y values.
        margin (float, optional): margin. Defaults to 0.1.
    """
    ymax = max(ys)
    ymin = min(ys)
    
    ax.yaxis.get_offset_text().set_visible(False)
    ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    # if the y values has a large range, then use scientific notation, and set log scale
    if ymax * ymin > 0 and ymax/ymin > log_scale_threshold:
        ax.set_yscale('log')
        ymax = ymax * 10
        ax.set_ylim(ymin,ymax)
        #ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        #ax.yaxis.set_minor_formatter(ScalarFormatter(useOffset=False))
        #ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    else:
        ymin = ymin - (ymax-ymin)*margin1
        ymax = ymax + (ymax-ymin)*margin2
        ax.set_ylim(ymin,ymax)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    return ymin, ymax
    
    