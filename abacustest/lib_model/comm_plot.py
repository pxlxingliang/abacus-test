from . import comm
import matplotlib.pyplot as plt
from typing import List

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
    