from . import comm
import matplotlib.pyplot as plt

def plot_line_point(x,ys, legends=None, title=None, xtitle=None, ytitle=None, xlabels=None,
                    fname=None,figsize=(8,6),fontsize=18,colors=None, markers=None,grid=True,legend_outbox=False):
    if colors == None:
        colors = ["b","g","r","c","m","y","k"]
    if markers == None:
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
        plt.xticks(x, xlabels, fontsize=fontsize-2)
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
    
