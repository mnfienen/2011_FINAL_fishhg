import pylab as plt
import numpy as np

def resid_plotter(x,y,DL,logFlag,commonname,cut,cspc):
    fig = plt.figure()
    fig.add_subplot(111,aspect='equal')
    print '...plotting...'
    if logFlag:
        x = np.log10(x)
        y = np.log10(y)
        xlabel = 'log10 Modeled Hg'
        ylabel = 'log10 Observed Hg'
    else:
        xlabel = 'Modeled Hg'
        ylabel = 'Observed Hg'

    minx = np.min(x)
    miny = np.min(y)
    minxy = np.min([minx,miny])
    maxx = np.max(x)
    maxy = np.max(y)
    maxxy = np.max([maxx,maxy])
    DLinds = np.nonzero(DL==1)
    Detectinds = np.nonzero(DL==0)
    plt.hold(True)
    # plot up the nondetects white
    plotx = x[DLinds]
    ploty = y[DLinds]
    plt.plot(plotx,ploty,'bo',markerfacecolor = 'None')
    # plot up the detects solid
    plotx = x[Detectinds]
    ploty = y[Detectinds]
    plt.plot(plotx,ploty,'bo',markerfacecolor = 'b')
    
    plt.plot([minxy,maxxy],[minxy,maxxy])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.title('Misfit plot for {0}: {1}'.format(commonname, cut))
    plt.savefig('spcfigures/spc_{0}.pdf'.format(cspc))
    plt.close('all')