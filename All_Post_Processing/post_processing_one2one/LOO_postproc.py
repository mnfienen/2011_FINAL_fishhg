import numpy as np
from matplotlib import pyplot as plt
import pickle
import os
import matplotlib as mpl
import shutil
from scipy import stats

mpl.rcParams['pdf.fonttype'] = 42

def killer(curr_plots):
    if os.path.exists(os.path.join(os.getcwd(),curr_plots)) == True:
        shutil.rmtree(os.path.join(os.getcwd(),curr_plots))

def resid_plotter(fig,cax,x,y,DL,logFlag,figtitle,figname,subf,crow,ccol,nrows):
    
    print '...plotting...'
#    if logFlag:
#        x = np.log10(x)
#        y = np.log10(y)
    xlabel = 'Modeled Hg'
    ylabel = 'Observed Hg'

    minx = np.min(x)
    miny = np.min(y)
    maxx = np.max(x)
    maxy = np.max(y)
    
    # override to specify axis limits
    minx = 0.001
    miny = minx
    maxx = 10.0
    maxy = maxx
    
    # limits
    minxy = np.min([minx,miny])
    maxxy = np.max([maxx,maxy])
    
    DLinds = np.nonzero(DL==1)
    Detectinds = np.nonzero(DL==0)
    cax.hold(True)
    mksize = 4.5
    # plot one to one line
    cax.plot([minxy,maxxy],[minxy,maxxy],'black')  

    # plot up the detects solid
    plotx = x[Detectinds]
    ploty = y[Detectinds]
    cax.plot(plotx,ploty,'bo',markerfacecolor = 'blue', markersize = mksize,markeredgecolor = 'black')
    # plot up the nondetects white
    plotx = x[DLinds]
    ploty = y[DLinds]
    cax.plot(plotx,ploty,'bo',markerfacecolor = 'white',markersize = mksize, markeredgecolor = 'black')
    cax.set_xlim([minxy, maxxy])
    cax.set_ylim([minxy, maxxy])

    if logFlag:
        cax.set_yscale('log')
        cax.set_xscale('log')
    

    ticknames = ['','0.01','0.1','1.0','10.0']
    if crow == nrows-1:
        cax.set_xticklabels(ticknames, size=4)
    else:
        cax.set_xticklabels([])
    if ccol == 0:
        cax.set_yticklabels(ticknames, size=4)
    else:
        cax.set_yticklabels([])

    cax.text(0.0015,2,'{0}'.format(figtitle),size=5)
    cax.xaxis.set_ticks_position('bottom')
    cax.yaxis.set_ticks_position('left')


# make a class for post-processing
class LOO_post:
    def __init__(self, infile,**args):
        self.infile = infile
        #Check to see if a logFlag has been specified in the function's arguments
        try: self.logFlag = args['logFlag']
        except KeyError: self.logFlag = True   # if there is no 'logFlag' key, then use the default (True)
    
    # method to load up the data
    def load_data(self):
        self.indat = np.genfromtxt(self.infile, delimiter = ',', names = True,dtype = None)

        
    def parse_input(self):           
        # now parse and plot
        self.SpC   = self.indat['SpC']
        self.eunich_spc = np.unique(self.indat['SpC'])
        self.Hgmod = self.indat['Hg_predicted']
        self.Hgobs = self.indat['Hg_observed']
        self.DL    = self.indat['DL']
        self.CommonName = self.indat['CommonName']
        self.Cut   = self.indat['SampleType2']

        
    def plot_each_series_misfit(self,logfile,**args):
        #Check to see if a cutoff has been specified in the function's arguments
        try: cutoff = args['cutoff']
        except KeyError: cutoff = 2   # if there is no 'cutoff' key, then use the default (2)
        try: cseries = args['series']
        except KeyError: os.error('Missing "series" arg\n')   # if there is no 'series' key, then error out
        try: subf = args['subf']
        except KeyError: os.error('Missing "subf" arg\n')   # if there is no 'subf' key, then error out
        logfile.write('\n############\n#Running {0} with cutoff of {1}\n############\n'.format(cseries,cutoff))
        logfile.write('Figures in folder --> {0}/\n'.format(subf))
        if os.path.exists(os.path.join(os.getcwd(),subf)) == False:
            os.mkdir(os.path.join(os.getcwd(),subf))
        allcat = self.indat[cseries]
        eunich = np.unique(allcat)
        i = 0
        bigfig = 0
        #set the parameters for the dimensions of subplot arrangement
        nrows = 8
        ncols = 6
        totsubplots = nrows*ncols
        w = 7 # width in inches
        h = nrows*w/ncols  #set length based on other dimensions      
        crow = 0
        ccol = -1
        hajime = True
        cplot= 0
        # make the figure panel
        cfig,ax = plt.subplots(nrows,ncols,figsize=(w,h)) 
        plt.subplots_adjust(wspace=0,hspace=0)
        for cs in eunich:
            i+=1
            inds = np.nonzero(allcat==cs)[0]
            print 'running {0} of {1} --> {2}: {3} samples included'.format(i,len(eunich),cs,len(inds))
                   
            if len(inds) > cutoff:
                ccol += 1
                cplot += 1
                currcol = np.mod(ccol,ncols)
                if ((currcol==0) & (hajime==False)):
                    crow+=1
                hajime = False
                figtitle = self.CommonName[inds][0] + ":\n" + self.Cut[inds][0]
                
                figname = figtitle.replace(' ','_')
                figname = figname.replace(":","")
                figname = figname.replace('/','_')
                if cplot == totsubplots:
                    resid_plotter(cfig,ax[crow,currcol],
                                  self.Hgmod[inds],
                                  self.Hgobs[inds],
                                  self.DL[inds],
                                  self.logFlag,
                                  figtitle,
                                  figname,
                                  subf,
                                  crow,
                                  currcol,
                                  nrows)                    
                    bigfig += 1
                    plt.savefig( subf + '/one2one_plots_%d.pdf' %(bigfig))
                    cfig,ax = plt.subplots(nrows,ncols,figsize=(w,h)) 
                    plt.subplots_adjust(wspace=0,hspace=0)
                    hajime=True
                    crow=0
                    ccol=-1
                    cplot = 0
                else:
                    resid_plotter(cfig,ax[crow,currcol],
                                  self.Hgmod[inds],
                                  self.Hgobs[inds],
                                  self.DL[inds],
                                  self.logFlag,
                                  figtitle,
                                  figname,
                                  subf,
                                  crow,
                                  currcol,
                                  nrows)
                logfile.write('%s.pdf:%d of %d -->  %d samples included\n' %(figname,i,len(eunich),len(inds)))
            else:
                logfile.write('   NOT PRINTED: %d of %d --> %s: N = %d samples --> N < cutoff\n' %(i,len(eunich),cs,len(inds)))