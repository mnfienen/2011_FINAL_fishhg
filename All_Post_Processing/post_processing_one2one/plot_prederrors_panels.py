import numpy as np
from matplotlib import pyplot as plt
import pickle
import os
import matplotlib as mpl
import shutil
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages

mpl.rcParams['pdf.fonttype'] = 42
def killer(curr_plots):
    if os.path.exists(os.path.join(os.getcwd(),curr_plots)) == True:
        shutil.rmtree(os.path.join(os.getcwd(),curr_plots))

def resid_plotter(cax,x,y,DL,logFlag,figtitle,figname,subf):
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
    plt.hold(True)
    mksize = 4.5
    # plot one to one line
    plt.plot([minxy,maxxy],[minxy,maxxy],'black')  

    # plot up the detects solid
    plotx = x[Detectinds]
    ploty = y[Detectinds]
    plt.plot(plotx,ploty,'bo',markerfacecolor = 'blue', markersize = mksize,markeredgecolor = 'black')
    # plot up the nondetects white
    plotx = x[DLinds]
    ploty = y[DLinds]
    plt.plot(plotx,ploty,'bo',markerfacecolor = 'white',markersize = mksize, markeredgecolor = 'black')
    plt.xlim([minxy, maxxy])
    plt.ylim([minxy, maxxy])

    if logFlag:
        plt.yscale('log')
        plt.xscale('log')
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    ticknames = ['0.001','0.01','0.1','1.0','10.0']
    plt.setp(ax1,xticklabels=ticknames)
    plt.setp(ax1,yticklabels=ticknames)

    plt.text(0.0015,5,'{0}'.format(figtitle))
    plt.savefig('{0}/{1}.pdf'.format(subf,figname))
    plt.close('all')
    if logFlag:
        slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(plotx),np.log10(ploty))
    else:
        slope, intercept, r_value, p_value, std_err = stats.linregress(plotx,ploty)
    return r_value**2

###################
# Special case for the ALL case
###################
def resid_plotter_all(x,y,DL,logFlag,figtitle,figname,subf):
    fig = plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111,aspect='equal')
    print '...plotting all-style...'
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
    minx = 0.05
    miny = minx
    maxx = 10**0.2
    maxy = maxx
    
    # limits
    minxy = np.min([minx,miny])
    maxxy = np.max([maxx,maxy])
    
    DLinds = np.nonzero(DL==1)
    Detectinds = np.nonzero(DL==0)
    plt.hold(True)
    mksize = 0.05
    # plot one to one line
    plt.plot([minxy,maxxy],[minxy,maxxy],'black')  

    # plot up the detects solid
    plotx = x[Detectinds]
    ploty = y[Detectinds]
    plt.plot(plotx,ploty,'bo',markerfacecolor = 'blue', markersize = mksize,markeredgecolor = 'blue')
    # plot up the nondetects white
    plotx = x[DLinds]
    ploty = y[DLinds]
    plt.plot(plotx,ploty,'bo',markerfacecolor = '0.6',markersize = mksize, markeredgecolor = '0.6')
    plt.xlim([minxy, maxxy])
    plt.ylim([minxy, maxxy])

    if logFlag:
        plt.yscale('log')
        plt.xscale('log')
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

#    ticknames = ['0.001','0.01','0.1','1.0','10.0']
    ticknames = ['0.01','0.1','1.0']
    plt.setp(ax1,xticklabels=ticknames)
    plt.setp(ax1,yticklabels=ticknames)

    plt.text(0.0015,5,'{0}'.format(figtitle))
    plt.savefig('{0}/{1}.pdf'.format(subf,figname))
    plt.savefig('{0}/{1}.eps'.format(subf,figname))
    plt.close('all')
    if logFlag:
        slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(plotx),np.log10(ploty))
    else:
        slope, intercept, r_value, p_value, std_err = stats.linregress(plotx,ploty)
    return r_value**2

# make a class for post-processing
class LOO_post:
    def __init__(self, infile,**args):
        self.infile = infile
        #Check to see if a logFlag has been specified in the function's arguments
        try: self.logFlag = args['logFlag']
        except KeyError: self.logFlag = True   # if there is no 'logFlag' key, then use the default (True)
    
    # method to load up the data
    def load_data(self):

        # first, see if a pickle file already exists, if not, make it from the CSV file
        pkl_exists = False
        if os.path.exists(os.path.join(os.getcwd(),self.infile+'.pkl')):
            if os.path.getsize(os.path.join(os.getcwd(),self.infile+'.pkl')) != 0:
                pkl_exists = True
        
        
        if pkl_exists == False:
            ifp = open(infile + '.pkl','wb')
            self.indat = np.genfromtxt(self.infile, delimiter = ',', names = True,dtype = None)
            pickle.dump(self.indat,ifp)
        else:
            ifp = open(self.infile+'.pkl','rb')
            self.indat = pickle.load(ifp)
        ifp.close()
        
    def parse_input(self):
            
        # now parse and plot
        self.SpC   = self.indat['SpC']
        self.eunich_spc = np.unique(self.indat['SpC'])
        self.Hgmod = self.indat['PredictedHg_Loo']
        self.Hgobs = self.indat['Hg']
        self.DL    = self.indat['DL']
        self.CommonName = self.indat['CommonName']
        self.Cut   = self.indat['SampleType2']
        

    def plot_consolidated_misfit(self,logfile,subf):
        logfile.write('Wrote all consolidated samples to --> "all.pdf"\n')
        if os.path.exists(os.path.join(os.getcwd(),subf)) == False:
            os.mkdir(os.path.join(os.getcwd(),subf))
        R2 = resid_plotter_all(self.Hgmod,self.Hgobs,self.DL,self.logFlag,'All','all',subf)
        logfile.write('R^2 = {0}\n'.format(R2))
        
    def plot_each_series_misfit(self,logfile,**args):
        #Check to see if a cutoff has been specified in the function's arguments
        try: cutoff = args['cutoff']
        except KeyError: cutoff = 2   # if there is no 'cutoff' key, then use the default (2)
        try: cseries = args['series']
        except KeyError: os.error('Missing "series" arg\n')   # if there is no 'seris' key, then error out
        try: subf = args['subf']
        except KeyError: os.error('Missing "subf" arg\n')   # if there is no 'subf' key, then error out
        logfile.write('\n############\n#Running {0} with cutoff of {1}\n############\n'.format(cseries,cutoff))
        logfile.write('Figures in folder --> {0}/\n'.format(subf))
        if os.path.exists(os.path.join(os.getcwd(),subf)) == False:
            os.mkdir(os.path.join(os.getcwd(),subf))
        allcat = self.indat[cseries]
        eunich = np.unique(allcat)
        i = 0
        for cs in eunich:
            i+=1
            inds = np.nonzero(allcat==cs)[0]
            print 'running {0} of {1} --> {2}: {3} samples included'.format(i,len(eunich),cs,len(inds))
            
            if len(inds) > cutoff:
                if allcat[0] == 'SpC':
                    figtitle = self.CommonName[inds][0] + ": " + self.Cut[inds][0]
                else:
                    figtitle = str(cs)
                figname = figtitle.replace(' ','_')
                figname = figname.replace(":","")
                figname = figname.replace('/','_')
                R2 = resid_plotter(self.Hgmod[inds],self.Hgobs[inds],self.DL[inds],self.logFlag,figtitle,figname,subf)
                logfile.write('%s.pdf:%d of %d -->  %d samples included:  R^2 = %f\n' %(figname,i,len(eunich),len(inds),R2))
            else:
                logfile.write('   NOT PRINTED: %d of %d --> %s: N = %d samples --> N < cutoff\n' %(i,len(eunich),cs,len(inds)))
        
    def boxplot(self,**args):
        #Check to see if a toprank has been specified in the function's arguments
        try: self.toprank = args['toprank']
        except KeyError: self.toprank = 5   # if there is no 'toprank' key, then use the default (5)