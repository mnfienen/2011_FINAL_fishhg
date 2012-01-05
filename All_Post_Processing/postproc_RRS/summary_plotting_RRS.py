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

# make a class for post-processing
class RRS_post:
    def __init__(self, infile,logfile):
        self.infile = infile
        self.logfile = logfile
    # method to load up the data
    def load_data(self):


        self.indat = np.genfromtxt(self.infile, delimiter = ',', names = True,dtype = None)
        
    def parse_input(self):
            
        # now parse and plot
        self.SpC   = self.indat['SpC']
        self.eunich_spc = np.unique(self.indat['SpC'])
        self.min_abs_res = self.indat['min_abs_res']
        self.mean_abs_res = self.indat['mean_res_abs']
        self.max_abs_res = self.indat['max_abs_res']
        self.min_res = self.indat['min_res']
        self.mean_res = self.indat['mean_res']
        self.max_res = self.indat['max_res']
        self.DataSource = self.indat['DataSource'] 
        self.CommonName = self.indat['CommonName']
        self.eunich_common = np.unique(self.CommonName)
        self.ID = self.indat['ID']
        self.N = self.indat['N']
        self.SD_abs = self.indat['SD_res_abs']
        self.Cut = self.indat['SampleType2']
       
        
    def barcharts(self):
        # plot all the residuals with error bars using 1SD +/-
        k = 0
        for cspc in self.eunich_spc:
            k+=1
            printst = 'Running Spc %d: is %d of %d' %(cspc,k,len(self.eunich_spc))
            print printst
            self.logfile.write(printst + '\n')
            currinds = np.nonzero(self.SpC==cspc)[0]
            cfig = plt.figure(figsize=(11,8.5))
            cax = cfig.add_subplot(111)
            x = np.arange(len(currinds)) + 1
            plt.hold = True
            #width = 1/3.0/len(currinds)
            #plt.bar(x,self.min_abs_res[currinds],width,color='r')
            plt.bar(x,self.mean_abs_res[currinds],color='g',yerr=self.SD_abs[currinds])
            #plt.bar(x+width*2,self.max_abs_res[currinds],width,color='b')
            plt.xticks(x,self.ID[currinds],rotation=90)
            plt.title(self.CommonName[currinds[0]] + ' ' + self.Cut[currinds[0]])
            plt.savefig('figures/bar_%d.pdf' %(cspc))
            plt.close(cfig)
            
    def boxplot_SpC(self,abs_flag = True):
        # plot all the residuals with error bars using 1SD +/-
        k = 0
        for cspc in self.eunich_spc:
            k+=1
            printst = 'Boxplots: Running Spc %d: is %d of %d' %(cspc,k,len(self.eunich_spc))
            print printst
            self.logfile.write(printst + '\n')
            currinds = np.nonzero(self.SpC==cspc)[0]
            if len(currinds) >= 10:
                cfig = plt.figure()
                cax = cfig.add_subplot(111)
                plt.hold = True
                if abs_flag:
                    data = [self.min_abs_res[currinds], self.mean_abs_res[currinds], self.max_abs_res[currinds], self.SD_abs[currinds]]
                else:                
                    data = [self.min_res[currinds], self.mean_res[currinds], self.max_res[currinds], self.SD_abs[currinds]]
                plt.boxplot(data)
                cax.grid(True)
                plt.title('n=%d, N=%d  ' %(len(currinds),np.sum(self.N[currinds])) + self.CommonName[currinds[0]] + ': ' + self.Cut[currinds[0]])            
                if abs_flag:
                    plt.xticks([1,2,3,4],('min(abs)','mean(abs)','max(abs)','SD'))                
                    plt.savefig('figures/abs_box_%d.pdf' %(cspc))
                else:
                    plt.xticks([1,2,3,4],('min','mean','max','SD'))
                    plt.savefig('figures/box_%d.pdf' %(cspc))
                plt.close(cfig)
            else:
                pass
    def boxplot_DataSource(self,abs_flag = True):
        # plot all the residuals with error bars using 1SD +/-
        k = 0
        self.eunich_DS = np.unique(self.DataSource)
        for cds in self.eunich_DS:
            k+=1
            printst = 'Boxplots: Running DataSource %s: is %d of %d' %(cds,k,len(self.eunich_DS))
            print printst
            self.logfile.write(printst + '\n')
            currinds = np.nonzero(self.DataSource==cds)[0]
            if len(currinds) >= 10:
                cfig = plt.figure()
                cax = cfig.add_subplot(111)
                plt.hold = True
                data = [self.min_res[currinds], self.mean_res[currinds], self.max_res[currinds], self.SD_abs[currinds]]
                plt.boxplot(data)
                cax.grid(True)
                plt.title('n=%d, N=%d  ' %(len(currinds),np.sum(self.N[currinds])) + self.DataSource[currinds[0]])            
            
                plt.xticks([1,2,3,4],('min','mean','max','SD'))
                plt.savefig('figures/box_%s.pdf' %(cds))
                plt.close(cfig)
            else:
                pass