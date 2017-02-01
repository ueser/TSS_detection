
from matplotlib import pylab as pl
import numpy as np
import genome.db

import seaborn as sns
import itertools
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
sns.set_style('white')


def main():
    savePath='/Users/umut/Projects/intragenicTranscription/analysis/'
    gdb = genome.db.GenomeDB(path='/Users/umut/Projects/genome/data/share/genome_db',assembly='sacCer3')

    with open('../analysis/TSSseq_replicates.txt','r') as f:
        datNames = [t.split('\n')[0].split('\t') for t in f.readlines()]

    dataDict = {}
    for datNam1,datNam2 in datNames:
        dataDict[datNam1+'_pos'] = gdb.open_track(datNam1+'_pos')
        dataDict[datNam1+'_neg'] = gdb.open_track(datNam1+'_neg')
        dataDict[datNam2+'_pos'] = gdb.open_track(datNam2+'_pos')
        dataDict[datNam2+'_neg'] = gdb.open_track(datNam2+'_neg')

    chname='chrIV'
    pos=0
    width=1500000
    xran  = np.array([1,2,5,10,15,20,100])

    logcorcoef={}
    corcoef={}
    for smps in tqdm(datNames):
        smp1 = smps[0]
        smp2 = smps[1]

        pdf = PdfPages(savePath+smp1+'_'+smp2+'_jointplots.pdf')
        smp1reads=np.r_[dataDict[smp1+'_pos'].get_nparray(chname,pos+1,(pos+width)),\
                dataDict[smp1+'_neg'].get_nparray(chname,pos+1,(pos+width))]
        smp2reads=np.r_[dataDict[smp2+'_pos'].get_nparray(chname,pos+1,(pos+width)),\
                dataDict[smp2+'_neg'].get_nparray(chname,pos+1,(pos+width))]

        qq=0
        logcorcoef[smp1+'_'+smp2] = np.zeros((len(xran)))
        corcoef[smp1+'_'+smp2] = np.zeros((len(xran)))
        for window_size in tqdm(xran):
            if window_size>1:
                smp1reads = running_mean(smp1reads,window_size)
                smp2reads = running_mean(smp2reads,window_size)

            logcorcoef[smp1+'_'+smp2][qq] = np.corrcoef(np.log10(smp1reads+1e-1),np.log10(smp2reads+1e-1))[1][0]
            corcoef[smp1+'_'+smp2][qq] = np.corrcoef(smp1reads,smp2reads)[1][0]
            pl.figure()
            pl.hexbin(np.log10(smp1reads+1),np.log10(smp2reads+1),bins='log', cmap=pl.cm.YlOrRd_r)
            pl.xlabel(smp1)
            pl.ylabel(smp2)
            pl.title('pearson log-corr: '+str(logcorcoef[smp1+'_'+smp2][qq])+'window '+str(window_size))
            pdf.savefig()
            pl.close()

            pl.figure()
            pl.hexbin(smp1reads,smp2reads,bins='log', cmap=pl.cm.YlOrRd_r)
            pl.xlabel(smp1)
            pl.ylabel(smp2)
            pl.title('pearson corr: '+str(corcoef[smp1+'_'+smp2][qq])+'window '+str(window_size))
            pdf.savefig()
            pl.close()

            qq+=1

        pl.figure()
        pl.plot(np.log10(xran),logcorcoef[smp1+'_'+smp2])
        pl.xlabel('log10 window size')
        pl.ylabel('pearson log-correlation score')
        pl.title(smp1 + ' vs. ' + smp2)
        pdf.savefig()
        pl.close()

        pl.figure()
        pl.plot(np.log10(xran),corcoef[smp1+'_'+smp2])
        pl.xlabel('log10 window size')
        pl.ylabel('pearson correlation score')
        pl.title(smp1 + ' vs. ' + smp2)
        pdf.savefig()
        pl.close()


        pdf.close()

    for datNam1,datNam2 in datNames:
        dataDict[datNam1+'_pos'].close()
        dataDict[datNam1+'_neg'].close()
        dataDict[datNam2+'_pos'].close()
        dataDict[datNam2+'_neg'].close()


    pdf = PdfPages(savePath+'All_corrScores.pdf')
    pl.figure()
    for smps in datNames:
        smp1 = smps[0]
        smp2 = smps[1]
        pl.plot(np.log10(xran),logcorcoef[smp1+'_'+smp2],label=smp1+'_'+smp2)
    pl.legend()
    pl.xlabel('log10 window size')
    pl.ylabel('pearson log-correlation score')
    pl.title('All')
    pdf.savefig()
    pl.close()

    pl.figure()
    for smps in datNames:
        smp1 = smps[0]
        smp2 = smps[1]
        pl.plot(np.log10(xran),corcoef[smp1+'_'+smp2],label=smp1+'_'+smp2)
    pl.legend()
    pl.xlabel('log10 window size')
    pl.ylabel('pearson correlation score')
    pl.title('All')
    pdf.savefig()
    pl.close()

    pdf.close()


def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / N

if __name__ == '__main__':
    main()
