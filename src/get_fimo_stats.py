
import pandas as pd
import numpy as np
from scipy.stats import poisson
import os
import sys
sys.path.append('/Users/umut/Projects/genome/python/lib')
sys.path.append('/Users/umut/Projects/genomNet/src')
import roman as rm
#
import matplotlib.pylab as pl
# from matplotlib.backends.backend_pdf import PdfPages as pp
import seaborn as sns
from tqdm import tqdm as tq
from optparse import OptionParser
import pdb


def main():

    usage = 'usage: %prog [options] '
    parser = OptionParser(usage)
    parser.add_option('-r', dest='runName', default='spt6', type=str, help='Name of the run')
    parser.add_option('-s', dest='save_path', type=str, default='/Users/umut/Projects/intragenicTranscription/results/')
    parser.add_option('-a', dest='annotation', type=str, default='/Users/umut/Projects/intragenicTranscription/data/annotation/ScerTSSannot.csv')
    parser.add_option("-o", action="store_true", dest="overwrite")
    (options,args) = parser.parse_args()

    # options.save_path = options.save_path+options.runName+'/peaks_motifs/'
    options.save_path = \
        options.save_path + options.runName + '/peaks_motifs_p_thresh_0.05/'
    # Make directory for the project
    if not os.path.exists(options.save_path):
        print('Directory does not exist, creating one...')
        os.makedirs(options.save_path)
    # Make directory for the project
    if not os.path.exists(options.save_path+'Figures'):
        os.makedirs(options.save_path+'Figures')

    yeastractSC = pd.read_csv('/Users/umut/Projects/ClassifySpecies/analysis/Fimo_Yeastract_SC.csv',sep='\t')

    # tf =  (yeastractSC['sequence name'] == 'chrII') | (yeastractSC['sequence name'] == 'chrIV')
    # yeastractSC =yeastractSC[tf]
    df = pd.read_csv(options.save_path+options.runName+'_intragenic_peaks.csv',index_col=0)

    for direction in ['both']:
        if not os.path.isfile(options.save_path+'intragenic_upsVSdws_TFbsFimo_'+direction+'.csv'):
            print('TF density file does not exist, calculating ...')

            print 'Calculating TF density for intragenic sense pseudo-promoters '
            senTFups, senLenups, senTFdws, senLendws = getTFdensity(df[df['orientation'] == 'sense'],
                                                                    direction=direction, yeastract=yeastractSC)
            print 'Calculating TF density for intragenic antisense pseudo-promoters'
            asenTFups, asenLenups, asenTFdws, asenLendws = getTFdensity(df[df['orientation'] == 'antisense'],
                                                                        direction=direction, yeastract=yeastractSC)
            pVals = pd.DataFrame({  'pVal_senUpsvsDws_'+direction: 1,
                                    'pVal_asenUpsvsDws_'+direction: 1,
                                    'sen_densityUps_'+direction: 0,
                                    'sen_densityDws_'+direction: 0,
                                    'asen_densityUps_'+direction: 0,
                                    'asen_densityDws_'+direction: 0},
                                 index=senTFups.index.values,
                                 dtype=float)
            for TFac in senTFups.index.values:
                # pdb.set_trace()
                pVals.loc[TFac, 'pVal_senUpsvsDws_'+direction] = getPvalue(senTFups[TFac],senLenups,senTFdws[TFac],senLendws)
                pVals.loc[TFac, 'pVal_asenUpsvsDws_'+direction] = getPvalue(asenTFups[TFac],asenLenups,asenTFdws[TFac],asenLendws)

                pVals.loc[TFac, 'sen_densityUps_'+direction] = np.float(senTFups[TFac])/np.float(senLenups+1)
                pVals.loc[TFac, 'sen_densityDws_'+direction] = np.float(senTFdws[TFac])/np.float(senLendws+1)
                pVals.loc[TFac, 'asen_densityUps_'+direction] = np.float(asenTFups[TFac])/np.float(asenLenups+1)
                pVals.loc[TFac, 'asen_densityDws_'+direction] = np.float(asenTFdws[TFac])/np.float(asenLendws+1)
            print pVals.head()
            pVals.to_csv(options.save_path+'intragenic_upsVSdws_TFbsFimo_'+direction+'.csv')

        sns.set_style('white')

        pVals= pd.read_csv(options.save_path+'intragenic_upsVSdws_TFbsFimo_'+direction+'.csv', index_col=0)

        print('Plotting for sense enriched ...')
        tf = (pVals['pVal_senUpsvsDws_'+direction]<0.05)
        pValsTmp = pVals[tf]
        idx = np.argsort(pValsTmp['pVal_senUpsvsDws_'+direction])
        pValsTmp = pValsTmp.iloc[idx]
        fig = pl.figure()
        pl.bar(np.arange(len(pValsTmp)),-np.log10(pValsTmp['pVal_senUpsvsDws_'+direction]),0.8,color='#d9544d',label = direction +' Strand')
        pl.xticks(np.arange(len(pValsTmp))+.4, pValsTmp.index.values,rotation='vertical')
        pl.ylabel('Log10 p-value')
        pl.legend()
        pl.savefig(options.save_path+'Figures/intragenic_sense_EnrichedTFs_'+direction+'Strand.pdf')
        pl.close()
        print('Plotting for sense depleted ...')

        pValsTmp['pVal_senUpsvsDws_'+direction] = 1-pVals['pVal_senUpsvsDws_'+direction]
        tf = (pValsTmp['pVal_senUpsvsDws_'+direction]<0.05)
        pValsTmpTmp = pValsTmp[tf]
        idx = np.argsort(pValsTmp['pVal_senUpsvsDws_'+direction])
        pValsTmp = pValsTmp.iloc[idx]

        fig = pl.figure()
        pl.bar(np.arange(len(pValsTmp)),-np.log10(pValsTmp['pVal_senUpsvsDws_'+direction]),0.8,color='#d9544d',label = direction+' Strands')
        pl.xticks(np.arange(len(pValsTmp))+.4, pValsTmp.index.values,rotation='vertical')
        pl.ylabel('Log10 p-value')
        pl.legend()
        pl.savefig(options.save_path+'Figures/intragenic_sense_DepletedTFs_'+direction+'Strand.pdf')
        pl.close()

        ### for antisense intragenic
        print('Plotting for antisense enriched ...')

        tf = (pVals['pVal_asenUpsvsDws_'+direction]<0.05)
        pValsTmp = pVals[tf]
        idx = np.argsort(pValsTmp['pVal_asenUpsvsDws_'+direction])
        pValsTmp = pValsTmp.iloc[idx]

        fig = pl.figure()
        pl.bar(np.arange(len(pValsTmp)),-np.log10(pValsTmp['pVal_asenUpsvsDws_'+direction]),0.8,color='#d9544d',label = direction +' Strand')
        pl.xticks(np.arange(len(pValsTmp))+.4, pValsTmp.index.values,rotation='vertical')
        pl.ylabel('Log10 p-value')
        pl.legend()
        pl.savefig(options.save_path+'Figures/intragenic_asense_EnrichedTFs_'+direction+'Strand.pdf')
        pl.close()
        print('Plotting for antisense depleted ...')

        pValsTmp['pVal_asenUpsvsDws_'+direction] = 1-pVals['pVal_asenUpsvsDws_'+direction]
        tf = (pValsTmp['pVal_asenUpsvsDws_'+direction]<0.05)#| (pValsTmp['pVal_se']<0.05) | (pValsTmp['pVal_as']<0.05)
        pValsTmpTmp = pValsTmp[tf]
        idx = np.argsort(pValsTmp['pVal_asenUpsvsDws_'+direction])
        pValsTmp = pValsTmp.iloc[idx]

        fig = pl.figure()
        pl.bar(np.arange(len(pValsTmp)),-np.log10(pValsTmp['pVal_asenUpsvsDws_'+direction]),0.8,color='#d9544d',label = direction+' Strands')
        pl.xticks(np.arange(len(pValsTmp))+.4, pValsTmp.index.values,rotation='vertical')
        pl.ylabel('Log10 p-value')
        pl.legend()
        pl.savefig(options.save_path+'Figures/intragenic_asense_DepletedTFs_'+direction+'Strand.pdf')

        pl.close()


def getTFdensity(annot,direction='sense', yeastract=None):

    tln=0
    tlnDws =0
    ttf = pd.Series(0, index=yeastract['TFlist'].unique())
    ttfDws = pd.Series(0, index=yeastract['TFlist'].unique())

    for i in tq(range(len(annot))):
        strn = annot['strand'].iloc[i]
        ch = annot['chr'].iloc[i]

        st = annot['peak_position'].iloc[i]-200 if strn == "+" else annot['peak_position'].iloc[i]
        en = annot['peak_position'].iloc[i] if strn == "+" else annot['peak_position'].iloc[i]+200


        tf1 = yeastract['start']>=st
        tf2 = yeastract['stop']<=en
        tf3 = yeastract['sequence name']==ch

        if direction == 'sense':
            tfa = yeastract['strand'] == strn
        elif direction == 'antisense':
            tfa = yeastract['strand'] == '-' if strn=='+' else yeastract['strand'] == '+'
        elif direction == 'both':
            tfa = np.ones((len(yeastract['strand'])),dtype=bool)

        df = yeastract[tf1&tf2&tf3&tfa]

        ttf = ttf.add(df['TFlist'].value_counts(),fill_value=0)
        tln += en-st

        st = annot['peak_position'].iloc[i] if strn == "+" else annot['peak_position'].iloc[i]-200
        en = annot['peak_position'].iloc[i]+200 if strn == "+" else annot['peak_position'].iloc[i]

        tf1 = yeastract['start']>=st
        tf2 = yeastract['stop']<=en

        df = yeastract[tf1&tf2&tf3&tfa]

        if any(df['TFlist'].value_counts()):
            if (en<st)| any(df['TFlist'].value_counts() < 0):
                print 'Alert!!! something is wrong!'
            ttfDws = ttfDws.add(df['TFlist'].value_counts(),fill_value=0)
            tlnDws += en-st
        # else:
        #     print df['TFlist'].value_counts()

    return ttf, tln, ttfDws, tlnDws


def getPvalue(x1,N1,x0,N0):
    mu = np.float(x0)*np.float(N1+1)/np.float(N0+1)
    return 1.0 - poisson.cdf(np.float(x1),mu)

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
