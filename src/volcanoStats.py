import numpy as np
import genome.db
import pandas as pd
import matplotlib.pylab as pl
import sys, os
sys.path.append('/Users/umut/Projects/Toolbox/')
import smooth as sm
from optparse import OptionParser
from tqdm import tqdm as tq

from scipy.stats import ttest_ind, ttest_ind_from_stats



def find_peaks(pntr_list,annot,ups=100,offset=100,smooth_param=50):
    """
        ups = 100 # how long to include upstream of TSS
        offset = 100 # how much to discard downstream of TSS for genic TSS detection
        smooth_param = 50 # approximate width of the TSS initiation region
        pntr_list: List of genome_db pointers formatted as [[pnt_pos,pnt_neg,Label], [pnt2_pos,pnt2_neg,Label2],...]
    """
    for pntr in pntr_list:
        df_peak= pd.DataFrame({ pntr[2]+'_Native_value':np.zeros((annot.shape[0])),
                            pntr[2]+'_Native_position':np.zeros((annot.shape[0])),
                            pntr[2]+'_Intragenic_value':np.zeros((annot.shape[0])),
                            pntr[2]+'_Intragenic_position':np.zeros((annot.shape[0])),
                            pntr[2]+'_Native_valueAS':np.zeros((annot.shape[0])),
                            pntr[2]+'_Native_positionAS':np.zeros((annot.shape[0])),
                            pntr[2]+'_Intragenic_valueAS':np.zeros((annot.shape[0])),
                            pntr[2]+'_Intragenic_positionAS':np.zeros((annot.shape[0]))
                          },index = annot['name'])
    for ix in tq(range(annot.shape[0])):
        gene = annot['name'].iloc[ix]
        chname = annot['chr'].iloc[ix]
        strand = annot['strand'].iloc[ix]
        start = annot['start'].iloc[ix]
        end = annot['end'].iloc[ix]
        for pntr in pntr_list:

            if strand=='+':
                tsP = pntr[0].get_nparray(chname,start-ups,end)
                tsN = pntr[1].get_nparray(chname,start-ups,end)
                vec = sm.smooth(tsP,smooth_param)[(smooth_param-1):-(smooth_param-1)]
                vecAS = sm.smooth(tsN,smooth_param)[(smooth_param-1):-(smooth_param-1)]
            elif strand=='-':
                tsP = pntr[0].get_nparray(chname,start,end+ups)
                tsN = pntr[1].get_nparray(chname,start,end+ups)
                vec = np.flipud(sm.smooth(tsN,smooth_param)[(smooth_param-1):-(smooth_param-1)])
                vecAS = np.flipud(sm.smooth(tsP,smooth_param)[(smooth_param-1):-(smooth_param-1)])

            df_peak.loc[gene,pntr[2]+'_Native_value'] = np.max(vec[:(ups+offset)])
            df_peak.loc[gene,pntr[2]+'_Native_position'] = np.argmax(vec[:(ups+offset)])-ups
            df_peak.loc[gene,pntr[2]+'_Intragenic_value'] = np.max(vec[(ups+offset):])
            df_peak.loc[gene,pntr[2]+'_Intragenic_position'] = np.argmax(vec[(ups+offset):])+offset

            df_peak.loc[gene,pntr[2]+'_Native_valueAS'] = np.max(vecAS[:(ups+offset)])
            df_peak.loc[gene,pntr[2]+'_Native_positionAS'] = np.argmax(vecAS[:(ups+offset)])-ups
            df_peak.loc[gene,pntr[2]+'_Intragenic_valueAS'] = np.max(vecAS[(ups+offset):])
            df_peak.loc[gene,pntr[2]+'_Intragenic_positionAS'] = np.argmax(vecAS[(ups+offset):])+offset

    return df_peak

def get_t_stat(xv,xv_s6):
    # Compute the descriptive statistics of a and b.
    abar = xv.mean(axis=1).values
    avar = np.sqrt(xv.var(ddof=1,axis=1)).values
    na = 2
    adof = 1

    bbar = xv_s6.mean(axis=1).values
    bvar = np.sqrt(xv_s6.var(ddof=1,axis=1)).values
    nb = 2
    bdof =1

    t2=np.zeros((len(abar)))
    p2=np.zeros((len(abar)))
    for ix in tq(range(len(abar))):
    # Use scipy.stats.ttest_ind_from_stats.
        t2[ix], p2[ix] = ttest_ind_from_stats(abar[ix], avar[ix], na,
                                  bbar[ix], bvar[ix], nb,
                                  equal_var=False)
    return t2,p2,(bbar-abar)

def filter_annotation(annot,min_gene_length=500,allow_overlap=False):
    tff = (annot['end']-annot['start'])>min_gene_length
    tff = tff&(annot['name'].str.find(',')==-1)
    return annot[tff]

def plot_volcano(logFC,p_val,sample_name,saveName,logFC_thresh):
    fig=pl.figure()
    ## To plot and save
    pl.scatter(logFC[(p_val>0.05)|(abs(logFC)<logFC_thresh)],-np.log10(p_val[(p_val>0.05)|(abs(logFC)<logFC_thresh)]),color='blue',alpha=0.5);
    pl.scatter(logFC[(p_val<0.05)&(abs(logFC)>logFC_thresh)],-np.log10(p_val[(p_val<0.05)&(abs(logFC)>logFC_thresh)]),color='red');
    pl.hlines(-np.log10(0.05),min(logFC),max(logFC))
    pl.vlines(-logFC_thresh,min(-np.log10(p_val)),max(-np.log10(p_val)))
    pl.vlines(logFC_thresh,min(-np.log10(p_val)),max(-np.log10(p_val)))
    pl.xlim(-3,3)
    pl.xlabel('Log Fold Change')
    pl.ylabel('-log10(p-value)')
    pl.savefig(saveName)
    pl.close(fig)


# def plot_histograms(df_peaks,pntr_list):
#
#     for pntr in pntr_list:
#         colName =pntr[2]+'_Intragenic_position'
#         pl.hist(df_peaks[colName])
#         pl.xlabel(colName)
#         pl.ylabel()
#         pl.show()



def main():

    usage = 'usage: %prog [options] '
    parser = OptionParser(usage)
    parser.add_option('-r', dest='runName', default='Intragenic', type=str, help='Align sizes with batch size')
    parser.add_option('-s', dest='save_path', type=str, default='/Users/umut/Projects/intragenicTranscription/results/')
    parser.add_option('-a', dest='annotation', type=str, default='/Users/umut/Projects/intragenicTranscription/data/annotation/ScerTSSannot.csv')
    parser.add_option('-t', dest='fold_change_threshold', type=float, default=10)

    (options,args) = parser.parse_args()



    options.save_path = options.save_path+options.runName+'/FC_threshold_'+str(options.fold_change_threshold)+'/'

    # Make directory for the project
    if not os.path.exists(options.save_path):
        print('Directory does not exist, creating one...')
        os.makedirs(options.save_path)
    # Make directory for the project
    if not os.path.exists(options.save_path+'Figures'):
        os.makedirs(options.save_path+'Figures')

    f_stats = open(options.save_path+'stats.txt','w')

    print('reading annotation...')
    # get the annotation specified form user
    annot = pd.read_csv(options.annotation)
    f_stats.write('number of genes initially: '+str(annot.shape[0])+'\n')
    annot = filter_annotation(annot,min_gene_length=500,allow_overlap=False)
    f_stats.write('number of genes after filtering : '+str(annot.shape[0])+'\n')



    gdb = genome.db.GenomeDB(path='/Users/umut/Projects/genome/data/share/genome_db',assembly='sacCer3')

    if options.runName=='spt6':
        WT_track_A = '81A'
        WT_track_B = '81B'
        condition_track_A = '80A'
        condition_track_B = '80B'
    elif options.runName=='diamide':
        WT_track_A = 'YPDA'
        WT_track_B = 'YPDB'
        condition_track_A = 'DiaA'
        condition_track_B = 'DiaB'


    pntr_list=[[gdb.open_track(WT_track_A+'_pos'),gdb.open_track(WT_track_A+'_neg'),'WT_A'],\
    [gdb.open_track(condition_track_A+'_pos'),gdb.open_track(condition_track_A+'_neg'),'condition_A'],\
    [gdb.open_track(WT_track_B+'_pos'),gdb.open_track(WT_track_B+'_neg'),'WT_B'],\
    [gdb.open_track(condition_track_B+'_pos'),gdb.open_track(condition_track_B+'_neg'),'condition_B']]

    print('finding peaks...')
    df_peaks = find_peaks(pntr_list,annot)
    df_peaks.to_csv(options.save_path+options.runName+'_peaks.csv')
    # plot_histograms()

    log10threshold = np.log10(options.fold_change_threshold)
    for cond in ['Native','Intragenic']:

        print('calculating statistics for '+cond+' ...')
        for direction in ['sense','antisense']:
            if direction=='sense':
                suffix=''
            else:
                suffix='AS'
            WT = df_peaks[['WT_A_'+cond+'_value'+suffix,'WT_B_'+cond+'_value'+suffix]]+1e-2
            condition = df_peaks[['condition_A_'+cond+'_value'+suffix,'condition_B_'+cond+'_value'+suffix]]+1e-2
            t_stat,p_val,logFC = get_t_stat(WT,condition)
            saveName= options.save_path+'Figures/'+options.runName+'_'+cond+suffix+'.pdf'

            print('saving the volcano plot  for '+cond+' TSS ...')
            plot_volcano(logFC,p_val,cond+'_TSS',saveName,log10threshold)

            repressed = df_peaks.index[(p_val<0.05)&(logFC<=-log10threshold)]
            activated = df_peaks.index[(p_val<0.05)&(logFC>=log10threshold)]

            f_stats.write('number of repressed txp in '+cond+'_'+direction+': '+str(len(repressed))+'\n')
            f_stats.write('number of activated txp in '+cond+'_'+direction+': '+str(len(activated))+'\n')


            print('saving the differential genes for '+cond+' TSS ...')

            with open(options.save_path+''+cond+'_TSS_'+options.runName+'_repressed'+suffix+'.txt','w') as f:
                for ln in repressed.tolist():
                    f.write(ln+'\n')
            with open(options.save_path+''+cond+'_TSS_'+options.runName+'_activated'+suffix+'.txt','w') as f:
                for ln in activated.tolist():
                    f.write(ln+'\n')

    for pntr in pntr_list:
        pntr[0].close()
        pntr[1].close()
    f_stats.close()

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
