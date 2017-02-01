import pandas as pd
import numpy as np
import matplotlib.pylab as pl
import h5py as h5
from matplotlib.backends.backend_pdf import PdfPages as pp
import sys,os
sys.path.append('/Users/umut/Projects/genome/python/lib')
sys.path.append('/Users/umut/Projects/deepTSS/analysis')
sys.path.append('/Users/umut/Projects/genomNet/src')
sys.path.append('/Users/umut/Projects/Toolbox/')

from scipy.stats import trim_mean as tm
import roman as rm
import genome.db
from scipy import stats as st
import smooth as sm
from tqdm import tqdm as tq
from optparse import OptionParser


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

    gdb = genome.db.GenomeDB(path='/Users/umut/Projects/genome/data/share/genome_db',assembly='sacCer3')

    # open data 'tracks' for DNase and MNase
    gerp = gdb.open_track('phyloP')



    pl.ioff()

    df= pd.read_csv(options.save_path+options.runName+'_intragenic_peaks.csv', index_col=0)

    consScoreSE = np.zeros((df[df['orientation']=='sense'].shape[0], 1000))

    print('getting conservation score for sense')
    for i in tq(range(len(df[df['orientation']=='sense']))):
        chname = df[df['orientation']=='sense']['chr'].iloc[i]
        st = df[df['orientation']=='sense']['peak_position'].iloc[i]
        strand = df[df['orientation']=='sense']['strand'].iloc[i]

        if strand =='+':
            consScoreSE[i,:]=gerp.get_nparray(chname,st-500+1,st+500)
        else:
            consScoreSE[i,:]=np.flipud(gerp.get_nparray(chname,st-500+1,st+500))

    consScoreAS = np.zeros((df[df['orientation']=='antisense'].shape[0],1000))
    print('getting conservation score for antisense')

    for i in tq(range(len(df[df['orientation']=='antisense']))):
        chname = df[df['orientation']=='antisense']['chr'].iloc[i]
        st = df[df['orientation']=='antisense']['peak_position'].iloc[i]
        strand = df[df['orientation']=='antisense']['strand'].iloc[i]

        if strand =='-':
            consScoreAS[i,:]=gerp.get_nparray(chname,st-500+1,st+500)
        else:
            consScoreAS[i,:]=np.flipud(gerp.get_nparray(chname,st-500+1,st+500))

    print('plotting...')
    pfile = pp(options.save_path+'Figures/intragenic_conservation_sense.pdf')
    xran=np.arange(-500,500)
    tmth = 0.1

    for wn in [3,50,75]:
        fig = pl.figure()
        pl.plot(xran,tm(np.apply_along_axis(sm.smooth,1,consScoreSE,wn)[:,(wn-1):-(wn-1)], tmth, axis=0), 'r', label='Sense')
        pl.xlabel('Position from intragenic TSS (bp)')
        pl.ylabel('Average GERP score (a.u.)')
        pl.title('Smoothing window: '+str(wn)+'bp')
        pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        pfile.savefig()
        pl.close(fig)

    pfile.close()

    pfile = pp(options.save_path+'Figures/intragenic_conservation_antisense.pdf')
    xran=np.arange(-500,500)
    tmth = 0.1

    for wn in [3,50,75]:
        fig = pl.figure()
        pl.plot(xran,tm(np.apply_along_axis(sm.smooth,1,consScoreAS,wn)[:,(wn-1):-(wn-1)],tmth,axis=0),'r',label='Antisense')
        pl.xlabel('Position from intragenic TSS (bp)')
        pl.ylabel('Average GERP score (a.u.)')
        pl.title('Smoothing window: '+str(wn)+'bp')
        pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        pfile.savefig()
        pl.close(fig)

    pfile.close()

    np.savetxt(options.save_path+'ConservationScoreGERP_sense.csv',consScoreSE)
    np.savetxt(options.save_path+'ConservationScoreGERP_antisense.csv',consScoreAS)


    gerp.close()

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
