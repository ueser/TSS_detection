from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pylab as pl
import pandas as pd
import numpy as np
import sys, os,subprocess
sys.path.append('/Users/umut/Projects/Toolbox/')
import smooth as sm
import detect_peaks as dp
from optparse import OptionParser
from tqdm import tqdm as tq
import genome.db


usage = 'usage: %prog [options] '
parser = OptionParser(usage)
parser.add_option('-r', dest='runName', default='spt6', type=str, help='Name of the run')
parser.add_option('-s', dest='save_path', type=str, default='/Users/umut/Projects/intragenicTranscription/results/')
parser.add_option('-t', dest='p_val_threshold', type=float, default=0.05)
parser.add_option('-d', dest='debug_idx', type=int, default=0)
parser.add_option('-a', dest='annotation', type=str, default='/Users/umut/Projects/intragenicTranscription/data/annotation/ScerTSSannot.csv')
parser.add_option("-o", action="store_true", dest="overwrite")
(options, args) = parser.parse_args()

options.save_path = options.save_path+options.runName+'/peaks_motifs_p_thresh_'+str(options.p_val_threshold)+'/'

gdb = genome.db.GenomeDB(path='/Users/umut/Projects/genome/data/share/genome_db',assembly='sacCer3')
pl.ioff()

# Track names as recorded in the database
if options.runName == 'spt6':
    WT_tssseq_pos = 'TSSseq_spt6Control_pos'
    WT_tssseq_neg = 'TSSseq_spt6Control_neg'
    condition_tssseq_pos = 'TSSseq_spt6_pos'
    condition_tssseq_neg = 'TSSseq_spt6_neg'

    WT_chipnexus_pos = 'ChIPnexus_spt6Control_pos'
    WT_chipnexus_neg = 'ChIPnexus_spt6Control_neg'
    condition_chipnexus_pos = 'ChIPnexus_spt6_pos'
    condition_chipnexus_neg = 'ChIPnexus_spt6_neg'

elif options.runName == 'diamide':
    WT_tssseq_pos = 'TSSseq_diamControl_pos'
    WT_tssseq_neg = 'TSSseq_diamControl_neg'
    condition_tssseq_pos = 'TSSseq_Dia_pos'
    condition_tssseq_neg = 'TSSseq_Dia_neg'

    WT_chipnexus_pos = 'ChIPnexus_diamControl_pos'
    WT_chipnexus_neg = 'ChIPnexus_diamControl_neg'
    condition_chipnexus_pos = 'ChIPnexus_Dia_pos'
    condition_chipnexus_neg = 'ChIPnexus_Dia_neg'

# Store the pointers into a dictionary
pntr_dict = {'WT_tssseq': [gdb.open_track(WT_tssseq_pos), gdb.open_track(WT_tssseq_neg)],
           options.runName+'_tssseq': [gdb.open_track(condition_tssseq_pos), gdb.open_track(condition_tssseq_neg)],
             'WT_chipnexus': [gdb.open_track(WT_chipnexus_pos), gdb.open_track(WT_chipnexus_neg)],
           options.runName+'_chipnexus': [gdb.open_track(condition_chipnexus_pos), gdb.open_track(condition_chipnexus_neg)]}


def plot_profiles_to_file(annot, pntr, ups=200, smooth_param=50):
    pp = PdfPages(options.save_path + 'Figures/individual_signals.pdf')
    clrs_ = ['red', 'blue', 'black', 'orange', 'magenta', 'cyan']
    vec_sense = {}
    vec_antisense = {}
    # for qq in tq(range(annot.shape[0])):
    for qq in tq(range(100)):

        chname = annot['chr'].iloc[qq]

        if annot['strand'].iloc[qq] == '+':
            start = annot['start'].iloc[qq] - ups
            stop = annot['end'].iloc[qq]
            for key in pntr.keys():
                vec_sense[key] = pntr[key][0].get_nparray(chname, start, stop - 1)
                vec_antisense[key] = pntr[key][1].get_nparray(chname, start, stop - 1)
            xran = np.arange(start, stop)
        else:
            start = annot['start'].iloc[qq]
            stop = annot['end'].iloc[qq] + ups
            for key in pntr.keys():
                vec_sense[key] = np.flipud(pntr[key][1].get_nparray(chname, start, stop))
                vec_antisense[key] = np.flipud(pntr[key][0].get_nparray(chname, start, stop))
            xran = np.arange(stop, start, -1)

        ax = {}
        fig = pl.figure()
        pl.title(annot['name'].iloc[qq])
        for i, key in enumerate(pntr.keys()):
            sm_vec_se = sm.smooth(vec_sense[key], smooth_param)[(smooth_param - 1):-(smooth_param - 1)]
            sm_vec_as = sm.smooth(vec_antisense[key], smooth_param)[(smooth_param - 1):-(smooth_param - 1)]
            ax[key] = pl.subplot(len(pntr), 1, i+1)
            ax[key].plot(xran, vec_sense[key], label=key, color=clrs_[i], alpha=0.5)
            ax[key].plot(xran, -vec_antisense[key], color=clrs_[i], alpha=0.5)
            ax[key].plot(xran, sm_vec_se,  color=clrs_[i], linewidth=2)
            ax[key].plot(xran, -sm_vec_as, color=clrs_[i], linewidth=2)
            ax[key].legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fontsize=6, ncol=1)
        pp.savefig()

        pl.close()
    pp.close()
    for pn in pntr.values():
        pn[0].close()
        pn[1].close()

def main():
    annot = pd.read_csv(options.annotation)
    peaks_df = pd.read_csv(options.save_path + options.runName + '_intragenic_peaks.csv', index_col=0)
    peaks_df = peaks_df[(peaks_df['orientation'] == 'sense') & (peaks_df['fold_change'] > 100) & (peaks_df['type'] == 'intragenic')]

    gene_names = peaks_df['gene'].unique()

    annot = annot[annot['name'].isin(gene_names)]

    print(len(pntr_dict))
    plot_profiles_to_file(annot, pntr_dict)


if __name__ == '__main__':
    main()