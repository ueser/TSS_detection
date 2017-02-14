import numpy as np
import genome.db
import pandas as pd
import matplotlib.pylab as pl
import sys, os,subprocess
sys.path.append('/Users/umut/Projects/Toolbox/')
import smooth as sm
import detect_peaks as dp
from optparse import OptionParser
from tqdm import tqdm as tq


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


def filter_annotation(annot, min_gene_length=500, allow_overlap=False):
    '''
    This function filters the annotation data frame

    :param annot: pandas data frame of the annotation file, contains the following columns: chr, start, end, strand, name
    :param min_gene_length: scalar, minimum gene length to pass the filter (default: 500)
    :param allow_overlap: boolean, if False, filters out the overlapping genes (default: False)
    :return: pandas data frame of filtered annotation
    '''

    if allow_overlap:
        raise NotImplementedError
    return annot[((annot['end']-annot['start']) > min_gene_length) & (annot['name'].str.find(',') == -1)]


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))


def get_peaks(vecs, offset, ups):

    filter_fun = np.vectorize(lambda pk: np.float(np.sum(vecs > vecs[pk])) / len(vecs))
    # detect all peaks in vecs
    peaks = dp.detect_peaks(vecs, mpd=25, show=False)
    if peaks.shape[0] == 0:
        return [], []
    else:
        try:
            return peaks[(filter_fun(peaks) < options.p_val_threshold) & (peaks > (ups+offset))], \
                   peaks[(filter_fun(peaks) < options.p_val_threshold) & (peaks < (ups+25))]
        except IndexError:
            print(peaks.shape, vecs.shape)
            raise


def detect_peaks_save_seq(debug_idx=0):

    f_stats = open(options.save_path+'stats.txt', 'w')

    print('reading annotation...')
    # get the annotation specified form user as csv format
    annot = pd.read_csv(options.annotation)
    f_stats.write('number of genes initially: '+str(annot.shape[0])+'\n')
    # filter the annotation
    annot = filter_annotation(annot, min_gene_length=500, allow_overlap=False)
    f_stats.write('number of genes after filtering : '+str(annot.shape[0])+'\n')

    # open a connection to TSS-seq database
    gdb = genome.db.GenomeDB(path='/Users/umut/Projects/genome/data/share/genome_db', assembly='sacCer3')

    # Track names as recorded in the database
    if options.runName == 'spt6':
        WT_track_A = '81A'
        WT_track_B = '81B'
        condition_track_A = '80A'
        condition_track_B = '80B'

    elif options.runName == 'diamide':
        WT_track_A = 'YPDA'
        WT_track_B = 'YPDB'
        condition_track_A = 'DiaA'
        condition_track_B = 'DiaB'

    elif options.runName == 'diamide_control':
        WT_track_A = 'YPDA'
        WT_track_B = 'YPDB'
        condition_track_A = 'YPDA'
        condition_track_B = 'YPDB'

    elif options.runName == 'spt6_control':
        WT_track_A = '81A'
        WT_track_B = '81B'
        condition_track_A = '81A'
        condition_track_B = '81B'

    elif options.runName == 'diamide_chipnexus':
        WT_track_A = 'ChIPnexus_DiamControl'
        WT_track_B = 'ChIPnexus_DiamControl'
        condition_track_A = 'ChIPnexus_Dia'
        condition_track_B = 'ChIPnexus_Dia'

    elif options.runName == 'nitrogen_chipnexus':
        WT_track_A = 'ChIPnexus_stressControl'
        WT_track_B = 'ChIPnexus_stressControl'
        condition_track_A = 'ChIPnexus_Nit'
        condition_track_B = 'ChIPnexus_Nit'

    elif options.runName == 'spt6_chipnexus':
        WT_track_A = 'ChIPnexus_spt6Control'
        WT_track_B = 'ChIPnexus_spt6Control'
        condition_track_A = 'ChIPnexus_spt6'
        condition_track_B = 'ChIPnexus_spt6'

    elif options.runName == 'aminoacid_chipnexus':
        WT_track_A = 'ChIPnexus_stressControl'
        WT_track_B = 'ChIPnexus_stressControl'
        condition_track_A = 'ChIPnexus_AA'
        condition_track_B = 'ChIPnexus_AA'

    # Store the pointers into a dictionary
    pntr_dict = {WT_track_A: [gdb.open_track(WT_track_A+'_pos'), gdb.open_track(WT_track_A+'_neg'), 'WT_A'],
               condition_track_A: [gdb.open_track(condition_track_A+'_pos'), gdb.open_track(condition_track_A+'_neg'), 'condition_A'],
               WT_track_B: [gdb.open_track(WT_track_B+'_pos'), gdb.open_track(WT_track_B+'_neg'), 'WT_B'],
               condition_track_B: [gdb.open_track(condition_track_B+'_pos'), gdb.open_track(condition_track_B+'_neg'), 'condition_B']}

    # pointer to DNA sequence
    seq = gdb.open_track("seq")

    ups = 100  # upstream of TSS
    offset = 100  # offset downstream from TSS to be considered as "intragenic"
    smooth_param = 50  # smoothing window
    buffSize = np.zeros((annot.shape[0]*10))  # Just to initialize the data frame that holds peaks (for speed)
    tt = 0

    # open the text i/o to record the "pseudo-promoters" for the intragenic TSSs
    f_sense_seq = open(options.save_path+options.runName+'_sense_seq.txt', 'w')
    f_antisense_seq = open(options.save_path+options.runName+'_antisense_seq.txt', 'w')

    # initialize the data frame for peaks
    df = pd.DataFrame({'chr': buffSize,
                       'peak_position': buffSize,
                       'relative_position': buffSize,
                       'strand': buffSize,
                       'fold_change': buffSize,
                       'peak_value': buffSize,
                       'gene': buffSize,
                       'orientation': buffSize,
                       'type': buffSize},
                      index=np.arange(annot.shape[0]*10))

    # loop over the genes in annotation data frame
    for ix in tq(range(debug_idx, annot.shape[0])):
        gene = annot['name'].iloc[ix]
        chname = annot['chr'].iloc[ix]
        strand = annot['strand'].iloc[ix]
        start = annot['start'].iloc[ix]
        end = annot['end'].iloc[ix]

        # loop over the pointers to the database than extract the signal arrays for the given gene from each experiment
        vec, vecAS = {}, {}
        for pntr_key, pntr_val in pntr_dict.items():
            if strand == '+':
                tsP = pntr_val[0].get_nparray(chname, start-ups, end)
                tsN = pntr_val[1].get_nparray(chname, start-ups, end)
                vec[pntr_key] = sm.smooth(tsP, smooth_param)[(smooth_param-1):-(smooth_param-1)]
                vecAS[pntr_key] = sm.smooth(tsN, smooth_param)[(smooth_param-1):-(smooth_param-1)]
            elif strand == '-':
                tsP = pntr_val[0].get_nparray(chname, start, end+ups)
                tsN = pntr_val[1].get_nparray(chname, start, end+ups)
                vec[pntr_key] = np.flipud(sm.smooth(tsN, smooth_param)[(smooth_param-1):-(smooth_param-1)])
                vecAS[pntr_key] = np.flipud(sm.smooth(tsP, smooth_param)[(smooth_param-1):-(smooth_param-1)])


        ### FOR SENSE DIRECTION ###
        # add up the replicates
        vecs = vec[condition_track_A]+vec[condition_track_B]
        # get the peaks located <offset> downstream of TSS
        peaks, tss = get_peaks(vecs,  offset, ups)
        if len(tss) > 0:
            xv = vec[WT_track_A][tss] + vec[WT_track_A][tss]
            xv_s6 = vecs[tss]
            ratios = (xv_s6 + 1e-5) / (xv + 1e-5)
            df.loc[tt:(tt + len(tss) - 1), 'chr'] = [chname] * len(tss)
            df.loc[tt:(tt + len(tss) - 1),
            'peak_position'] = tss + start - ups if strand == '+' else -tss + end + ups
            df.loc[tt:(tt + len(tss) - 1), 'relative_position'] = tss - ups
            df.loc[tt:(tt + len(tss) - 1), 'strand'] = [strand] * len(tss)
            df.loc[tt:(tt + len(tss) - 1), 'fold_change'] = ratios
            df.loc[tt:(tt + len(tss) - 1), 'peak_value'] = xv_s6
            df.loc[tt:(tt + len(tss) - 1), 'gene'] = [gene] * len(tss)
            df.loc[tt:(tt + len(tss) - 1), 'orientation'] = ['sense'] * len(tss)
            df.loc[tt:(tt + len(tss) - 1), 'type'] = ['native'] * len(tss)
            # for each peak position, record the DNA sequence of pseudo-promoter (>200bp upstream) to the txt files
            for ww in range(len(tss)):
                ps = tss[ww] + start - ups - 150 if strand == '+' else -tss[ww] + end + ups - 50
                seqVec = seq.get_seq_str(chname, ps + 1, ps + 200)
                if strand == '-':
                    seqVec = reverse_complement(seqVec)
                f_sense_seq.write('>' + gene + '_pk:' + str(tt + ww) + '\n')
                f_sense_seq.write(seqVec + '\n')
            tt += len(tss)
            
        # if there are some peaks, record the peak attributes
        if len(peaks) > 0:
            xv = vec[WT_track_A][peaks]+vec[WT_track_A][peaks]
            xv_s6 = vecs[peaks]
            ratios = (xv_s6+1e-5)/(xv+1e-5)
            df.loc[tt:(tt+len(peaks)-1), 'chr'] = [chname]*len(peaks)
            df.loc[tt:(tt+len(peaks)-1), 'peak_position'] = peaks+start-ups if strand == '+' else -peaks+end+ups
            df.loc[tt:(tt+len(peaks)-1), 'relative_position'] = peaks-ups
            df.loc[tt:(tt+len(peaks)-1), 'strand']= [strand]*len(peaks)
            df.loc[tt:(tt+len(peaks)-1), 'fold_change'] = ratios
            df.loc[tt:(tt+len(peaks)-1), 'peak_value'] = xv_s6
            df.loc[tt:(tt+len(peaks)-1), 'gene'] = [gene]*len(peaks)
            df.loc[tt:(tt+len(peaks)-1), 'orientation'] = ['sense']*len(peaks)
            df.loc[tt:(tt + len(peaks) - 1), 'type'] = ['intragenic']*len(peaks)
            # for each peak position, record the DNA sequence of pseudo-promoter (>200bp upstream) to the txt files
            for ww in range(len(peaks)):
                ps = peaks[ww]+start-ups-150 if strand == '+' else -peaks[ww]+end+ups-50
                seqVec = seq.get_seq_str(chname, ps+1, ps+200)
                if strand == '-':
                    seqVec = reverse_complement(seqVec)
                f_sense_seq.write('>'+gene+'_pk:'+str(tt+ww)+'\n')
                f_sense_seq.write(seqVec+'\n')
            tt += len(peaks)

        ### FOR ANTISENSE DIRECTION ###
        vecs = vecAS[condition_track_A]+vecAS[condition_track_B]
        peaks, _ = get_peaks(vecs, offset, ups)

        if len(peaks) > 0:
            xv = vecAS[WT_track_A][peaks]+vecAS[WT_track_A][peaks]
            xv_s6 = vecs[peaks]
            ratios = (xv_s6+1e-5)/(xv+1e-5)
            df.loc[tt:(tt+len(peaks)-1), 'chr'] = [chname]*len(peaks)
            df.loc[tt:(tt+len(peaks)-1), 'peak_position'] = peaks+start-ups if strand == '+' else -peaks+end+ups
            df.loc[tt:(tt+len(peaks)-1), 'relative_position'] = peaks-ups
            df.loc[tt:(tt+len(peaks)-1), 'strand'] = ['+']*len(peaks) if strand == '-' else ['-']*len(peaks)
            df.loc[tt:(tt+len(peaks)-1), 'fold_change'] = ratios
            df.loc[tt:(tt+len(peaks)-1), 'peak_value'] = xv_s6
            df.loc[tt:(tt+len(peaks)-1), 'gene'] = [gene]*len(peaks)
            df.loc[tt:(tt+len(peaks)-1), 'orientation'] = ['antisense']*len(peaks)
            df.loc[tt:(tt + len(peaks) - 1), 'type'] = ['intragenic']*len(peaks)
            for ww in range(len(peaks)):
                ps = peaks[ww]+start-ups-50 if strand=='+' else -peaks[ww]+end+ups-150
                seqVec = seq.get_seq_str(chname, ps+1, ps+200)
                if strand == '+':
                    seqVec = reverse_complement(seqVec)
                f_antisense_seq.write('>'+gene+'_pk:'+str(tt+ww)+'\n')
                f_antisense_seq.write(seqVec+'\n')
            tt += len(peaks)
    # close txt i/o connections
    f_sense_seq.close()
    f_antisense_seq.close()
    # trim off the unrecorded part of the data frame
    df = df.iloc[:tt]

    # save the data frame to csv table with column names as headers
    df.to_csv(options.save_path+options.runName+'_intragenic_peaks.csv')
    print('Total ' + str(df.shape[0]) + ' peaks detected...')
    # close the database connections
    for pntr in pntr_dict.values():
        pntr[0].close()
        pntr[1].close()
    seq.close()
    f_stats.close()


def main():

    # Make directory for the project
    if not os.path.exists(options.save_path):
        print('Directory does not exist, creating one...')
        os.makedirs(options.save_path)
    # Make Figure directory for the project if it does not exist
    if not os.path.exists(options.save_path+'Figures'):
        os.makedirs(options.save_path+'Figures')

    print(os.path.isfile(options.save_path+options.runName+'_sense_seq.txt'))
    print(options.overwrite)
    # if "pseudo-promoter" sequence file does not exist, run peak detection, else, continue to MEME suit
    if (not os.path.isfile(options.save_path+options.runName+'_sense_seq.txt')) | options.overwrite:
        print('Sequence file does not exist, detecting peaks and saving sequences ...')
        detect_peaks_save_seq(options.debug_idx)

    meme_db = '/Users/umut/Projects/DataManagement/data/motif_databases/YEAST/YEASTRACT_20130918.meme'

    print('Running meme suit for sense intragenic txp...')
    seq_path = options.save_path+options.runName+'_sense_seq.txt'

    print('     Running dreme for sense intragenic txp...')

    subprocess.call('dreme -v 1 -oc %s/dreme_%s -dna -p %s -norc -t 18000 -e 0.05'
                    % (options.save_path, options.runName + '_sense', seq_path), shell=True)
    print('     Running tomtom for sense intragenic txp...')
    subprocess.call('tomtom -no-ssc -min-overlap 5 -dist pearson -evalue -thresh 10.0 -oc %s/tomtom_%s %s/dreme_%s/dreme.txt %s'
                    % (options.save_path, options.runName+'_sense', options.save_path, options.runName+'_sense', meme_db), shell=True)
    print('     Running centrimo for sense intragenic txp...')
    subprocess.call('centrimo --oc %s/centrimo_%s --verbosity 1 --score 5.0 --ethresh 10.0 %s %s'
                    % (options.save_path,options.runName+'_sense', seq_path, meme_db), shell=True)

    # subprocess.call('fimo -v 1 -oc %s/dreme_%s %s %s' % (options.save_path,options.runName+'_antisense',meme_db,seq_path), shell=True)

    print('Running meme suit for antisense intragenic txp...')
    seq_path = options.save_path+options.runName+'_antisense_seq.txt'
    print('     Running dreme for antisense intragenic txp...')
    subprocess.call('dreme -v 1 -oc %s/dreme_%s -dna -p %s -norc -t 18000 -e 0.05'
                    % (options.save_path,options.runName+'_antisense', seq_path), shell=True)
    print('     Running tomtom for antisense intragenic txp...')
    subprocess.call('tomtom -no-ssc -min-overlap 5 -dist pearson -evalue -thresh 10.0 -oc %s/tomtom_%s %s/dreme_%s/dreme.txt %s'
                    % (options.save_path,options.runName+'_antisense', options.save_path, options.runName+'_antisense', meme_db), shell=True)
    print('     Running centrimo for antisense intragenic txp...')
    subprocess.call('centrimo --oc %s/centrimo_%s --verbosity 1 --score 5.0 --ethresh 10.0 %s %s'
                    % (options.save_path,options.runName+'_antisense', seq_path, meme_db), shell=True)

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
