import pandas as pd
from optparse import OptionParser
from intragenic_peaks_motifs import filter_annotation
import numpy as np
import subprocess

usage = 'usage: %prog [options] '
parser = OptionParser(usage)
parser.add_option('-r', dest='runName', default='spt6', type=str, help='Name of the run')
parser.add_option('-s', dest='save_path', type=str, default='/Users/umut/Projects/intragenicTranscription/results/')
parser.add_option('-t', dest='p_val_threshold', type=float, default=0.05)
parser.add_option('-d', dest='debug_idx', type=int, default=0)
parser.add_option('-a', dest='annotation', type=str, default='/Users/umut/Projects/intragenicTranscription/data/annotation/ScerTSSannot.csv')
parser.add_option("-o", action="store_true", dest="overwrite")
(options, args) = parser.parse_args()

options.save_path = options.save_path + options.runName+'/peaks_motifs_p_thresh_'+str(options.p_val_threshold)+'/'


def main():

    f1 = np.vectorize(lambda x, y: int(x) - 200 if y == '+' else int(x))
    f2 = np.vectorize(lambda x, y: int(x) if y == '+' else int(x) + 200)

    print('processing intragenic peaks')
    peaks_df = pd.read_csv(options.save_path + options.runName+'_intragenic_peaks.csv', index_col=0)
    peaks_df = peaks_df[peaks_df['orientation'] == 'sense']
    peaks_df['start'] = f1(peaks_df['peak_position'], peaks_df['strand'])
    peaks_df['stop'] = f2(peaks_df['peak_position'], peaks_df['strand'])
    peaks_df[['chr', 'start', 'stop', 'strand', 'gene', 'peak_value', 'fold_change', 'type']].to_csv(
        options.save_path + 'all_promoter_regions.bed', sep='\t', header=None, index=None)

    del peaks_df

    # print('processing intergenic peaks')
    # # get the annotation specified form user as csv format
    # annot = pd.read_csv(options.annotation)
    # # filter the annotation
    # annot = filter_annotation(annot, min_gene_length=500, allow_overlap=False)
    #
    # annot['start'] = f1(annot['tss'], annot['strand'])
    # annot['stop'] = f2(annot['tss'], annot['strand'])
    # annot['name'] = annot['name'].map(lambda x: 'i'+x)
    # annot[['chr', 'start', 'stop', 'strand', 'name']].to_csv(options.save_path + 'intergenic_promoter_regions.bed',
    #                                                          sep='\t', header=None, index=None)
    # del annot

    print('processing fimo hits')
    # fimo_df = pd.read_csv('/Users/umut/Projects/ClassifySpecies/analysis/Fimo_Yeastract_SC.csv', sep='\t')
    fimo_df = pd.read_csv('/Users/umut/Projects/intragenicTranscription/analysis/fimo_mcisaac.tsv', sep='\t')

    fimo_df[['sequence name', 'start', 'stop', 'strand', 'TFlist', 'q-value', 'p-value', 'score']].to_csv(
        options.save_path + 'fimo_regions.bed', sep='\t', header=None, index=None)

    toRun = 'bedtools intersect -wa -wb -a {} -b {} >> {}'.format(options.save_path + 'fimo_regions.bed',
                                                                       options.save_path + 'all_promoter_regions.bed',
                                                                   options.save_path + 'TF_binding.txt')
    header_list = ['sequence name', 'start', 'stop', 'strand', 'TFlist', 'q-value', 'p-value', 'score'] + \
                  ['chr', 'start', 'stop', 'strand', 'gene', 'peak_value', 'fold_change', 'type']

    with open(options.runName + '_tf_binding_run.sh', 'w') as f:
        f.write('echo {} > {}\n\n'.format('\t'.join(header_list), options.save_path + 'TF_binding.txt'))
        f.write(toRun)

    print('Run the following command: \n $ bash {}_tf_binding_run.sh'.format(options.runName))


if __name__ == '__main__':
    main()
