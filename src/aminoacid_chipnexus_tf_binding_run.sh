echo sequence name	start	stop	strand	TFlist	q-value	p-value	score	chr	start	stop	strand	gene	peak_value	fold_change	type > /Users/umut/Projects/intragenicTranscription/results/aminoacid_chipnexus/peaks_motifs_p_thresh_0.05/TF_binding.txt

bedtools intersect -wa -wb -a /Users/umut/Projects/intragenicTranscription/results/aminoacid_chipnexus/peaks_motifs_p_thresh_0.05/fimo_regions.bed -b /Users/umut/Projects/intragenicTranscription/results/aminoacid_chipnexus/peaks_motifs_p_thresh_0.05/all_promoter_regions.bed >> /Users/umut/Projects/intragenicTranscription/results/aminoacid_chipnexus/peaks_motifs_p_thresh_0.05/TF_binding.txt