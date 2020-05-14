# do all the heavy lifting here
# gotta get all the imports for different classes
"""
### THE NEXT THING TO FIGURE OUT HERE IS HOW TO HANDLE AND STORE INDIVIDUAL OLIGOS
#TODO: oligos
1. get all the rRNA sequence from the genbank file and store it in fasta file.
1. We need to get all the fasta files and get the rRNA sequences
2. Use the rRNA sequences to generate a rRNA object
3. Take the old oligos and align them
4. If those oligos align, then make an object? Maybe make object in 3. Add that object to rRNA attr.
5. Then calculate the gaps and add that as attributes
6. Maybe all this should be done on initialization


1. Do all of those of each of the three types of rRNA
2. After you run it for one of the rRNA, add those oligos into the od oligos df so you can see if you can reuse them
3. 
"""
	#TODO: if preRNA, we need to add preRNA for each of the three rRNA type


def design_oligos(gbk_files, oligo_fa=None):

	exp = Experiment()

	#TODO: this should not be hardcoded
	#might be ok for now since bacteria mostly have these
	r23 = RRNA(gbk_files)
	r16 = RRNA(gbk_files)  
	r5 =  RRNA(gbk_files)

	for r in [r23, r16, r5]:

		if oligo_fa:
			exp.align_oligos(r, oligo_fa)
			r.oligos_df = exp.find_old_oligos(r, oligo_fa)

		new_oligos = exp.gapfill(r, oligo_fa)
		if r.oligos_df:
			r.oligos_df = pd.concat([r.oligos_df, new_oligos])
		else:
			r.oligos_df = new_oligos

	return r23, r16, r5	