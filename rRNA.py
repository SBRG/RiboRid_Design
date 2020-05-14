class RRNA():

	import os
	from Bio import SeqIO, AlignIO


	#TODO: make it so we can pass a dictonary for the riborid conditions e.g. [Na] and temp.
	def __init__(self, name, gbk_files, rtype, outdir='rrd', pre=0, oligos_df):
		"""
		Parameters
		------------
		name: string
			name of the rRNA, usually 23S, 16S or 5S
		rRNA: string
			path to fasta file containing sequences for 23S, 16S or 5S rRNA
		oligo_db: string
			path to fasta file containing sequences for old oligos that you want to reuse
		"""
  		if rtpye not in ['23S', '16S', '5S']:
  			raise AttributeError('rtype must be 23S, 16S or 5S')

		self.name = name
		self.outdir = outdir
		if not os.path.isdir(self.outdir):
			os.mkdir(self.outdir)

	    self.pre = pre
		oname = '{}/{}_'.format(self.outdir,self.name)

		self.gbk_files = gbk_files
		if self.gbk_files and type(self.gbk_files) == str:
	        self.gbk_files = [self.gbk_files]

	    self.rRNA_fa = rRNA_fa
	    if not rRNA_fa:
	    	self.rRNA_fa = os.path.join(self.outdir, '{}_{}.fa'.format(self.name, self.rtype))
			get_rRNA()

		self.consensus = get_consensus(self.rRNA, clst=oname + 'clustal.clst',
			cs=oname + 'consensus_seq', cs_name=self.name + '_consensus_seq')

		self.oligos_df = oligos_df

		
	def get_rRNA(self):
	    """ Generate fasta files of all rRNA in genbank file
	        Parameters
	        ----------
	        gbk: string or list
	        	location of the genbank file, or a list containing locations
	        of multiple genbank file
	        fasta_out: string, Default 'rRNA'
	        	name of output fasta file
	        pre: int, default None
	        	number of preRNA base-pairs to include for removal
	    """	    
        for gbk in self.gbk_files:
            rRNA_count = 0
            for ref_seq in SeqIO.parse(gbk,'genbank'):
                for feat in [f for f in ref_seq.features if f.type == 'rRNA']:
                    product = feat.qualifiers['product'][0]
                    #distinguish between 23S, 16S and 5S
                    if product.startswith(self.rtype):
                    	write_fasta(feat, ref_seq)
                    	rRNA_count += 1
            if rRNA_count == 0:
                raise ValueError('No {} rRNA found in {}. Please reannotate with prokka.'.format(self.rtype, gbk))
            else:
                print((str(rRNA_count) + ' rRNA sequences found in ' + gbk))

    def write_fasta(self, feat, ref_seq):
    	"""
    	Write the feature sequence into the fasta file
    	"""
    	with open(self.rRNA, 'w+') as fa_out:
    		start, end = feat.location.start - pre, feat.location.end
    		seq = ref_seq.seq[start: end]
	    	if feat.strand == -1:
	    		seq = seq.reverse_complement()
    		fa_out.write('>{}|{}\n{}\n'.format(feat.qualifiers['locus_tag'][0],
                                               feat.qualifiers['product'][0],
                                               str(seq)))

	def get_consensus(self, fasta, clst='clustal.clst', cs_file='consensus.fa', cs_name='consensus_seq'):
		""" 
		Generates consensus sequence for all the fasta sequences in the
	    in the fasta file. Internally runs MUSCLE multiple sequence alignment
	    with default parameters. Usually not necessary if only working with one
	    strain as all copies of their rRNA tend to be identical(or close to it).
	    Parameters
	    ----------
	    fasta: fasta files containing sequences to be aligned
	    clst: name of output clustal file
	    cs_file: fasta file where consensus sequence will be written
	    cs_name: name for consensus sequence
	    """
	    
	    if not cs_file.endswith('.fa'):
	        cs_file = cs_file + '.fa'

	    if os.path.isfile(cs_file) and cs_name in [r.id for r in SeqIO.parse(cs_file, 'fasta')]:
	        raise Exception('Consensus sequence for {} already exists in {}'.format(cs_name, cs_file))

	    #get multiple sequence alignment with clstalo
	    call = ['muscle', '-in', fasta, '-out', clst, '-clw']
	    subprocess.call(call)
	    
	    #get consensus sequence from genetics file
	    align = AlignIO.read(open(clst), "clustal")
	    summary_align = AlignInfo.SummaryInfo(align)

	    with open(cs_file, 'a+') as cs_out:
	        cs_out.write('>{}\n{}\n'.format(cs_name, str(consensus)))

	    return cs_file
