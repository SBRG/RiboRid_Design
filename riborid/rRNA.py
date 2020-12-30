from os import path, mkdir
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo


class RRNA:
    # TODO: what to do if rrna_fa already exists.
    # TODO: Figure out how to implement oligo_df and check if it messes with gap filling
    def __init__(self, name, rtype, infile, ftype, outdir='rrd', pre=0, oligos_df=None):
        """
        Parameters
        ------------
        name: string
            name of the rRNA, usually 23S, 16S or 5S
        rtype: str
            type of rRNA e.g. 16S, 23S, 5S
        infile: string, list
            path to input file or list of paths to multiple input files; input files must be
            either genbank or fasta.
        ftype: string
            type of input file; 'genbank' or 'fasta'
        outdir: str, default rrd
            path to the output directory
        pre: int, default 0
            number of base-pairs upstream of rRNA to generate oligos for. Must provide genbank
            files to use this option.
        oligos_df: string
            path to csv file containing sequences for old oligos that you want to reuse
        """
        # run param checks
        if ftype not in ['fasta', 'genbank']:
            raise AttributeError('ftype must be either \'fasta\' or \'genbank\' ')
        if pre !=0 and ftype=='fasta':
            raise AttributeError('Cannot pass both rrna fasta and pre rRNA')

        self.name = name
        self.rtype = rtype
        self.outdir = outdir
        self.pre = pre
        self.infile = infile
        self.oname = path.abspath(f'{self.outdir}/{self.name}_{self.rtype}')
        self.ftype = ftype
        self.consensus = self.get_consensus(clst=self.oname + '_clustal.clst',
                                            cs_file=self.oname + '_consensus',
                                            cs_name=self.name + '_consensus')

        self.oligos_df = oligos_df

    @property
    def ftype(self):
        """Get loc of rrna fasta file"""
        return self.__ftype

    @ftype.setter
    def ftype(self, ftype):
        if ftype == 'fasta':
            self.rrna_fa = path.abspath(self.infile)
        else:  #generate fasta file from infile infile
            self.rrna_fa = self.oname + '.fa'
            self.get_rRNA()
        self.__ftype = ftype

    @property
    def outdir(self):
        """ Get dir for the output files"""
        return self.__outdir

    @outdir.setter
    def outdir(self, outdir):
        if not path.isdir(outdir):
            mkdir(outdir)
        self.__outdir = outdir

        #update full outdir name
        self.__oname = f'{self.outdir}/{self.name}_{self.rtype}_'

    @property
    def infile(self):
        """ Get the list of gbk_files"""
        return self.__infile

    @infile.setter
    def infile(self, infile):
        if type(infile) == str:
            infile = [infile]
        self.__infile = infile



    def get_rRNA(self):
        """ Generate fasta files of all rRNA in genbank file"""
        for gbk in self.infile:
            rrna_count = False

            for ref_seq in SeqIO.parse(gbk, 'genbank'):
                for feat in [f for f in ref_seq.features if f.type == 'rRNA']:
                    product = feat.qualifiers['product'][0]
                    # distinguish between 23S, 16S and 5S
                    if product.startswith(self.rtype):
                        self.write_fasta(feat, ref_seq)
                        rrna_count = True
            if not rrna_count:
                raise ValueError(f'No {self.rtype} rRNA found in {gbk}. Please reannotate with prokka.')
            else:
                print((str(rrna_count) + ' rRNA sequences found in ' + gbk))

    def write_fasta(self, feat, ref_seq):
        """
        Write the feature sequence into the fasta file
        """
        with open(self.rrna_fa, 'a+') as fa_out:
            start, end = feat.location.start - self.pre, feat.location.end

            if feat.strand == -1:
                seq = ref_seq.seq[start: end + self.pre]
                seq = seq.reverse_complement()
            else:
                seq = ref_seq.seq[start - self.pre: end]

            fa_out.write('>{}|{}\n{}\n'.format(feat.qualifiers['locus_tag'][0],
                                               feat.qualifiers['product'][0],
                                               str(seq)))

    def get_consensus(self, clst='clustal.clst', cs_file='consensus.fa', cs_name='consensus'):
        """ 
        Generates consensus sequence for all the fasta sequences in the
        in the fasta file. Internally runs MUSCLE multiple sequence alignment
        with default parameters. Usually not necessary if only working with one
        strain as all copies of their rRNA tend to be identical(or close to it).
        Parameters
        ----------
        clst: name of output clustal file
        cs_file: fasta file where consensus sequence will be written
        cs_name: name for consensus sequence
        """

        if not cs_file.endswith('.fa'):
            cs_file = cs_file + '.fa'

        if path.isfile(cs_file) and cs_name in [r.id for r in SeqIO.parse(cs_file, 'fasta')]:
            Warning('Consensus sequence for {} already exists in {}.'
                          'File will not be overwritten'.format(cs_name, cs_file))
            return cs_file

        # get multiple sequence alignment with clstalo
        call = ['muscle', '-in', self.rrna_fa, '-out', clst, '-clw']
        print('Running MUSCLE with command: ' + ' '.join(call))
        subprocess.call(call)

        # get consensus sequence from genetics file
        align = AlignIO.read(open(clst), "clustal")
        summary_align = AlignInfo.SummaryInfo(align)  
        consensus = summary_align.dumb_consensus()
        
        with open(cs_file, 'w') as cs_out:
            cs_out.write('>{}_{}\n{}\n'.format(cs_name, self.rtype, str(consensus)))

        return cs_file
