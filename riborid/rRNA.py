from os import path, mkdir
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo


class RRNA:
    # TODO: what to do if rrna_fa already exists.
    # TODO: Figure out how to implement oligo_df and check if it messes with gap filling
    def __init__(self, rrna_fa, name, rtype, outdir='rrd', pre=0, oligos_df=None):
        """
        Parameters
        ------------
        name: str
            name of the rRNA, usually 23S, 16S or 5S
        rtype: str
            type of rRNA e.g. 16S, 23S, 5S
        rrna_fa: str, list
            path to input file or list of paths to multiple input rrna file; each fasta file
            must contain sequences from single rtype e.g. '23S', '5S' etc.
        outdir: str, default rrd
            path to the output directory
        pre: int, default 0
            number of base-pairs upstream of rRNA to generate oligos for. Must provide genbank
            files to use this option.
        oligos_df: str, default None
            path to csv file containing sequences for old oligos that you want to reuse
        """

        self.name = name
        self.rtype = rtype
        self.outdir = outdir
        self.pre = pre
        self.rrna_fa = rrna_fa
        self.oname = path.abspath(f'{self.outdir}/{self.name}_{self.rtype}')
        self.consensus = self.get_consensus(clst=self.oname + '_clustal.clst',
                                            cs_file=self.oname + '_consensus',
                                            cs_name=f'{self.name}_{self.rtype}_consensus')

        self.oligos_df = oligos_df

    @property
    def outdir(self):
        """ Get dir for the output files"""
        return self.__outdir

    @outdir.setter
    def outdir(self, outdir):
        if not path.isdir(outdir):
            mkdir(outdir)
        self.__outdir = outdir

        # update full outdir name
        self.__oname = f'{self.outdir}/{self.name}_{self.rtype}_'

    @property
    def rrna_fa(self):
        """ Get the list of gbk_files"""
        return self.__rrna_fa

    @rrna_fa.setter
    def rrna_fa(self, rrna_fa):
        self.__rrna_fa = rrna_fa

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
            Warning(f'Consensus sequence for {cs_name} already exists in {cs_file}.'
                    'File will not be overwritten')
            return cs_file

        # get multiple sequence alignment with clstalo
        call = ['muscle', '-in', self.rrna_fa, '-out', clst, '-clw']
        subprocess.call(call)

        # get consensus sequence from genetics file
        align = AlignIO.read(open(clst), "clustal")
        summary_align = AlignInfo.SummaryInfo(align)  
        consensus = summary_align.dumb_consensus()
        
        with open(cs_file, 'w') as cs_out:
            cs_out.write(f'>{cs_name}\n{str(consensus)}\n')

        return cs_file
