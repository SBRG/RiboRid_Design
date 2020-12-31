#!/usr/bin/env python 

"""
THE NEXT THING TO FIGURE OUT HERE IS HOW TO HANDLE AND STORE INDIVIDUAL OLIGOS
"""
from rRNA import RRNA
from experiment import Experiment

from Bio import SeqIO
import pandas as pd
from os import path


def design_oligos(infile, ftype, name, pre=0, rrna_type=None, oligos_fa=None, max_gap=50, max_shift=10, oligo_len=32, mt_thresh=65,
                  mt_err=3, na=100, mg=4, oligoc=150, outdir='rrd', oligo_df=None):
    """
    Main function used to run riborid oligo design. This will generate the desired oligos given
    the sequence information from the organism and the experimental conditions.
    """

    # initialize the experiment
    exp = Experiment(max_gap, max_shift, oligo_len, mt_thresh, mt_err, na, mg, oligoc)

    if not rrna_type:
        rrna_type = ['16S', '23S']

    if type(infile) == str:
        infile = [infile]

    # generate rRNA from genbank
    if ftype == 'genbank':
        rrna_fa = path.join(outdir, 'combined_rRNA.fa')
        for gbk in infile:
            parse_genbank(gbk, pre, rrna_type, rrna_fa)
    else:
        rrna_fa = infile
    split_names = split_rrna(rrna_fa, rrna_type, outdir)

    rrna_dict = {rt: RRNA(rrna_fa=split_names[rt], name=name,
                          rtype=rt, outdir=outdir, pre=pre,
                          oligos_df=oligo_df) for rt in rrna_type}

    # generate oligos for each rRNA type
    for rname, rna_obj in rrna_dict.items():
        if oligos_fa:
            exp.align_oligos(rna_obj, oligos_fa)
            rna_obj.oligos_df = exp.find_old_oligos(rna_obj, oligos_fa)
        new_oligos = exp.gapfill(rna_obj)
        if rna_obj.oligos_df:
            rna_obj.oligos_df = pd.concat([rna_obj.oligos_df, new_oligos])
        else:
            rna_obj.oligos_df = new_oligos

    oligos = pd.concat([i.oligos_df for i in rrna_dict.values]).reset_index()
    oligos.to_csv('test.csv')


def split_rrna(rrna_fa, rrna_type, outdir):
    """ Splits a full rRNA file into individual files """
    split_names = {}
    for rt in rrna_type:
        rt_found = False
        rt_name = path.join(outdir, str(rt) + '_rrna.fa')
        split_names.update({rt: rt_name})
        with open(rt_name, 'w') as rt_out:
            for refseq in SeqIO.parse(rrna_fa, 'fasta'):
                if refseq.id.split('|')[1].startswith(rt):
                    rt_out.write(f'>{refseq.id}\n{refseq.seq}\n')
            if not rt_found:
                raise ValueError(f'No {rt} rRNA found in {rrna_fa}')
    return split_names


def parse_genbank(gbk, pre, rrna_type, rrna_fa):
    rrna_found = False
    for ref_seq in SeqIO.parse(gbk, 'genbank'):
        for feat in [f for f in ref_seq.features if f.type == 'rRNA']:
            product = feat.qualifiers['product'][0]
            # distinguish between 23S, 16S and 5S

            if product.startswith(tuple(rrna_type)):
                write_fasta(feat, ref_seq, pre, rrna_fa)
                rrna_found = True

    if rrna_found:
        return
    raise ValueError(f'No {rrna_type} rRNA found in {gbk}. Please reannotate with prokka.')


def write_fasta(feat, ref_seq, pre, rrna_fa):
    """
    Write the feature sequence into the fasta file
    """
    with open(rrna_fa, 'a+') as fa_out:
        start, end = feat.location.start - pre, feat.location.end
        if feat.strand == -1:
            seq = ref_seq.seq[start: end + pre]
            seq = seq.reverse_complement()
        else:
            seq = ref_seq.seq[start - pre: end]

        fa_out.write('>{}|{}\n{}\n'.format(feat.qualifiers['locus_tag'][0],
                                           feat.qualifiers['product'][0],
                                           str(seq)))


if __name__ == '__main__':

    """
    (rrna_fa, ftype, name, pre=0, rrna_type=None, oligos_fa=None, max_gap=50, max_shift=10, oligo_len=32, mt_thresh=65,
    mt_err=3, na=100, mg=4, oligoc=150, outdir='rrd', oligo_df=None):
    """


    # TODO: add force argument on whether to overwrite any of the files in there.
    import argparse
    p = argparse.ArgumentParser(description='Design oligos for rRNA removal using RiboRid protocol.')
    p.add_argument('rrna_fa', help='Paths to input files of the target organism. Must have rRNA annotated for genbank file.',
                   nargs='*')
    p.add_argument('--ftype', help='type of input file; fasta or genbank', choices=['genbank', 'fasta'],
                   type=str)
    p.add_argument('-n', '--name', help='Name of the experiment design. Makes it easier to track multiple'
                   'designs; default rrd', type=str, default='rrd')
    p.add_argument('-p', '--pre', help='Number of bp upstream of rRNA start site to include as oligo targets.'
                   'This deals with some organisms that have pre-rRNA; default 0', type=int, default=0)
    p.add_argument('--rrna_type', help='Type of rRNAs to design oligos against; Default 16s and 23S.',
                   nargs='*')
    p.add_argument('--oligos_fa', help='Path to fasta file with sequences of old oligos to be reused.',
                  type=str)
    p.add_argument('--max_gap', help='Maximum gap allowed between oligos; default 50',
                  type=int, default=50)
    p.add_argument('--max_shift', help='Maximum number of bp oligos can be shifted to find the one with'
                   'highest melting temperature; default 10', type=int, default=10)
    p.add_argument('--oligo_len', help='Length of the oligos design. All oligos (old and new) must have'
                   'same oligo length; default 32', type=int, default=32)
    p.add_argument('--mt_thresh', help='Minimum melting temperature allowed for any given oligo; default 65',
                  type=int, default=65)
    p.add_argument('--mt_err', help='Error on melting temp estimator. This value is added to \'mt_thresh\''
                   'to ensure that oligo is always above the provided melting threshold; default 3',
                   type=int, default=3)
    p.add_argument('--na', help='Sodium concentration (uM) in the riboRid reaction; default 100',
                   type=int, default=100)
    p.add_argument('--mg', help='Magnesium concentration (uM) in the riboRid reaction; default 4',
                   type=int, default=4)
    p.add_argument('--oligoc', help='Oligos concentration (nM) in the riboRid reaction; default 150',
                   type=int, default=150)
    p.add_argument('-o', '--outdir', help='Output directory. WARNING: Will overwrite'
                   'files(fix coming soon); default rrd', type=str, default='rrd')
    p.add_argument('--oligo_df', help='Path to csv file containing previously designed'
                   'oligos for target organism', type=str)
    
    params = vars(p.parse_args())
    design_oligos(**params)