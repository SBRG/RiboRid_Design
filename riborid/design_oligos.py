#!/usr/bin/env python 

from rRNA import RRNA
from experiment import Experiment
from get_auth import get_authentication
import getpass

from Bio import SeqIO
import pandas as pd
from os import path, mkdir


def design_oligos(infile, ftype, name='rrd', pre=0, rrna_type=['16S', '23S'], oligos_fa=None, max_gap=50, max_shift=10,
                  oligo_len=32, mt_thresh=65, mt_err=3, na=100, mg=4, oligoc=150, outdir='rrd_res',
                  oligos_df=None, log_exp=False, idt_calc=False, id_file=None):
    """
    Main function used to run riborid oligo design. This will generate the desired oligos given
    the sequence information from the organism and the experimental conditions.
    """

    # initialize the experiment
    exp = Experiment(max_gap, max_shift, oligo_len, mt_thresh, mt_err, na, mg, oligoc, idt_calc)

    _make_checks(rrna_type, ftype, pre, outdir, idt_calc, id_file)  # check input params

    # generate rRNA from genbank
    if ftype == 'genbank':
        if type(infile) == str:
            infile = [infile]
        rrna_fa = path.join(outdir, 'combined_rRNA.fa')
        for gbk in infile:
            parse_genbank(gbk, pre, rrna_type, rrna_fa)
    # rrna fasta passed
    else:
        if type(infile) == list:
            if len(infile) > 1:
                raise ValueError('More than one rrna fasta passed. Concat all sequences into single fasta file.')
            else:
                rrna_fa = infile[0]
        else:  # str or file type passed for fasta file <- Maybe should double check?
            rrna_fa = infile

    split_names = split_rrna(rrna_fa, rrna_type, outdir)

    rrna_dict = {rt: RRNA(rrna_fa=split_names[rt], name=name,
                          rtype=rt, outdir=outdir, pre=pre,
                          oligos_df=oligos_df) for rt in rrna_type}
    # generate IDT token
    if idt_calc:
        usr = input('IDT username: ')
        pwd = getpass.getpass()
        with open(id_file, 'r') as f:
            client_info = f.readline().strip()
        get_authentication(client_info, usr, pwd)

    # generate oligos for each rRNA type
    for rname, rna_obj in rrna_dict.items():
        # if old oligos are provided, use them
        if oligos_fa:
            exp.align_oligos(rna_obj, oligos_fa)
            rna_obj.oligos_df = exp.find_old_oligos(rna_obj, oligos_fa)

        new_oligos = exp.gapfill(rna_obj)
        if rna_obj.oligos_df:
            rna_obj.oligos_df = pd.concat([rna_obj.oligos_df, new_oligos])
        else:
            rna_obj.oligos_df = new_oligos

    oligos = pd.concat([i.oligos_df for i in rrna_dict.values()]).reset_index()
    oligo_name = path.join(outdir, name + '_oligosdf.csv')

    # write the output files
    oligos.to_csv(oligo_name, index=False)

    if log_exp:
        logfile = path.join(outdir,  name + '_exp.log')
        exp.log_exp(logfile)

    # TODO: Check if you can still write primers if you have old and new mixed
    write_primers(name, oligos, outdir)


def _make_checks(rrna_type, ftype, pre, outdir, idt_calc, id_file):
    """
    run sanity checks to make sure correct params were passed
    """
    if not rrna_type:
        rrna_type = ['16S', '23S', '5S']

    if ftype not in ['genbank', 'fasta']:
        raise ValueError(f"ftype must be either 'genbank' or 'fasta.' {ftype} passed instead.")

    if ftype == 'fasta' and pre != 0:
        raise ValueError("Cannot generate oligos for pre-rRNA with rRNA fasta file. Need to pass genbank instead.")
    if not path.isdir(outdir):
        mkdir(outdir)

    if idt_calc:
        if id_file is None:
            raise ValueError('Must pass id_file for temp check with IDT Analyzer.')
        elif not path.isfile(id_file):
            raise FileNotFoundError(f"{id_file} not found.")

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
                    rt_found = True
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


def write_primers(name, oligos_df, outdir):
    primerfa = path.join(outdir, name + '_primers.fa')
    with open(primerfa, 'w') as pout:
        for rname, grp in oligos_df.groupby('target_rRNA'):
            fa_in = path.join(outdir, rname + '.fa')
            seq = next(SeqIO.parse(fa_in, 'fasta'))
            for idx, row in grp.iterrows():
                pseq = seq.seq[row.rRNA_start: row.rRNA_end + 1]
                pseq = pseq.reverse_complement()
                pout.write(f'>{row.oligo_name}|{rname}\n{pseq}\n')


if __name__ == '__main__':

    # TODO: add force argument on whether to overwrite any of the files in there.
    import argparse
    p = argparse.ArgumentParser(description='Design oligos for rRNA removal using RiboRid protocol.')
    p.add_argument('infile', help='Paths to input files of the target organism. Must have '
                                  'rRNA annotated for genbank file.', nargs='*')
    p.add_argument('--ftype', help='type of input file; fasta or genbank', choices=['genbank', 'fasta'],
                   type=str, required=True)
    p.add_argument('-n', '--name', help='Name of the experiment design. Makes it easier to track multiple'
                   'designs; default rrd', type=str, default='rrd')
    p.add_argument('-p', '--pre', help='Number of bp upstream of rRNA start site to include as oligo targets.'
                   'This deals with some organisms that have pre-rRNA; default 0', type=int, default=0)
    p.add_argument('--rrna_type', help='Type of rRNAs to design oligos against; Default 16S and 23S.',
                   nargs='*', default=['16S', '23S'])
    p.add_argument('--oligos_fa', help='Path to fasta file with sequences of old oligos to be reused.',
                   type=str)
    p.add_argument('--max_gap', help='Maximum gap allowed between oligos; default 50',
                   type=int, default=50)
    p.add_argument('--max_shift', help='Maximum number of bp oligos can be shifted to find the one with'
                   'highest melting temperature; default 10', type=int, default=10)
    p.add_argument('--oligo_len', help='Length of the oligos design. All oligos (old and new) must have'
                   'same oligo length; default 32', type=int, default=32)
    p.add_argument('--mt_thresh', help='Minimum melting temperature allowed for any given oligo; default 65',
                   type=float, default=65)
    p.add_argument('--mt_err', help='Error on melting temp estimator. This value is added to \'mt_thresh\''
                   'to ensure that oligo is always above the provided melting threshold; default 3',
                   type=float, default=3)
    p.add_argument('--na', help='Sodium concentration (uM) in the riboRid reaction; default 100',
                   type=int, default=100)
    p.add_argument('--mg', help='Magnesium concentration (uM) in the riboRid reaction; default 4',
                   type=int, default=4)
    p.add_argument('--oligoc', help='Oligos concentration (nM) in the riboRid reaction; default 150',
                   type=int, default=150)
    p.add_argument('-o', '--outdir', help='Output directory. WARNING: Will overwrite '
                   'files(fix coming soon); default rrd_res', type=str, default='rrd_res')
    p.add_argument('--oligos_df', help='Path to csv file containing previously designed'
                   'oligos for target organism', type=str)
    p.add_argument('--log_exp', help='Writes the experimental setup in a log file when flag is present',
                   action='store_true')
    p.add_argument('--idt_calc', help='Whether to use IDT API to calculate the melting temp instead of local calculator',
                    action='store_true')
    p.add_argument('--id_file', help='Path to file containing the client_id and client_secret for IDT API', type=str)
    
    params = vars(p.parse_args())
    design_oligos(**params)
