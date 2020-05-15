#!/usr/bin/env python 

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
from rRNA import RRNA
from experiment import Experiment
import pandas as pd

def design_oligos(gbk, name, pre, oligos_fa, max_gap, max_shift, oligo_len, mt_thresh,
                 mt_err, Na, Mg, oligoc, outdir, oligo_df, rRNA_fa):

    exp = Experiment(max_gap, max_shift, oligo_len, mt_thresh, mt_err, Na, Mg, oligoc)
    print(exp.Na, exp.Mg, exp.oligoc)
    #TODO: this should not be hardcoded
    #might be ok for now since bacteria mostly have these
    #TODO: parse oligos_df and assign it to rRNA specific
    r23 = RRNA(name, gbk, rtype='23S', rRNA_fa=rRNA_fa, outdir=outdir, pre=pre, oligos_df=oligo_df)
    r16 = RRNA(name, gbk, rtype='16S', rRNA_fa=rRNA_fa, outdir=outdir, pre=pre, oligos_df=oligo_df)  
    r5 =  RRNA(name, gbk, rtype='5S', rRNA_fa=rRNA_fa, outdir=outdir, pre=pre, oligos_df=oligo_df)

    for r in [r23, r16, r5]:
        if oligos_fa:
            exp.align_oligos(r, oligos_fa)
            r.oligos_df = exp.find_old_oligos(r, oligos_fa)  
        new_oligos = exp.gapfill(r)
        if r.oligos_df:
            r.oligos_df = pd.concat([r.oligos_df, new_oligos])
        else:
            r.oligos_df = new_oligos

    oligos = pd.concat([r23.oligos_df, r16.oligos_df, r5.oligos_df]).reset_index()
    oligos.to_csv('test.csv')

if __name__== '__main__':
    #TODO: accept rRNA file instead of gbk as positional argument
    #TODO: add force argument on whether to overwrite any of the files in there.
    import argparse
    p = argparse.ArgumentParser(description='Design oligos for rRNA removal using riboRid protocol.')
    p.add_argument('gbk', help='Paths to genbank files of the target organism. Must have rRNA annotated.',
                  nargs='*')
    p.add_argument('-n', '--name', help='Name of the experiment design. Makes it easier to track multiple'
                   'designs; default rrd', type=str, default='rrd')
    p.add_argument('-p', '--pre', help='Number of bp upstream of rRNA start site to include as oligo targets.'
                   'This deals with some organisms that have pre-rRNA; default 0', type=int, default=0)
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
    p.add_argument('--Na', help='Sodium concentration (uM) in the riboRid reaction; default 100',
                   type=int, default=100)
    p.add_argument('--Mg', help='Magnesium concentration (uM) in the riboRid reaction; default 4',
                   type=int, default=4)
    p.add_argument('--oligoc', help='Oligos concentration (nM) in the riboRid reaction; default 150',
                   type=int, default=150)
    p.add_argument('-o', '--outdir', help='Output directory. WARNING: Will overwrite'
                   'files(fix coming soon); default rrd',
                    type=str, default='rrd')
    p.add_argument('--oligo_df', help='Path to csv file containing previously designed'
                   'oligos for target organism', type=str)
    p.add_argument('--rRNA_fa', help='Path to fasta file containing rRNA sequence. Still under construction,'
                   'doesn\'t do anything yet', type=str)
    
    params = vars(p.parse_args())
    if len(params['gbk']) == 0:
        p.error('Must provide at least one genbank file.')
    design_oligos(**params)