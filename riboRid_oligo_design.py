from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Blast import NCBIXML

import subprocess
import os
import itertools
import pandas as pd

# class riboRid():
#     """ This class is used to track oligo design for sepcific organism 
#     ...
    
#     Attributes
#     ----------
#     rRNA: Biopython.Seq (or it might be a generator, double check). 
#         target rRNA sequence that the oligos will be designed for
#     mt_thresh: int
#         melting temperature threshold thresh celcius. Should be the annealing temp. of
#         reaction, default=65.
#     mt_err: int
#         temperature above the mt_thresh used for cutoff, needed because most melting
#         temp. calculators have error of +/- 3C, default=3.
#     Na: int 
#         concentration of Na in mM, default=100.
#     Mg: int 
#         Concetration of Mg in mM, default=4.
#     oligoc: int
#         concentration of oligos in nM, default=150
    
#     Methods
#     -------
#     TO BE COMPLETED AT THE END
#     """ 
    
#     def __init__(self):
        
def get_rRNA(gbk_files, fasta_out='rRNA_seq'):
    """ Generate fasta files of all rRNA in genbank file
        Parameters
        ----------
        gbk: location of the genbank file, or a list-like object containing locations
        of multiple genbank file
        fasta_out: name of output fasta file
    """
    if type(gbk_files) != list:
        gbk_files = [gbk_files]
    with open(fasta_out + '_23S.fa', 'w') as fa_23, open(fasta_out + '_16S.fa', 'w') as fa_16,\
    open(fasta_out + '_5S.fa', 'w') as fa_5:
        for gbk in gbk_files:
            rRNA_count = 0
            for ref_seq in SeqIO.parse(gbk,'genbank'):
                for feat in [f for f in ref_seq.features if f.type == 'rRNA']:
                    product = feat.qualifiers['product'][0]
                    if product.startswith('23S'):
                        fa_23.write('>{}|{}\n{}\n'.format(feat.qualifiers['locus_tag'][0],
                                                       feat.qualifiers['product'][0],
                                                       str(feat.extract(ref_seq.seq))))
                    elif product.startswith('16S'):
                        fa_16.write('>{}|{}\n{}\n'.format(feat.qualifiers['locus_tag'][0],
                                                       feat.qualifiers['product'][0],
                                                       str(feat.extract(ref_seq.seq))))
                    elif product.startswith('5S'):
                        fa_5.write('>{}|{}\n{}\n'.format(feat.qualifiers['locus_tag'][0],
                                                         feat.qualifiers['product'][0],
                                                       str(feat.extract(ref_seq.seq))))
                    else:
                        print(('No product name found for' + feat.qualifiers['locus_tag'][0]))
                        continue
                    rRNA_count += 1
            if rRNA_count == 0:
                raise ValueError('No rRNA found in {}. Please reannotate with prokka.'.format(gbk))
            else:
                print((str(rRNA_count) + ' rRNA sequences found in ' + gbk))
                

#For each of the fasta file, we need to find a consensus sequence
def get_consensus(fasta, clst='clustal.clst', cs_file='consensus.fa', cs_name='consensus_seq'):
    """ Generates consensus sequence for all the fasta sequences in the
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
    consensus = summary_align.dumb_consensus()
    with open(cs_file, 'a+') as cs_out:
        cs_out.write('>{}\n{}\n'.format(cs_name, str(consensus)))
        
        
def align_oligos(rRNA_fa, oligos_fa, outdir='blast'):
    """ align current oligos library to the rRNA sequence file.
    Parameters:
    ----------
    rRNA: fasta file with consensus rRNA sequence, only one sequence allowed per file
    oligos_fa: fasta containing the library of oligos
    """
    
    #TODO: Figure out why the database is not maade into output dir
    #create blastDB
    oligos_db = os.path.join(outdir,os.path.split(oligos_fa)[-1])
    subprocess.call(['cp', oligos_fa, oligos_db])
    subprocess.call(['makeblastdb', '-dbtype', 'nucl', '-in', oligos_fa, '-out',
                     os.path.join(outdir, oligos_fa)])
    #blast-off!
    blast = ['blastn', '-db', oligos_db, '-query', rRNA_fa, '-word_size', 
             '20', '-outfmt', '5', '-out', os.path.join(outdir,'blast_out.xml')]
    print((' '.join(blast)))
    subprocess.call(blast)
        
def find_old_oligos(bl_res, mt_thresh=65, mt_err=3,
                 Na=100, Mg=4, oligoc=150, max_gap=9):
    """ Find oligos from the library that can be reused. Oligos with length of exact sequence
        match >= 20 (blast word_size) and melting temp > mt_thresh - mt_err can be reused. For 
        oligos with mismatches the melting temperature is calculated for the longest stretch
        of matching subsequence.
        
        Parameters:
        ----------
        bl_res: File containg blast results from aligning oligos to rRNA. Must be XML format
        generated with '-outfmt 6' setting in blast.
        mt_thresh: melting temperature threshold thresh celcius. Should be the annealing temp. of
        reaction, default=65.
        mt_err: temperature above the mt_thresh used for cutoff, needed because most melting
        temp. calculators have error of +/- 3C, default=3.
        Na: concentration of Na in mM, default=100.
        Mg: Concetration of Mg in mM, default=4.
        oligoc: concentration of oligos in nM, default=150
        Returns:
        --------
        old_oligos_df: Pandas dataframe containing aligned oligos and their alignment properties.
    """
    
    #open and parse the blast output XML file
    with open(bl_res, 'r') as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))

    oligo_list = []
    all_aln = [(rec.query, aln) for rec in blast_records for aln in rec.alignments]
    #go through each alignment and remove those with low mt
    for aln_info in all_aln:
        align = aln_info[1]
        hsp = align.hsps[0]
        if hsp.gaps == 0:
            matchseq = Seq(hsp.query)
        else:
            maxidx, maxlen = longest_stretch(hsp.match) #grab idx of longest match strectch
            matchseq = Seq(hsp.query[maxidx: maxidx + maxlen])
        mt_tm = mt.Tm_NN(matchseq.transcribe(),nn_table=mt.R_DNA_NN1, 
                         Na=Na, Mg=Mg, dnac1=oligoc)
        if mt_tm >= mt_thresh - mt_err:
            oligo_list.append([align.title.split(' ')[-1], aln_info[0],  float(hsp.identities)/len(hsp.query),
                     hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, align.length,
                     mt_tm, 'old_oligo'])
    
    old_oligos_df = pd.DataFrame(oligo_list,columns=['oligo_name', 'target_rRNA', 'pident','rRNA_start','rRNA_end', 
						'oligo_start', 'oligo_end', 'oligo_len','melting_temp','set'])
    #filter out overlapping oligos, drop the one with lower mt
#     drop_idx = []
#     for name, grp in old_oligos_df.groupby(['target_rRNA']):
#         idx = 0
#         grp = grp.sort_values(['rRNA_start']).reset_index()
#         while idx < len(grp) - 1:
#             oend  = grp.loc[idx, 'rRNA_end']
#             next_ostart = grp.loc[idx+1, 'rRNA_start']
#             if  oend - next_ostart > 9: #oligos overlap by greater than 9 bp
#                 lower_mt = grp.loc[idx: idx+1, 'melting_temp'].idxmin()
#                 ovl_oligo = old_oligos_df[old_oligos_df.oligo_name == grp.loc[lower_mt,
#                                                                           'oligo_name']].index
#                 drop_idx.extend(ovl_oligo)
#             idx += 1
    drop_idx = []
    kept = True
    for name, grp in old_oligos_df.groupby(['target_rRNA']):
        idx = 0
        grp = grp.sort_values(['rRNA_start']).reset_index()
        while idx < len(grp) - 1:
            if kept: #if oligos wasn't dropped use this to measure distance
                oend  = grp.loc[idx, 'rRNA_end']
            next_oend = grp.loc[idx+1, 'rRNA_end']
            apart = next_oend - oend
            if apart <= max_gap: #ends of oligos are closer than max gap
                lower_mt = grp.loc[idx: idx+1, 'melting_temp'].idxmin()
    #             print lower_mt
                ovl_oligo = old_oligos_df[old_oligos_df.oligo_name == grp.loc[lower_mt,
                                                                          'oligo_name']].index
                drop_idx.extend(ovl_oligo)
                kept = False
            else: 
                kept = True
            idx += 1
    return old_oligos_df.drop(drop_idx)

def longest_stretch(matches): #TODO: probably should make this a private method
    """
    Find the longest stretch of repeated '|' in a string. '|' represents
    match in a BLAST outcome
    Parameters:
    -----------
    matches: output from hsp.match, where '|' represents a match beteween sbj. and
    query and ' ' represents a mismatch
    returns 
    -------
    maxidx: idx where the longest stretch of '|' starts
    maxlen: length of the longest stretch of '|'
    """
    
    idx = 0
    maxidx, maxlen = 0, 0
    for _, group in itertools.groupby(matches):
        if _ != '|':
            idx += grouplen
            continue
        grouplen = sum(1 for match in group)
        if grouplen > maxlen:
            maxidx, maxlen = idx, grouplen
        idx += grouplen
    return maxidx, maxlen

def gapfill(rRNA_fa, old_oligos_df=None,  max_gap=9, oligo_len=32, mt_thresh=65, mt_err=3, Na=100,
           Mg=4, oligoc=150):
    """
    Fill the gaps between the old oligos and design new oligos that to fill it in.
    Parameters:
    ----------
    old_oligos_df: DataFrame containing old oligos alignment information. Generated from 
    find_old_oligos. If no DataFrame is passed, only new oligos will be used. Default=None.
    rRNA_fa: rRNA fasta to align to. Only the first sequence in the file is considered
    max_gap: allowed gap between oligos, default=9.
    oligo_len: length of newly designed oligos, default=32.
    mt_thresh: melting temperature threshold celcius. Should be the annealing temp. of
    reaction, default=65.
    mt_err: temperature above the mt_thresh used for cutoff, needed because most melting
    temp. calculators have error of +/- 2.7C, default=3.
    Na: concentration of Na in mM, default=100.
    Mg: Concetration of Mg in mM, default=4.
    oligoc: concentration of oligos in nM, default=150
    Returns:
    --------
    new_oligos_df:Pandas dataframe containing aligned oligos and their alignment properties.
    """
    
    all_rRNA = {seq.id: seq for seq in SeqIO.parse(rRNA_fa, 'fasta')}
    all_gaps = {}
    
    if old_oligos_df is None:
        gaps = [[(0, len(all_rRNA[rRNA].seq))] for rRNA in all_rRNA]
        all_gaps.update(dict(zip([s for s in all_rRNA], gaps)))
    else:
        for seq_id in all_rRNA:
            rRNA_seq = all_rRNA[seq_id].seq
            odf = old_oligos_df[old_oligos_df.target_rRNA == seq_id]
            gaps = find_gaps(odf, rRNA_seq, max_gap=max_gap)
            all_gaps.update({seq_id:gaps})
    new_oligos = []
    oligo_n = 0 #track number of new oligos added
    for seq_id, gaps in list(all_gaps.items()):
        if len(gaps) == 0: #no gaps to fill
            continue
        rRNA_seq = all_rRNA[seq_id].seq
        for gap in sorted(gaps):
            gap_str = (rRNA_seq[gap[0]: gap[1]]) # get the sequence of gap
            gap_len = gap[1] - gap[0]
            gpos = gap[0]
            
            itr = max(max_gap, oligo_len)
            while gpos + itr <= gap[1]:
                if gpos + max_gap + oligo_len < gap[1]:
                    temp= gpos
                    gpos = gpos + max_gap
                #adding oligo will lead to overlap with existing one
                if gpos + oligo_len >= gap_len + gap[0]:
                    gap_mid = round((gap[1] + gpos) / 2)
                    gpos = int(gap_mid - round(float(oligo_len)/2))
                    
                oseq = rRNA_seq[gpos: gpos + oligo_len]
                mt_tm = mt.Tm_NN(oseq.transcribe(),nn_table=mt.R_DNA_NN1, 
                         Na=Na, Mg=Mg, dnac1=oligoc)

                #if mt_tm is too low, shift the frame to right and try again
                if mt_tm < mt_thresh:
                    gpos, mt_tm = __find_maxtm(gpos, rRNA_seq, oseq, oligo_len, max_shift=10)
                    
                name = 'New_Oligo_' + str(oligo_n)
                oligo_n += 1
                new_oligos.append([name, seq_id, 1.0, gpos, gpos + oligo_len, 
                                   1, oligo_len,oligo_len, mt_tm, 'new_oligo'])
                gpos = gpos + oligo_len
   
    return pd.DataFrame(new_oligos,columns=['oligo_name', 'target_rRNA','pident', 'rRNA_start',
                                            'rRNA_end', 'oligo_start', 'oligo_end', 'oligo_len', 
                                            'melting_temp','set'])

#TODO: Make a class that holds all this so you dont have to pass everything down. BUT For now do the dirty way. WILL FIX LATER.
def __find_maxtm(gpos, rRNA_seq, oseq, oligo_len, max_shift=10, Na=100,
           Mg=4, oligoc=150):
    """
    scans the region to find the location with highest melting temp.
    Parameters:
    -----------
    max_shift: maximum number of shifts allowed
    gpos: original starting position of the oligo
    
    Returns:
    --------
    spos: starting position of the region with the highest melting temp
    max_tm: melting temp of oligo starting at the spos
    """
    
    max_tm = 0
    shift_n = 0
    while shift_n < max_shift:
        gpos = gpos  - 1 #only shift "left", shifting right would increase gap length
        oseq = rRNA_seq[gpos: gpos + oligo_len]
        mt_tm = mt.Tm_NN(oseq.transcribe(),nn_table=mt.R_DNA_NN1, 
             Na=Na, Mg=Mg, dnac1=oligoc)
        if max(max_tm, mt_tm) == mt_tm:
            max_tm = mt_tm
            spos = gpos
        shift_n += 1
    return spos, max_tm

def find_gaps(oligo_df, rRNA_seq, max_gap=9): #dont make private, useful for troubleshooting
    """
    Find all the gaps between oligos in oligo_df that are greater than max_gap.
    Parameters:
    -----------
    oligo_df: pandas dataframe with oligo info. generated from find_old_oligos or gapfill method.
    rRNA_seq: Bio Seq object of the rRNA sequence.
    max_gap: maximum allowed gap between oligos.
    
    Returns:
    --------
    gaps: list of tuples containing the starting and ending position of gaps.
    """
    
    gap_start = 0
    gaps = []
    for i, data in oligo_df.sort_values(['rRNA_start']).iterrows():
        gap_end = data.rRNA_start
        if gap_end - gap_start > max_gap:
            gaps.append((gap_start, gap_end))
        gap_start = data.rRNA_end
    #calc gap from last oligo to end of rRNA
    gap_end = len(rRNA_seq)
    if gap_end - gap_start > max_gap:
        gaps.append((gap_start, gap_end))
    return gaps

def fasta_out(fasta_name, rRNA_fa, oligo_df):
    """
    Write fasta file for the NEW oligos in the df. Do NOT use this to write the fasta file for
    old oligos.
    Parameters:
    -----------
    fasta_name: name of the fasta file.
    rRNA_fa: fasta file with rRNA_seq.
    oligo_df: pandas dataframe with oligo information; generated from gapfill method.
    """
    
    if fasta_name[-3:] != '.fa':
        fast_name += '.fa'
    
    seq_iter = list(SeqIO.parse(rRNA_fa, "fasta"))
    with open(fasta_name, 'w') as fa:
        for ref_seq in seq_iter:
            df = oligo_df[oligo_df['target_rRNA'] == ref_seq.id]
            for idx, data in df.iterrows():
                end = data.rRNA_start + data.oligo_len
                fa.write('>{}\n{}\n'.format(data.oligo_name, ref_seq.seq[data.rRNA_start:end]))

if __name__ == '__main__':
    print('Writing into main: Feature is under construction.')
    
        
