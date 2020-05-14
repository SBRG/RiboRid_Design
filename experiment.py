class Experiment():
    """
    This class holds all the information about the experimental setup
    e.g. oligco concentration, melting temp etc.

    """
    def __init__(self, max_gap=9, oligo_len=32, mt_thresh=65, mt_err=3, Na=100, Mg=4, oligoc=150, max_shift=10):

        self.max_gap = max_gap
        self.oligo_len = oligo_len
        self.mt_thresh = mt_thresh
        self.mt_err = mt_err
        self.Na = Na
        self.Mg = Mg
        self.oligoc = oligoc
        self.max_shift = max_shift


    def align_oligos(self, rrna, oligos_fa):
        """ 
        Align current oligos library to the rRNA sequence file.
        Parameters:
        ----------
        rRNA: RRNA object. Contains pertinent information about the rRNA. 
        oligos_fa: fasta containing the library of oligos
        """

        #TODO: Figure out why the database is not made into output dir
        #create blastDB
        oligos_db = os.path.join(rrna.outdir, os.path.split(rrna.fa)[-1])
        subprocess.call(['cp', oligos_fa, oligos_db])
        subprocess.call(['makeblastdb', '-dbtype', 'nucl', '-in', rrna.rRNA_fa, '-out',
                         rrna.name])
        #blast-off!
        blast = ['blastn', '-db', oligos_db, '-query', rRNA_fa, '-word_size', 
                 '20', '-outfmt', '5', '-out', os.path.join(outdir,'blast_out.xml')]
        print((' '.join(blast)))
        subprocess.call(blast)

        #save the blast res to rrna


    def find_old_oligos(self, rrna):
        """ 
        Find oligos from the library that can be reused. Oligos with length of exact sequence
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
        with open(rrna.bl_res, 'r') as result_handle:
            blast_records = list(NCBIXML.parse(result_handle))

        oligo_list = []
        all_aln = [(rec.query, aln) for rec in blast_records for aln in rec.alignments]
        #go through each alignment and remove those with low mt
        for aln_info in all_aln:
            align = aln_info[1]
            hsp = align.hsps[0]

            # oligo matches rRNA with no gaps
            if hsp.gaps == 0: 
                matchseq = Seq(hsp.query)
            #if gaps in match, find the positions of longest stretch of ungapped sub-sequence
            else: 
                maxidx, maxlen = longest_stretch(hsp.match)
                matchseq = Seq(hsp.query[maxidx: maxidx + maxlen])
            mt_tm = mt.Tm_NN(matchseq.transcribe(),nn_table=mt.R_DNA_NN1, 
                             Na=self.Na, Mg=self.Mg, dnac1=self.oligoc)
            if mt_tm >= self.mt_thresh - self.mt_err:
                oligo_list.append([align.title.split(' ')[-1], aln_info[0], float(hsp.identities)/len(hsp.query),hsp.query_start,
             hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, align.length, mt_tm, 'old_oligo'])
            
        old_oligos_df = pd.DataFrame(oligo_list,columns=['oligo_name', 'target_rRNA','pident','rRNA_start',
                                                         'rRNA_end','oligo_start', 'oligo_end',                                                                        'oligo_len','melting_temp','set'])

        #find matching oligos with overlaps, drop one with lower melting temp.
        drop_idx = []
        kept = True

        #TODO: remove groupby if doing for only one rRNA e.g. 23S
        for name, grp in old_oligos_df.groupby(['target_rRNA']):
            print('Name of the group in each BLAST res: ', name)
            idx = 0
            grp = grp.sort_values(['rRNA_start']).reset_index()
            while idx < len(grp) - 1:
                if kept: #if oligos wasn't dropped use this to measure distance
                    oend  = grp.loc[idx, 'rRNA_end']
                next_oend = grp.loc[idx+1, 'rRNA_end']
                apart = next_oend - oend
                if apart <= max_gap: #ends of oligos are closer than max gap
                    lower_mt = grp.loc[idx: idx+1, 'melting_temp'].idxmin()
                    ovl_oligo = old_oligos_df[old_oligos_df.oligo_name == grp.loc[lower_mt,
                                                                              'oligo_name']].index
                    drop_idx.extend(ovl_oligo)
                    kept = False
                else: 
                    kept = True
                idx += 1
        return old_oligos_df.drop(drop_idx)

    def longest_stretch(matches):
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

    def gapfill(self,  rrna):
        """
        Fill the gaps between the old oligos and design new oligos that to fill it in.
        Parameters:
        ----------
        old_oligos_df: DataFrame containing old oligos alignment information. Generated from 
        find_old_oligos. If no DataFrame is passed, only new oligos will be used. Default=None.
        rRNA_fa: rRNA fasta to align to
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

        #TODO: Do we really need to look at all rRNA? Shouldn't this contain one consensus sequence?
        all_rRNA = {seq.id: seq for seq in SeqIO.parse(rrna.rRNA_fa, 'fasta')}
        all_gaps = {}

        if rrna.old_oligos_df is None:
            all_gaps.update(dict(zip([s for s in all_rRNA], [(50, len(seqs)) for seqs in all_rRNA])))
            return all_gaps #TODO: why are we returning all_gaps?
        else:
            for seq_id in all_rRNA:
                rRNA_seq = all_rRNA[seq_id].seq
                odf = rrna.old_oligos_df[rrna.old_oligos_df.target_rRNA == seq_id]
                gaps = find_gaps(odf, rRNA_seq) #TODO fix parameters, might not need to pass both
                all_gaps.update({seq_id:gaps})
        new_oligos = []
        oligo_n = 0 #track number of new oligos added
        for seq_id, gaps in list(all_gaps.items()):
            rRNA_seq = all_rRNA[seq_id].seq
            for gap in sorted(gaps):
                gap_str = (rRNA_seq[gap[0]: gap[1]]) # get the sequence of gap
                gap_len = gap[1] - gap[0]
                gpos = gap[0]

                itr = max(self.max_gap, oligo_len)
                while gpos + itr <= gap[1]:
                    if gpos + self.max_gap + self.oligo_len < gap[1]:
                        temp= gpos
                        gpos = gpos + self.max_gap
                    #adding oligo will lead to overlap with existing one
                    if gpos + self.oligo_len >= self.gap_len + gap[0]:
                        gap_mid = round((gap[1] + gpos) / 2)
                        gpos = int(gap_mid - round(float(self.oligo_len)/2))

                    oseq = rRNA_seq[gpos: gpos + self.oligo_len]
                    mt_tm = mt.Tm_NN(oseq.transcribe(),nn_table=mt.R_DNA_NN1, 
                             Na=self.Na, Mg=self.Mg, dnac1=self.oligoc)

                    #if mt_tm is too low, shift the frame to right and try again
                    #TODO: consider moving both sides, but that might create complications
                    if mt_tm < self.mt_thresh:
                        gpos, mt_tm = __find_maxtm(gpos, rRNA_seq, oseq, self.oligo_len)

                    name = 'New_Oligo_' + str(oligo_n)
                    oligo_n += 1
                    new_oligos.append([name, seq_id, 1.0, gpos, gpos + self.oligo_len, 
                                       1, self.oligo_len, self.oligo_len, mt_tm, 'new_oligo'])
                    gpos = gpos + self.oligo_len

        return pd.DataFrame(new_oligos,columns=['oligo_name', 'target_rRNA','pident', 'rRNA_start',
                                                'rRNA_end', 'oligo_start', 'oligo_end', 'oligo_len', 
                                                'melting_temp','set'])

    def __find_maxtm(self, gpos, rRNA_seq, oseq):
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
        while shift_n < self.max_shift:
            gpos = gpos  - 1 #only shift "left", shifting right would increase gap length
            oseq = rRNA_seq[gpos: gpos + self.oligo_len]
            mt_tm = mt.Tm_NN(oseq.transcribe(),nn_table=mt.R_DNA_NN1, 
                 Na=self.Na, Mg=self.Mg, dnac1=self.oligoc)
            if max(max_tm, mt_tm) == mt_tm:
                max_tm = mt_tm
                spos = gpos
            shift_n += 1
        return spos, max_tm

    def find_gaps(self, oligo_df, rRNA_seq): #dont make private, useful for troubleshooting
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

        #TODO: should gap start be 50? 
        gap_start = 0
        gaps = []
        for i, data in oligo_df.sort_values(['rRNA_start']).iterrows():
            gap_end = data.rRNA_start
            if gap_end - gap_start > self.max_gap:
                gaps.append((gap_start, gap_end))
            gap_start = data.rRNA_end
        #calc gap from last oligo to end of rRNA
        gap_end = len(rRNA_seq)
        if gap_end - gap_start > self.max_gap:
            gaps.append((gap_start, gap_end))
            return gaps