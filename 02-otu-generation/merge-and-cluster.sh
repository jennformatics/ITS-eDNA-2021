SP=`basename $PWD`

. ~/init/underc $SP


#####################################
###  Read merging and clustering  ###
#####################################

# Merging, filtering, and fq/fa conversion (keeping the filt.fq for later reference).
# This is the crazy, hairy way to process stderr on its way into the log:   >> $MAINLOG 2> >(sed -n '/Totals/,$p')

mkdir -p logs; rm $MAINLOG
SEP="\n\n===================================\n"
echo -e "$SEP`date`  Merging with usearch fastq_mergepairs.\n" > $MAINLOG
usearch -fastq_mergepairs $DEPRXDIR/3*_R1.fq -relabel @ -fastqout ${PRIMER}_merged.fq >> $MAINLOG 2> >(sed -n '/Totals/,$p')
echo -e "$SEP`date`  Filtering with usearch fastq_filter.\n" >> $MAINLOG
usearch -fastq_filter ${PRIMER}_merged.fq -fastq_maxee 0.5 -fastq_maxns 1 -fastqout ${PRIMER}_merged.filt.fq >> $MAINLOG 2> >(sed -n '/Reads/,$p')
echo -e "$SEP`date`  Converting to fasta with usearch fastq_filter.\n" >> $MAINLOG
usearch -fastq_filter ${PRIMER}_merged.filt.fq -fastaout ${PRIMER}_merged.filt.fa >> $MAINLOG 2> >(sed -n '/Reads/,$p')
echo -e "$SEP`date`  Converting to single-line fasta, manually, no thanks to usearch.\n" >> $MAINLOG
$SCRIPTSDIR/unwrap-fasta.sh ${PRIMER}_merged.filt.fa > tmp && mv tmp ${PRIMER}_merged.filt.fa
echo -e "$SEP`date`  Uniquing with YY's script; generates 'cast' file.\n" >> $MAINLOG
perl $EDNA/unique.pl ${PRIMER}_merged.filt.fa ${PRIMER}_merged.filt.uniq.fa
echo -e "$SEP`date`  Removing singletons with usearch sortbysize.\n" >> $MAINLOG
usearch -sortbysize ${PRIMER}_merged.filt.uniq.fa -fastaout ${PRIMER}_merged.filt.uniq.min2.fa -minsize 2 >> $MAINLOG 2> >(grep 'Sorting')
echo -e "$SEP`date`  Clustering with usearch cluster_otus.\n" >> $MAINLOG
usearch -cluster_otus ${PRIMER}_merged.filt.uniq.min2.fa -otus otus.fa -uparseout uparse.txt -relabel OTU >> $MAINLOG 2> >(tr '\r' '\n' | grep '100.0%')
echo -e "$SEP`date`  Rewriting headers with YY's rename_4.0.pl.\n" >> $MAINLOG
perl $EDNA/rename_4.0.pl otus.fa ${PRIMER}_OTU otus.id.fa
sed -i "s/^>_/>/" otus.id.fa  # Underscore fix
echo -e "$SEP`date`  Assigning original reads to clusters with usearch_global. Generates 'uc' file. \n" >> $MAINLOG
$EDNA/usearch -usearch_global ${PRIMER}_merged.filt.fa -db otus.id.fa -strand plus -id 0.97 -uc otus.id.uc # JD
echo -e "$SEP`date`  Summarizing 'uc' file with YY's one-liner. \n" >> $MAINLOG
UC=otus.id.uc; less -S $UC | perl -e 'while(<>){@s=split; if($s[-1] eq "*"){}else{$hash{$s[-1]} ++;} } foreach my $k (sort {$hash{$b}<=>$hash{$a}} keys %hash){ print "$k\t$hash{$k}\n"; }' > $UC.summary
echo -e "$SEP`date`  Adding abundance (sequence count) info to OTU headers." >> $MAINLOG
perl $EDNA/add_abundance.pl otus.id.fa otus.id.uc.summary NA otus.id.abund.fa
echo -e "$SEP`date`  Final size sort and unwrapping, generating final OTU FASTA ($SP.otus.fa)." >> $MAINLOG
usearch -sortbysize otus.id.abund.fa -fastaout $SP.otus.fa -minsize 2
$SCRIPTSDIR/unwrap-fasta.sh $SP.otus.fa > foo; mv foo $SP.otus.fa >> $MAINLOG 2>&1



#############################################################
###  VSEARCH VERSION (mutually exclusive with the above)  ###
#############################################################

echo -e "$LOGSEP`date`  Merging with vsearch fastq_mergepairs.\n" > $MAINLOG
for F in $DEPRXDIR/3*_R1.fq; do 
    R=`echo $F | sed 's|R1|R2|'`; LABEL=`basename $F | cut -d '_' -f 1 | tr '-' '_'`;
    (echo; echo -e "Merging pair $F.\n") >> $MAINLOG;
    vsearch --fastq_mergepairs $F --reverse $R --fastqout merged/$LABEL.mrgflt.fq \
            --fastq_maxee 0.5 --fastq_maxns 0 --fastq_allowmergestagger --relabel "$LABEL.";
done >> $MAINLOG 2>&1 &

echo -e "$LOGSEP`date`  Concatenating sample-level merge files into ${PRIMER}.mergfilt.fq.\n" >> $MAINLOG 
cat merged/*.fq > ${PRIMER}.mergfilt.fq

echo -e "$LOGSEP`date`  Converting to fasta with vsearch fastq_filter; discarding sequences < 150 bp.\n" >> $MAINLOG
vsearch --fastq_filter ${PRIMER}.mergfilt.fq --fastaout ${PRIMER}.mergfilt.fa --fasta_width 0 --fastq_minlen 150 >> $MAINLOG 2>&1

echo -e "$LOGSEP`date`  Uniquing with vsearch derep_fulllength." >> $MAINLOG
vsearch --derep_fulllength ${PRIMER}.mergfilt.fa --output ${PRIMER}.mergfilt.uniq.fa --fasta_width 0 --sizeout >> $MAINLOG 2>&1

echo -e "$LOGSEP`date`  Compressing source FQ and FA." >> $MAINLOG
gzip ${PRIMER}.mergfilt.fq ${PRIMER}.mergfilt.fa &

N=8
echo -e "$LOGSEP`date`  Denoising with vsearch cluster_unoise and keeping clusters of at least $N (default 8)." >> $MAINLOG
vsearch --cluster_unoise ${PRIMER}.mergfilt.uniq.fa --centroids asvs.withchimeras.fa --fasta_width 0 --minsize $N --relabel preASV_ --sizein --sizeout >> $MAINLOG 2>&1

# Chimera removal with uchime2 (ratio of 2, not 16); does sorting by sizein first
echo -e "$LOGSEP`date`  Chimera removal with uchime2_denovo." >> $MAINLOG
vsearch --uchime2_denovo asvs.withchimeras.fa --nonchimeras $SP.asvs.fa --fasta_width 0 --sizein --sizeout --relabel ${PRIMER}_ASV_ >> $MAINLOG 2>&1

echo -e "$LOGSEP`date`  Searching for original reads with usearch_global using ASVs as database.\n" >> $MAINLOG
vsearch --usearch_global ${PRIMER}.mergfilt.fa.gz --db $ALLASVFASTA --id 0.97 --otutabout $OTUTABLE >> $MAINLOG 2>&1 &
# Optional at this point: sed -i '1s/#OTU ID/Primer_OTU/' $OTUTABLE
