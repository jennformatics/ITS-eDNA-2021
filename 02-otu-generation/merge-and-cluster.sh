. ~/init/underc

EDNA=$HOME/spe-software/notre-dame/eDNA_pipeline_v3.7
D=/storage/hpc/work/spe1/data/edna
ACDIR=/projects/ac53/underc
PROJDIR=/projects/spe1/edna
METADIR=/projects/spe1/edna/metadata

SCRATCHDIR=/scratch/jenn/underc
MAINDIR=$SCRATCHDIR/for-realz
RSRCDIR=$MAINDIR/rsrc
SCRIPTSDIR=$MAINDIR/scripts
LOGDIR=$MAINDIR/logs
BLASTDIR=$MAINDIR/blast
DEPRXDATA=$MAINDIR/plantITS_FR_fq

mkdir -p $MAINDIR/blast
mkdir -p $LOGDIR/deprx
cd $MAINDIR
ln -s /home/jenn/scripts/underc scripts
mkdir -p rsrc; cd rsrc
for R in $PROJDIR/rsrc/*.{dmp,fas}; do ln -s $R; done
ln -s $METADIR/sample-names-and-order.txt
# Keep this as a backup somewhere: tar -hczf underc-basic-rsrc.tgz rsrc

# OR:
cd $MAINDIR
rmdir rsrc
tar -xzf $PROJDIR/metadata/underc-basic-rsrc.tgz


#####################################
###  Read merging and clustering  ###
#####################################

# Check out difference between usearch -fastq_mergepairs and...what else do we use? DADA? Pear? BBDuk?

# Merging, filtering, and fq/fa conversion (keeping the filt.fq for later reference).
# This is the crazy, hairy way to process stderr on its way into the log:   >> $MAINLOG 2> >(sed -n '/Totals/,$p')

### ========================= ###

MAINLOG=$LOGDIR/otu-generation.log
DEPRXDIR=/scratch/jenn/underc/deprx-then-trim/plantITS_FR_fq

SEP="\n\n===================================\n"
echo -e "$SEP`date`  Merging with usearch fastq_mergepairs.\n" > $MAINLOG
usearch -fastq_mergepairs $DEPRXDIR/3*_R1.fq -relabel @ -fastqout plantITS_merged.fq >> $MAINLOG 2> >(sed -n '/Totals/,$p')
echo -e "$SEP`date`  Filtering with usearch fastq_filter.\n" >> $MAINLOG
usearch -fastq_filter plantITS_merged.fq -fastq_maxee 0.5 -fastq_maxns 1 -fastqout plantITS_merged.filt.fq >> $MAINLOG 2> >(sed -n '/Reads/,$p')
echo -e "$SEP`date`  Converting to fasta with usearch fastq_filter.\n" >> $MAINLOG
usearch -fastq_filter plantITS_merged.filt.fq -fastaout plantITS_merged.filt.fa >> $MAINLOG 2> >(sed -n '/Reads/,$p')
echo -e "$SEP`date`  Converting to single-line fasta, manually, no thanks to usearch.\n" >> $MAINLOG
$SCRIPTSDIR/unwrap-fasta.sh plantITS_merged.filt.fa > tmp; mv tmp plantITS_merged.filt.fa

# Comparison of *uniquing* (technically dereplication, but not yet clustering)

echo -e "$SEP`date`  Uniquing with YY's script; generates 'cast' file.\n" >> $MAINLOG
perl $EDNA/unique.pl plantITS_merged.filt.fa plantITS_merged.filt.uniq.fa                      # YY
# usearch -fastx_uniques plantITS_merged.filt.fa -sizeout -fastaout plantITS_merged.filt.uniq.fa    # JD

# YY's script creates the cast file (list of reads aggregated) and sometimes picks different (later) representative sequences, but the number and size distros of reads are the same. // Actually, no, on the full set, his misses 17 sequences, including one with size=2280.

# Now do different min-size 2 operations on each, since they had different steps.
echo -e "$SEP`date`  Removing singletons with usearch sortbysize.\n" >> $MAINLOG
usearch -sortbysize plantITS_merged.filt.uniq.fa -fastaout plantITS_merged.filt.uniq.min2.fa -minsize 2 >> $MAINLOG 2> >(grep 'Sorting')         # YY
# usearch -fastx_uniques plantITS_merged.filt.fa -minuniquesize 2 -sizeout -fastaout plantITS_merged.filt.uniq.min2.fa  # JD

# These should now both be counted and sorted, and nearly identical. So, cluster them.
echo -e "$SEP`date`  Clustering with usearch cluster_otus.\n" >> $MAINLOG
usearch -cluster_otus plantITS_merged.filt.uniq.min2.fa -otus otus.fa -uparseout uparse.txt -relabel YYOtu >> $MAINLOG 2> >(tr '\r' '\n' | grep '100.0%') # YY
# Note that my version is going straight into a file called otus.id.fa, since I'm relabeling inside usearch instead of using rename_4.0.pl later.
# usearch -cluster_otus plantITS_merged.filt.uniq.min2.fa -otus otus.fa -uparseout uparse.txt -relabel plantITS_OTU -sizeout -width 0 &   # JD

# Regardless, we lose the size info. Grrr, usearch! -cluster_otus just doesn't include a -sizeout option.
echo -e "$SEP`date`  Adding OTU cluster size back into headers with YY's rename_4.0.pl.\n" >> $MAINLOG
perl $EDNA/rename_4.0.pl otus.fa OTU otus.id.fa   # YY
sed -i 's/^>_/>plantITS_/' otus.id.fa  # Underscore fix since we're prepending the primer name differently.

echo -e "$SEP`date`  Doing something mysterious with usearch_global. Generates 'uc' file. \n" >> $MAINLOG
$EDNA/usearch -usearch_global plantITS_merged.filt.fa -db otus.id.fa -strand plus -id 0.97 -uc otus.id.uc >> $MAINLOG 2> >(tr '\r' '\n' | grep '100.0%') # YY

# Skipping the renaming Perl script on my version, because I did it directly in the -cluster_otus command.
# Honestly, I'm not quite sure what this is doing, since it's not -cluster_otus. Where *is* the 0.97 in cluster_otus, above?
# $EDNA/usearch -usearch_global plantITS_merged.filt.fa -db otus.id.fa -strand plus -id 0.97 -uc otus.id.uc # JD

echo -e "$SEP`date`  Summarizing 'uc' file with YY's one-liner. \n" >> $MAINLOG
UC=otus.id.uc
less -S $UC | perl -e 'while(<>){@s=split; if($s[-1] eq "*"){}else{$hash{$s[-1]} ++;} } foreach my $k (sort {$hash{$b}<=>$hash{$a}} keys %hash){ print "$k\t$hash{$k}\n"; }' > $UC.summary

echo -e "$SEP`date`  Adding abundance (sequence count) info to OTU headers. Again." >> $MAINLOG
perl $EDNA/add_abundance.pl otus.id.fa otus.id.uc.summary NA otus.id.abund.fa  # YY
# perl $EDNA/add_abundance.pl otus.id.fa otus.id.uc.summary NA otus.id.abund.fa           # JD

# And why aren't we just letting them do the OTU table? Oh, because it takes forever.
# usearch -otutab plantITS_merged.fq -otus otus.fa -otutabout otutab.txt -notmatched unmapped.fa -dbmatched otus_with_sizes.fa -sizeout
# or: J=cluster; sbatch -J $J -o $J.log -e $J.log --mem-per-cpu=16000m -n 1 -c 1 -t 24:00:00 <(echo '#!/bin/sh'; echo "usearch -otutab plantITS_merged.fq -otus otus.fa -otutabout otutab.txt -notmatched unmapped.fa -dbmatched otus_with_sizes.fa -sizeout -threads 1")
# And that ended up doing it by RUNS, not SAMPLES! Bleh. Forget it.


#############################################
###  Comparative counts up to this point  ###
#############################################

# Counting stats for all steps
grep -c '^@M' $DEPRXDATA/3*plantITS_R1.fq | sed 's|_plantITS_R1.fq:|\t|'  | tr '-' '\t' | sort -k1,1n -k2,2n | sed 's/\t/-/' | tr ' ' '\t' > $LOGDIR/01-deprimerplexed.counts
grep '^@3' plantITS_merged.fq | cut -d '.' -f 1 | cut -c 2- | $SCRIPTSDIR/prettify.sh > $LOGDIR/02-merged.counts
grep '^@3' plantITS_merged.filt.fq | cut -d '.' -f 1 | cut -c 2- | $SCRIPTSDIR/prettify.sh > $LOGDIR/03-filt.counts
grep '^>3' plantITS_merged.filt.uniq.fa | cut -d '.' -f 1 | cut -c 2- | $SCRIPTSDIR/prettify.sh > $LOGDIR/04-uniq.counts
grep '^>3' plantITS_merged.filt.uniq.min2.fa | cut -d '.' -f 1 | cut -c 2- | $SCRIPTSDIR/prettify.sh > $LOGDIR/05-min2.counts
cp otus.id.uc.summary $LOGDIR/06-otus.counts

gzip plantITS_merged*f{a,q}  # Should be done with all the uncompressed fastx's now, so zip them up.

