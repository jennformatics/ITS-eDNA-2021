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


#####################
###  Trimmomatic  ###
#####################

# Oh, crud; did I not trim the first time when I restarted?
# So, trim the already-deprx'd ones (deprx-then-trim),
#   and then trim a copy of the originals and deprx them afterwards (trim-then-deprx).

cat ~/software/Trimmomatic-0.38/adapters/NexteraPE-PE.fa  ~/spe-software/notre-dame/eDNA_pipeline_v3.7/primers/MiSeq.adapter.fas > ~/scripts/underc/all-adapters-ever.fa

mkdir -p /scratch/jenn/underc/ugh-trim/{32,33,37,38}

ARRAYNUMS=`ls *.fastq.gz | cut -d '_' -f 1 | sort -u | sed "s/^$RUNNUM-//" | sort -n | tr '\n' ',' | sed 's/,$/\n/'`;
sbatch -a "${ARRAYNUMS}%20" -n 1 -c 1 --mem-per-cpu=2000m -t 01:00:00 ~/scripts/underc/just-run-trimmomatic.sh;
done

cd /scratch/jenn/underc/ugh-trim/trimmed
for RUNNUM in 32 33 37 38; do (for S in {1..50}; do
    SN=$RUNNUM-$S; echo -ne "$SN";
    for F in 1.pe.fq 2.pe.fq 1.se.fq 2.se.fq; do echo -ne "\t`grep -c '^@' $SN.$F`"; done
    echo; done); done > trimmed-stats.log &

cd /scratch/jenn/underc/ugh-trim/deprx
ls ../plant | while read F; do ln -sf ../plant/$F; done
for RUNNUM in 32 33 37 38; do (for S in {1..50}; do
    SN=$RUNNUM-$S; ~/scripts/underc/just-run-trimmomatic.sh $SN; done); done > logs/trim.log 2>&1

grep -v 'Using Long' logs/trim.log > foo; mv foo logs/trim.log


##########################################################
###  Primer separation: same as in recapitulation.txt  ###
##########################################################

# Either do the below,
#  or link to the previous set: ln -s $SCRATCHDIR/plant/plantITS_R1R2_fq plantITS_FR_fq
#  or unzip the stored copies in $D/deprimerplexed/plant/*.gz into $DEPRXDATA

# ln -s $SCRATCHDIR/plant/plantITS_R1R2_fq plantITS_FR_fq

# Deprimerplex runs 32 and 33 with the original primers, and 37-38 with the later primers.
LOGDIR=./logs
ls 3{2,3,7,8}/3*.1.pe.fq | while read F; do
    # F is the forward-reads file path, and R is the corresponding reverse-reads file path
    R=`echo $F | sed 's/\.1\./.2./'`  # Corresponding reverse file name
    SAMPLEID=`basename $F | cut -d '.' -f 1`;
    RUN=`echo $SAMPLEID | cut -d '-' -f 1`
    if [[ "$RUN" == "32" || "$RUN" == "33" ]]; then PRIMERFILE=$RSRCDIR/primer_degen.fas;
      else PRIMERFILE=$RSRCDIR/37-38-primer_degen.fas; fi
    perl $EDNA/Demultiplex_primer_v1.3.pl $F $R $PRIMERFILE $SAMPLEID $LOGDIR/deprx/$SAMPLEID.deprimerplex.log;
    mv *plantITS* plant; mv *.* other;
done

cd /scratch/jenn/underc/ugh-trim/plant
ls *.F.* | while read F; do NEWNAME=`echo $F | sed 's/\.F\./_R1./'`; mv $F $NEWNAME; done
ls *.R.* | while read R; do NEWNAME=`echo $R | sed 's/\.R\./_R2./'`; mv $R $NEWNAME; done
ls *.1.* | while read F; do NEWNAME=`echo $F | sed 's/\.1\./plantITS_R1./'`; mv $F $NEWNAME; done
ls *.2.* | while read R; do NEWNAME=`echo $R | sed 's/\.2\./plantITS_R2./'`; mv $R $NEWNAME; done

ls 3* | while read R; do NEWNAME=`echo $R | sed 's/plantITS_R/_plantITS_R/'`; mv $R $NEWNAME; done

for RUNNUM in 32 33 37 38; do (for S in {1..50}; do
    SN=$RUNNUM-$S; echo -ne "$SN";
    for F in _plantITS_R1.fq _plantITS_R2.fq; do echo -ne "\t`grep -c '^@' $SN$F`"; done
    echo; done); done > plant-deprx-stats.log &


##################################
###  Read prep and clustering  ###
##################################

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


#####################################
###  BLAST lookup and annotation  ###
#####################################

# blastn -task blastn -max_target_seqs 20 -remote -db nt -outfmt '6 qseqid sseqid sacc staxids sscinames scomnames sblastnames pident nident evalue score bitscore'

cd $MAINDIR
ln -s otus.id.fa asvs.fasta
# or temporarily,  head -n 200 otus.id.abund.fa > asvs.fasta

cat $MAINDIR/asvs.fasta | tr '\n' ' ' | sed 's/>/\n>/g' | awk '{FS=" "} {printf("%s\t%s\n", $1, length($2))}' | sed 's/;size=[0-9]*//g' | sed '1d' | cut -c 2- | sed 's/;//' > $RSRCDIR/asv-sequence-lengths.txt
cd $BLASTDIR

# Do this only once, or you'll overwrite your results.
echo "ASV|GenBankGI|GenBankNucl|Taxon Acc|Sci Name|Common Name|Taxon Category|Pct Ident|Len Ident|e-Value|Score|Bitscore" | tr '|' '\t' > all-asvs.blast

mkdir -p $BLASTDIR/split; cd $BLASTDIR/split
split -l 200 --numeric-suffixes=1 $MAINDIR/asvs.fasta OTUs.subset.
for F in OTUs.subset.*; do NEWNAME=`echo $F | sed 's/\.0*/./g'`; if [ "$F" != "$NEWNAME" ]; then mv $F $NEWNAME; fi; done

for F in /projects/ac53/underc/analysis/rsrc/databases/ncbi-taxdb/taxdb*; do ln -s `readlink -m $F`; done
# cp /storage/hpc/work/spe1/blast/taxdb.tar.gz .; tar xzf taxdb.tar.gz; rm taxdb.tar.gz
ARRAYNUMS=`ls OTUs.subset.* | cut -d '.' -f 3 | sort -n | tr '\n' ',' | sed 's/,$//'`
sbatch -a "$ARRAYNUMS%5" -n 1 -c 1 --mem-per-cpu=1000m -t 01:00:00 $SCRIPTSDIR/remote-blast.sh

###  To repeat until done (because of inevitable errors)  ###

# After BLAST is done, append the results, in order, to the master file. Verify the results, then remove the split directory.
cd $BLASTDIR/split
for N in `ls OTUs.subset.* | cut -d '.' -f 3 | sort -n`; do cat OTUs.subset.$N.blast; done >> ../all-asvs.blast

# See how many BLASTs succeeded, and which failed.
diff <(cat $BLASTDIR/all-asvs.blast | cut -f 1 | sort -u) <(grep '^>' $MAINDIR/asvs.fasta | cut -c 2- | sort) | grep '^> ' | cut -c 3- > failures.txt

# Cleanup: make sure you've copied the results into all-asvs.blast before doing the following!
mkdir -p logs; mv *.out logs;
N=2
mkdir blastruns/$N; mv *.blast blastruns/$N; find blastruns/$N -size 0 -exec rm {} \;
rm *.subset.* *.fasta

# Make a new FASTA with the failures, split it, and run again.
grep --no-group-separator -w -A 1 -f failures.txt $MAINDIR/asvs.fasta > OTUs.fasta
split -l 200 --numeric-suffixes=1 OTUs.fasta OTUs.subset.
for F in OTUs.subset.*; do NEWNAME=`echo $F | sed 's/\.0*/./g'`; if [ "$F" != "$NEWNAME" ]; then mv $F $NEWNAME; fi; done
ARRAYNUMS=`ls *.subset.* | cut -d '.' -f 3 | sort -n | tr '\n' ',' | sed 's/,$//'`
sbatch -a "$ARRAYNUMS%5" -n 1 -c 1 --mem-per-cpu=4000m -t 04:00:00 ~/scripts/underc/remote-blast.sh

# Repeat as necessary until only unrecoverable errors show up in the *.out files.

###  End repeat  ###

cd $MAINDIR
cat otus.id.abund.fa | tr '\n' ' ' | sed 's/>/\n>/g' | awk '{FS=" "} {printf("%s\t%s\n", $1, length($2))}' | sed 's/;size=[0-9]*;*//g' | sed '1d' | cut -c 2- > $RSRCDIR/otu-sequence-lengths.txt
cd $BLASTDIR; rm -rf split; ln -s $RSRCDIR
cut -f 4 all-asvs.blast | sed '1d' | tr ';' '\n' | sort -u | $SCRIPTSDIR/get-taxonomy.py > $RSRCDIR/taxonomy-trees.txt

$SCRIPTSDIR/evaluate-blast.py --rsrcdir ../rsrc --blastdir .   # Generates all-asvs.blast.evaluated

# evaluate-blast.py    (jd)
#  * Collects summary data for the up to 20 returned BLAST hits for each OTU
#  * Splits multiple-result BLAST lines (with semicolon-separated taxids) into separate items so they can be annotated properly
#  * Checks for ties in the top hits (taxon ambiguity)
#  * Writes an ambiguity report listing all the tied BLAST hits for each ambiguous OTU
#  * Distills BLAST stats into a few summary flags (like bitscore < 100)
#  * Categorizes each OTU as pass, near-pass, ambiguous, or failing based on combinations of flags
#  * Calculates various summary stats across OTUs and organisms
#  * Adds full tree of canonical taxon levels (K/P/C/O/F/G/S) from NCBI taxdmp


######################################################
###  Manual creation of OTU table plus taxon info  ###
######################################################

cd $MAINDIR

cat otus.id.abund.fa | tr '\n' ' ' | sed 's/>/\n>/g' | awk '{FS=" "} {printf("%s\t%s\n", $1, length($2))}' | sed 's/;size=[0-9]*;*//g' | sed '1d' | cut -c 2- > $RSRCDIR/otu-sequence-lengths.txt

# Generate raw data for OTU table by comparing the very verbose ".uc" file with the ordered list of samples.
rm progress.log fine-scaled-read-quantities.txt
cut -f 1 $RSRCDIR/sample-names-and-order.txt | sed '1d' | sort -t'-' -k1,1n -k2,2n | \
    while read SAMPLENAME; do echo $SAMPLENAME >> progress.log; \
        cut -f 9,10 otus.id.uc | grep -w "$SAMPLENAME" | tr '.' '\t' | cut -f 1,3 | \
        ~/scripts/pretty-count.sh | sed 's/  */\t/g' | awk '{print $2, $3, $1}' | sort -t'_' -k2,2n | \
        grep -v '\*' >> fine-scaled-read-quantities.txt; done &
tail -f progress.log

# Rearrange that "long" file into "wide" format.
~/scripts/melt.py 2 1 fine-scaled-read-quantities.txt > tmp; head -n 1 tmp | sed 's/Index/Primer_OTU/' > otu-table.txt; sed '1d' tmp | sort -t'_' -k3,3n >> otu-table.txt; rm tmp

# Create source files that will allow integration scrpt to spread known taxon categories to related taxa that appear elsewhere in the file (maybe should write this functionality in to evaluate-blast.py).
BLASTEVAL=blast/AllOTUs.blast.evaluated; cat $BLASTEVAL | grep 'needs lookup' | cut -f 5 | cut -d ' ' -f 1 | sort -u | while read TAXON; do grep -w "$TAXON" $BLASTEVAL | cut -f 7 | sort -u | grep -v 'needs lookup' | sed "s/^/$TAXON\t/"; done > taxcats.txt
cut -f 1 taxcats.txt | sort | uniq -c | grep -v ' 1 ' | cut -c 9- > multcats.txt
grep -wvf multcats.txt taxcats.txt > foo; mv foo taxcats.txt
cat multcats.txt | sed 's/$/\tmultiple/' >> taxcats.txt

$SCRIPTSDIR/integrate-otus-and-samples.py

cd $MAINDIR
gzip *.fq *.uc *.cast uparse.txt &
rm progress.log multcats.txt taxcats.txt

cp blast/ambiguity-report.txt blast/AllOTUs.blast.evaluated .
tar czf output-bundle.tgz cumulative* otu-table.txt otus.id.abund.fa ambiguity-report.txt ambiguity-report.txt AllOTUs.blast.evaluated logs
rm $MAINDIR/ambiguity-report.txt $MAINDIR/AllOTUs.blast.evaluated


#integrate-orgs-and-samples.py   (jd -- combines the BLAST results with the OTU-by-sample info and NCBI taxdmp)
#  * Merges sequencing sample names (like 32-1) with lake sample names (like Red1S)
#  * Matches OTU IDs from the annotated BLAST results with numbers of reads per sample from the clustering output
#  * Counts up total numbers of lakes and samples in which each OTU is present
#  * Counts total p/np/a/f OTUs for each organism, and puts that summary info on each of that organism's OTUs
#  ** Generates "otus" file with the above information
#  * Totals reads for all the OTUs in each primer/organism combo (can leave out certain categories like near-pass at this step, if requested)
#  ** Generates "primers.organisms" file with the above information
#  * Combines all primer records for a given organism, choosing the read counts from the one with the most reads in each sample
#  * Diverts "uncultured eukaryote" and "environmental sample" OTUs to a separate file (taxon-exclusions), probably to be ignored forever
#  * Retains stats on total organism reads, total included reads (if some categories were left out), and total maxed reads (added up from the samples)
#  ** Generates "organisms" file with the above information
#
#cat organisms | awk '{OFS="\t"; FS="\t"} {if($7>0) print $0}' > passing.organisms
#  * Generates "passing.organisms" file from all the organisms that contained passing reads.


###############################################
###  Auto-curation, to the extent possible  ###
###############################################

cat cumulative.4.organisms.passing