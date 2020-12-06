#!/bin/bash

SP=`basename $PWD`

. ~/init/underc $SP


#####################################
###  BLAST lookup and annotation  ###
#####################################

# blastn -task blastn -max_target_seqs 20 -remote -db nt -outfmt '6 qseqid sseqid sacc staxids sscinames scomnames sblastnames pident nident evalue score bitscore'

ALLOTUFASTA=$MAINDIR/$SP.otus.fa
BLASTFILE=$BLASTDIR/$SP.new.otus.blast

mkdir -p $BLASTDIR; cd $BLASTDIR

# Do this only once, or you'll overwrite your results.
echo "Primer_OTU_Size|GenBankGI|GenBankNucl|Taxon Acc|Sci Name|Common Name|Taxon Category|Pct Ident|Len Ident|e-Value|Score|Bitscore" | tr '|' '\t' > $BLASTFILE

# First, see if some of the OTUs have previously been BLASTed and evaluated.

OTUDBDIR=/scratch/jenn/underc/universal-otus
DBFILE=$OTUDBDIR/universal-otu-db.txt
OCCFILE=$OTUDBDIR/occurrences-db.txt

ALLOTUFASTA=$MAINDIR/$SP.otus.fa
NEWOTUFASTA=$BLASTDIR/$SP.new.otus.fa
BLASTFILE=$BLASTDIR/$SP.new.otus.blast
BLASTEVALFILE=$BLASTDIR/$SP.new.otus.blast.evaluated
SEQFILENAME=$SP.otus.seq

sed -n '2~2p' $ALLOTUFASTA > $SEQFILENAME
grep -F -wf <(cut -f 29 $DBFILE) $SEQFILENAME > already-in-db.seq
grep -F -wvf <(cut -f 29 $DBFILE) $SEQFILENAME > not-in-db.seq

# Put a second record in OCCFILE for the ones already in DBFILE.
grep -Fwf already-in-db.seq -B 1 --no-group-separator $ALLOTUFASTA \
    | tr '\n' '|' | sed 's/|>/\n/g' | sed 's/>//' | tr '|' '\t' | sed 's/\t$//' \
    | sed "s|^|$ALLOTUFASTA\t|" >> $OCCFILE   # Be careful not to do this twice; if you do, clean up with sort -u.
echo >> $OCCFILE  # For some reason there's no trailing newline.

# Create the FASTA of new sequences for BLASTing
grep -B 1 --no-group-separator -Fwf not-in-db.seq $ALLOTUFASTA > $NEWOTUFASTA

mkdir -p $BLASTDIR/split; cd $BLASTDIR/split
split -l 100 --numeric-suffixes=1 $NEWOTUFASTA OTUs.subset.
for F in OTUs.subset.*; do NEWNAME=`echo $F | sed 's/\.0*/./g'`; if [ "$F" != "$NEWNAME" ]; then mv $F $NEWNAME; fi; done

# BLAST the new OTUs. Will generate OTUs.subset.XX.blast.
cp /storage/hpc/work/spe1/blast/taxdb.tar.gz .; tar xzf taxdb.tar.gz; rm taxdb.tar.gz
ARRAYNUMS=`ls OTUs.subset.* | cut -d '.' -f 3 | sort -n | tr '\n' ',' | sed 's/,$//'`
sbatch -a "$ARRAYNUMS%10" -n 1 -c 1 --mem-per-cpu=2000m -t 04:00:00 $SCRIPTSDIR/remote-blast.sh

###  If partial repeat is needed due to errors, repeat starts here  ###

# After BLAST is done, append the results, in order, to the master file. Verify the results, then remove the split directory.
cd $BLASTDIR/split   # In case you went somewhere else while waiting for the BLAST.
for N in `ls OTUs.subset.* | cut -d '.' -f 3 | sort -n`; do cat OTUs.subset.$N.blast; done >> $BLASTFILE

# See how many BLASTs succeeded, and which failed.
diff <(cat $BLASTFILE | cut -f 1 | sort -u) <(grep '^>' $NEWOTUFASTA | cut -c 2- | sort) | grep '^> ' | cut -c 3- > failures.txt

# Cleanup: make sure you've copied the results into $BLASTFILE before doing the following!
N=1  # Increment with every rerun.
mkdir -p logs; mv *.out logs;
mkdir -p blastruns/$N; mv *.blast blastruns/$N; find blastruns/$N -size 0 -exec rm {} \;
rm *.subset.*

# IF THERE WERE FAILURES: Make a new FASTA with the failures, split it, and run again.
grep --no-group-separator -w -A 1 -f failures.txt $ALLOTUFASTA > tmp.reblast.fasta
split -l 200 --numeric-suffixes=1 tmp.reblast.fasta OTUs.subset.
for F in OTUs.subset.*; do NEWNAME=`echo $F | sed 's/\.0*/./g'`; if [ "$F" != "$NEWNAME" ]; then mv $F $NEWNAME; fi; done
ARRAYNUMS=`ls *.subset.* | cut -d '.' -f 3 | sort -n | tr '\n' ',' | sed 's/,$//'`
sbatch -a "$ARRAYNUMS%5" -n 1 -c 1 --mem-per-cpu=4000m -t 04:00:00 ~/scripts/underc/remote-blast.sh

# Repeat as necessary until only unrecoverable errors show up in the *.out files.

###  End repeat  ###

cd $MAINDIR
wc -l $RSRCDIR/taxonomy-trees.txt
cut -f 4 $BLASTFILE | sed '1d' | tr ';' '\n' | sort -u | grep -v 'None' | $SCRIPTSDIR/get-taxonomy.py > local-taxonomy-trees.txt
cat local-taxonomy-trees.txt >> $RSRCDIR/taxonomy-trees.txt; sort -u $RSRCDIR/taxonomy-trees.txt > foo; mv foo $RSRCDIR/taxonomy-trees.txt
wc -l $RSRCDIR/taxonomy-trees.txt

cd $BLASTDIR
$SCRIPTSDIR/evaluate-blast.py -b $BLASTFILE -o $ALLOTUFASTA -e $BLASTFILE.evaluated   # Generates $BLASTFILE.evaluated

# evaluate-blast.py    (jd)
#  * Collects summary data for the up to 20 returned BLAST hits for each OTU
#  * Splits multiple-result BLAST lines (with semicolon-separated taxids) into separate items so they can be annotated properly
#  * Checks for ties in the top hits (taxon ambiguity)
#  * Writes an ambiguity report listing all the tied BLAST hits for each ambiguous OTU
#  * Distills BLAST stats into a few summary flags (like bitscore < 100)
#  * Categorizes each OTU as pass, near-pass, ambiguous, or failing based on combinations of flags
#  * Calculates various summary stats across OTUs and organisms
#  * Adds full tree of canonical taxon levels (K/P/C/O/F/G/S) from NCBI taxdmp

# Now put those new evaluated BLAST results in the universal database.
OTUDB=/scratch/jenn/underc/universal-otus
DBFILE=$OTUDB/universal-otu-db.txt
OCCFILE=$OTUDB/occurrences-db.txt
ALLOTUFASTA=`readlink -m $SP.otus.fa`
NEWOTUFASTA=`readlink -m $SP.new.otus.fa`
BLASTFILE=`readlink -m $BLASTDIR/$SP.new.otus.blast`
BLASTEVALFILE=`readlink -m $BLASTFILE.evaluated`

# Put NEW OTUs ONLY in the database.
# Put OTU and its ID on one line, then prepend OTU file source path, then the evaluated BLAST data which now includes the sequence.
cat $NEWOTUFASTA | tr '\n' '|' | sed 's/|>/\n/g' | sed 's/|/\t/g' | sed 's/>//g' | sed "s|^|$ALLOTUFASTA\t|" \
    | while read P OTU SEQ; do echo -e "`date '+%Y-%m-%d %H:%m'`\t$P\t`grep -w $OTU $BLASTEVALFILE`"; done | grep -vP '\t$' >> $DBFILE &

# Put ALL OTUs in the occurrences list.
# Same for OCCFILE, though we're using the whole OTU file rather than just the newly-BLASTED ones.
(cat $ALLOTUFASTA | tr '\n' '|' | sed 's/|>/\n/g' | sed 's/>//' | tr '|' '\t' | sed 's/\t$//' | sed "s|^|$ALLOTUFASTA\t|"; echo) >> $OCCFILE

# Now USE the database to generate $SP.prev.otus.blast.evaluated for the previous set.
cd $BLASTDIR
cat already-in-db.seq | while read SEQ; do LOCALOTU=`grep -Fw $SEQ $OCCFILE | grep "$SP" | cut -f 2 | head -n 1`; grep -Fw $SEQ $DBFILE | cut -f 4-29 | sed "s/^/$LOCALOTU\t/"; done > $SP.prev.otus.blast.evaluated &

# $SP.new.otus.blast has a header, from the fresh BLAST we just did, and ...evaluated doesn't, which is perfect.
cat $BLASTEVALFILE $SP.prev.otus.blast.evaluated > $BLASTDIR/$SP.blast.evaluated


######################################################
###  Manual creation of OTU table plus taxon info  ###
######################################################

cd $MAINDIR

# Generate raw data for OTU table by summarizing the very verbose ".uc" file.
cut -f 9,10 otus.id.uc | tr '.' '\t' | cut -f 1,3 | ~/scripts/underc/pretty-count.sh | sed 's/  */\t/g' | awk '{print $2, $3, $1}' > fine-scaled-read-quantities.txt

OTUTABLE=$SP.otu-table.txt

# Rearrange that "long" file into "wide" format.
$HOME/scripts/melt.py 2 1 fine-scaled-read-quantities.txt > tmp; head -n 1 tmp | sed 's/Index/Primer_OTU/' > $OTUTABLE; sed '1d' tmp | sed "s/^_/${PRIMER}_/" | sort -t'_' -k3,3n >> $OTUTABLE; rm tmp

# Create source files that will allow integration script to spread known taxon categories to related taxa that appear elsewhere in the file (should write this functionality in to evaluate-blast.py, along with handling of taxa that aren't actually very similar OTUs even though they share a high-level taxid).
BLASTEVAL=$BLASTDIR/$SP.blast.evaluated
cat $BLASTEVAL | grep 'needs lookup' | cut -f 5 | cut -d ' ' -f 1 | sort -u | while read TAXON; do grep -w "$TAXON" $BLASTEVAL | cut -f 7 | sort -u | grep -v 'needs lookup' | sed "s/^/$TAXON\t/"; done > taxcats.txt
cut -f 1 taxcats.txt | sort | uniq -c | grep -v ' 1 ' | cut -c 9- > multcats.txt
grep -wvf multcats.txt taxcats.txt > foo; mv foo taxcats.txt
cat multcats.txt | sed 's/$/\tmultiple/' >> taxcats.txt

$SCRIPTSDIR/integrate-otus-and-samples.py -s $OTUTABLE -b $BLASTEVAL -p $SP -o $ALLOTUFASTA

cd $MAINDIR
gzip *.fq *.uc *.cast uparse.txt &
rm fine-scaled-read-quantities.txt local-* progress.log multcats.txt taxcats.txt
mkdir intermediates
mv otus.fa otus.id.* ${PRIMER}_* uparse* intermediates
tar -czf intermediates.tgz intermediates && rm -rf intermediates

cp blast/ambiguity-report.txt $BLASTEVAL .
tar czf output-bundle.tgz $SP.* $OTUFASTA ambiguity-report.txt $SP.blast.evaluated logs
rm $MAINDIR/ambiguity-report.txt $MAINDIR/$SP.blast.evaluated 


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

cd $MAINDIR

for SP in amph; do 
(cd $SP; echo $SP;
    PASSING=$SP.4.organisms.passing
    CURATED=$SP.5.orgs.curated
    grep -f $RSRCDIR/autocurate-yes.txt $PASSING > tmp1
    grep -vf $RSRCDIR/autocurate-no.txt $PASSING | grep -wvf $RSRCDIR/taxon-exclusions.txt > tmp2
    cat tmp1 tmp2 | sort -u > $CURATED
    # sed -i "s/Richardson's pondweed/\"Richardson's pondweed\"/" $CURATED); done

rm tmp1 tmp2
