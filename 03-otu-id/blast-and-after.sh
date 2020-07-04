. ~/init/underc

EDNA=$HOME/spe-software/notre-dame/eDNA_pipeline_v3.7
D=/storage/hpc/work/spe1/data/edna
ACDIR=/projects/ac53/underc
PROJDIR=/projects/spe1/edna
METADIR=/projects/spe1/edna/metadata

SCRATCHDIR=/scratch/jenn/underc
MAINDIR=$SCRATCHDIR/sixlakes  # wkg-td
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

# for F in /projects/ac53/underc/analysis/rsrc/databases/ncbi-taxdb/taxdb*; do ln -s `readlink -m $F`; done
cp /storage/hpc/work/spe1/blast/taxdb.tar.gz .; tar xzf taxdb.tar.gz; rm taxdb.tar.gz
ARRAYNUMS=`ls OTUs.subset.* | cut -d '.' -f 3 | sort -n | tr '\n' ',' | sed 's/,$//'`
sbatch -a "$ARRAYNUMS%10" -n 1 -c 1 --mem-per-cpu=2000m -t 04:00:00 $SCRIPTSDIR/remote-blast.sh

###  To repeat until done (because of inevitable errors)  ###

# After BLAST is done, append the results, in order, to the master file. Verify the results, then remove the split directory.
cd $BLASTDIR/split
for N in `ls OTUs.subset.* | cut -d '.' -f 3 | sort -n`; do cat OTUs.subset.$N.blast; done >> ../all-asvs.blast

# See how many BLASTs succeeded, and which failed.
diff <(cat $BLASTDIR/all-asvs.blast | cut -f 1 | sort -u) <(grep '^>' $MAINDIR/asvs.fasta | cut -c 2- | sort) | grep '^> ' | cut -c 3- > failures.txt

# Cleanup: make sure you've copied the results into all-asvs.blast before doing the following!
mkdir -p logs; mv *.out logs;
N=2
mkdir -p blastruns/$N; mv *.blast blastruns/$N; find blastruns/$N -size 0 -exec rm {} \;
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
cut -f 4 all-asvs.blast | sed '1d' | tr ';' '\n' | sort -u | grep -v 'None' | $SCRIPTSDIR/get-taxonomy.py > $RSRCDIR/taxonomy-trees.txt

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

# Create source files that will allow integration script to spread known taxon categories to related taxa that appear elsewhere in the file (maybe should write this functionality in to evaluate-blast.py).
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