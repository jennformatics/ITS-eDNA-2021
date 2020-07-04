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
