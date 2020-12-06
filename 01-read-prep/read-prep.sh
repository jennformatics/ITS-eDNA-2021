SP=sixlakes  # subproject name

. ~/init/underc $SP

# EDNA=$SOFTWARE/notre-dame/eDNA_pipeline_v3.7
# D=/storage/hpc/work/spe1/data/edna
# PROJDIR=/projects/spe1/edna
# METADIR=/projects/spe1/edna/metadata
# 
# SCRATCHDIR=/scratch/jenn/underc
# MAINDIR=$SCRATCHDIR/sixlakes
# RSRCDIR=$MAINDIR/rsrc
# SCRIPTSDIR=$MAINDIR/scripts
# LOGDIR=$MAINDIR/logs
# BLASTDIR=$MAINDIR/blast
# DEPRXDATA=$MAINDIR/plantITS_FR_fq

# mkdir -p $MAINDIR/blast
mkdir -p $LOGDIR/deprx
cd $MAINDIR
ln -s /home/jenn/scripts/underc scripts
mkdir -p rsrc; cd rsrc
for R in $PROJDIR/rsrc/*.{dmp,fas}; do ln -s $R; done
ln -s $METADIR/sample-names-and-order.txt
# Keep this as a backup somewhere: tar -hczf $PROJDIR/metadata/underc-basic-rsrc.tgz rsrc
# To retrieve: cd $MAINDIR; rmdir rsrc; tar -xzf $PROJDIR/metadata/underc-basic-rsrc.tgz


#####################
###  Trimmomatic  ###
#####################

# Trim raw Illumina read files.
# If not already created: cat $SOFTWARE/Trimmomatic-0.38/adapters/NexteraPE-PE.fa  ~/spe-software/notre-dame/eDNA_pipeline_v3.7/primers/MiSeq.adapter.fas > ~/scripts/underc/all-adapters-ever.fa

mkdir -p $MAINDIR/{37,38}; cd $MAINDIR

for RUNNUM in 37 38; do cd $MAINDIR/$RUNNUM;
    cp -va $D/miseq_run_$RUNNUM/$RUNNUM*.fastq.gz .;
    ARRAYNUMS=`ls *.fastq.gz | cut -d '_' -f 1 | sort -u | sed "s/^$RUNNUM-//" | sort -n | tr '\n' ',' | sed 's/,$/\n/'`;
    sbatch -a "${ARRAYNUMS}%20" -n 1 -c 1 --mem-per-cpu=2000m -t 01:00:00 ~/scripts/underc/just-run-trimmomatic.sh;
    cd ..
done

mkdir trimmed; mv 37 38 trimmed

# Generate four-file stats for each original Illumina pair.
cd trimmed
for RUNNUM in 37 38; do (for S in {1..50}; do
    SN=$RUNNUM-$S; echo -ne "$SN";
    for F in 1.pe.fq 2.pe.fq 1.se.fq 2.se.fq; do echo -ne "\t`grep -c '^@' $RUNNUM/$SN.$F`"; done
    echo; done); done > logs/trimmed-stats.log &


###########################
###  Primer separation  ###
###########################

# Separate primers, which I'm calling "deprimerplexing", or "deprx", by analogy to "demux".

# Deprimerplex runs 32 and 33 with the original primers, and 37-38 with the later primers.
cd $MAINDIR/trimmed
LOGDIR=$PWD/logs; mkdir -p $LOGDIR/deprx plant other
for RUNNUM in 37 38; do
  ls $RUNNUM/3*.1.pe.fq | while read F; do
    # F is the forward-reads file path, and R is the corresponding reverse-reads file path
    R=`echo $F | sed 's/\.1\./.2./'`  # Corresponding reverse file name
    SAMPLEID=`basename $F | cut -d '.' -f 1`;
    RUN=`echo $SAMPLEID | cut -d '-' -f 1`
    if [[ "$RUN" == "32" || "$RUN" == "33" ]]; then PRIMERFILE=$RSRCDIR/primer_degen.fas;
        else PRIMERFILE=$RSRCDIR/37-38-primer_degen.fas; fi
    perl $EDNA/Demultiplex_primer_v1.3.pl $F $R $PRIMERFILE $SAMPLEID $LOGDIR/deprx/$SAMPLEID.deprimerplex.log;
    mv *plantITS* plant; mv *.* other;
  done
done

mv other/trimmed-stats.log .

cd $MAINDIR/plant
ls *.F.* | while read F; do NEWNAME=`echo $F | sed 's/\.F\./_R1./'`; mv $F $NEWNAME; done
ls *.R.* | while read R; do NEWNAME=`echo $R | sed 's/\.R\./_R2./'`; mv $R $NEWNAME; done
for RUNNUM in 37 38; do (for S in {1..50}; do
    SN=$RUNNUM-$S; echo -ne "$SN";
    for F in _plantITS_R1.fq _plantITS_R2.fq; do echo -ne "\t`grep -c '^@' $SN$F`"; done
    echo; done); done > ../logs/plant-deprx-stats.log

mv plant $MAINDIR/trim-deprx
cd $MAINDIR
mkdir data
mv trimmed trim-deprx data

# Other renaming bits, if needed
ls *.1.* | while read F; do NEWNAME=`echo $F | sed 's/\.1\./plantITS_R1./'`; mv $F $NEWNAME; done
ls *.2.* | while read R; do NEWNAME=`echo $R | sed 's/\.2\./plantITS_R2./'`; mv $R $NEWNAME; done

ls 3* | while read R; do NEWNAME=`echo $R | sed 's/plantITS_R/_plantITS_R/'`; mv $R $NEWNAME; done

