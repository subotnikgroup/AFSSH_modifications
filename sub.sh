#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V

SCRATCHDIR=/scratch/ambjain/$JOB_ID

if [[ ! -d “$SCRATCHDIR” ]] 
  then
    mkdir -p $SCRATCHDIR 
fi

cp *.inp $SCRATCHDIR
cp aout $SCRATCHDIR

cd $SCRATCHDIR
./aout

cp output* $SGE_O_WORKDIR
cp fort.* $SGE_O_WORKDIR
cp *.out $SGE_O_WORKDIR
cd ..
rm -r $SCRATCHDIR
￼
