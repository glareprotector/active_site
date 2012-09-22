#!/bin/bash

BSUBOUTFOLDER="../bsub_outs/"

WIFS=( 0 1 )

WHICHFOLDS=( 0 1 )

WREGS=( 0 1 )

REGS=( 0.0 25.0 100.0 500.0 )

# for all obj fxns, vary reg type and reg constant
# for for nodewise, vary the smooth_f tapering constant, and vary the "true value" which tapers are you go further away from actual active_site

for WIF in ${WIFS[@]}
do
    for WHICHFOLD in ${WHICHFOLDS[@]}
    do
	for REG in ${REGS[@]}
	do

	    if [ $REG == "0.0" ]
	    then
#		echo OOHHHH
		WREGS_TO_USE=( 0 )
	    else
		export WREGS_TO_USE=${WREGS[@]}
	    fi

#	    echo REG$REG
#	    echo WREGS_TO_USE$WREGS_TO_USE
#	    echo WREGS$WREGS

	    for WREG in ${WREGS_TO_USE[@]}
	    do
		OUTFILE=$BSUBOUTFOLDER"WIF-"$WIF"WHICHFOLD-"$WHICHFOLD"WREG-"$WREG"REG-"$REG
		CMD="bsub -R "rusage[mem=3000]" -o $OUTFILE -a openmpi -n 5 -q shared_2h mpirun.lsf ./prog wif $WIF wfld $WHICHFOLD wreg $WREG reg $REG wob 0"
		echo $CMD
		eval "$CMD"
	    done
	done
    done
done

NWCS=( -0.5 -1.0 -10.0 -100.0 )

SFCS=( 15.0 50.0 200.0 )

for WIF in ${WIFS[@]}
do
    for WHICHFOLD in ${WHICHFOLDS[@]}
    do
	for REG in ${REGS[@]}
	do
	    if [ "$REG" == 0.0 ]
	    then
		WREGS_TO_USE=( 0 )
	    else
		WREGS_TO_USE=$WREGS
	    fi

	    for WREG in ${WREGS_TO_USE[@]}
	    do
		for NWC in ${NWCS[@]}
		do
		    for SFC in ${SFCS[@]}
		    do
			OUTFILE=$BSUBOUTFOLDER"WIF-"$WIF"WHICHFOLD-"$WHICHFOLD"WREG-"$WREG"REG-"$REG"NWC-"$NWC"SFC-"$SFC
			CMD="bsub -R "rusage[mem=3000]" -o $OUTFILE -a openmpi -n 5 -q shared_2h mpirun.lsf ./prog wif $WIF wfld $WHICHFOLD wreg $WREG reg $REG nwc $NWC sfc $SFC wob 2"
			echo $CMD
			eval "$CMD"
		    done
		done
	    done
	done
    done
done