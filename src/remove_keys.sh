BIN_FOLDER=../test8/

WRAPS=( bwW caW cbW ccW cdW ahW qW pW mW nW mfdW bhW buW brW ceW btW)

for WRAP in ${WRAPS[@]}
do
    FILES=$BIN_FOLDER*$WRAP*
    CMD="rm "$FILES
    echo $CMD
    eval "$CMD"
done