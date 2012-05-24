# load shared functions and variables
. `dirname $0`/run-tests.functions 

# forget to give 2 files
$VH $DIR/fas $K -shortPaired -fasta.gz -separate $FAL 1> /dev/null 2> /dev/null
if [ $? -eq 0 ]; then
  problem "$VH did not fail when given only one file for -shortPaired -separate"
else
  inform "ok"
fi

# files have diff number of seqs in them
$VH $DIR/fas $K -shortPaired -fasta.gz -separate $FAL /dev/zero 1> /dev/null 2> /dev/null  
if [ $? -eq 0 ]; then
  problem "$VH did not fail when given mis-matched reads files in -shortPaired -separate" 
else
  inform "ok"
fi

