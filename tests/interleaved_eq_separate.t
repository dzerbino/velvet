# load shared functions and variables
. `dirname $0`/run-tests.functions 

$VH $DIR/fai $K -shortPaired -fasta.gz $FAI                > /dev/null
$VH $DIR/fas $K -shortPaired -fasta.gz -separate $FAL $FAR > /dev/null
$VH $DIR/fqi $K -shortPaired -fastq.gz $FQI                > /dev/null
$VH $DIR/fqs $K -shortPaired -fastq.gz -separate $FQL $FQR > /dev/null

cmp -s $DIR/fai/Sequences $DIR/fas/Sequences
if [ $? -ne 0 ]; then
  problem "$FAI and $FAL+$FAR produced different Sequences file"
else
  inform "ok"
fi

cmp -s $DIR/fai/Sequences $DIR/fqs/Sequences
if [ $? -ne 0 ]; then
  problem "$FAI and $FQL+$FQR produced different Sequences file"
else
  inform "ok"
fi

