# load shared functions and variables
. `dirname $0`/run-tests.functions

# velveth ascii
$VH $DIR/fai $K -shortPaired -fasta.gz $FAI > /dev/null
$VH $DIR/fqi $K -shortPaired -fastq.gz $FQI > /dev/null

# check fai seqs
cmp -s $SEQ $DIR/fai/Sequences
if [ $? -ne 0 ]; then
  problem "$FAI Sequences different to $SEQ"
else
  inform "ok"
fi

# check fqi seqs
cmp -s $SEQ $DIR/fqi/Sequences
if [ $? -ne 0 ]; then
  problem "$FQI Sequences different to $SEQ"
else
  inform "ok"
fi

# check fai roadmap
cmp -s $ROADMAP $DIR/fai/Roadmaps
if [ $? -ne 0 ]; then
  problem "$FAI Roadmap differs to $ROADMAP"
else
  inform "ok"
fi
  
# check fqi roadmap
cmp -s $ROADMAP $DIR/fqi/Roadmaps
if [ $? -ne 0 ]; then
  problem "$FQI Roadmap differs to $ROADMAP"
else
  inform "ok"
fi
