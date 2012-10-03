# load shared functions and variables
. `dirname $0`/run-tests.functions

# velveth binary
$VH $DIR/fai.bin $K -create_binary -shortPaired -fasta.gz $FAI > /dev/null
$VH $DIR/fqi.bin $K -create_binary -shortPaired -fastq.gz $FQI > /dev/null

# check fai
cmp -s $ROADMAP $DIR/fai.bin/Roadmaps
if [ $? -ne 0 ]; then
  problem "$FAI binary mode Roadmap differs to $ROADMAP"
else
  inform "ok"
fi

# check fqi
cmp -s $ROADMAP $DIR/fqi.bin/Roadmaps
if [ $? -ne 0 ]; then
  problem "$FQI binary mode Roadmap differs to $ROADMAP"
else
  inform "ok"
fi

