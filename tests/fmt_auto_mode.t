# load shared functions and variables
. `dirname $0`/run-tests.functions

# velveth binary
$VH $DIR/fmtAuto_fai.bin $K -create_binary -shortPaired -fmtAuto $FAI > /dev/null
$VH $DIR/fmtAuto_fqi.bin $K -create_binary -shortPaired -fmtAuto $FQI > /dev/null

# check fai
cmp -s $ROADMAP $DIR/fmtAuto_fai.bin/Roadmaps
if [ $? -ne 0 ]; then
  problem "$FAI binary mode Roadmap differs to $ROADMAP"
else
  inform "ok"
fi

# check fqi
cmp -s $ROADMAP $DIR/fmtAuto_fqi.bin/Roadmaps
if [ $? -ne 0 ]; then
  problem "$FQI binary mode Roadmap differs to $ROADMAP"
else
  inform "ok"
fi

