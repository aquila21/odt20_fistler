#!/bin/bash

# Single realization run
# run as "./runOne.sh" or as "./runOne.sh -r"
# See the line near the bottom: runcase "someCase". You should change "someCase"

###############################################################################
echo "the start time is"
date
###############################################################################

inputDir="../input/shearFlow"
#inputDir="../input/particleJet"

###############################################################################

runCase () {

    caseName=$1

    rm -rf "../data/$caseName" > /dev/null 2>&1
    mkdir  "../data/$caseName"
    mkdir  "../data/$caseName/data"
    mkdir  "../data/$caseName/input"
    mkdir  "../data/$caseName/runtime"
        mkdir  "../data/$caseName/particle_data"
    cp     "$inputDir/"*        "../data/$caseName/input/" > /dev/null 2>&1
    cp -r  "$inputDir/restart"* "../data/$caseName/input/" > /dev/null 2>&1

    #--------------------------------------------------------------------------

    echo '*** RUNNING ***'
    ./odt.x $caseName 0         # 0 is the shift (realization # here)

}

###############################################################################

rebuild () {
  echo '*** REBUILDING ***'
  cd ../build
  make -j8
  if [ $? -ne 0 ] ; then
    echo ; echo 'FATAL: error in the build' ; echo
    exit 0
  fi
  echo '*** DONE REBUILDING ***'
  cd ../run
}

###############################################################################

if [ "$1" == "-r" ]; then rebuild; fi

runCase "test"

###############################################################################
echo
echo "the end simulation time is"
date
###############################################################################

exit 0

