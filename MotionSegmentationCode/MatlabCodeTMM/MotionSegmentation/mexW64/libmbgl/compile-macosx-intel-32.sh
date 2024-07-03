#!/bin/bash -e

BOOST_DIR=/opt/homebrew/Cellar/boost/1.85.0/include
YASMIC_DIR=.

source ccfiles.sh
OFILES=`echo ${CCFILES} | sed -e 's/\.cc/\.o/g'`

CFLAGS="-O2 -c -I${BOOST_DIR} -I${YASMIC_DIR}"

function echocmd {
	echo $@
	$@
}

for file in ${CCFILES}; do
	echocmd g++ $CFLAGS $file
done

echocmd ar rc libmbgl-macosx-intel-32.a ${OFILES} 
	
echocmd rm ${OFILES}	
