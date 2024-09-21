#!/bin/sh
mkdir -p ../shl/squeeze
cat ../inc/squeeze/squeeze.h > ../shl/squeeze/squeeze.h
echo '' >> ../shl/squeeze/squeeze.h
echo '#ifdef squeeze_implementation' >> ../shl/squeeze/squeeze.h
cat ../src/squeeze.c >> ../shl/squeeze/squeeze.h
echo '' >> ../shl/squeeze/squeeze.h
echo '#endif // squeeze_implementation' >> ../shl/squeeze/squeeze.h
echo '' >> ../shl/squeeze/squeeze.h
