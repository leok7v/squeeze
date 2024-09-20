#!/bin/sh
mkdir -p ../shl/squeeze
sed '/#endif.*squeeze_h/,$d' squeeze.h > ../shl/squeeze/squeeze.h
echo '#ifdef squeeze_implementation' >> ../shl/squeeze/squeeze.h
cat squeeze.c >> ../shl/squeeze/squeeze.h
echo '' >> ../shl/squeeze/squeeze.h
echo '#endif // squeeze_implementation' >> ../shl/squeeze/squeeze.h
echo '' >> ../shl/squeeze/squeeze.h
echo '#endif // squeeze_h' >> ../shl/squeeze/squeeze.h
