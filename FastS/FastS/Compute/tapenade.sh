#!/bin/sh
#echo $1
#echo $2
#echo $3
#echo $4
tapenade -tangent  -inputlanguage F77 -O $4 -output $3 -head $2  $1
#tapenade -tangent -vars $3 -outvars $2 -inputlanguage F77 -output $4 -O $5  $1
#tapenade -tangent -inputlanguage F77 -head "spsource_sa(drodm)/(rop)" spsource_SA.f
