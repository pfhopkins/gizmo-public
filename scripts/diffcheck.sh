#!/bin/sh
## This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO ##
for i in *.c *.h
do
echo $i
diff $i $1/$i
done
for i in */*.c */*.h
do
echo $i
diff $i $1/$i
done
for i in */*/*.c */*/*.h
do
echo $i
diff $i $1/$i
done
