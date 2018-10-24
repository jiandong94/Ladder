#!/bin/bash

# auto qsub 
PBS=PBSrun.sh
APP=ladder_d5
IF=IF_U2U4_1

# generate second inputfile
NAME=`echo "$IF" | grep -o "[^_]\+[0-9]"`
IF2=IF_${NAME}_2
cp $IF $IF2
echo $IF2

sed -i "s/ReadPsi.*no/ReadPsi = yes/g" $IF2
sed -i "s/ReadNum.*0/ReadNum = 1/g" $IF2
sed -i "s/WriteNum.*1/WriteNum = 2/g" $IF2
sed -i "s/sweeps$/1sweeps/g" $IF2
sed -i "s/21sweeps$/sweeps/g" $IF2
sed -i "s/nsweeps.*=.*[0-9]*/nsweeps = 20/g" $IF2


# change PBSrun.sh
sed -i "s/PBS -N .*/PBS -N $NAME/g" $PBS
sed -i "s/INPUTFILE=.*/INPUTFILE=$IF/g" $PBS
sed -i "s/INPUTFILE2=.*/INPUTFILE2=IF_${NAME}_2/g" $PBS
sed -i "s/APP=.*/APP=$APP/g" $PBS

# qsub 
qsub $PBS


