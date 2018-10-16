#!/bin/bash

LD=ladder_divide
LDN=ladder_d
# clean up 
echo "Enter ./$LD"
cd ./$LD
make clean
if [ -d PH* ]; then
    #echo "Please clean up PH_* files"
    make cleanPH
fi
cd ..
echo "Leave ./$LD"


# generate
for((i=2;i<=5;i++))
do
    if [ ! -d "ladder_d$i" ]
    then
        mkdir ladder_d$i

        echo "Enter ladder_d$i"
        cd $LDN$i
        echo "Copy file"
        cp -rf ../$LD/* ./
        ./ChangeLocalDimension.sh $i
        mv ladder_divide.cc ladder_d$i.cc
        sed -i "s/^APP=.*/APP=ladder_d$i/g" Makefile
        echo "Leave ladder_d$i"
        cd ..
    else
        echo "./ladder_d$i exist"
    fi
done
