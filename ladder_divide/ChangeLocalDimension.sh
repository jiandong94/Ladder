#!/bin/bash
if [ $# != 1 ]
then
    echo "Pls Enter ./ChangeLocalDim [dim]"
    exit 1
fi

if [[ $1 -gt 6 || $1 -lt 2 ]]
then
    echo "Local dimension should >=2 and <= 5"
    exit 1
fi

sed -i "s/hubbard_d[0-9]_divide.h/hubbard_d$1_divide.h/g" ladder_divide.cc
sed -i "s/HubbardD[0-9]Divide/HubbardD$1Divide/g" ladder_divide.cc 
sed -i "s/hubbard_d[0-9]_divide.h/hubbard_d$1_divide.h/g" observer.cc
sed -i "s/HubbardD[0-9]Divide/HubbardD$1Divide/g" observer.cc
