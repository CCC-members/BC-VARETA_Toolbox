#!/bin/bash 

matlab_workspace=$1
function=$2
idnode=$3
N=$4

echo ""
echo "-->> Running process in node: "$idnode" of "$N" nodes."
echo ""
cd $matlab_workspace
matlab -nodisplay -nodesktop -r "$function($idnode,$N);quit"

echo ""
#echo "-->> Process completed in node: "$idnode" name: "$node    
sleep 10
exit
