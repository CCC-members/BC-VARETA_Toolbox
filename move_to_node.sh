#!/bin/bash 

script_path=$1
function_workspace=$2
function=$3
node=$4
idnode=$5
N=$6
password=$7


#echo "Login done in node: "$node
#sshpass -p $password ssh -tt $node << EOF
echo "-->> Input params"
echo "----------------------------------------------------------"
echo "Path: "$function_workspace
echo "Function: "$function
echo "Node: "$node
echo "Id node: "$idnode
echo "Total nodes: "$N
#echo "Password: "$password 
if test -f "$script_path/run_process.sh"; then
  cd $script_path
   ./run_process.sh $function_workspace $function $idnode $N
else
  echo "The file: 'run_process.sh' cannot be found in the current folder"
fi
exit
>/dev/null
EOF
exit
