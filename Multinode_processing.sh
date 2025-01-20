#!/bin/bash 
tittle="------Executing multi process on Linux------"
echo ""
echo $tittle
echo ""
echo "-->> Starting process..."
username="username"
password="userpass"
function_workspace="/mnt/Develop/BC-VARETA_Toolbox"
script_path="`pwd`"
function_name="bcvareta"

if test -d "$function_workspace"; then
  declare -a nodes=("node1" "node2" "node3" "node4" "node5" "node6" "node7" "node8")
  
  N=${#nodes[@]}
  echo "Total of instances: "$N
  
  script_path="`pwd`"
  function_name="bcvareta"
  idnode=1
  if test -f "$script_path/move_to_node.sh"; then
    for node in "${nodes[@]}"; do
      (   
        echo "Open instance :"$idnode
        #if [[ $idnode -eq "7" ]] ; then
        #gnome-terminal -- "bash -c '$script_path/move_to_node.sh $script_path $function_workspace $function_name $node $idnode $N $password'"
        #gnome-terminal --geometry=25x2 --tab --title="$node" -- "bash -c '$script_path/move_to_node.sh $script_path $function_workspace $function_name $node $idnode $N $password'"
      gnome-terminal --geometry=25x2 --tab --title="$node" --command="bash -c '$script_path/move_to_node.sh $script_path $function_workspace $function_name $node $idnode $N $password'"
        #xterm -e "/home/$username/Actual_Process/move_to_node.sh" $workspace_path $function_name $node $idnode $N $password
        #fi
      ) &
      ((idnode=idnode+1))
      echo ""
      
     # allow to execute up to $N jobs in parallel
      if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
          # now there are $N jobs already running, so wait here for any job
          # to be finished so there is a place to start next one.
          wait
      fi
    done
  else
    echo "The file: 'move_to_node.sh' cannot be found in the current folder"
  fi
  
  # no more jobs to be started but wait for pending jobs
  # (all need to be finished)
  wait
else
  echo "The address: function_workspace='$function_workspace' is not a directory"
fi

echo "Process finished!!!"
