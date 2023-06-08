function subjects = multinode_subjects(subjects,idnode,total_node)
if(length(subjects)< idnode)
    warning("BC-V-->> There is not enough data to run in this node");
    disp(strcat("BC-V-->> Closing process in nodo: ",num2str(idnode),"."));
    return;
end
sub_count                       = fix(length(subjects)/total_node);
rest_sub                        = mod(length(subjects),total_node);
start_ind                       = idnode * sub_count - sub_count + 1;
% end_ind = idnode * sub_count;
if(start_ind>length(subjects))
    fprintf(2,strcat('\nBC-V-->> Error: The follow folder: \n'));
    disp(root_path);
    fprintf(2,strcat('BC-V-->> Error: Does not contain enough data to process on this node.\n'));
    fprintf(2,strcat('BC-V-->> Error: Or Does not contain any subject information file.\n'));
    disp("Please verify the configuration of the input data and start the process again.");
    return;
end
if(~isequal(rest_sub,0))
    if(idnode<=rest_sub)
        start_ind                                               = start_ind + idnode - 1;
        end_ind                                                 = start_ind + sub_count;
    else
        start_ind                                               = start_ind + rest_sub;
        end_ind                                                 = start_ind + sub_count - 1;
    end
else
    end_ind                                                     = start_ind + sub_count - 1;
end
if(sub_count == 0)
    start_ind                                                   = idnode;
    end_ind                                                     = idnode;
end
subjects                                                        = subjects(start_ind:end_ind);

end