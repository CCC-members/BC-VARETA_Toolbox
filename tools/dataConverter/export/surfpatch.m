function [neig_indexes]=surfpatch(source_ind,sub_faces)
%% Initialization
index           = source_ind; %Indices of patch's elements
findex          = source_ind; %Indices of patch's frontier elements

findex_new  = [];
for cont = 1:length(findex)
    fpoint        = findex(cont);             %Pick point at the frontier 'fpoint'
    %% Search the neighbors of 'fpoint' out of the patche 'nfpoint'
    [row,col]      = find(sub_faces==fpoint);
    neig_fpoint    = sub_faces(row,:);
    neig_fpoint    = neig_fpoint(:);
    neig_fpoint    = setdiff(neig_fpoint,index);
    findex_new     = [findex_new; setdiff(neig_fpoint,findex_new)];
end
findex = findex_new;
neig_indexes  = [index; findex];

end