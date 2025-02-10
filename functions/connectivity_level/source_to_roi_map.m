function [ThetaSS] = source_to_roi_map(Thetajj,Scouts,Nsources)
        
    % Initialize R with zeros
    R = zeros(length(Scouts), Nsources);    
    for i = 1:length(Scouts)
        R(i, Scouts(i).Vertices) = 1;        
    end    
    % Normalize the rows of R
    for i = 1:size(R, 1)
        norm_i = norm(R(i, :));  % Calculate the norm of the row
        if norm_i > 0
            R(i, :) = R(i, :) / norm_i;  % Normalize the row
        end
    end
    % Get Connectivity in the rois
    ThetaSS = R*Thetajj*R';
end



