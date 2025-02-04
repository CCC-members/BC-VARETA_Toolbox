function [ThetaRR] = source_to_roi_map(Thetajj,Atlas,Nsources)
        
    % Initialize R with zeros
    R = zeros(length(Atlas.Scouts), Nsources);
    
    for i = 1:length(Atlas.Scouts)
        scout = Atlas.Scouts(i).Vertices;        
        % Calculate the mapping values for the current ROI
        for j = 1:Nsources
            R(i, j) = sum(scout == j);
        end
    end    
    % Normalize the rows of R
    for i = 1:size(R, 1)
        norm_i = norm(R(i, :));  % Calculate the norm of the row
        if norm_i > 0
            R(i, :) = R(i, :) / norm_i;  % Normalize the row
        end
    end

    % Get Connectivity in the rois
    ThetaRR = R*Thetajj*R';
end



